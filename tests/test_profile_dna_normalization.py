"""Native/public allele identities require a verified biological scheme and locus."""

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import sqlite3
import sys
from types import SimpleNamespace

import pytest

from chronoclade import pasteur_profile_normalization as binding
from chronoclade.datasets import DatasetError
from chronoclade.pathogenwatch import content_hash
from chronoclade.profile_dna_normalization import COMPARISON_VERSION, ProfileDNANormalizer
from chronoclade.query_typing import QueryTypingError
from test_query_typing import indexed_scheme


@pytest.fixture
def caller_reference(tmp_path):
    loci = ['a', 'b', 'c']
    version = 'snapshot:' + 'a' * 64
    index = tmp_path / 'index'
    indexed_scheme(index, scheme='saureus_1', version=version, genes=loci)
    metadata = index / 'saureus_1_metadata.json'
    details = json.loads(metadata.read_text())
    details.update(host='pubmlst', hostPath='pubmlst_saureus_seqdef', schemeId='2')
    metadata.write_text(json.dumps(details))
    # The SQLite exact-hit index is authoritative; FASTAs may be representative subsets.
    digests = {locus: hashlib.sha1((locus + 'DNA').encode()).hexdigest() for locus in loci}
    with sqlite3.connect(index / 'allele.db') as connection:
        connection.execute('DELETE FROM alleles')
        connection.executemany('INSERT INTO alleles VALUES (?,?,?,0)',
                               [(digest, locus, 1) for locus, digest in digests.items()])
        # Reverse-complement lookup rows do not define the oriented allele identity.
        connection.execute('INSERT INTO alleles VALUES (?, ?, 1, 1)', ('f' * 40, 'a'))
    config = tmp_path / 'typing.json'
    config.write_text(json.dumps({'ready': True, 'organisms': {'Staphylococcus aureus': {
        'mode': 'profiles-only', 'cgmlst': {
            'command': [sys.executable], 'scheme': 'saureus_1', 'scheme_version': version,
            'tool_version': 'fixture-1', 'index_dir': str(index),
        },
    }}}))
    return config, index, digests, loci


def row(ident='query', **changes):
    value = {'sample_id': ident, 'species': 'Staphylococcus aureus',
             'cgmlst_scheme': 'pathogenwatch:cgmlst-1280-2',
             'cgmlst_scheme_version': 'cgmlst-1280-2',
             'cgmlst_profile': {'a': '1'}, 'cgmlst_loci': ['a'],
             'cgmlst_locus_universe_complete': False}
    return dict(value, **changes)


def test_known_numbers_and_hashes_share_locus_specific_dna_without_assigner(caller_reference, tmp_path):
    config, _, digests, loci = caller_reference
    original = [row('native', cgmlst_scheme='saureus_1', cgmlst_scheme_version='snapshot:new'),
                row('public', cgmlst_profile={'a': digests['a'].upper()},
                    source_genome_id='exact-provider-identity')]
    before = deepcopy(original)
    normalizer = ProfileDNANormalizer(config, tmp_path / 'normalized')
    output = normalizer.normalize(original)
    assert original == before
    assert [r['cgmlst_profile'] for r in output] == [{'a': digests['a']}] * 2
    assert all(r['cgmlst_scheme'] == 'pubmlst:pubmlst_saureus_seqdef:2' for r in output)
    assert all(r['cgmlst_scheme_version'] == COMPARISON_VERSION for r in output)
    assert all(r['cgmlst_loci'] == loci and r['cgmlst_locus_universe_complete'] for r in output)
    assert output[1]['source_genome_id'] == 'exact-provider-identity'
    assert all('hiercc_codes' not in r and 'cglin_raw' not in r for r in output)
    evidence = json.loads(Path(normalizer.audits[0]['path']).read_text())
    assert evidence['reference']['typing_mode'] == 'profiles-only'
    assert evidence['records'][0]['original_typing'] == {
        k: v for k, v in original[0].items() if k.startswith('cgmlst_')}
    assert normalizer.audits[0]['unsupported_calls'] == 0


def test_disjoint_observed_loci_keep_complete_panel_and_no_intersection(caller_reference, tmp_path):
    config, _, digests, loci = caller_reference
    output = ProfileDNANormalizer(config, tmp_path / 'normalized').normalize([
        row('q', cgmlst_profile={'a': '1'}), row('context', cgmlst_profile={'b': '1'}),
    ])
    assert output[0]['cgmlst_profile'] == {'a': digests['a']}
    assert output[1]['cgmlst_profile'] == {'b': digests['b']}
    assert all(r['cgmlst_loci'] == loci for r in output)
    assert output[0]['cgmlst_missing_loci'] == ['b', 'c']
    assert output[1]['cgmlst_missing_loci'] == ['a', 'c']
    assert output[0]['cgmlst_database_sha256'] == output[1]['cgmlst_database_sha256']


def test_unavailable_numbers_become_auditable_missing_and_novel_hash_remains(caller_reference, tmp_path):
    normalizer = ProfileDNANormalizer(caller_reference[0], tmp_path / 'normalized')
    output = normalizer.normalize([row(cgmlst_profile={'a': '999'},
                                       cgmlst_novel_alleles={'b': 'B' * 40},
                                       cgmlst_ambiguous_loci=['c'])])[0]
    assert output['cgmlst_profile'] == {'b': 'b' * 40}
    assert output['cgmlst_novel_alleles'] == {}
    assert output['cgmlst_missing_loci'] == ['a']
    assert output['cgmlst_ambiguous_loci'] == ['c']
    evidence = json.loads(Path(normalizer.audits[0]['path']).read_text())
    assert evidence['records'][0]['unsupported_calls'] == [{
        'locus': 'a', 'allele': '999',
        'reason': 'numbered_allele_unavailable_in_verified_reference',
    }]
    assert normalizer.audits[0]['unsupported_calls'] == 1


@pytest.mark.parametrize('changes', [
    {'species': 'Escherichia coli'},
    {'cgmlst_scheme': 'pathogenwatch:cgmlst-1280-99'},
    {'cgmlst_scheme': 'pathogenwatch:cgmlst-573-2'},
    {'cgmlst_profile': {'outside_panel': '1'}},
    {'cgmlst_ambiguous_loci': ['outside_panel']},
    {'cgmlst_loci': ['outside_panel']},
    {'cgmlst_locus_universe_complete': True},
    {'cgmlst_novel_alleles': {'a': 'f' * 40}},
    {'cgmlst_profile': {'a': 'partial~1'}},
])
def test_incompatible_source_species_panel_and_calls_rejected(caller_reference, tmp_path, changes):
    normalizer = ProfileDNANormalizer(caller_reference[0], tmp_path / 'normalized')
    with pytest.raises(DatasetError):
        normalizer.normalize([row(**changes)])


@pytest.mark.parametrize('field,value', [
    ('host', 'pasteur'), ('hostPath', 'different_seqdef'), ('schemeId', '1'),
])
def test_unverified_biological_scheme_metadata_rejected(caller_reference, tmp_path, field, value):
    config, index, _, _ = caller_reference
    metadata = index / 'saureus_1_metadata.json'
    changed = json.loads(metadata.read_text())
    changed[field] = value
    metadata.write_text(json.dumps(changed))
    normalizer = ProfileDNANormalizer(config, tmp_path / 'normalized')
    with pytest.raises(DatasetError, match='verified public biological scheme'):
        normalizer.normalize([row()])


@pytest.mark.parametrize('change', ['conflict', 'missing_locus', 'invalid_hash', 'invalid_number'])
def test_allele_index_definitions_must_be_complete_and_unambiguous(caller_reference, tmp_path, change):
    config, index, _, _ = caller_reference
    with sqlite3.connect(index / 'allele.db') as connection:
        if change == 'conflict':
            connection.execute('INSERT INTO alleles VALUES (?, ?, 1, 0)', ('f' * 40, 'a'))
        elif change == 'missing_locus':
            connection.execute("DELETE FROM alleles WHERE gene='c'")
        elif change == 'invalid_hash':
            connection.execute("UPDATE alleles SET hash='invalid' WHERE gene='a'")
        else:
            connection.execute("UPDATE alleles SET st=0 WHERE gene='a'")
    normalizer = ProfileDNANormalizer(config, tmp_path / 'normalized')
    with pytest.raises(DatasetError):
        normalizer.normalize([row()])


def test_caller_only_descriptor_auto_binding_and_resolution(caller_reference):
    config, index, _, loci = caller_reference
    descriptor = binding.reference_descriptor(config)
    assert descriptor == binding.reference_descriptor(config, profiles_only=True)
    assert descriptor['typing_mode'] == 'profiles-only'
    species = descriptor['resources']['Staphylococcus aureus']
    assert set(species) == {'caller'}
    assert species['caller']['canonical_loci_sha256'] == content_hash(loci)
    assert species['caller']['indexed_loci_sha256'] == content_hash(loci)
    assert Path(species['caller']['resources']['allele.db']['path']) == index / 'allele.db'
    dataset = SimpleNamespace(parameters={'native_typing_reference': descriptor})
    assert binding.resolve_reference(dataset) == config
    with pytest.raises(DatasetError, match='requires full caller and lineage readiness'):
        binding.reference_descriptor(config, profiles_only=False)


@pytest.mark.parametrize('resource', ['config', 'allele.db', 'metadata', 'details'])
def test_saved_caller_only_reference_rejects_changed_content(caller_reference, resource):
    config, index, _, _ = caller_reference
    descriptor = binding.reference_descriptor(config)
    dataset = SimpleNamespace(parameters={'native_typing_reference': descriptor})
    if resource == 'config':
        config.write_text(config.read_text() + '\n')
    elif resource == 'allele.db':
        with sqlite3.connect(index / 'allele.db') as connection:
            connection.execute('INSERT INTO alleles VALUES (?, ?, 2, 0)', ('e' * 40, 'a'))
    elif resource == 'metadata':
        p = index / 'saureus_1_metadata.json'
        changed = json.loads(p.read_text())
        changed['evidence_note'] = 'changed bytes with unchanged biological panel'
        p.write_text(json.dumps(changed))
    else:
        p = index / 'metadata.json.gz'
        p.write_bytes(p.read_bytes() + b'changed')
    with pytest.raises(DatasetError, match='changed|unavailable'):
        binding.resolve_reference(dataset)


@pytest.mark.parametrize('suffix', ['-wal', '-journal'])
def test_caller_only_descriptor_rejects_uncheckpointed_sqlite_writes(caller_reference, suffix):
    config, index, _, _ = caller_reference
    (index / ('allele.db' + suffix)).write_bytes(b'uncheckpointed write')
    with pytest.raises((DatasetError, QueryTypingError), match='uncheckpointed SQLite writes|unreadable'):
        binding.reference_descriptor(config)


def test_caller_only_descriptor_hash_tampering_is_rejected(caller_reference):
    descriptor = binding.reference_descriptor(caller_reference[0])
    descriptor['typing_mode'] = 'full'
    with pytest.raises(DatasetError, match='legacy/unverified'):
        binding.resolve_reference(SimpleNamespace(parameters={'native_typing_reference': descriptor}))


def test_all_numeric_forty_character_dna_hash_is_preserved(caller_reference, tmp_path):
    token = '1' * 40
    normalizer = ProfileDNANormalizer(caller_reference[0], tmp_path / 'normalized')
    output = normalizer.normalize([row(cgmlst_profile={'a': token})])[0]
    assert output['cgmlst_profile'] == {'a': token}
    assert normalizer.audits[0]['unsupported_calls'] == 0
