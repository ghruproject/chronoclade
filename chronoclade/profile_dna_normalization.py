"""Compare native and provider profiles through verified locus-specific DNA IDs.

Numbered alleles are translated only through an installed, content-bound caller
index for the exact advertised biological scheme. No lineage assigner is used.
"""

from copy import deepcopy
import json
from pathlib import Path
import re
import sqlite3

from chronoclade.artifacts import file_sha256, write_json
from chronoclade.datasets import DatasetError
from chronoclade.pathogenwatch import content_hash
from chronoclade.query_typing import load_typing_config, normalise_typing_species
from chronoclade.pasteur_profile_normalization import reference_descriptor, _assert_checkpointed_file

COMPARISON_VERSION = 'locus-oriented-dna-sha1-v1'
ROUTES = {
    'klebsiella pneumoniae': ('klebsiella_1', 'pathogenwatch:cgmlst-573-2',
                            'pasteur', 'pubmlst_klebsiella_seqdef', '18'),
    'staphylococcus aureus': ('saureus_1', 'pathogenwatch:cgmlst-1280-2',
                             'pubmlst', 'pubmlst_saureus_seqdef', '2'),
    'escherichia coli': ('ecoli_1', 'pathogenwatch:cgmlst-562-1',
                         'enterobase', 'Escherichia.cgMLSTv1', 'Escherichia.cgMLSTv1'),
}


class ProfileDNANormalizer:
    """Preserve complete reference panels and expose every unsupported call."""

    def __init__(self, config_path, output):
        self.config_path = Path(config_path).expanduser().resolve()
        self.reference = reference_descriptor(self.config_path, profiles_only=True)
        config = load_typing_config(self.config_path, profiles_only=True)
        self.callers = {normalise_typing_species(k): v['cgmlst'] for k, v in config['organisms'].items()}
        self.panels, self.audits = {}, []
        self.output = Path(output)
        self.output.mkdir(parents=True, exist_ok=True)

    def _panel(self, species):
        if species in self.panels:
            return self.panels[species]
        caller, route = self.callers.get(species), ROUTES.get(species)
        if not caller or not route:
            raise DatasetError(f'No verified native/public cgMLST scheme bridge for {species}')
        native, public, host, host_path, scheme_id = route
        metadata = json.loads((caller['index_dir'] / f'{native}_metadata.json').read_text())
        if caller['scheme'] != native or (
                metadata.get('host'), metadata.get('hostPath'), str(metadata.get('schemeId'))
        ) != (host, host_path, scheme_id):
            raise DatasetError(f'Native index is not the verified public biological scheme for {species}')
        loci = caller['index_loci']
        locus_set = set(loci)
        canonical = f'{host}:{host_path}:{scheme_id}'
        database = _assert_checkpointed_file(caller['index_dir'] / 'allele.db')
        sequences = {gene: {} for gene in loci}
        try:
            with sqlite3.connect(f'file:{database}?mode=ro&immutable=1', uri=True) as connection:
                for digest, gene, number in connection.execute('SELECT hash, gene, st FROM alleles WHERE reverse=0'):
                    if gene not in locus_set or not re.fullmatch(r'[0-9a-f]{40}', str(digest)):
                        raise DatasetError('Native allele index has an invalid locus or DNA hash')
                    token = str(number)
                    if not re.fullmatch(r'[1-9][0-9]*', token):
                        raise DatasetError('Native allele index has an invalid allele number')
                    previous = sequences[gene].get(token)
                    if previous is not None and previous != digest:
                        raise DatasetError('Conflicting locus-specific allele DNA definitions')
                    sequences[gene][token] = digest
        except sqlite3.Error as error:
            raise DatasetError('Native allele DNA index is unreadable') from error
        if any(not values for values in sequences.values()):
            raise DatasetError('Native allele index does not cover the complete reference locus panel')
        panel_digest = content_hash({'scheme': canonical, 'policy': COMPARISON_VERSION, 'loci': sorted(loci)})
        panel = {'loci': loci, 'locus_set': locus_set, 'sequences': sequences,
                 'scheme': canonical, 'digest': panel_digest,
                 'sources': {native, public, canonical}, 'allele_index_sha256': file_sha256(database)}
        self.panels[species] = panel
        return panel

    def normalize(self, records):
        adapted, evidence = [], []
        for original in records:
            row = deepcopy(original)
            species = normalise_typing_species(row.get('species'))
            panel = self._panel(species)
            if row.get('cgmlst_scheme') not in panel['sources']:
                raise DatasetError(f"Unverified cgMLST source namespace for {row['sample_id']}")
            regular, novel = row.get('cgmlst_profile') or {}, row.get('cgmlst_novel_alleles') or {}
            if set(regular) & set(novel):
                raise DatasetError('Known and novel calls overlap during DNA normalization')
            calls = dict(regular, **novel)
            ambiguous = set(row.get('cgmlst_ambiguous_loci') or [])
            if (set(calls) | ambiguous | set(row.get('cgmlst_loci') or [])) - panel['locus_set']:
                raise DatasetError('Source profile contains loci outside the verified biological scheme')
            if (row.get('cgmlst_locus_universe_complete') is True
                    and set(row.get('cgmlst_loci') or []) != panel['locus_set']):
                raise DatasetError('Declared complete source panel differs from the native reference')
            normalized, unsupported = {}, []
            for locus, value in sorted(calls.items()):
                if locus in ambiguous or value in (None, '', 0, '0', '-', '?'):
                    continue
                token = str(value)
                if re.fullmatch(r'[0-9a-fA-F]{40}', token):
                    digest = token.lower()
                elif re.fullmatch(r'[1-9][0-9]*', token):
                    digest = panel['sequences'][locus].get(token)
                    if digest is None:
                        unsupported.append({'locus': locus, 'allele': token,
                                            'reason': 'numbered_allele_unavailable_in_verified_reference'})
                        continue
                else:
                    raise DatasetError(f"Invalid allele token for {row['sample_id']} at {locus}")
                normalized[locus] = digest
            evidence.append({'sample_id': row['sample_id'],
                             'original_typing': {k: v for k, v in original.items() if k.startswith('cgmlst_')},
                             'unsupported_calls': unsupported, 'normalized_called_loci': len(normalized),
                             'reference_panel_sha256': panel['digest'],
                             'allele_index_sha256': panel['allele_index_sha256']})
            row.update(cgmlst_scheme=panel['scheme'], cgmlst_scheme_version=COMPARISON_VERSION,
                       cgmlst_database_version=None, cgmlst_database_sha256=panel['digest'],
                       cgmlst_profile=normalized, cgmlst_novel_alleles={},
                       cgmlst_loci=panel['loci'], cgmlst_locus_universe_complete=True,
                       cgmlst_missing_loci=sorted(panel['locus_set'] - normalized.keys() - ambiguous),
                       cgmlst_status='resolved' if normalized else 'unassigned')
            adapted.append(row)
        if evidence:
            path = self.output / (content_hash(evidence) + '-normalization.json')
            write_json(path, {'reference': self.reference, 'records': evidence,
                              'policy': 'Verified locus-specific DNA identities; full reference panel; '
                                        'unavailable numbered calls remain missing; lineage labels optional'})
            self.audits.append({'path': str(path), 'records': len(adapted),
                                'unsupported_calls': sum(len(r['unsupported_calls']) for r in evidence)})
        return adapted
