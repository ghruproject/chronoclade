"""Reference-wise Verticall pairwise/mask workflow with exact cohort retention."""

import csv
import json
from pathlib import Path
import re

from Bio import SeqIO
import numpy as np

from chronoclade.artifacts import file_sha256
from chronoclade.errors import WorkflowError


def repair_assembly(source, target):
    """Snapshot called segments, splitting at ambiguity without inventing bases."""
    original_bases = kept_bases = segments = 0
    with Path(target).open('w') as stream:
        for index, record in enumerate(SeqIO.parse(source, 'fasta')):
            sequence = str(record.seq).upper()
            original_bases += len(sequence)
            for part in re.finditer('[ACGT]+', sequence):
                stream.write(f'>contig_{index}_segment_{segments}\n{part[0]}\n')
                kept_bases += len(part[0])
                segments += 1
    if not segments:
        raise WorkflowError(f'Assembly has no called segments for Verticall: {source}')
    return {'original_sha256': file_sha256(Path(source)), 'snapshot_sha256': file_sha256(Path(target)),
            'original_bases': original_bases, 'called_bases': kept_bases,
            'ambiguous_bases_removed': original_bases - kept_bases, 'segments': segments}


def read_regions(path, aliases, length):
    """Validate the first native result per reference/sample, including full coverage."""
    selected = {}
    secondary = 0
    with Path(path).open() as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        required = {'assembly_a', 'assembly_b', 'assembly_a_vertical_regions',
                    'assembly_a_horizontal_regions', 'assembly_a_unaligned_regions'}
        if not required <= set(reader.fieldnames or []):
            raise WorkflowError('Verticall pairwise TSV lacks required reference-region columns')
        for row in reader:
            if row['assembly_a'] != 'reference':
                raise WorkflowError('Unexpected reference in Verticall pairwise result')
            alias = row['assembly_b']
            if alias not in aliases:
                raise WorkflowError('Unexpected sample in Verticall pairwise result')
            if alias in selected:
                secondary += 1
                continue
            regions = {}
            for kind in ['vertical', 'horizontal', 'unaligned']:
                regions[kind] = []
                for entry in filter(None, row[f'assembly_a_{kind}_regions'].split(',')):
                    match = re.fullmatch(r'reference:(\d+)-(\d+)', entry)
                    if match is None:
                        raise WorkflowError('Malformed Verticall reference-region coordinates')
                    begin, end = map(int, match.groups())
                    if not 0 <= begin < end <= length:
                        raise WorkflowError('Verticall region is outside its reference sequence')
                    regions[kind].append((begin, end))
            position = 0
            for begin, end in sorted(r for values in regions.values() for r in values):
                if begin != position:
                    raise WorkflowError('Verticall regions overlap or leave unclassified reference bases')
                position = end
            if position != length:
                raise WorkflowError('Verticall regions do not cover the complete reference sequence')
            selected[alias] = regions
    if set(selected) != set(aliases):
        missing = ', '.join(aliases[a] for a in set(aliases) - set(selected))
        raise WorkflowError('Verticall omitted reference comparisons; no genomes dropped: ' + missing)
    return selected, secondary


def run_verticall_mask(files, members, reference, *, threads, force):
    from chronoclade.native_pipeline import _run_command

    directory = files.directory / 'verticall'
    assemblies = directory / 'assemblies'
    assemblies.mkdir(parents=True, exist_ok=True)
    aliases = {f'sample{i:06d}': member.sample_id for i, member in enumerate(members)}
    raw = {r.id: str(r.seq).upper() for r in SeqIO.parse(files.alignment, 'fasta')}
    if len(raw) != len(members) or set(raw) != {m.sample_id for m in members} or len({len(s) for s in raw.values()}) != 1:
        raise WorkflowError('Verticall input alignment must contain the exact selected cohort')
    records = list(SeqIO.parse(reference, 'fasta'))
    length = len(next(iter(raw.values())))
    if not records or sum(len(r) for r in records) != length:
        raise WorkflowError('Verticall reference boundaries disagree with the mapped alignment')
    if any(set(str(r.seq).upper()) - set('ACGT') for r in records):
        raise WorkflowError('Verticall needs an unambiguous reference; choose another --reference-sample-id')
    repairs = []
    paths = []
    for alias, member in zip(aliases, members, strict=True):
        target = assemblies / (alias + '.fasta')
        receipt = repair_assembly(member.assembly, target)
        repairs.append({'sample_id': member.sample_id, 'alias': alias, **receipt})
        paths.append(target)
    (directory / 'assembly_repair.json').write_text(json.dumps(repairs, indent=2) + '\n')
    output = {ident: [] for ident in raw}
    horizontal = np.zeros(length, dtype=bool)
    unaligned = np.zeros(length, dtype=bool)
    region_rows, record_summaries = [], []
    offset = secondary = 0
    for index, record in enumerate(records):
        work = directory / f'record-{index:04d}'
        work.mkdir(exist_ok=True)
        size = len(record)
        ref_path = work / 'reference.fasta'
        ref_path.write_text('>reference\n' + str(record.seq).upper() + '\n')
        alignment = work / 'alignment.fasta'
        with alignment.open('w') as stream:
            stream.write(ref_path.read_text())
            for alias, ident in aliases.items():
                stream.write(f'>{alias}\n{raw[ident][offset:offset+size]}\n')
        pairwise = work / 'pairwise.tsv'
        _run_command(['verticall', 'pairwise', '-i', str(assemblies), '-o', str(pairwise),
                      '-r', str(ref_path), '-t', str(threads)],
                     log=work / 'pairwise.log', expected=pairwise, force=force, cwd=work,
                     inputs=(ref_path, *paths))
        regions, extra = read_regions(pairwise, aliases, size)
        secondary += extra
        masked = work / 'masked.fasta'
        _run_command(['verticall', 'mask', '-i', str(pairwise), '-a', str(alignment), '-o', str(masked),
                      '--reference', 'reference', '--exclude_reference', '--multi', 'first'],
                     log=work / 'mask.log', expected=masked, force=force, cwd=work,
                     inputs=(pairwise, alignment))
        results = list(SeqIO.parse(masked, 'fasta'))
        if len(results) != len(aliases) or {r.id for r in results} != set(aliases):
            raise WorkflowError('Verticall mask omitted or duplicated selected genomes')
        for result in results:
            alias = result.id
            expected = bytearray(raw[aliases[alias]][offset:offset+size].encode('ascii'))
            for kind, values in regions[alias].items():
                for begin, end in values:
                    region_rows.append({'sample_id': aliases[alias], 'reference_record': record.id,
                                        'classification': kind, 'record_start_1based': begin + 1,
                                        'record_end_1based': end, 'alignment_start_1based': offset + begin + 1,
                                        'alignment_end_1based': offset + end})
                    if kind == 'horizontal':
                        expected[begin:end] = b'N' * (end - begin)
                        horizontal[offset+begin:offset+end] = True
                    elif kind == 'unaligned':
                        expected[begin:end] = b'-' * (end - begin)
                        unaligned[offset+begin:offset+end] = True
            if bytes(result.seq).upper() != bytes(expected):
                raise WorkflowError('Verticall masked bases disagree with its saved reference-region calls')
            output[aliases[alias]].append(str(result.seq).upper())
        record_summaries.append({'record': record.id, 'bases': size, 'samples': len(results),
                                 'secondary_results': extra, 'alignment_offset': offset})
        offset += size
    combined = {ident: ''.join(parts) for ident, parts in output.items()}
    incomplete = np.zeros(length, dtype=bool)
    retain = np.ones(length, dtype=bool)
    for ident, seq in combined.items():
        a = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)
        retain &= (a == 65) | (a == 67) | (a == 71) | (a == 84)
        original = np.frombuffer(raw[ident].encode('ascii'), dtype=np.uint8)
        incomplete |= ~((original == 65) | (original == 67) | (original == 71) | (original == 84))
    if not retain.any():
        raise WorkflowError('Verticall masking leaves no complete shared alignment sites')
    masked_full = files.directory / 'verticall.masked.fasta'
    filtered = files.directory / 'verticall.filtered.fasta'
    with masked_full.open('w') as full, filtered.open('w') as clean:
        for member in members:
            seq = combined[member.sample_id]
            full.write(f'>{member.sample_id}\n{seq}\n')
            kept = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)[retain].tobytes().decode('ascii')
            clean.write(f'>{member.sample_id}\n{kept}\n')
    with (directory / 'regions.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(region_rows[0]))
        writer.writeheader()
        writer.writerows(region_rows)
    summary = {'method': 'verticall', 'correction_applied': True,
               'workflow': 'reference-to-each-assembly pairwise; native mask; complete shared columns; fresh IQ-TREE',
               'alignment_policy': 'Native sample-specific mask retained separately; downstream alignment contains only columns called A/C/G/T in every selected genome',
               'reference_policy': 'Each reference FASTA record processed independently, without artificial joins',
               'multi_result_policy': 'first native result', 'secondary_results': secondary,
               'reference_records': record_summaries, 'samples': len(members),
               'repaired_assemblies': sum(r['ambiguous_bases_removed'] > 0 for r in repairs),
               'ambiguous_assembly_bases_removed': sum(r['ambiguous_bases_removed'] for r in repairs),
               'alignment_sites': length, 'recombination_removed_sites': int(horizontal.sum()),
               'unaligned_only_removed_sites': int((unaligned & ~horizontal).sum()),
               'incomplete_only_removed_sites': int((incomplete & ~(horizontal | unaligned)).sum()),
               'retained_clonal_sites': int(retain.sum()),
               'limitation': 'Reference-relative horizontal/unaligned classifications are not branch-specific import calls. The common complete-site step is an additional downstream policy; native masked alignments are also retained.'}
    assert summary['retained_clonal_sites'] + summary['recombination_removed_sites'] + summary['unaligned_only_removed_sites'] + summary['incomplete_only_removed_sites'] == length
    (files.directory / 'verticall_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    return filtered, summary


def report_html(summary):
    from html import escape

    return ('<section><h2>Recombination filtering: Verticall</h2><p>'
            f"{summary['samples']:,} genomes retained; {len(summary['reference_records'])} reference records processed independently. "
            f"Horizontal masks cover {summary['recombination_removed_sites']:,} mapped columns; "
            f"{summary['retained_clonal_sites']:,} fully called columns remain for a rebuilt IQ-TREE phylogeny.</p><p>"
            + escape(summary['alignment_policy']) + '</p><p>' + escape(summary['limitation']) + '</p></section>')
