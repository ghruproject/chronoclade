"""Fast, coordinate-preserving PHIPack profiling and exploratory window masking."""

from concurrent.futures import ThreadPoolExecutor
import csv
import json
import math
from pathlib import Path
import re

from Bio import SeqIO
import numpy as np

from chronoclade.errors import WorkflowError


def validate_options(method, window, step, threshold):
    if method not in {'phipack', 'clonalframeml', 'verticall'}:
        raise WorkflowError('Recombination method must be phipack, clonalframeml or verticall')
    if type(window) is not int or type(step) is not int or not 200 <= window or not 1 <= step <= window:
        raise WorkflowError('PHIPack window must be at least 200 bases; step must be between 1 and window')
    if not math.isfinite(threshold) or not 0 < threshold <= 1:
        raise WorkflowError('PHIPack FDR threshold must be greater than zero and at most one')


def read_profile(csv_path, log_path, *, starts, window):
    """Never trust CSV values for windows whose normal approximation was unavailable.

    Profile leaves its CSV value uninitialised when there are too few informative
    sites. Its console record is the authority for whether a test was performed;
    numeric CSV values provide precision only after that check.
    """
    log = Path(log_path).read_text()
    blocks = re.findall(r'Checking for recombination at (\d+) between (\d+) and (\d+)\s+(.*?)(?=Checking for recombination at|Number of tests performed|\Z)', log, re.S)
    with Path(csv_path).open() as stream:
        values = list(csv.reader(stream))
    if len(blocks) != len(starts) or len(values) != len(starts):
        raise WorkflowError('PHIPack output does not cover the exact requested windows')
    result = []
    for start, block, value in zip(starts, blocks, values, strict=True):
        center, begin, end, text = block
        if (int(begin), int(end), int(center)) != (start, start + window, start + window // 2):
            raise WorkflowError('PHIPack window coordinates disagree with the planned alignment')
        if len(value) != 2 or int(value[0]) != int(center):
            raise WorkflowError('Malformed PHIPack profile coordinates')
        status = re.search(r'PHI \(Normal\):\s*(\S+)', text)
        if status is None:
            raise WorkflowError('PHIPack did not report a test status for a requested window')
        reported = status[1].lower()
        if reported == '--' or 'nan' in reported or 'inf' in reported:
            p_value = None
        else:
            try:
                p_value = float(value[1])
                rounded = float(reported)
            except ValueError as error:
                raise WorkflowError('Malformed PHIPack probability') from error
            if not math.isfinite(p_value) or not 0 <= p_value <= 1 or not math.isclose(p_value, rounded, rel_tol=0.006, abs_tol=1e-300):
                raise WorkflowError('PHIPack CSV and console probabilities disagree')
        result.append({'start': start, 'end': start + window, 'p_value': p_value})
    return result


def adjust_fdr(windows):
    """Benjamini–Hochberg adjustment over all finite tests in this alignment."""
    tested = sorted((r for r in windows if r['p_value'] is not None), key=lambda r: r['p_value'])
    upper = 1.0
    for rank in range(len(tested), 0, -1):
        row = tested[rank - 1]
        upper = min(upper, row['p_value'] * len(tested) / rank)
        row['q_value'] = upper
    for row in windows:
        row.setdefault('q_value', None)


def run_profile_mask(files, reference, *, threads, force, window=10000, step=5000, threshold=0.05):
    from chronoclade.native_pipeline import _run_command

    validate_options('phipack', window, step, threshold)
    directory = files.directory / 'phipack'
    directory.mkdir(exist_ok=True)
    sequences = {r.id: str(r.seq).upper() for r in SeqIO.parse(files.alignment, 'fasta')}
    if not sequences or len({len(s) for s in sequences.values()}) != 1:
        raise WorkflowError('PHIPack requires an equal-length selected alignment')
    length = len(next(iter(sequences.values())))
    records = list(SeqIO.parse(reference, 'fasta'))
    if sum(len(r) for r in records) != length:
        raise WorkflowError('PHIPack reference boundaries disagree with mapped alignment length')
    specs = []
    offset = 0
    short_bases = 0
    for record in records:
        size = len(record)
        actual_window = min(window, size)
        if actual_window < 200:
            short_bases += size
            offset += size
            continue
        starts = list(range(0, size - actual_window + 1, step))
        # Batch regular overlapping windows, retaining physical distances and
        # never crossing a reference-record boundary. Include the final tail.
        groups = [starts[i:i + 32] for i in range(0, len(starts), 32)]
        if starts[-1] != size - actual_window:
            groups.append([size - actual_window])
        for group in groups:
            specs.append((offset + group[0], offset + group[-1] + actual_window,
                          actual_window, record.id, group[0], len(group)))
        offset += size

    def profile(item):
        index, (begin, end, width, record, local_start, count) = item
        work = directory / f'block-{index:05d}'
        work.mkdir(exist_ok=True)
        source = work / 'alignment.fasta'
        with source.open('w') as stream:
            for ident, seq in sequences.items():
                stream.write(f'>{ident}\n{seq[begin:end]}\n')
        log = work / 'profile.log'
        _run_command(['Profile', '-f', str(source), '-n', str(width), '-m', str(step), '-w', '100'],
                     log=log, expected=work / 'Profile.csv', force=force, cwd=work, inputs=(source,))
        found = read_profile(work / 'Profile.csv', log, starts=[i * step for i in range(count)], window=width)
        for row in found:
            row.update(start=row['start'] + begin, end=row['end'] + begin,
                       reference_record=record, record_offset=begin - local_start)
        return found

    with ThreadPoolExecutor(max_workers=threads) as executor:
        windows = [row for block in executor.map(profile, enumerate(specs)) for row in block]
    adjust_fdr(windows)
    mask = np.zeros(length, dtype=bool)
    intervals = []
    for row in windows:
        row['masked'] = row['q_value'] is not None and row['q_value'] <= threshold
        if row['masked']:
            mask[row['start']:row['end']] = True
            intervals.append((row['start'], row['end'], row['reference_record']))
    merged = []
    for begin, end, record in sorted(intervals):
        if merged and merged[-1][2] == record and begin <= merged[-1][1]:
            merged[-1][1] = max(end, merged[-1][1])
        else:
            merged.append([begin, end, record])
    incomplete = np.zeros(length, dtype=bool)
    for seq in sequences.values():
        bases = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)
        incomplete |= ~((bases == 65) | (bases == 67) | (bases == 71) | (bases == 84))
    retain = ~(mask | incomplete)
    if not retain.any():
        raise WorkflowError('PHIPack filtering leaves no complete alignment sites; inspect the profile evidence')
    output = files.directory / 'phipack.filtered.fasta'
    with output.open('w') as stream:
        for ident, seq in sequences.items():
            filtered = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)[retain].tobytes().decode('ascii')
            stream.write(f'>{ident}\n{filtered}\n')
    with (directory / 'windows.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['reference_record', 'alignment_start_1based', 'alignment_end_1based', 'record_start_1based', 'record_end_1based', 'p_value', 'q_value', 'masked'])
        for row in windows:
            writer.writerow([row['reference_record'], row['start'] + 1, row['end'],
                             row['start'] - row['record_offset'] + 1, row['end'] - row['record_offset'],
                             row['p_value'], row['q_value'], row['masked']])
    summary = {'method': 'phipack', 'policy': 'Profile windows; Benjamini-Hochberg FDR; union of significant whole windows',
               'correction_applied': True, 'window_bases': window, 'step_bases': step,
               'phi_window_bases': 100, 'fdr_threshold': threshold,
               'coordinate_system': '1-based closed mapped alignment; windows do not cross reference records',
               'tested_windows': sum(r['p_value'] is not None for r in windows),
               'untested_windows': sum(r['p_value'] is None for r in windows),
               'unscanned_short_record_bases': short_bases,
               'significant_windows': sum(r['masked'] for r in windows),
               'alignment_sites': length, 'recombination_removed_sites': int(mask.sum()),
               'incomplete_only_removed_sites': int((incomplete & ~mask).sum()),
               'retained_clonal_sites': int(retain.sum()),
               'intervals': [{'start': a + 1, 'end': b, 'reference_record': r} for a, b, r in merged],
               'limitation': 'Exploratory window masking, not branch-specific import inference. Untested windows are retained; no significant window does not establish absence of recombination.'}
    (files.directory / 'phipack_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    return output, summary


def report_html(summary):
    from html import escape

    return ('<section><h2>Recombination filtering: PHIPack Profile</h2><p>'
            f"{summary['significant_windows']:,} significant windows; {summary['tested_windows']:,} tested; "
            f"{summary['untested_windows']:,} untested. Masked {summary['recombination_removed_sites']:,} "
            f"of {summary['alignment_sites']:,} mapped columns at FDR ≤ {summary['fdr_threshold']}. "
            'The tree was rebuilt on the retained complete alignment.</p><p>'
            + escape(summary['limitation']) + '</p></section>')
