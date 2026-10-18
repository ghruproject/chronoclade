"""Exact classical MLST calls for the provider's primary species/ST search.

This is independent of cgMLST and of LIN/HierCC. Unknown or partial profiles
never acquire an inferred ST from a similar allele or neighbouring genome.
"""

import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

from chronoclade.artifacts import file_sha256, write_json
from chronoclade.progress import progress_task

# Explicitly matched to Pathogenwatch's primary `mlst`, not `mlst2`, namespace.
# In particular, tseemann/mlst's scheme named `ecoli` is NOT the Achtman scheme.
PRIMARY_SCHEMES = {
    'klebsiella pneumoniae': ('klebsiella', {'gapA', 'infB', 'mdh', 'pgi', 'phoE', 'rpoB', 'tonB'}),
    'escherichia coli': ('ecoli_achtman_4', {'adk', 'fumC', 'gyrB', 'icd', 'mdh', 'purA', 'recA'}),
    'staphylococcus aureus': ('saureus', {'arcC', 'aroE', 'glpF', 'gmk', 'pta', 'tpi', 'yqiL'}),
}


def numeric_st(value):
    return bool(re.fullmatch(r'[1-9][0-9]*', str(value or '').strip()))


def _resources():
    executable = shutil.which('mlst')
    if not executable:
        candidate = Path(sys.executable).parent / 'mlst'
        executable = str(candidate) if candidate.is_file() else None
    if not executable:
        # Keep classical MLST's Perl dependencies separate from tree/time tools.
        candidate = Path(sys.executable).parent.parent.parent / 'mlst' / 'bin' / 'mlst'
        executable = str(candidate) if candidate.is_file() else None
    if not executable:
        return None
    executable = Path(executable).resolve()
    root = Path(os.environ.get('MLST_DBDIR') or executable.parent.parent / 'db').resolve()
    environment = dict(os.environ, PATH=str(executable.parent) + os.pathsep + os.environ.get('PATH', ''))
    version = subprocess.run([str(executable), '--version'], capture_output=True,
                             text=True, env=environment, timeout=30, check=True).stdout.strip()
    return executable, root, environment, version


def fill_missing_mlst(samples, assembly_root, evidence_directory):
    """Reuse supplied STs; call missing input STs when a local assembly exists.

    Returns sample rows plus portable per-sample audit records. Missing resources
    and unsuccessful calls remain explicit and do not block local cgMLST work.
    """
    rows = [dict(row) for row in samples]
    audits, resources, resource_error, databases = [], None, None, {}
    attempted_resources = False
    evidence_directory = Path(evidence_directory)
    for row in rows:
        if row.get('role') != 'input':
            continue
        if numeric_st(row.get('mlst_st')):
            row['mlst_st'] = str(row['mlst_st']).strip()
            audits.append({'sample_id': row['sample_id'], 'status': 'reused',
                           'st': str(row['mlst_st']).strip(), 'scheme': row.get('mlst_scheme'),
                           'source': row.get('mlst_source') or 'supplied_metadata'})
            continue
        audit = {'sample_id': row['sample_id'], 'status': 'unavailable',
                 'input_st': row.get('mlst_st'), 'source': 'classical_mlst'}
        audits.append(audit)
        species = ' '.join(str(row.get('species') or '').replace('_', ' ').casefold().split())
        spec = PRIMARY_SCHEMES.get(species)
        reference = row.get('assembly_reference')
        if not spec:
            audit['reason'] = 'primary_provider_mlst_scheme_not_configured_for_species'
        elif not reference or '://' in reference or not (Path(assembly_root) / reference).is_file():
            audit['reason'] = 'local_assembly_required_for_missing_st'
        else:
            scheme, loci = spec
            audit['scheme'] = scheme
            if not attempted_resources:
                attempted_resources = True
                try:
                    resources = _resources()
                except (OSError, subprocess.SubprocessError) as error:
                    resource_error = type(error).__name__
            if resources is None:
                audit['reason'] = ('mlst_tool_failed_preflight' if resource_error
                                   else 'mlst_tool_unavailable_run_pixi_install_--frozen_-e_mlst')
            else:
                executable, root, environment, version = resources
                directory = root / 'pubmlst' / scheme
                files = [directory / (scheme + '.txt'), *(directory / (gene + '.tfa') for gene in sorted(loci))]
                blast = root / 'blast' / 'mlst.fa'
                if any(not path.is_file() for path in files) or not list(blast.parent.glob(blast.name + '.*')):
                    audit['reason'] = 'mlst_scheme_or_blast_database_unavailable'
                else:
                    if scheme not in databases:
                        hashes = {path.name: file_sha256(path) for path in files}
                        databases[scheme] = {'files': hashes, 'sha256': hashlib.sha256(
                            json.dumps(hashes, sort_keys=True).encode()).hexdigest()}
                    audit.update(tool_version=version, tool_sha256=file_sha256(executable),
                                 database=databases[scheme], search_field='mlst')
                    assembly = (Path(assembly_root) / reference).resolve()
                    audit['assembly_sha256'] = file_sha256(assembly)
                    evidence_directory.mkdir(parents=True, exist_ok=True)
                    prefix = hashlib.sha256(row['sample_id'].encode()).hexdigest()[:24]
                    raw = evidence_directory / (prefix + '.json')
                    command = [str(executable), '--quiet', '--threads', '1', '--scheme', scheme,
                               '--datadir', str(root / 'pubmlst'), '--blastdb', str(blast),
                               '--json', str(raw), str(assembly)]
                    audit['command'] = command
                    try:
                        with progress_task('mlst', f"Classical MLST for {row['sample_id']} ({scheme})"):
                            result = subprocess.run(command, env=environment, capture_output=True,
                                                    text=True, timeout=600, check=False)
                        (evidence_directory / (prefix + '.stdout.raw')).write_text(result.stdout)
                        (evidence_directory / (prefix + '.stderr.raw')).write_text(result.stderr)
                        if result.returncode:
                            audit.update(status='failed', reason='mlst_caller_nonzero_exit',
                                         exit_code=result.returncode)
                        else:
                            parsed = json.loads(raw.read_text())
                            if (not isinstance(parsed, list) or len(parsed) != 1
                                    or Path(parsed[0]['filename']).resolve() != assembly
                                    or parsed[0]['scheme'] != scheme):
                                raise ValueError('MLST output identity or scheme differs from input')
                            call = parsed[0]
                            alleles = call.get('alleles') or {}
                            audit['alleles'] = alleles
                            audit['raw_result'] = raw.name
                            audit['raw_result_sha256'] = file_sha256(raw)
                            if (numeric_st(call.get('sequence_type')) and set(alleles) == loci
                                    and all(numeric_st(value) for value in alleles.values())):
                                st = str(call['sequence_type'])
                                row.update(mlst_st=st, mlst_scheme=scheme, mlst_status='resolved',
                                           mlst_source='classical_mlst', mlst_alleles=alleles,
                                           mlst_database_sha256=databases[scheme]['sha256'],
                                           mlst_tool_version=version,
                                           mlst_assembly_sha256=audit['assembly_sha256'])
                                audit.update(status='resolved', st=st)
                            else:
                                audit.update(status='unresolved', reason='novel_partial_missing_or_ambiguous_mlst_profile')
                    except (OSError, ValueError, KeyError, TypeError, subprocess.SubprocessError) as error:
                        audit.update(status='failed', reason='mlst_call_or_output_validation_failed',
                                     error_type=type(error).__name__)
        if audit['status'] != 'resolved':
            row.update(mlst_status=audit['status'], mlst_source='classical_mlst',
                       mlst_reason=audit['reason'])
    if audits:
        evidence_directory.mkdir(parents=True, exist_ok=True)
        write_json(evidence_directory / 'audit.json', {
            'schema': 'chronoclade.classical-mlst', 'schema_version': 1,
            'policy': 'reuse supplied ST; exact complete primary-scheme calls for missing ST',
            'samples': audits,
        })
    return tuple(rows), audits
