"""Public EnteroBase cgMLST references with frozen, resumable file evidence.

Availability of allele definitions does not establish native or HierCC support.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import gzip
import hashlib
from html import unescape
import json
from pathlib import Path
import re
import tempfile
from urllib.request import Request, urlopen

from chronoclade.artifacts import file_sha256, write_json
from chronoclade.progress import progress

BASE_URL = 'https://enterobase.warwick.ac.uk/schemes/'
SCHEMES = {
    'Escherichia.cgMLSTv1': 'Escherichia/Shigella',
    'Photorhabdus.cgMLSTv1': 'Photorhabdus',
    'Salmonella.cgMLSTv2': 'Salmonella',
    'Streptococcus.cgMLSTv1': 'Streptococcus',
    'VIBwgMLST.cgMLSTv1': 'Vibrio',
    'Yersinia.cgMLSTv1': 'Yersinia',
    'clostridium.cgMLSTv1': 'Clostridioides (historical directory name)',
    'klebsiella.pasteurcgmlst': 'Klebsiella pneumoniae species complex',
}
_ALLELE_FILE = re.compile(r'([A-Za-z0-9][A-Za-z0-9_.-]*)\.fasta\.gz')


class EnterobaseSchemeError(ValueError):
    """Public scheme retrieval or integrity validation failed."""


def _open_url(url):
    return urlopen(Request(url, headers={'User-Agent': 'ChronoClade-public-cgMLST-reference'}), timeout=60)


def parse_index(body, scheme_id):
    """Accept only local, explicitly listed allele files and their byte lengths."""
    if scheme_id not in SCHEMES:
        raise EnterobaseSchemeError('Unregistered EnteroBase cgMLST scheme')
    files = []
    for line in body.decode('utf-8').splitlines():
        href = re.search(r'<a\s+href="([^"]+)"', line)
        if not href:
            continue
        name = unescape(href.group(1))
        match = _ALLELE_FILE.fullmatch(name)
        if not match:
            continue
        text = unescape(re.sub(r'<[^>]+>', '', line))
        size = re.search(r'\s(\d+)\s*$', text)
        if not size or int(size.group(1)) <= 0:
            raise EnterobaseSchemeError(f'Missing allele-file size: {name}')
        files.append({'locus': match.group(1), 'filename': name, 'compressed_bytes': int(size.group(1))})
    files.sort(key=lambda item: item['locus'])
    loci = [item['locus'] for item in files]
    if not loci or len(loci) != len(set(loci)) or len(loci) > 10000:
        raise EnterobaseSchemeError('Empty, duplicate or unexpectedly large cgMLST locus list')
    return {
        'schema': 'chronoclade.enterobase-scheme-index', 'schema_version': 1,
        'scheme_id': scheme_id, 'organism_group': SCHEMES[scheme_id],
        'url': BASE_URL+scheme_id+'/', 'index_sha256': hashlib.sha256(body).hexdigest(),
        'loci': loci, 'allele_files': files,
        'compressed_bytes': sum(item['compressed_bytes'] for item in files),
        'fetched_utc': datetime.now(timezone.utc).isoformat(),
        'support_scope': 'Public allele definitions; native calling and HierCC require separate validation',
    }


def _validate_snapshot(snapshot, scheme_id):
    """Validate frozen local metadata before using it as download paths or URLs."""
    if not isinstance(snapshot, dict):
        raise EnterobaseSchemeError('Invalid reference scheme snapshot')
    if (snapshot.get('schema') != 'chronoclade.enterobase-scheme-index'
            or snapshot.get('schema_version') != 1
            or snapshot.get('scheme_id') != scheme_id
            or snapshot.get('url') != BASE_URL+scheme_id+'/'):
        raise EnterobaseSchemeError('Invalid or different reference scheme snapshot')
    entries = snapshot.get('allele_files', [])
    if not isinstance(entries, list) or not 1 <= len(entries) <= 10000:
        raise EnterobaseSchemeError('Invalid frozen locus list')
    loci = []
    total = 0
    for entry in entries:
        if not isinstance(entry, dict):
            raise EnterobaseSchemeError('Invalid frozen allele file entry')
        match = _ALLELE_FILE.fullmatch(str(entry.get('filename', '')))
        size = entry.get('compressed_bytes')
        if (not match or match.group(1) != entry.get('locus')
                or type(size) is not int or size <= 0):
            raise EnterobaseSchemeError('Invalid frozen allele file entry')
        loci.append(entry['locus'])
        total += size
    if (loci != sorted(set(loci)) or loci != snapshot.get('loci')
            or total != snapshot.get('compressed_bytes')):
        raise EnterobaseSchemeError('Inconsistent frozen locus metadata')


def inspect_scheme(scheme_id):
    if scheme_id not in SCHEMES:
        raise EnterobaseSchemeError('Unregistered EnteroBase cgMLST scheme')
    with _open_url(BASE_URL+scheme_id+'/') as response:
        body = response.read(3_000_001)
    if len(body) > 3_000_000:
        raise EnterobaseSchemeError('Scheme index exceeds the metadata size limit')
    return parse_index(body, scheme_id)


def validate_alleles(path, *, max_uncompressed_bytes=256*1024**2):
    """Check gzip FASTA structure, unique allele IDs and locus-oriented DNA."""
    ids = set()
    length = 0
    total = 0
    try:
        with gzip.open(path, 'rt', encoding='ascii') as handle:
            for line in handle:
                total += len(line)
                if total > max_uncompressed_bytes:
                    raise EnterobaseSchemeError('Allele file exceeds decompressed size limit')
                text = line.strip()
                if not text:
                    continue
                if text.startswith('>'):
                    if ids and length == 0:
                        raise EnterobaseSchemeError('Empty allele sequence')
                    ident = text[1:].split()[0] if text[1:].strip() else ''
                    if not ident or ident in ids:
                        raise EnterobaseSchemeError('Missing or duplicate allele identifier')
                    ids.add(ident)
                    length = 0
                else:
                    if not ids or set(text.upper()) - set('ACGTRYSWKMBDHVN'):
                        raise EnterobaseSchemeError('Invalid allele FASTA or DNA symbols')
                    length += len(text)
        if not ids or not length:
            raise EnterobaseSchemeError('Empty allele FASTA')
    except (OSError, EOFError, UnicodeError) as error:
        raise EnterobaseSchemeError(f'Cannot read allele FASTA: {Path(path).name}') from error
    return {'allele_count': len(ids), 'uncompressed_bytes': total}


@contextmanager
def _locked(directory):
    with (directory/'.reference.lock').open('a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        yield


def _retrieve(entry, snapshot, destination):
    target = destination/'alleles'/entry['filename']
    evidence = target.with_name(target.name+'.json')
    if target.is_symlink() or evidence.is_symlink():
        raise EnterobaseSchemeError('Reference allele paths must not be symlinks')
    if target.exists() and evidence.exists():
        previous = json.loads(evidence.read_text())
        if (isinstance(previous, dict) and all(previous.get(key) == value for key, value in entry.items())
                and previous.get('path') == 'alleles/'+entry['filename']
                and previous.get('url') == snapshot['url']+entry['filename']
                and target.stat().st_size == entry['compressed_bytes']
                and previous.get('sha256') == file_sha256(target)):
            return previous
        raise EnterobaseSchemeError(f'Cached allele file changed: {target.name}')
    temporary = None
    url = snapshot['url']+entry['filename']
    try:
        with _open_url(url) as response, tempfile.NamedTemporaryFile(dir=target.parent, delete=False) as handle:
            temporary = Path(handle.name)
            size = 0
            while chunk := response.read(1024*1024):
                size += len(chunk)
                if size > entry['compressed_bytes']:
                    raise EnterobaseSchemeError(f'Allele file changed after index freeze: {target.name}')
                handle.write(chunk)
            modified = response.headers.get('Last-Modified')
        if size != entry['compressed_bytes']:
            raise EnterobaseSchemeError(f'Incomplete allele file: {target.name}')
        statistics = validate_alleles(temporary)
        record = {**entry, **statistics, 'url': url, 'last_modified': modified,
                  'path': 'alleles/'+entry['filename'], 'sha256': file_sha256(temporary)}
        temporary.replace(target)
        write_json(evidence, record)
        return record
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def install_scheme(scheme_id, destination, *, workers=4, max_compressed_bytes=2*1024**3):
    """Install one selected scheme; resume only its pinned index and checked files."""
    if scheme_id not in SCHEMES or not isinstance(workers, int) or not 1 <= workers <= 8:
        raise EnterobaseSchemeError('Select a registered scheme and 1–8 download workers')
    destination = Path(destination).expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    for name in ('index.json', 'reference.json', 'README.md', 'alleles', '.reference.lock'):
        if (destination/name).is_symlink():
            raise EnterobaseSchemeError('Reference paths must not be symlinks')
    with _locked(destination):
        frozen = destination/'index.json'
        if frozen.exists():
            snapshot = json.loads(frozen.read_text())
            if not isinstance(snapshot, dict) or snapshot.get('scheme_id') != scheme_id:
                raise EnterobaseSchemeError('Destination belongs to a different reference scheme')
        else:
            if any(path.name != '.reference.lock' for path in destination.iterdir()):
                raise EnterobaseSchemeError('Reference destination is not empty or a resumable scheme')
            snapshot = inspect_scheme(scheme_id)
            write_json(frozen, snapshot)
        _validate_snapshot(snapshot, scheme_id)
        if snapshot['compressed_bytes'] > max_compressed_bytes:
            raise EnterobaseSchemeError('Selected scheme exceeds the compressed download budget')
        (destination/'alleles').mkdir(exist_ok=True)
        progress('setup', f"{scheme_id}: {len(snapshot['loci'])} loci; {snapshot['compressed_bytes']:,} compressed bytes")
        records = []
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(_retrieve, entry, snapshot, destination) for entry in snapshot['allele_files']]
            for future in as_completed(futures):
                records.append(future.result())
                if len(records) % 50 == 0:
                    progress('setup', f"{scheme_id}: {len(records)}/{len(futures)} locus files verified")
        records.sort(key=lambda item: item['locus'])
        identity = {'scheme_id': scheme_id, 'loci': snapshot['loci'],
                    'files': [(item['locus'], item['sha256']) for item in records]}
        manifest = {**snapshot, 'schema': 'chronoclade.enterobase-allele-reference',
                    'files': records, 'status': 'alleles_ready',
                    'reference_sha256': hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest(),
                    'hiercc_status': 'not_installed'}
        write_json(destination/'reference.json', manifest)
        (destination/'README.md').write_text(
            f"# {scheme_id} allele reference\n\n"
            f"Verified catalogue: {len(records):,} locus files.\n"
            f"Reference SHA256: `{manifest['reference_sha256']}`.\n\n"
            "- `index.json`: frozen public index, ordered loci and listed byte lengths.\n"
            "- `reference.json`: content identity and per-file verification records.\n"
            "- `alleles/`: original compressed FASTA files and checksum sidecars.\n\n"
            "This is an allele catalogue, not a caller index or HierCC reference.\n"
            "The source is a live directory, not an upstream atomic release; content\n"
            "hashes identify the exact files acquired. Repeating setup in this folder\n"
            "verifies/resumes this snapshot. Use a new directory for newer data.\n",
            encoding='utf-8',
        )
        return manifest
