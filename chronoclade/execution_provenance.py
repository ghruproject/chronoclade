"""Content-addressed software identity for numerical reuse, separate from imports.

Historical results remain readable. Callers explicitly bind this identity to a
computation's cache contract; merely rendering an old result must not replace it.
"""
from importlib import metadata
from pathlib import Path
import platform

from chronoclade.artifacts import file_sha256
from chronoclade.pathogenwatch import content_hash

_LIBRARIES = ('numpy', 'scipy', 'biopython', 'pandas', 'dendropy', 'phylo-treetime')


def _runtime():
    libraries = {}
    for name in _LIBRARIES:
        try:
            libraries[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            libraries[name] = None
    return {'python': platform.python_version(), 'implementation': platform.python_implementation(),
            'system': platform.system(), 'machine': platform.machine(), 'libraries': libraries}


def execution_provenance(component, method_revision, *, source_root=None, software_version=None):
    """Hash installed Python source, method revision and relevant runtime versions.

    Source location and modification times are deliberately excluded. All package
    Python files count, including numerical routines under report_components.
    """
    import chronoclade

    if not all(isinstance(value, str) and value.strip() for value in (component, method_revision)):
        raise ValueError('Execution identity requires a component and explicit method revision')
    version = software_version if software_version is not None else chronoclade.__version__
    if not isinstance(version, str) or not version.strip():
        raise ValueError('Execution identity requires a software version')
    root = Path(source_root) if source_root is not None else Path(chronoclade.__file__).parent
    files = sorted(root.rglob('*.py'))
    if not files or any(path.is_symlink() for path in files):
        raise ValueError('Execution identity requires non-symlink package Python sources')
    source = [(path.relative_to(root).as_posix(), file_sha256(path)) for path in files]
    payload = {
        'schema': 'chronoclade.execution-provenance', 'schema_version': 1,
        'component': component, 'method_revision': method_revision,
        'software': {'name': 'ChronoClade', 'version': version,
                     'source_sha256': content_hash(source), 'source_file_count': len(source)},
        'runtime': _runtime(),
    }
    return {**payload, 'identity_sha256': content_hash(payload)}


def require_execution_match(saved, current):
    """Refuse automatic numerical reuse when provenance is absent or different."""
    if not isinstance(saved, dict) or saved != current:
        raise ValueError('Saved numerical execution identity differs or is absent; use a new output directory')
