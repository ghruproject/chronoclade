from pathlib import Path

import pytest

import chronoclade
from chronoclade import execution_provenance as provenance


def package(root: Path):
    root.mkdir()
    (root/'__init__.py').write_text('version = 1\n')
    (root/'method.py').write_text('result = 2\n')
    return root


def test_identity_binds_source_version_method_runtime_not_location(monkeypatch, tmp_path):
    first, second = package(tmp_path/'first'), package(tmp_path/'second')
    monkeypatch.setattr(provenance, '_runtime', lambda: {'numpy':'1'})
    identity = provenance.execution_provenance('cgmlst', '1', source_root=first)
    assert identity == provenance.execution_provenance('cgmlst', '1', source_root=second)
    assert identity != provenance.execution_provenance('cgmlst', '2', source_root=first)
    (second/'method.py').write_text('result = 3\n')
    assert identity != provenance.execution_provenance('cgmlst', '1', source_root=second)
    monkeypatch.setattr(chronoclade, '__version__', 'changed')
    assert identity != provenance.execution_provenance('cgmlst', '1', source_root=first)
    changed = provenance.execution_provenance('cgmlst', '1', source_root=first)
    monkeypatch.setattr(provenance, '_runtime', lambda: {'numpy':'2'})
    assert changed != provenance.execution_provenance('cgmlst', '1', source_root=first)


def test_legacy_load_is_separate_from_strict_reuse(tmp_path):
    current = provenance.execution_provenance('tree', '1', source_root=package(tmp_path/'source'))
    provenance.require_execution_match(current.copy(), current)
    for legacy in (None, {}, {'identity_sha256':current['identity_sha256']}):
        with pytest.raises(ValueError, match='new output'):
            provenance.require_execution_match(legacy, current)


def test_invalid_sources_and_methods(tmp_path):
    with pytest.raises(ValueError):
        provenance.execution_provenance('tree', '1', source_root=tmp_path)
    with pytest.raises(ValueError):
        provenance.execution_provenance('tree', '')
    (tmp_path/'module.py').symlink_to(tmp_path/'elsewhere')
    with pytest.raises(ValueError):
        provenance.execution_provenance('tree', '1', source_root=tmp_path)
