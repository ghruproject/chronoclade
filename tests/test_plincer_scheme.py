import json
from pathlib import Path
import pytest

from chronoclade import plincer_scheme as m

ORIGINAL = (Path(__file__).parent / "fixtures/plincer_scheme_upstream.toml").read_bytes()


@pytest.mark.parametrize("policy", ["level6-only", "published-kpsc-629"])
def test_exact_correction_preserves_original_and_receipt(tmp_path, policy):
    source = tmp_path / "downloaded.toml"
    source.write_bytes(ORIGINAL)
    working = tmp_path / "db/scheme.published.toml"
    receipt = m.prepare_scheme(source, working, policy=policy)
    assert source.read_bytes() == ORIGINAL
    assert working.with_name("scheme.upstream.toml").read_bytes() == ORIGINAL
    if policy == "published-kpsc-629":
        assert m.validate_scheme(working) == receipt
    else:
        with pytest.raises(m.SchemeValidationError, match="published KpSC"):
            m.validate_scheme(working)
    assert receipt["original_sha256"] == m.ORIGINAL_SHA256
    assert receipt["working_sha256"] != receipt["original_sha256"]
    assert len(receipt["transformations"]) == (1 if policy == "level6-only" else 3)
    assert m.prepare_scheme(source, working, policy=policy) == receipt


@pytest.mark.parametrize(
    "mutation",
    [
        lambda b: b.replace(b"798.8871", b"798.8872"),
        lambda b: b.replace(b"798.8871", b"98.8871").replace(b"min = 98.8871", b"min = 98.9"),
        lambda b: b.replace(b"798.8871", b"98.8871").replace(b"min = 98.8871", b"min = 98.8"),
        lambda b: b.replace(b"798.8871", b"inf"),
        lambda b: b.replace(b"798.8871", b"nan"),
        lambda b: b.replace(b"798.8871", b"true"),
        lambda b: b.replace(b"798.8871", b"98.8871").replace(b"max = 100", b"max = 99.99"),
        lambda b: b.replace(b"798.8871", b"98.8871").replace(b"scheme_id = 18", b"scheme_id = 19"),
        lambda b: b.replace(b"798.8871", b"98.8871").replace(b'name = "7"', b'name = "eight"'),
        lambda b: b"not toml",
    ],
)
def test_unknown_malformed_never_corrected(tmp_path, mutation):
    source = tmp_path / "original"
    source.write_bytes(mutation(ORIGINAL))
    working = tmp_path / "db/scheme.published.toml"
    with pytest.raises(m.SchemeValidationError):
        m.prepare_scheme(source, working)
    assert not working.exists()


@pytest.mark.parametrize(
    "target", ["original", "working", "receipt", "missing_original", "missing_receipt"]
)
def test_corrected_reference_tamper_fails_closed(tmp_path, target):
    source = tmp_path / "original"
    source.write_bytes(ORIGINAL)
    working = tmp_path / "db/scheme.published.toml"
    m.prepare_scheme(source, working)
    if target == "original":
        working.with_name("scheme.upstream.toml").write_bytes(ORIGINAL + b"\n")
    if target == "working":
        working.write_bytes(working.read_bytes() + b"\n")
    if target == "receipt":
        p = working.with_name("scheme.published.correction.json")
        v = json.loads(p.read_text())
        v["transformations"][0]["to"] = 98.9
        p.write_text(json.dumps(v))
    if target == "missing_original":
        working.with_name("scheme.upstream.toml").unlink()
    if target == "missing_receipt":
        working.with_name("scheme.published.correction.json").unlink()
    with pytest.raises(m.SchemeValidationError):
        m.validate_scheme(working)


def test_original_bad_reference_rejected_without_mutation(tmp_path):
    path = tmp_path / "scheme.toml"
    path.write_bytes(ORIGINAL)
    with pytest.raises(m.SchemeValidationError):
        m.validate_scheme(path)
    assert path.read_bytes() == ORIGINAL


def test_valid_future_scheme_allowed_without_automatic_edit(tmp_path):
    data = (
        ORIGINAL.replace(b"798.8871", b"98.8871").replace(b"99.3461", b"99.3641")
        + b"\n# future upstream valid source\n"
    )
    source = tmp_path / "source"
    source.write_bytes(data)
    dest = tmp_path / "db/scheme.published.toml"
    assert m.prepare_scheme(source, dest) is None
    assert dest.read_bytes() == data
    assert m.validate_scheme(dest) is None


def test_derivation_preserves_old_paths_and_reuses_all_large_resources(tmp_path):
    db = tmp_path / "db"
    db.mkdir()
    old = db / "scheme.toml"
    old.write_bytes(ORIGINAL)
    profiles = db / "profiles.json.xz"
    profiles.write_bytes(b"unchanged profiles")
    alleles = db / "alleles.sqlite"
    alleles.write_bytes(b"unchanged alleles")
    index = tmp_path / "index"
    index.mkdir()
    config = {
        "ready": True,
        "organisms": {
            "Klebsiella pneumoniae": {
                "cgmlst": {"index_dir": "index", "command": ["python"]},
                "assignment": {
                    "kind": "plincer",
                    "scheme_file": "db/scheme.toml",
                    "profiles_file": "db/profiles.json.xz",
                    "alleles_db": "db/alleles.sqlite",
                    "command": ["python"],
                },
            }
        },
    }
    source = tmp_path / "query.json"
    source.write_text(json.dumps(config))
    old_bytes = {p: p.read_bytes() for p in (source, old, profiles, alleles)}
    old_mtime = {p: p.stat().st_mtime_ns for p in old_bytes}
    target = tmp_path / "elsewhere/query.published.json"
    assert m.derive_config(source, target) == target.resolve()
    result = json.loads(target.read_text())
    assign = result["organisms"]["Klebsiella pneumoniae"]["assignment"]
    assert assign["scheme_file"] == str(db / "scheme.published.toml")
    assert assign["profiles_file"] == str(profiles)
    assert assign["alleles_db"] == str(alleles)
    assert result["organisms"]["Klebsiella pneumoniae"]["cgmlst"]["index_dir"] == str(index)
    assert m.validate_scheme(Path(assign["scheme_file"]))["policy"] == "published-kpsc-629"
    for p, expected in old_bytes.items():
        assert p.read_bytes() == expected
        assert p.stat().st_mtime_ns == old_mtime[p]
    with pytest.raises(m.SchemeValidationError, match="already exists"):
        m.derive_config(source, target)
    with pytest.raises(m.SchemeValidationError, match="new path"):
        m.derive_config(source, source)


def test_never_overwrites_historical_scheme_path(tmp_path):
    old = tmp_path / "scheme.toml"
    old.write_bytes(ORIGINAL)
    for target in (old, tmp_path / "another-scheme.toml"):
        with pytest.raises(m.SchemeValidationError, match="separate scheme.published.toml"):
            m.prepare_scheme(old, target)
    assert old.read_bytes() == ORIGINAL


def test_structurally_valid_noncanonical_future_scheme_not_ready(tmp_path):
    data = ORIGINAL.replace(b"798.8871", b"98.8871") + b"\n# not the exact known source\n"
    source = tmp_path / "unknown"
    source.write_bytes(data)
    with pytest.raises(m.SchemeValidationError, match="custom schemes are unverified"):
        m.prepare_scheme(source, tmp_path / "db/scheme.published.toml")


def test_canonical_boundaries_match_published_integer_loci():
    assert m.PUBLISHED_BOUNDS == [
        round(i / 629 * 100, 4) for i in [0, 19, 44, 439, 586, 619, 622, 625, 627, 628, 629]
    ]
