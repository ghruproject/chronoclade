"""Checked receipts for interrupted public catalogue captures.

This captures a paginated public query, not a transactional upstream snapshot.
Only the caller holding the catalogue lock may update these receipts.
"""

from pathlib import Path
from datetime import datetime, timezone

from chronoclade.artifacts import write_json
from chronoclade.prepare_recovery import (
    RecoveryError, canonical_digest, load_object, owned_directory,
)


class CatalogueCapture:
    def __init__(self, root: Path, binding: dict, *, resume: bool):
        self.root = Path(root)
        self.reused = {"pages": 0, "details": 0}
        if self.root.is_symlink():
            raise RecoveryError("Catalogue capture directory is a symlink")
        manifest_path = self.root / "capture.json"
        if self.root.exists():
            if not resume:
                raise RecoveryError(
                    f"Catalogue capture exists at {self.root}; resume the interrupted context retrieval"
                )
            if manifest_path.is_symlink():
                raise RecoveryError("Catalogue capture manifest is a symlink")
            manifest = load_object(manifest_path, "catalogue capture")
            if (manifest.get("schema") != "chronoclade.catalogue-capture.v1"
                    or manifest.get("binding") != binding
                    or not isinstance(manifest.get("started_at"), str)):
                raise RecoveryError("Catalogue capture binding mismatch or invalid manifest")
            for child in self.root.rglob("*"):
                if child.is_symlink():
                    raise RecoveryError("Catalogue capture contains a symlink")
        else:
            self.root.mkdir(parents=True)
            manifest = {"schema": "chronoclade.catalogue-capture.v1", "binding": binding,
                        "started_at": datetime.now(timezone.utc).isoformat()}
            write_json(manifest_path, manifest)
        self.started_at = manifest["started_at"]
        for kind in self.reused:
            owned_directory(self.root, kind, create=True)

    def _path(self, kind, request):
        return self.root / kind / (canonical_digest(request) + ".json")

    def read(self, kind, request):
        path = self._path(kind, request)
        if not path.exists():
            return None
        receipt = load_object(path, "catalogue response receipt")
        if (receipt.get("request") != request
                or receipt.get("sha256") != canonical_digest(receipt.get("response"))):
            raise RecoveryError(f"Catalogue {kind} receipt failed identity/checksum validation")
        self.reused[kind] += 1
        return receipt["response"]

    def write(self, kind, request, response):
        write_json(self._path(kind, request), {
            "request": request, "response": response, "sha256": canonical_digest(response),
        })
