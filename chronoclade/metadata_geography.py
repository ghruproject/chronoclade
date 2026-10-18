"""Aggregate explicitly reported sample geography for profile summaries."""

from collections import Counter
from collections.abc import Iterable, Mapping


def aggregate_metadata_geography(rows: Iterable[Mapping]) -> list[dict]:
    """Count country, region, city and legacy NUTS2 metadata without inference."""
    counts = Counter(
        (
            row.get("origin", "context"),
            row.get("country") or "Unknown",
            row.get("region") or "Unknown",
            row.get("city") or "Unknown",
            row.get("nuts2") or "Unknown",
        )
        for row in rows
    )
    return [
        dict(origin=key[0], country=key[1], region=key[2], city=key[3], nuts2=key[4], count=count)
        for key, count in sorted(counts.items())
    ]
