"""Small SVG forests that distinguish within-run intervals from selection spread."""

from html import escape


def date_forest(rows, sensitivity):
    """Render explicit saved node-date intervals; absent bounds are never inferred."""
    values = [
        value
        for row in rows
        for value in (row.get("estimate"), row.get("lower_bound"), row.get("upper_bound"))
        if value is not None
    ]
    if not values:
        return None
    low, high = min(values), max(values)
    padding = max((high - low) * 0.06, 0.05)
    low, high = low - padding, high + padding
    left, right, width = 190, 740, 940
    height = 110 + 34 * (len(rows) + 1)

    def x(value):
        return left + (value - low) / (high - low) * (right - left)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        'role="img" aria-labelledby="forest-title forest-description">',
        '<title id="forest-title">Shared-anchor MRCA dates</title>',
        '<desc id="forest-description">Teal bars show saved within-run 90% max-posterior '
        "regions. Orange shows between-selection minimum, median and maximum of point "
        "estimates, which is selection sensitivity and not a confidence interval. "
        "Unavailable runs remain listed.</desc>",
        '<rect width="100%" height="100%" fill="white"/>',
        '<g font-family="Arial,sans-serif" font-size="13" fill="#252525">',
        '<text x="190" y="22">Within-run TreeTime 90% max-posterior regions</text>',
    ]
    for index, row in enumerate(rows):
        y = 52 + index * 34
        parts.append(f'<text x="12" y="{y + 4}">{escape(row["label"])}</text>')
        point, lower, upper = (row.get(key) for key in ("estimate", "lower_bound", "upper_bound"))
        if point is None:
            parts.append(f'<text x="{left}" y="{y + 4}" fill="#666">Unavailable</text>')
            continue
        parts.append(
            "<g><title>" + escape(f"{row['label']}: {point}; bounds {lower}, {upper}") + "</title>"
        )
        if lower is not None and upper is not None:
            parts.append(
                f'<path d="M{x(lower):.3f},{y}H{x(upper):.3f}" stroke="#167d8d" stroke-width="3"/>'
            )
        parts.append(f'<circle cx="{x(point):.3f}" cy="{y}" r="4" fill="#167d8d"/></g>')
        text = f"{point:.4f}" + (
            f" [{lower:.4f}, {upper:.4f}]"
            if lower is not None and upper is not None
            else " (bounds unavailable)"
        )
        parts.append(f'<text x="760" y="{y + 4}" font-size="11">{escape(text)}</text>')
    y = 52 + len(rows) * 34
    parts.append(f'<text x="12" y="{y + 4}">Selection sensitivity</text>')
    if sensitivity.get("median") is not None:
        minimum, median, maximum = (sensitivity[key] for key in ("minimum", "median", "maximum"))
        parts.append(
            f'<path d="M{x(minimum):.3f},{y}H{x(maximum):.3f}" stroke="#ad631b" stroke-width="5"/>'
        )
        parts.append(f'<circle cx="{x(median):.3f}" cy="{y}" r="5" fill="#ad631b"/>')
        parts.append(f'<text x="760" y="{y + 4}" font-size="11">Median {median:.4f}</text>')
    axis = y + 30
    parts.append(f'<path d="M{left},{axis}H{right}" stroke="#555"/>')
    for tick in range(5):
        value = low + (high - low) * tick / 4
        parts.append(
            f'<text x="{x(value):.3f}" y="{axis + 19}" text-anchor="middle">{value:.3f}</text>'
        )
    parts.append(
        f'<text x="190" y="{axis + 43}">Calendar year · orange range is not a confidence interval</text></g></svg>'
    )
    return "".join(parts)
