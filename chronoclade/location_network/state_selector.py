"""Small shared selector for reports with more than one reconstructed state."""

from html import escape
from uuid import uuid4


def state_network_panel(network, *, basis="", widget_id="state-network"):
    """Present an observed field consistently in every analysis module."""
    from .widget import interactive_network_html

    if not network:
        return None
    field = network.get("state_label", "country")
    plural = network.get("state_plural", f"{field} values")
    if network.get("status") == "insufficient_metadata":
        total = network.get("sample_count", 0)
        missing = network.get("unknown_state_count", total)
        known = network.get("known_state_count", 0)
        return (
            f'<p class="empty">{escape(field.capitalize())} metadata is recorded for '
            f'{max(0, total - missing)} of {total} genomes. This analysis has {known} '
            f'distinct known {escape(field)} state(s). At least two distinct known '
            f'{escape(plural)} are required to reconstruct {escape(field)} changes.</p>'
        )
    heading = f"{basis} {field} history".strip()
    note = (
        '<p>Place labels include known parent locations to distinguish places with the same name.</p>'
        if network.get("state_field") in {"region", "city"} else ""
    )
    return f'<h3>{escape(heading)}</h3>' + note + interactive_network_html(network, widget_id)


def state_network_selector(country_panel: str | None, host_panel: str | None = None, *, panels=None) -> str:
    """Switch observed metadata fields; omit optional fields without evidence."""
    fields = {"country": country_panel} if country_panel else {}
    if host_panel:
        fields["host"] = host_panel
    fields.update({key: value for key, value in (panels or {}).items() if value})
    if not fields:
        return '<p class="empty">No observed host or location metadata is available for a state-change network.</p>'
    if len(fields) == 1:
        return next(iter(fields.values()))
    selected_field = next(iter(fields))
    identifier = "cc-state-network-field-" + uuid4().hex[:12]
    labels = {"country": "Location (country)", "host": "Host", "region": "Region", "city": "City"}
    options = "".join(
        f'<option value="{escape(field, quote=True)}">{escape(labels.get(field, field.replace("_", " ").title()))}</option>'
        for field in fields
    )
    content = "".join(
        f'<div data-state-panel="{escape(field, quote=True)}"'
        + ("" if field == selected_field else " hidden") + ">" + panel + "</div>"
        for field, panel in fields.items()
    )
    return (
        '<div class="cc-state-network-selector">'
        f'<label for="{identifier}">Reconstruct state </label>'
        f'<select id="{identifier}" data-state-network-field>'
        + options + '</select>' + content +
        '<script data-state-network-switch>(function(){'
        'document.querySelectorAll(".cc-state-network-selector").forEach(root=>{'
        'if(root.dataset.stateSwitchReady)return;root.dataset.stateSwitchReady="true";'
        'const select=root.querySelector("[data-state-network-field]");'
        'select.addEventListener("change",()=>{'
        'root.querySelectorAll("[data-state-panel]").forEach(panel=>{'
        'panel.hidden=panel.dataset.statePanel!==select.value;'
        '});'
        'if(window.ChronoCladeNetworks)window.ChronoCladeNetworks.mount(root);'
        '});'
        '});})();</script></div>'
    )
