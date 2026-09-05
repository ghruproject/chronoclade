from __future__ import annotations

import pytest

from validation.komori2024_st131.prepare_dataset import build_provenance


def paper_row(*, strain: str = "131", accession: str = "NZ_OX030701") -> dict[str, str]:
    return {
        "Strain name": strain,
        "GenBank accession no.": accession,
        "Year": "2018",
        "Country": "Spain",
        "Source": "Human",
        "Clade": "C1",
        "C1-subclade": "8",
    }


def ncbi_record(*, organism: str = "Escherichia coli") -> dict[str, str]:
    return {
        "accession_version": "NZ_OX030701.1",
        "organism": organism,
        "definition": "Escherichia coli isolate 131 chromosome, complete sequence",
        "sequence_length": "5070006",
        "strain": "131",
        "collection_date": "2018",
        "country": "Spain",
        "biosample": "SAMEA8065784",
        "assembly": "GCF_905330895.2",
    }


def test_provenance_accepts_valid_infraspecific_ncbi_name() -> None:
    rows = build_provenance(
        [paper_row()],
        {"NZ_OX030701": ncbi_record(organism="Escherichia coli O25b:H4-ST131")},
    )

    assert rows[0]["species_check"] == "pass"
    assert rows[0]["strain_name_check"] == "pass"
    assert rows[0]["year_check"] == "pass"
    assert rows[0]["country_check"] == "pass"


def test_provenance_stops_on_a_different_organism() -> None:
    with pytest.raises(RuntimeError, match="non-E. coli records: 131"):
        build_provenance(
            [paper_row()],
            {"NZ_OX030701": ncbi_record(organism="Klebsiella pneumoniae")},
        )
