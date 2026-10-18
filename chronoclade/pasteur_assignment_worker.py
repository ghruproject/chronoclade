"""Run inside installed plincer's environment; reuse its published classifier.

No extra models or scientific methods are implemented here. Profiles are loaded
once, then the same upstream functions and conservative result parser are used.
"""
from __future__ import annotations

import json
import lzma
from pathlib import Path
import sqlite3
import sys


def main():
    import toml
    from plincer.allelestore import lookup_st
    from plincer.classify import assign_bin, build_result, closest_profiles

    requests, scheme_file, profiles_file, alleles_file = map(Path, sys.argv[1:])
    with lzma.open(profiles_file, "rt") as stream:
        profiles = json.load(stream)
    with scheme_file.open() as stream:
        scheme = toml.load(stream)
    outputs = []
    cached = {}
    encoded_cache = {}
    allele_cache = {}
    with sqlite3.connect(f"file:{alleles_file}?mode=ro", uri=True) as connection:
        cursor = connection.cursor()
        for raw in json.loads(requests.read_text()):
            encoded = raw["code"]
            if encoded in encoded_cache:
                outputs.append(encoded_cache[encoded])
                continue
            profile = []
            for index, token in enumerate(raw["code"].split("_")):
                if token.isdigit() or token == "":
                    profile.append(token)
                elif len(token) == 40:
                    lookup_key = (index, token[:15])
                    if lookup_key not in allele_cache:
                        allele_cache[lookup_key] = str(lookup_st(cursor, token[:15], index) or sys.maxsize)
                    profile.append(allele_cache[lookup_key])
                else:
                    raise ValueError("Invalid cgMLST token passed to plincer")
            key = tuple(profile)
            if key in cached:
                outputs.append(cached[key])
                encoded_cache[encoded] = outputs[-1]
                continue
            matches = closest_profiles(profile, profiles)
            if not matches:
                outputs.append(build_result(len(profile), "unassigned"))
                cached[key] = outputs[-1]
                encoded_cache[encoded] = outputs[-1]
                continue
            identity = matches[0]["identity"]
            selected_bin = assign_bin(identity, scheme["levels"])
            code = ["*"] * 10
            for index in range(selected_bin + 1):
                code[index] = matches[0]["LINcode"][index]
            outputs.append(build_result(
                len(profile), matches[0]["st"] if selected_bin == 9 else "unassigned",
                matches=matches, identical=matches[0]["identical"], identity=identity,
                closest_cgst="/".join(sorted(m["st"] for m in matches)), lincode=code,
                clonal_group=matches[0]["Clonal Group"] if selected_bin > 2 else "",
                sublineage=matches[0]["Sublineage"] if selected_bin > 1 else "",
            ))
            cached[key] = outputs[-1]
            encoded_cache[encoded] = outputs[-1]
    print(json.dumps({"assignments": outputs}))


if __name__ == "__main__":
    main()
