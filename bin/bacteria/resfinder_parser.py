#!/usr/bin/env python3
"""
Parse ResFinder / PointFinder tabular outputs and emit a de-duplicated JSON
summary compatible with the rest of the pzh pipeline.
"""

import sys
import json
from collections import defaultdict
import click


@click.command()
@click.option(
    "-i", "--input_file_resfinder",
    type=click.Path(exists=True, readable=True),
    default="ResFinder_results_tab.txt",
    required=True,
    help="[INPUT] path to the ResFinder tab file",
)
@click.option(
    "-j", "--input_file_pointfinder",
    type=click.Path(exists=True, readable=True),
    default="PointFinder_results.txt",
    required=True,
    help="[INPUT] path to the PointFinder tab file",
)
@click.option(
    "-s", "--status",
    type=click.Choice(["tak", "nie", "blad"], case_sensitive=False),
    required=True,
    help="[INPUT] PRE-DEFINED status transferred to the output JSON",
)
@click.option(
    "-r", "--error",
    type=str,
    default="",
    help="[INPUT] Optional error message (used when status ≠ 'tak')",
)
@click.option(
    "-o", "--output",
    type=click.Path(writable=True),
    required=True,
    help="[OUTPUT] filename for the JSON report",
)
def main_program(status, input_file_resfinder, input_file_pointfinder, output, error=""):
    """
    Main entry point – either emit a stub JSON (status ≠ tak) or parse both
    files and write the de-duplicated report.
    """
    if status.lower() != "tak":
        json_output = {
            "program_name": "ResFinder/PointFinder",
            "status": status,
            "error_message": error,
        }
    else:
        # ------------------------------------------------------------------ #
        # 1. gather raw data in 'store', using sets to avoid duplicates       #
        # ------------------------------------------------------------------ #
        store: dict[str, set[tuple]] = defaultdict(set)

        # ----------- ResFinder genes ------------------------------------- #
        with open(input_file_resfinder, encoding="utf-8") as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if cols[0] == "Resistance gene":
                    continue  # skip header

                (
                    resistance_gene,
                    seq_identity,
                    _alignment_len,            # unused
                    coverage,
                    _pos_ref,                  # unused
                    contig_name,
                    _pos_contig,               # unused
                    antibiotic,
                    reference_name,
                ) = cols

                # normalise & deduplicate antibiotic names on the line
                # in somce cases in antibiotic the same antibiotic can appear multiple times (e.g with extra " ")
                # we consider such entries as technical issues with Resfinder so we ignore these "duplicates"
                antibiotic_set = {ab.strip() for ab in antibiotic.split(",") if ab.strip()}

                entry = (
                    resistance_gene,
                    contig_name,
                    int(float(seq_identity)),
                    int(float(coverage)),
                    reference_name,
                    "gen",
                )

                for ab in antibiotic_set:
                    store[ab].add(entry)

        # ----------- PointFinder point mutations ------------------------- #
        with open(input_file_pointfinder, encoding="utf-8") as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if cols[0] == "Mutation":
                    continue  # skip header

                mutation, *_unused, antibiotic = cols[0:4]
                gene_name, mutation_name = mutation.split(" ", maxsplit=1)

                antibiotic_set = {ab.strip() for ab in antibiotic.split(",") if ab.strip()}

                entry = (
                    gene_name,
                    mutation_name,
                    "mutacja_punktowa",
                )

                for ab in antibiotic_set:
                    store[ab].add(entry)

        # ------------------------------------------------------------------ #
        # 2. build JSON, converting each per-antibiotic *set* back to a list  #
        # ------------------------------------------------------------------ #
        json_output = {
            "program_name": "ResFinder/PointFinder",
            "status": status,
            "program_data": [],
        }

        for antibiotic, factors in sorted(store.items()):
            antibiotic_dict = {
                "antibiotic_name": antibiotic,
                "antibiotic_status": "oporny",
                "antibiotic_resistance_data": [],
            }

            for factor in sorted(factors):
                if factor[-1] == "gen":
                    (
                        gene,
                        contig,
                        identity,
                        overlap,
                        ref_name,
                        _marker,
                    ) = factor
                    antibiotic_dict["antibiotic_resistance_data"].append(
                        {
                            "factor_name": gene,
                            "factor_contig_name": contig,
                            "factor_sequence_similarity_to_reference_value": identity,
                            "factor_degree_of_overlap_with_reference_value": overlap,
                            "factor_reference_name": ref_name,
                            "factor_type_name": "gen",
                        }
                    )
                else:  # point mutation
                    gene, mutation_name, _marker = factor
                    antibiotic_dict["antibiotic_resistance_data"].append(
                        {
                            "factor_name": gene,
                            "factor_mutation": mutation_name,
                            "factor_type_name": "mutacja_punktowa",
                        }
                    )

            json_output["program_data"].append(antibiotic_dict)

    # ---------------------------------------------------------------------- #
    # 3. write out the finished JSON                                         #
    # ---------------------------------------------------------------------- #
    with open(output, "w", encoding="utf-8") as f_out:
        json.dump(json_output, f_out, indent=4, ensure_ascii=False)


if __name__ == "__main__":
    # allow the usual '--help' behaviour when no args are passed
    if len(sys.argv) == 1:
        main_program(["--help"])
    else:
        main_program()

