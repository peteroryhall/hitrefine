#!/usr/bin/env python3
"""
HitRefine - post-docking hit refinement tool.
Filters docking hits by physicochemical properties, structural alerts,
drug-likeness rules, and similarity to reference compounds.
"""

import argparse
import csv
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors

from filters.physicochemical import calc_properties, apply_cutoffs
from filters.rules import lipinski, veber, egan
from filters.scores import calc_qed, calc_sa_score
from filters.alerts import build_pains_catalog, build_brenk_catalog, check_pains, check_brenk
from filters.similarity import tanimoto, murcko_scaffold, load_reference_fps


def parse_args():
    p = argparse.ArgumentParser(
        description="HitRefine: filter docking hits by physicochemical and structural properties."
    )
    p.add_argument("-i", "--input",    required=True,  help="Input SDF file")
    p.add_argument("-o", "--output",   required=True,  help="Output filtered SDF file")
    p.add_argument("--csv",            default=None,   help="Output CSV with all properties (optional)")
    p.add_argument("--plots",          action="store_true", help="Generate property distribution plots")

    # Physicochemical cutoffs
    g = p.add_argument_group("Physicochemical filters (upper limits unless noted)")
    g.add_argument("--logp",     type=float, default=None, help="Max logP (e.g. 5)")
    g.add_argument("--tpsa",     type=float, default=None, help="Max TPSA in Å² (e.g. 90)")
    g.add_argument("--hba",      type=int,   default=None, help="Max H-bond acceptors (e.g. 7)")
    g.add_argument("--hbd",      type=int,   default=None, help="Max H-bond donors (e.g. 5)")
    g.add_argument("--rotbonds", type=int,   default=None, help="Max rotatable bonds (e.g. 10)")
    g.add_argument("--fsp3-min", type=float, default=None, help="Min Fsp3 fraction (e.g. 0.2)")

    # Rule-based filters
    g2 = p.add_argument_group("Rule-based filters")
    g2.add_argument("--lipinski", action="store_true", help="Apply Lipinski Ro5")
    g2.add_argument("--veber",    action="store_true", help="Apply Veber rules")
    g2.add_argument("--egan",     action="store_true", help="Apply Egan rules")

    # Score-based filters
    g3 = p.add_argument_group("Score-based filters")
    g3.add_argument("--qed",     type=float, default=None, help="Min QED score (e.g. 0.3)")
    g3.add_argument("--sa-max",  type=float, default=None, help="Max SA score (e.g. 6.0, lower=easier to synthesise)")

    # Structural alerts
    g4 = p.add_argument_group("Structural alert filters")
    g4.add_argument("--pains", action="store_true", help="Remove PAINS compounds")
    g4.add_argument("--brenk", action="store_true", help="Remove Brenk alert compounds")

    # Similarity
    g5 = p.add_argument_group("Similarity filters")
    g5.add_argument("--ref",         default=None,   help="Reference SDF for Tanimoto similarity")
    g5.add_argument("--sim-min",     type=float,     default=None, help="Min Tanimoto similarity to reference (e.g. 0.3)")
    g5.add_argument("--sim-max",     type=float,     default=None, help="Max Tanimoto similarity to reference (e.g. 0.9)")

    return p.parse_args()


def make_plots(records, output_dir):
    import matplotlib.pyplot as plt
    props = ["logP", "TPSA", "HBA", "HBD", "RotBonds", "Fsp3", "QED", "SA_score"]
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.flatten()
    for i, prop in enumerate(props):
        vals_all = [r[prop] for r in records if r.get(prop) is not None]
        vals_pass = [r[prop] for r in records if r.get(prop) is not None and r["passes"]]
        axes[i].hist(vals_all, bins=30, alpha=0.5, color="gray", label="all")
        axes[i].hist(vals_pass, bins=30, alpha=0.7, color="#2196F3", label="passing")
        axes[i].set_title(prop)
        axes[i].legend(fontsize=7)
    plt.suptitle("HitRefine property distributions", fontweight="bold")
    plt.tight_layout()
    out = Path(output_dir) / "hitrefine_plots.png"
    plt.savefig(str(out), dpi=150, bbox_inches="tight")
    print(f"Plots saved: {out}")


def main():
    args = parse_args()

    suppl = Chem.SDMolSupplier(args.input, removeHs=True)
    if suppl is None:
        sys.exit(f"Could not read input: {args.input}")

    # Build catalogs once
    pains_cat = build_pains_catalog() if args.pains else None
    brenk_cat = build_brenk_catalog() if args.brenk else None
    ref_fps = load_reference_fps(args.ref) if args.ref else []

    writer = Chem.SDWriter(args.output)
    records = []
    n_total = n_pass = 0

    for mol in suppl:
        if mol is None:
            continue
        n_total += 1
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"mol_{n_total}"
        props = calc_properties(mol)
        qed_val = calc_qed(mol)
        sa_val = calc_sa_score(mol)
        scaffold = murcko_scaffold(mol)
        sim = tanimoto(mol, ref_fps) if ref_fps else None

        fail_reasons = []

        # Physicochemical cutoffs
        fail_reasons += apply_cutoffs(props, args)

        # Rule filters
        if args.lipinski:
            violations, ok = lipinski(props)
            if not ok:
                fail_reasons.append(f"Lipinski({violations} violations)")
        if args.veber and not veber(props):
            fail_reasons.append("Veber")
        if args.egan and not egan(props):
            fail_reasons.append("Egan")

        # Score filters
        if args.qed is not None and qed_val < args.qed:
            fail_reasons.append(f"QED={qed_val:.2f}<{args.qed}")
        if args.sa_max is not None and sa_val > args.sa_max:
            fail_reasons.append(f"SA={sa_val:.2f}>{args.sa_max}")

        # Structural alerts
        if pains_cat:
            ok, desc = check_pains(mol, pains_cat)
            if not ok:
                fail_reasons.append(f"PAINS:{desc}")
        if brenk_cat:
            ok, desc = check_brenk(mol, brenk_cat)
            if not ok:
                fail_reasons.append(f"Brenk:{desc}")

        # Similarity
        if sim is not None:
            if args.sim_min is not None and sim < args.sim_min:
                fail_reasons.append(f"Sim={sim:.2f}<{args.sim_min}")
            if args.sim_max is not None and sim > args.sim_max:
                fail_reasons.append(f"Sim={sim:.2f}>{args.sim_max}")

        passes = len(fail_reasons) == 0
        if passes:
            n_pass += 1
            writer.write(mol)

        records.append({
            "name": name,
            "passes": passes,
            "fail_reasons": "; ".join(fail_reasons),
            "scaffold": scaffold,
            "tanimoto": sim,
            "QED": qed_val,
            "SA_score": sa_val,
            **props,
        })

    writer.close()

    print(f"\nResults: {n_pass} / {n_total} compounds passed")
    print(f"Output:  {args.output}")

    # CSV output
    csv_path = args.csv or str(Path(args.output).with_suffix(".csv"))
    fieldnames = ["name", "passes", "fail_reasons", "MW", "logP", "TPSA", "HBA", "HBD",
                  "RotBonds", "Fsp3", "QED", "SA_score", "tanimoto", "scaffold"]
    with open(csv_path, "w", newline="") as f:
        writer_csv = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer_csv.writeheader()
        writer_csv.writerows(records)
    print(f"CSV:     {csv_path}")

    if args.plots:
        make_plots(records, Path(args.output).parent)


if __name__ == "__main__":
    main()
