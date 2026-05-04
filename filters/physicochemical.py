from rdkit.Chem import Descriptors, rdMolDescriptors


def calc_properties(mol):
    return {
        "MW":              Descriptors.MolWt(mol),
        "logP":            Descriptors.MolLogP(mol),
        "TPSA":            Descriptors.TPSA(mol),
        "HBA":             rdMolDescriptors.CalcNumHBA(mol),
        "HBD":             rdMolDescriptors.CalcNumHBD(mol),
        "RotBonds":        rdMolDescriptors.CalcNumRotatableBonds(mol),
        "Fsp3":            rdMolDescriptors.CalcFractionCSP3(mol),
        "HeavyAtoms":      mol.GetNumHeavyAtoms(),
        "RingCount":       rdMolDescriptors.CalcNumRings(mol),
        "AromaticRings":   rdMolDescriptors.CalcNumAromaticRings(mol),
    }


def apply_cutoffs(props, args):
    reasons = []
    checks = [
        ("logP",     args.logp,     "logP"),
        ("TPSA",     args.tpsa,     "TPSA"),
        ("HBA",      args.hba,      "HBA"),
        ("HBD",      args.hbd,      "HBD"),
        ("RotBonds", args.rotbonds, "RotBonds"),
        ("Fsp3",     args.fsp3_min, "Fsp3_min", "min"),
    ]
    for item in checks:
        if len(item) == 3:
            key, cutoff, label = item
            if cutoff is not None and props[key] > cutoff:
                reasons.append(f"{label}={props[key]:.2f}>{cutoff}")
        else:
            key, cutoff, label, mode = item
            if cutoff is not None and props[key] < cutoff:
                reasons.append(f"{label}={props[key]:.2f}<{cutoff}")
    return reasons
