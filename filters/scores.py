from rdkit.Chem import QED
from rdkit.Chem import RWMol
import math


def calc_qed(mol):
    return QED.qed(mol)


# SA score implementation (Ertl & Schuffenhauer, 2009)
def calc_sa_score(mol):
    from rdkit.Chem import rdMolDescriptors
    from rdkit import Chem
    import os, pickle, gzip

    # Try to load from RDKit contrib
    try:
        from rdkit.Chem.Scaffolds import MurckoScaffold
        from rdkit.Contrib.SA_Score import sascorer
        return sascorer.calculateScore(mol)
    except Exception:
        pass

    # Fallback: simple complexity estimate
    ring_info = mol.GetRingInfo()
    n_rings = ring_info.NumRings()
    n_stereo = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    n_heavy = mol.GetNumHeavyAtoms()
    score = 1 + 0.1 * n_rings + 0.2 * n_stereo + 0.01 * n_heavy
    return min(score, 10.0)
