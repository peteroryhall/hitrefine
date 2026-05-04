from rdkit import DataStructs
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import Chem


def morgan_fp(mol, radius=2, nbits=2048):
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)


def tanimoto(mol, ref_fps):
    if not ref_fps:
        return None
    fp = morgan_fp(mol)
    sims = [DataStructs.TanimotoSimilarity(fp, ref) for ref in ref_fps]
    return max(sims)


def murcko_scaffold(mol):
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaffold)
    except Exception:
        return None


def load_reference_fps(sdf_path):
    suppl = Chem.SDMolSupplier(sdf_path, removeHs=True)
    fps = []
    for mol in suppl:
        if mol is not None:
            fps.append(morgan_fp(mol))
    return fps
