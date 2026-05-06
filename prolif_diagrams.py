#!/usr/bin/env python3
"""
Generate 2D protein-ligand interaction diagrams for docking hits using ProLIF.
Outputs one PNG per ligand showing H-bonds and VdW contacts to receptor residues.
"""

import argparse
import warnings
warnings.filterwarnings('ignore')

import prolif as plf
import MDAnalysis as mda
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from PIL import Image
from pathlib import Path
import io


hbond_types = {'HBDonor', 'HBAcceptor'}
color_map = {
    'HBDonor':    '#2196F3',
    'HBAcceptor': '#2196F3',
    'VdWContact': '#9E9E9E',
    'Hydrophobic':'#9E9E9E',
    'PiStacking': '#4CAF50',
    'CationPi':   '#FF9800',
    'Anionic':    '#F44336',
    'Cationic':   '#9C27B0',
}


def make_diagram(mol, prot_mol, out_path):
    orig_to_2d = {}
    new_idx = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 1:
            orig_to_2d[atom.GetIdx()] = new_idx
            new_idx += 1

    mol_2d = Chem.RemoveHs(Chem.RWMol(mol))
    AllChem.Compute2DCoords(mol_2d)
    lig_mol = plf.Molecule(mol)

    fp = plf.Fingerprint()
    fp.run_from_iterable([lig_mol], prot_mol)
    ifp = fp.ifp[0]
    if not ifp:
        return False

    canvas_w, canvas_h = 500, 450
    drawer = rdMolDraw2D.MolDraw2DCairo(canvas_w, canvas_h)
    drawer.drawOptions().addStereoAnnotation = False
    drawer.DrawMolecule(mol_2d)
    drawer.FinishDrawing()
    mol_img = Image.open(io.BytesIO(drawer.GetDrawingText()))

    atom_coords = {}
    for i in range(mol_2d.GetNumAtoms()):
        pt = drawer.GetDrawCoords(i)
        atom_coords[i] = (pt.x / canvas_w, 1.0 - pt.y / canvas_h)

    img_x0, img_x1 = 0.1, 0.9
    img_y0, img_y1 = 0.1, 0.9

    def to_ax(nx, ny):
        return img_x0 + nx * (img_x1 - img_x0), img_y0 + ny * (img_y1 - img_y0)

    all_coords = [to_ax(*atom_coords[i]) for i in range(mol_2d.GetNumAtoms())]
    cx = np.mean([c[0] for c in all_coords])
    cy = np.mean([c[1] for c in all_coords])

    entries = []
    for (lig_res, prot_res), interactions in ifp.items():
        res_label = f'{prot_res.name}{prot_res.number}'
        for itype, idata in interactions.items():
            orig_idx = idata[0]['indices']['ligand'][0]
            mapped_idx = orig_to_2d.get(orig_idx)
            if mapped_idx is None:
                atom = mol.GetAtomWithIdx(orig_idx)
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() != 1:
                        mapped_idx = orig_to_2d.get(nbr.GetIdx())
                        break
            if mapped_idx is None:
                continue
            ax_x, ax_y = to_ax(*atom_coords[mapped_idx])
            color = color_map.get(itype, '#9E9E9E')
            entries.append((res_label, ax_x, ax_y, color, itype in hbond_types))

    fig, ax = plt.subplots(figsize=(10, 9))
    ax.imshow(mol_img, extent=[img_x0, img_x1, img_y0, img_y1],
              transform=ax.transAxes, aspect='auto', zorder=1)

    placed = []
    for res_label, ax_x, ax_y, color, is_hbond in entries:
        dx, dy = ax_x - cx, ax_y - cy
        norm = max(np.sqrt(dx**2 + dy**2), 0.01)
        lx = ax_x + 0.18 * dx / norm
        ly = ax_y + 0.18 * dy / norm
        for px, py in placed:
            while np.sqrt((lx - px)**2 + (ly - py)**2) < 0.08:
                lx += 0.03 * dx / norm
                ly += 0.03 * dy / norm
        lx, ly = np.clip(lx, 0.03, 0.97), np.clip(ly, 0.03, 0.97)
        placed.append((lx, ly))
        ax.annotate('', xy=(ax_x, ax_y), xytext=(lx, ly),
                    xycoords='axes fraction', textcoords='axes fraction',
                    arrowprops=dict(arrowstyle='-', color=color, lw=1.8,
                                    linestyle='-' if is_hbond else '--'))
        ax.text(lx, ly, res_label, transform=ax.transAxes,
                ha='center', va='center', fontsize=8, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.25', facecolor='white',
                          edgecolor=color, linewidth=1.5))

    patches = [mpatches.Patch(color='#2196F3', label='H-bond'),
               mpatches.Patch(color='#9E9E9E', label='VdW/Hydrophobic'),
               mpatches.Patch(color='#4CAF50', label='Pi-stacking')]
    ax.legend(handles=patches, loc='upper right', fontsize=9)
    ax.axis('off')
    mol_name = mol.GetProp('_Name') if mol.HasProp('_Name') else 'ligand'
    ax.set_title(mol_name, fontsize=11, fontweight='bold')
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    return True


def parse_args():
    p = argparse.ArgumentParser(
        description="Generate 2D interaction diagrams for docking hits using ProLIF."
    )
    p.add_argument('-i', '--input',    required=True, help='Input SDF file of docking hits')
    p.add_argument('-r', '--receptor', required=True, help='Receptor PDB file')
    p.add_argument('-o', '--output',   required=True, help='Output directory for PNG files')
    return p.parse_args()


def main():
    args = parse_args()
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    u = mda.Universe(args.receptor)
    prot_mol = plf.Molecule.from_mda(u.select_atoms('protein'))

    suppl = Chem.SDMolSupplier(args.input, removeHs=False)
    total = sum(1 for m in suppl if m is not None)
    suppl = Chem.SDMolSupplier(args.input, removeHs=False)

    for i, mol in enumerate(suppl):
        if mol is None:
            continue
        mol_name = mol.GetProp('_Name') if mol.HasProp('_Name') else f'lig_{i}'
        safe_name = mol_name.replace(':', '_').replace('/', '_')[:50]
        out_path = out_dir / f'{safe_name}.png'
        if out_path.exists():
            continue
        try:
            success = make_diagram(mol, prot_mol, str(out_path))
            status = 'OK' if success else 'SKIP'
        except Exception as e:
            status = f'FAIL ({e})'
        print(f'[{status}] {i+1}/{total} {safe_name}')

    print(f'Done. PNGs saved to {out_dir}')


if __name__ == '__main__':
    main()
