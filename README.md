# HitRefine

Post-docking hit refinement tool. Filters docking hits by physicochemical properties, structural alerts, drug-likeness rules, and similarity to reference compounds - extending what tools like Ringtail cannot do.

## Installation

```bash
git clone https://github.com/peteroryhall/hitrefine.git
cd hitrefine
micromamba create -n hitrefine python=3.10 -y
micromamba activate hitrefine
micromamba install -n hitrefine -c conda-forge rdkit matplotlib pandas -y
```

## Usage

```bash
python hitrefine.py -i hits.sdf -o filtered.sdf [options]
```

## Examples

Filter by logP and TPSA only:
```bash
python hitrefine.py -i hits.sdf -o filtered.sdf --logp 5 --tpsa 90
```

Full drug-likeness filter with structural alerts and plots:
```bash
python hitrefine.py -i hits.sdf -o filtered.sdf \
  --logp 5 --tpsa 90 --hba 7 --hbd 5 --rotbonds 10 \
  --lipinski --veber \
  --qed 0.3 --sa-max 6 \
  --pains --brenk \
  --plots
```

Filter by similarity to known reference hits:
```bash
python hitrefine.py -i hits.sdf -o filtered.sdf \
  --ref known_hits.sdf --sim-min 0.3 --sim-max 0.9
```

## Options

| Flag | Description |
|------|-------------|
| `--logp` | Max logP - lipophilicity; higher values reduce solubility and increase off-target binding |
| `--tpsa` | Max TPSA (Å²) - polar surface area; predicts membrane permeability |
| `--hba` | Max H-bond acceptors - too many reduces cell permeability |
| `--hbd` | Max H-bond donors - too many reduces membrane permeability |
| `--rotbonds` | Max rotatable bonds - high flexibility reduces oral bioavailability |
| `--fsp3-min` | Min Fsp3 - fraction of sp3 carbons; higher values indicate more 3D character and better developability |
| `--lipinski` | Apply Lipinski Ro5 - classic oral drug-likeness rule set |
| `--veber` | Apply Veber rules - oral bioavailability filter (TPSA + rotatable bonds) |
| `--egan` | Apply Egan rules - passive intestinal absorption filter |
| `--qed` | Min QED score - overall drug-likeness estimate (0=bad, 1=ideal) |
| `--sa-max` | Max SA score - synthetic accessibility (1=easy, 10=hard to synthesise) |
| `--pains` | Remove PAINS - pan-assay interference compounds that give false positives in screens |
| `--brenk` | Remove Brenk alerts - unstable or potentially toxic functional groups |
| `--ref` | Reference SDF for Tanimoto similarity |
| `--sim-min` | Min Tanimoto similarity to reference - keep compounds similar to known hits |
| `--sim-max` | Max Tanimoto similarity to reference - remove near-duplicates of known hits |
| `--plots` | Generate property distribution plots |
| `--csv` | Custom path for CSV output |
