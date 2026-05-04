# HitRefine

Post-docking hit refinement tool. Filters docking hits by physicochemical properties, structural alerts, drug-likeness rules, and similarity to reference compounds - extending what tools like Ringtail cannot do.

## Installation

```bash
git clone https://github.com/peteroryhall/hitrefine.git
cd hitrefine
micromamba env create -f environment.yml
micromamba activate hitrefine
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
| `--logp` | Max logP |
| `--tpsa` | Max TPSA (Å²) |
| `--hba` | Max H-bond acceptors |
| `--hbd` | Max H-bond donors |
| `--rotbonds` | Max rotatable bonds |
| `--fsp3-min` | Min Fsp3 fraction |
| `--lipinski` | Apply Lipinski Ro5 |
| `--veber` | Apply Veber rules |
| `--egan` | Apply Egan rules |
| `--qed` | Min QED score |
| `--sa-max` | Max SA score (1=easy, 10=hard) |
| `--pains` | Remove PAINS compounds |
| `--brenk` | Remove Brenk alert compounds |
| `--ref` | Reference SDF for Tanimoto similarity |
| `--sim-min` | Min Tanimoto similarity to reference |
| `--sim-max` | Max Tanimoto similarity to reference |
| `--plots` | Generate property distribution plots |
| `--csv` | Custom path for CSV output |
