# HighRes_Builder

A small toolkit to build **high-quality monomer inputs** for **AlphaFold3** and **PDB CCD (Chemical Component Dictionary)** workflows.

This repository covers a practical pipeline:

1. **AA_Monomeric**: screen *amino-acid-like monomers* from **ChEMBL chemreps** (SMILES), audit/filter them, and cross-reference them against **PDB CCD**.
2. **Smile_to_SDF**: for screened monomers **without a CCD entry**, generate a reasonable **3D conformer** (SDF) using **RDKit**.
3. **conformation_to_ccd**: convert the **SDF conformer** into a **CCD-style component description** (mmCIF-like fields exported as CSV), for downstream structure prediction.

---

## Repository layout

```
HighRes_Builder/
├─ AA_Monomeric/
│  ├─ aa_from_chemreps_audit_strict.py
│  └─ AA_monomers_audited.xlsx               # example output
├─ Smile_to_SDF/
│  └─ excel_smiles_to_individual_sdf_3d.py
├─ conformation_to_ccd/
│  └─ CCD_generation.py
└─ README.md
```

---

## Dependencies

- Python 3.9+ (recommended: 3.10/3.11)
- RDKit
- pandas
- openpyxl (Excel I/O)

Example conda environment:

```bash
conda create -n highres_builder python=3.11 -y
conda activate highres_builder
conda install -c conda-forge rdkit pandas openpyxl -y
```

---

## End-to-end workflow

```
ChEMBL chemreps (SMILES)
   │
   ▼
AA_Monomeric: screening + auditing + CCD mapping
   │  output: AA_monomers_audited.xlsx
   │
   ├─ has CCD → can use CCD directly
   └─ no CCD  → proceed
          │
          ▼
Smile_to_SDF: RDKit 3D conformer generation (SDF)
          │  output: one SDF per CHEMBL_ID
          ▼
conformation_to_ccd: SDF → CCD-like fields (CSV)
          │  output: *_ccd.csv (mmCIF-like columns)
          ▼
Downstream structure prediction / modeling
```

---

## Step 1 — AA_Monomeric

**Goal:** screen *monomeric amino-acid-like* molecules (α/β/γ; discard δ+) from ChEMBL chemreps, apply strict filters, and map to CCD when possible.

### Run

```bash
cd AA_Monomeric

python aa_from_chemreps_audit_strict.py \
  --chemreps chembl_36_chemreps.txt.gz \
  --ccd components-pub.sdf.gz \
  --out AA_monomers_audited.xlsx \
  --jobs 20 \
  --chunk-size 50000 \
  --verbose \
  --acid-only \
  --max-amide-bonds 1 \
  --max-esters 0 \
  --max-carboxyl 2 \
  --max-mw 350 \
  --max-heavy 30 \
  --max-rings 2 \
  --max-hetero 10 \
  --max-oxy 6 \
  --max-chiral 4
```

### Output

The script writes `AA_monomers_audited.xlsx`. Typical sheets include:

- `Monomers_OK`: passed candidates
- `Rejected_with_Reasons`: rejected entries + explicit rejection reasons
- `All_Annotated`: merged annotations
- `Summary` / `Summary_by_Class`: screening statistics

### Algorithm overview

This script is designed to be **auditable** and **streaming-friendly**:

1. **Build CCD lookup tables**
   - Stream-read `components-pub.sdf.gz`.
   - Extract or compute **InChIKey** for each CCD component.
   - Build maps:
     - `InChIKey (27 chars) → [CCD_IDs]` (exact match)
     - optional `InChIKey[:14] → [CCD_IDs]` (connectivity-layer fallback)

2. **Stream ChEMBL chemreps**
   - Read the gzipped chemreps text in chunks (`--chunk-size`).
   - Parallelize chunk processing with `--jobs`.

3. **Strict backbone detection (α/β/γ amino-acid-like monomers)**
   - By default uses strict SMARTS patterns for:
     - **free amine** (not amide nitrogen)
     - **terminal carboxyl group** (optionally acid-only)
     - carbon-chain length consistent with α/β/γ

4. **Hard filters (functional groups + complexity)**
   - Count functional motifs (e.g., amide bonds, esters, carboxyl groups).
   - Reject molecules exceeding thresholds for: MW, heavy atoms, rings, hetero atoms, oxygens, chiral centers.

5. **Robust InChIKey generation + CCD mapping**
   - Optional salt stripping and standardization (when RDKit MolStandardize is available).
   - Prefer exact InChIKey match; optionally fall back to connectivity-layer match.

6. **Excel export for auditability**
   - Every molecule is annotated with class (alpha/beta/gamma), metrics, InChIKey, CCD hit status/type, and rejection reasons.

**Tip:** Step 2 typically consumes only `Monomers_OK` rows where `CCD` is empty (or the CCD match column indicates no hit).

---

## Step 2 — Smile_to_SDF

**Goal:** generate 3D conformers (SDF) for screened monomers that **do not have a corresponding CCD entry**.

### Run (defaults)

```bash
cd Smile_to_SDF
python excel_smiles_to_individual_sdf_3d.py
```

### Default input contract

- Reads `smile.xlsx` from the current directory (first sheet by default).
- Auto-detects columns:
  - ID: `CHEMBL_ID` / `molecule_chembl_id` / `id`
  - SMILES: `SMILES` / `canonical_smiles` / `standard_smiles`

### Output

- Writes one SDF per molecule to `sdf_out_3d/` (default).
- Failed items are recorded into `sdf_failed.xlsx`.

### Algorithm overview

For each SMILES:

1. Parse molecule: `Chem.MolFromSmiles`
2. Add explicit H: `Chem.AddHs`
3. Embed a 3D conformer with **ETKDGv3**
   - fixed random seed (default `2025`) for reproducibility
   - `enforceChirality=True` to respect chirality encoded in SMILES
4. Optional geometry optimization
   - UFF or MMFF (off by default)
5. Write individual SDF with traceable properties (`CHEMBL_ID`, `SMILES`)

> Note on “L amino acids”: RDKit will enforce **the stereochemistry present in the SMILES**. If you require strictly L-forms, ensure your SMILES encodes the desired stereocenter configuration.

---

## Step 3 — conformation_to_ccd

**Goal:** convert an SDF conformer to **CCD-style component information** for downstream structure prediction.

### Run (defaults)

```bash
cd conformation_to_ccd

# Put your SDF files under ./sdf
# Output will be written to ./output
python CCD_generation.py
```

### Output

For each `*.sdf` in `conformation_to_ccd/sdf/`, the script writes into `conformation_to_ccd/output/`:

- `<name>_atom_names.csv` — atom index → renamed atom name mapping
- `<name>_ccd.csv` — a single-row CSV containing mmCIF-like CCD fields:
  - `_chem_comp.*`
  - `_chem_comp_atom.*`
  - `_chem_comp_bond.*`

### Algorithm overview (what the script does)

This module follows a practical **"amino-acid-like naming + CCD field extraction"** strategy:

1. **Backbone detection (N, CA, C, O, OXT)**
   - First locate a **carboxyl-like carbon** (C) by oxygen bonding pattern.
   - Then BFS along **C–C bonds** up to 3 steps to find a **CA candidate** connected to an N.
     - distance 1/2/3 corresponds to α/β/γ backbones.
   - Score candidates (favor NH2-like N, sp3 CA, shorter distance) and pick the best.

2. **Side-chain naming (starting from CB)**
   - Choose CB among carbon neighbors of CA not in the backbone.
   - BFS outward purely by bond connectivity.
   - Use Greek-letter hierarchy for naming levels: `CB, CG, CD, CE, CZ, ...` (extended).
   - For hetero atoms, use `Element + level` patterns (e.g., `OG`, `OD1`, `NZ`).
   - Resolve branches with numeric suffixes based on subtree “weight”.

3. **AlphaFold3 compatibility checks**
   - Validate names against an ATOM37-like set.
   - Apply a conservative rule to label terminal `CH2` only when carbon has exactly 2 H and is terminal.
   - Globally uniquify atom names to avoid duplicates.

4. **CCD field export (mmCIF-like columns, written as CSV)**
   - Derive formula, molecular weight.
   - Export per-atom fields: charges, aromatic flags, atom IDs, and ideal Cartesian coordinates from the conformer.
   - Export per-bond fields: atom pairs, bond order (with kekulization fallback for aromatics), aromatic flags.
   - Order atoms “backbone-first” and bonds “backbone-bond-first” to resemble protein conventions.

---

## Practical notes

- **Input data files** (Step 1) are not included:
  - `chembl_36_chemreps.txt.gz` (ChEMBL chemreps)
  - `components-pub.sdf.gz` (PDB CCD SDF dump)

- **Reproducibility**: Step 2 uses a fixed ETKDG seed by default.

- **Quality caveat**: A single embedded conformer is often enough for CCD-style “ideal coordinates” in modeling pipelines, but you may want to generate multiple conformers and pick the lowest-energy one for hard cases.

---

## Troubleshooting

- **RDKit installation issues**: prefer conda-forge builds.
- **Embedding failures**: try enabling `--opt uff` or increasing RDKit embedding attempts (requires script changes).
- **Chirality ambiguities**: ensure your SMILES encodes stereochemistry; otherwise RDKit may produce arbitrary chiral assignments.

