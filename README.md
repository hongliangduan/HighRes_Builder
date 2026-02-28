# HighRes_Builder

HighRes_Builder is a small toolkit that **builds AlphaFold3-compatible inputs** for **noncanonical / CCD-absent residue-like molecules** by generating:

- **CCD-compatible residue topology** (atoms + bonds) and **residue-style atom naming**
- A reasonable set of **“ideal” Cartesian coordinates** (from RDKit conformer generation)
- **CSV outputs** that can be loaded as custom chemical components in a local AlphaFold3 setup

> **Important clarification (method scope):** HighRes_Builder **does not transplant Cartesian coordinates from CCD templates into the generated residue**. CCD is used for **identity/matching** and **conventions** (topology/atom naming); coordinates for generated residue definitions come from this pipeline (RDKit conformer generation and optional refinement).

---

## Quick start

```bash
# 0) Clone
git clone https://github.com/hongliangduan/HighRes_Builder.git
cd HighRes_Builder

# 1) Create env
conda create -n highres_builder python=3.11 -y
conda activate highres_builder
conda install -c conda-forge rdkit pandas openpyxl -y

# 2) Run screening + CCD mapping (ChEMBL 36)
cd AA_Monomeric
python aa_from_chemreps_audit_strict.py --help
```

---

## Repository layout

```text
HighRes_Builder/
├─ AA_Monomeric/
│  ├─ aa_from_chemreps_audit_strict.py     # ChEMBL screening + CCD mapping (audited Excel output)
│  └─ AA_monomers_audited.xlsx             # example output
├─ Smile_to_SDF/
│  └─ excel_smiles_to_individual_sdf_3d.py # RDKit 3D conformer generation (SDF)
├─ conformation_to_ccd/
│  └─ CCD_generation.py                    # SDF → CCD-like fields (CSV)
└─ README.md
```

---

## What “residue-like” means (chemical scope)

This repo targets **amino-acid-like monomers** that can be embedded into peptides/proteins and represented as AlphaFold3 chemical components.

Operational scope (Step 1):

- **Backbone class:** α / β / γ only (discard δ+)
- **Required motifs:** a **free amine** (non-amide N with ≥1 H) connected by an aliphatic chain to a **terminal carboxyl group** (–C(=O)O(H/–))
- **Hard exclusions (by default):**
  - peptides/oligomers (excess amide bonds)
  - esters/carbonates (masked acids / prodrug-like forms)
  - multi-acid species beyond the configured limit
  - metal/coordination complexes and other chemotypes outside residue-definition scope
  - overly large/complex molecules (configurable thresholds)

These rules are intentionally strict to keep the dataset in a regime where residue-definition construction and downstream modeling are stable and interpretable. See `AA_Monomeric/aa_from_chemreps_audit_strict.py --help` for all thresholds.

## User guidance: applicability and tuning

HighRes_Builder is intentionally scoped to **peptide-embeddable, amino-acid-like monomers**. If you are unsure whether your chemistry is in scope, use this checklist.

### Intended targets (good fit)
- **α/β/γ amino acids** with side-chain substitutions (halogenation, methylation, hydroxylation, thioether/selenide, etc.)
- Common **residue-like modifications** that still form standard peptide bonds (free amine + terminal carboxyl group)
- Molecules with **explicit stereochemistry** in SMILES when chirality matters (recommended)
- We curated residue-like molecules from ChEMBL 36, using the chemreps bulk file chembl_36_chemreps.txt.gz (canonical SMILES). The CCD reference for mapping is the public wwPDB Chemical Component Dictionary SDF dump (components-pub.sdf.gz).
Reproducible inclusion/exclusion checklist:
Parsing & normalization: RDKit parsing of canonical SMILES; salt stripping; optional RDKit standardization/tautomer canonicalization prior to InChIKey generation.
Residue-like scaffold detection: strict SMARTS identification of an amino-acid-like scaffold with a free amine (non-amide N with ≥1 H) connected via an aliphatic chain to a terminal carboxyl group; classification restricted to α/β/γ (δ+ discarded).
Hard functional-group filters: exclude molecules with C(=O)N (amide) > 1, C(=O)O–C (ester/carbonate) > 0, or carboxyl groups > 2; ambiguous multi-class hits are excluded and logged.
CCD mapping (InChIKey logic): map each candidate to CCD by (i) exact full InChIKey, then (ii) optional connectivity-layer fallback (first 14 characters) to tolerate protonation/tautomer differences.

### When to expect exclusions
- **Peptides/oligomers** (multiple amide bonds) — not a single-residue definition
- **Masked acids** (esters/carbonates) — often prodrug/protected intermediates
- **Multi-acid** species (> configured carboxyl limit) and highly charged complexes
- **Metal/coordination** complexes and inorganic/organometallic species
- **Large / highly complex** scaffolds (high MW, many rings/hetero atoms), which are more likely to fail conformer generation or produce ambiguous residue definitions
- **Unspecified stereochemistry** in SMILES — RDKit may choose an arbitrary assignment

### How to broaden chemical coverage (tradeoffs)
You can relax filters in `AA_Monomeric/aa_from_chemreps_audit_strict.py` via flags. Broader coverage typically increases:
- RDKit embedding/optimization failures,
- ambiguous “residue-like” matches,
- and noisy chemotypes that are less meaningful as peptide residues.

Practical approach:
1) Run with the **default strict** settings.
2) Inspect `Rejected_with_Reasons` to see dominant exclusion tags (e.g., `MW>`, `RINGS>`, `HETERO>`).
3) Relax **one threshold at a time**, re-run, and record the exact command/config.

Example: relax size/complexity limits (accept larger chemotypes)

```bash
cd AA_Monomeric
python aa_from_chemreps_audit_strict.py   --chemreps chembl_36_chemreps.txt.gz   --ccd components-pub.sdf.gz   --out AA_monomers_audited_relaxed.xlsx   --max-mw 500 --max-heavy 45 --max-rings 4 --max-hetero 15
```

You can also tighten/loosen chemistry constraints:
- `--max-amide-bonds`, `--max-esters`, `--max-carboxyl`
- `--no-conn14` (stricter CCD mapping; fewer ambiguous hits)
- `--keep-all-ccd` (retain all CCD candidates for manual review)

---

---

## End-to-end workflow

```text
ChEMBL chemreps (SMILES; e.g., chembl_36_chemreps.txt.gz)
   │
   ▼
AA_Monomeric: screening + auditing + CCD mapping
   │  output: AA_monomers_audited.xlsx
   │
   ├─ CCD hit → recommend using the official CCD component directly
   └─ no CCD  → proceed
          │
          ▼
Smile_to_SDF: RDKit 3D conformer generation (SDF)
          │  output: one SDF per CHEMBL_ID
          ▼
conformation_to_ccd: SDF → CCD-like fields (CSV)
          │  output: *_ccd.csv (mmCIF-like columns) + *_atom_names.csv
          ▼
Downstream AlphaFold3 modeling (local) / other workflows
```

---

## Step 1 — AA_Monomeric (screen + audit + CCD mapping)

**Goal:** screen monomeric residue-like molecules from **ChEMBL 36** chemreps, apply strict filters, and map to CCD when possible.

### Required inputs (not included)
- `chembl_36_chemreps.txt.gz` (ChEMBL 36 chemreps)
- `components-pub.sdf.gz` (PDB CCD SDF dump)

### Run (example)

```bash
cd AA_Monomeric

python aa_from_chemreps_audit_strict.py   --chemreps chembl_36_chemreps.txt.gz   --ccd components-pub.sdf.gz   --out AA_monomers_audited.xlsx   --jobs 20   --chunk-size 50000   --verbose   --acid-only   --max-amide-bonds 1   --max-esters 0   --max-carboxyl 2   --max-mw 350   --max-heavy 30   --max-rings 2   --max-hetero 10   --max-oxy 6   --max-chiral 4
```

### Output

`AA_monomers_audited.xlsx` (auditable spreadsheet):

- `Monomers_OK`: passed candidates
- `Rejected_with_Reasons`: rejected entries + machine-readable reasons
- `All_Annotated`: merged annotations
- `Summary` / `Summary_by_Class`: screening statistics

### CCD matching and ambiguity handling

CCD mapping uses a two-stage InChIKey strategy:

1) **Exact match** on the full 27-character InChIKey  
2) Optional **connectivity-layer fallback** using `InChIKey[:14]` (tolerates protonation/tautomer changes while preserving connectivity)

**Ambiguity:** multiple CCD components may share the same connectivity layer. By default, the script records all candidates in the audit table and selects one deterministically unless configured otherwise. Useful options:

- `--no-conn14`: disable connectivity-layer fallback (stricter; fewer hits)
- `--keep-all-ccd`: keep all CCD candidates instead of selecting only one

---

## Step 2 — Smile_to_SDF (RDKit conformer generation)

**Goal:** generate 3D conformers (SDF) for screened monomers that **do not have a corresponding CCD entry**.

```bash
cd Smile_to_SDF
python excel_smiles_to_individual_sdf_3d.py
```

Default behavior:
- reads `smile.xlsx` in the current directory
- writes one SDF per molecule to `sdf_out_3d/`
- failures are recorded to `sdf_failed.xlsx`

Embedding:
- ETKDGv3 with a fixed seed (default `2025`)
- `enforceChirality=True` to respect chirality encoded in SMILES  
  (If your SMILES does not specify stereochemistry, RDKit may choose an arbitrary assignment.)

---

## Step 3 — conformation_to_ccd (SDF → CCD-like CSV)

**Goal:** convert an SDF conformer into a **CCD-style component description** for downstream structure prediction.

```bash
cd conformation_to_ccd

# Put your SDF files under ./sdf
# Output will be written to ./output
python CCD_generation.py
```

Outputs per molecule:
- `<name>_atom_names.csv` — atom index → renamed atom name mapping
- `<name>_ccd.csv` — a single-row CSV containing mmCIF-like fields:
  - `_chem_comp.*`, `_chem_comp_atom.*`, `_chem_comp_bond.*`

---

## Using the generated CSV with local AlphaFold3

If you run AlphaFold3 locally and want it to load custom chemical components from CSV, you may need to patch AlphaFold3’s chemical-component loader (depending on your fork/version). This repo includes patched versions of:

- `./src/alphafold3/common/folding_input.py`
- `./src/alphafold3/constants/chemical_components.py`
- `./src/alphafold3/constants/residue_names.py`

Copy these files into the corresponding locations of your local AlphaFold3 checkout **only if needed**.

---

## Reproducibility and limitations

- **Determinism:** Step 2 uses a fixed RDKit embedding seed by default.
- **Chirality note:** HighRes_Builder preserves stereochemistry as encoded in the input where available, but final stereochemical correctness in predicted complexes can still depend on the downstream predictor. We recommend post-hoc stereochemical validation for chirality-critical applications.
- **Quality caveat:** a single conformer is often sufficient for “ideal coordinates” in CCD-like workflows, but hard chemotypes may benefit from multiple conformers and energy-based selection.

---

## Troubleshooting

- RDKit installation issues: prefer conda-forge builds
- Embedding failures: try enabling force-field optimization (UFF/MMFF) or increasing embedding attempts (script changes)
- Chirality ambiguities: ensure SMILES encodes stereochemistry; otherwise RDKit may assign chirality arbitrarily
