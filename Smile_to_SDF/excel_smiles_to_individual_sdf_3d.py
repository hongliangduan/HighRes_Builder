import argparse, os, sys, re
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

def sanitize_name(s: str) -> str:
    return re.sub(r'[\\/:*?"<>|]+', "_", str(s).strip()) or "NONAME"

def build_3d(smi: str, seed: int, do_opt: str|None):
    """Return a 3D molecule (with explicit hydrogens). do_opt: None / 'uff' / 'mmff'"""
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None, "parse_failed"
    m = Chem.AddHs(m)  # 3D: keep explicit hydrogens
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    params.useRandomCoords = False
    params.enforceChirality = True  # preserve chirality in SMILES
    if AllChem.EmbedMolecule(m, params) != 0:
        return None, "embed_failed"

    if do_opt == "uff":
        try:
            AllChem.UFFOptimizeMolecule(m, maxIters=200)
        except Exception:
            pass
    elif do_opt == "mmff":
        try:
            if AllChem.MMFFHasAllMoleculeParams(m):
                AllChem.MMFFOptimizeMolecule(m, maxIters=200)
        except Exception:
            pass

    return m, ""

def main():
    ap = argparse.ArgumentParser(description="Excel(CHEMBL_ID/SMILES) -> individual RDKit 3D SDF (V2000, with H)")
    ap.add_argument("--in", dest="inp", default="smile.xlsx", help="Input Excel (default: smile.xlsx)")
    ap.add_argument("--sheet", default=None, help="Worksheet name; if not provided, read the first sheet")
    ap.add_argument("--id-col", default=None, help="ID column name (default: auto-detect CHEMBL_ID/molecule_chembl_id/id)")
    ap.add_argument("--smiles-col", default=None, help="SMILES column name (default: auto-detect SMILES/canonical_smiles/standard_smiles)")
    ap.add_argument("--out-dir", default="sdf_out_3d", help="SDF output directory (default: sdf_out_3d)")
    ap.add_argument("--seed", type=int, default=2025, help="ETKDG random seed (default: 2025, for reproducibility)")
    ap.add_argument("--opt", choices=["none","uff","mmff"], default="none", help="Geometry optimization: none/uff/mmff (default: none)")
    ap.add_argument("--skip-existing", action="store_true", help="Skip if the output file already exists")
    ap.add_argument("--dedup", action="store_true", help="Deduplicate by CHEMBL_ID and keep the first entry")
    ap.add_argument("--failed-xlsx", default="sdf_failed.xlsx", help="Failed items log (default: sdf_failed.xlsx)")
    args = ap.parse_args()


    if not os.path.exists(args.inp):
        print(f"[ERROR] Input file not found: {os.path.abspath(args.inp)}", file=sys.stderr); sys.exit(1)
    df = pd.read_excel(args.inp, sheet_name=(args.sheet if args.sheet else 0))


    low2orig = {str(c).strip().lower(): c for c in df.columns}
    id_col = args.id_col or low2orig.get("chembl_id") or low2orig.get("molecule_chembl_id") or low2orig.get("id")
    sm_col = args.smiles_col or low2orig.get("smiles") or low2orig.get("canonical_smiles") or low2orig.get("standard_smiles")
    if id_col is None or sm_col is None:
        print(f"[ERROR] CHEMBL_ID/SMILES columns not found. Actual columns: {', '.join(map(str, df.columns))}", file=sys.stderr); sys.exit(1)
    df = df.rename(columns={id_col:"CHEMBL_ID", sm_col:"SMILES"})

    if args.dedup:
        df = df.drop_duplicates(subset=["CHEMBL_ID"], keep="first").copy()

    os.makedirs(args.out_dir, exist_ok=True)

    failures = []
    n_ok = 0
    do_opt = None if args.opt=="none" else args.opt

    for _, r in df.iterrows():
        cid_raw = str(r.get("CHEMBL_ID","")).strip()
        smi = str(r.get("SMILES","")).strip()
        if not smi:
            failures.append((cid_raw, smi, "empty_smiles")); continue

        fname = os.path.join(args.out_dir, f"{sanitize_name(cid_raw)}.sdf")
        if args.skip_existing and os.path.exists(fname):
            continue

        m, why = build_3d(smi, args.seed, do_opt)
        if m is None:
            failures.append((cid_raw, smi, why)); continue

        m.SetProp("_Name", cid_raw if cid_raw else sanitize_name(cid_raw))
        m.SetProp("CHEMBL_ID", cid_raw)
        m.SetProp("SMILES", smi)

        w = Chem.SDWriter(fname)  
        w.write(m)
        w.close()
        n_ok += 1

    print(f"[DONE] Successfully wrote {n_ok} 3D SDF files -> {os.path.abspath(args.out_dir)}")
    if failures:
        pd.DataFrame(failures, columns=["CHEMBL_ID","SMILES","reason"]).to_excel(args.failed_xlsx, index=False)
        print(f"[INFO] Failed/skipped: {len(failures)} entries, logged to {os.path.abspath(args.failed_xlsx)}")

if __name__ == "__main__":
    main()