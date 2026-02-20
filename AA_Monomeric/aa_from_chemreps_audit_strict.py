# -*- coding: utf-8 -*-
# aa_from_chemreps_audit_strict.py
#
# Parallel screening of “amino-acid monomers” (α/β/γ; discard δ+ ) from ChEMBL chemreps,
# auditing each entry and cross-referencing with PDB CCD.
# Export:
#   Monomers_OK (passed set) / Rejected_with_Reasons (rejected set + reasons)
#   All_Annotated (merged details) / Summary / Summary_by_Class
#
# Dependencies: rdkit, pandas, openpyxl

import argparse, gzip, sys, os, time, math
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdmolops, Descriptors, rdMolDescriptors
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem.inchi import MolToInchiKey
from rdkit import RDLogger; RDLogger.DisableLog("rdApp.*")

try:
    from rdkit.Chem.MolStandardize import rdMolStandardize
    HAS_STD = True
except Exception:
    HAS_STD = False

VERBOSE = False
def vlog(msg: str, *, force=False):
    if VERBOSE or force:
        print(msg, flush=True)

# ----------------- strict SMARTS (monomer backbone) -----------------
def build_strict_smarts(acid_only=False):
    # Free amine (N has ≥1 H, non-amide N) + aliphatic chain + terminal carboxylic acid/carboxylate
    # (acid_only: only OH)
    o_tail = "[O;H1]" if acid_only else "[O;H1,-1]"
    alpha = f"[N;X3;H1,H2,H3;!$([N]-C(=O)-*)]-[C;X4;H1]-C(=O){o_tail}"
    beta  = f"[N;X3;H1,H2,H3;!$([N]-C(=O)-*)]-[CH2]-[C;X4;H1]-C(=O){o_tail}"
    gamma = f"[N;X3;H1,H2,H3;!$([N]-C(=O)-*)]-[CH2]-[CH2]-[C;X4;H1]-C(=O){o_tail}"
    return (Chem.MolFromSmarts(alpha),
            Chem.MolFromSmarts(beta),
            Chem.MolFromSmarts(gamma))

# Fallback “path-based method” (only when --no-strict-smarts)
_AMIDE = Chem.MolFromSmarts("C(=O)N")
def count_amide_bonds(m: Chem.Mol) -> int:
    return len(m.GetSubstructMatches(_AMIDE))

def is_carboxyl_carbon(a: Chem.Atom) -> bool:
    if a.GetAtomicNum()!=6: return False
    m = a.GetOwningMol()
    has_d = has_s = False
    for nb in a.GetNeighbors():
        if nb.GetAtomicNum()!=8: continue
        b = m.GetBondBetweenAtoms(a.GetIdx(), nb.GetIdx())
        if not b: continue
        if b.GetBondType()==Chem.rdchem.BondType.DOUBLE: has_d=True
        elif b.GetBondType()==Chem.rdchem.BondType.SINGLE: has_s=True
    return has_d and has_s

def classify_path_based(mol: Chem.Mol, free_only=True) -> str|None:
    m = SaltRemover().StripMol(mol, dontRemoveEverything=True)
    m = Chem.RemoveHs(m)
    if count_amide_bonds(m) >= 2:  # avoid dipeptides/polypeptides
        return None
    cs = [a for a in m.GetAtoms() if is_carboxyl_carbon(a)]
    ns = [a for a in m.GetAtoms() if a.GetAtomicNum()==7]
    if free_only:
        ns = [n for n in ns if n.GetTotalNumHs()>0]
    best=None
    for c in cs:
        for n in ns:
            path = rdmolops.GetShortestPath(m, c.GetIdx(), n.GetIdx())
            Lb = len(path)-1
            if Lb==2:
                aidx=path[1]
                a=m.GetAtomWithIdx(aidx)
                if not (a.GetAtomicNum()==6 and
                        a.GetHybridization()==Chem.rdchem.HybridizationType.SP3 and
                        a.GetTotalNumHs()>=1 and not a.GetIsAromatic()):
                    continue
            if best is None or Lb<best: best=Lb
    if best==2: return "alpha"
    if best==3: return "beta"
    if best==4: return "gamma"
    return None

# ----------------- global hard filters (functional groups + size/complexity) -----------------
_SMARTS_ESTER   = Chem.MolFromSmarts("C(=O)O[#6]")       # ester/carbonate
_SMARTS_CARBOXY = Chem.MolFromSmarts("C(=O)[O;H1,-1]")   # carboxylic acid/carboxylate

# These thresholds are broadcast via global variables in workers
_G_MAX_AMIDE=1; _G_MAX_ESTER=0; _G_MAX_CARBOXY=2
_G_MAX_MW=-1.0; _G_MAX_HEAVY=-1; _G_MAX_RINGS=-1; _G_MAX_HETERO=-1; _G_MAX_OXY=-1; _G_MAX_CHIRAL=-1

def count_smarts(m: Chem.Mol, patt: Chem.Mol) -> int:
    try:
        return len(m.GetSubstructMatches(patt))
    except Exception:
        return 0

def calc_complexity_metrics(m: Chem.Mol):
    # Compute after salt stripping and removing explicit H
    mm = SaltRemover().StripMol(m, dontRemoveEverything=True)
    mm = Chem.RemoveHs(mm)
    heavy  = Descriptors.HeavyAtomCount(mm)
    mw     = Descriptors.ExactMolWt(mm)
    rings  = rdMolDescriptors.CalcNumRings(mm)
    hetero = sum(1 for a in mm.GetAtoms() if a.GetAtomicNum() not in (1,6))
    oxy    = sum(1 for a in mm.GetAtoms() if a.GetAtomicNum()==8)
    chiral = len(Chem.FindMolChiralCenters(mm, includeUnassigned=True))
    return mm, heavy, mw, rings, hetero, oxy, chiral

def passes_global_filters(m: Chem.Mol):
    reasons=[]
    mm, heavy, mw, rings, hetero, oxy, chiral = calc_complexity_metrics(m)
    n_amide = count_amide_bonds(mm)                # C(=O)N
    n_ester = count_smarts(mm, _SMARTS_ESTER)      # C(=O)O-C
    n_carb  = count_smarts(mm, _SMARTS_CARBOXY)    # C(=O)O(H/-)

    if _G_MAX_AMIDE  >=0 and n_amide > _G_MAX_AMIDE:    reasons.append(f"AMIDE>{_G_MAX_AMIDE}({n_amide})")
    if _G_MAX_ESTER  >=0 and n_ester > _G_MAX_ESTER:    reasons.append(f"ESTER>{_G_MAX_ESTER}({n_ester})")
    if _G_MAX_CARBOXY>=0 and n_carb  > _G_MAX_CARBOXY:  reasons.append(f"CARBOXY>{_G_MAX_CARBOXY}({n_carb})")
    if _G_MAX_MW     > 0 and mw     > _G_MAX_MW:        reasons.append(f"MW>{_G_MAX_MW:.1f}({mw:.1f})")
    if _G_MAX_HEAVY  > 0 and heavy  > _G_MAX_HEAVY:     reasons.append(f"HEAVY>{_G_MAX_HEAVY}({heavy})")
    if _G_MAX_RINGS  > 0 and rings  > _G_MAX_RINGS:     reasons.append(f"RINGS>{_G_MAX_RINGS}({rings})")
    if _G_MAX_HETERO > 0 and hetero > _G_MAX_HETERO:    reasons.append(f"HETERO>{_G_MAX_HETERO}({hetero})")
    if _G_MAX_OXY    > 0 and oxy    > _G_MAX_OXY:       reasons.append(f"OXY>{_G_MAX_OXY}({oxy})")
    if _G_MAX_CHIRAL > 0 and chiral > _G_MAX_CHIRAL:    reasons.append(f"CHIRAL>{_G_MAX_CHIRAL}({chiral})")

    ok = (len(reasons)==0)
    return ok, reasons, n_amide, n_ester, n_carb, heavy, mw, rings, hetero, oxy, chiral

# ----------------- CCD -----------------
def _get_prop_ci(mol: Chem.Mol, keys: list[str]) -> str|None:
    for k in keys:
        if mol.HasProp(k):
            v=mol.GetProp(k).strip()
            if v: return v
    return None

def load_ccd_maps(sdf_gz: str, build_conn14=True):
    ids_by_ik, ids_by_ik14 = {}, {}
    vlog(f"[CCD] Reading {sdf_gz} ...", force=True)
    with gzip.open(sdf_gz, "rb") as fh:
        suppl = Chem.ForwardSDMolSupplier(fh, sanitize=False, removeHs=False)
        n=0
        for mol in suppl:
            n+=1
            if mol is None: continue
            comp_id = _get_prop_ci(mol, ["id","ID","pdbx_id","chem_comp.id"]) \
                      or (mol.GetProp("_Name") if mol.HasProp("_Name") else None)
            if not comp_id: continue
            ik = _get_prop_ci(mol, ["InChIKey","INCHIKEY","InchiKey","INCHI_KEY"])
            if not ik:
                try: Chem.SanitizeMol(mol)
                except Exception: pass
                try: ik = MolToInchiKey(mol)
                except Exception: continue
            if not ik or len(ik)<14: continue
            ids_by_ik.setdefault(ik, set()).add(comp_id)
            if build_conn14:
                ids_by_ik14.setdefault(ik[:14], set()).add(comp_id)
        vlog(f"[CCD] Done. Entries traversed: {n}", force=True)
    ids_by_ik = {k: sorted(v) for k,v in ids_by_ik.items()}
    ids_by_ik14 = {k: sorted(v) for k,v in ids_by_ik14.items()} if build_conn14 else {}
    return ids_by_ik, ids_by_ik14

def mol_to_inchikey_robust(m: Chem.Mol, remove_salts=True, standardize=True) -> str:
    if m is None: return ""
    mm = m
    if remove_salts:
        mm = SaltRemover().StripMol(mm, dontRemoveEverything=True)
    if standardize and HAS_STD:
        try:
            mm = rdMolStandardize.Cleanup(mm)
            te = rdMolStandardize.TautomerEnumerator()
            mm = te.Canonicalize(mm)
        except Exception:
            pass
    try: Chem.SanitizeMol(mm)
    except Exception: pass
    try: return MolToInchiKey(mm)
    except Exception: return ""

# ----------------- parallel (shared objects) -----------------
_IK_MAP=None; _IK14_MAP=None
_STRICT_ALPHA=None; _STRICT_BETA=None; _STRICT_GAMMA=None
_USE_STRICT=True

def _init_worker(ik_map, ik14_map, alpha_smarts, beta_smarts, gamma_smarts, use_strict,
                 g_max_amide, g_max_ester, g_max_carboxy,
                 g_max_mw, g_max_heavy, g_max_rings, g_max_hetero, g_max_oxy, g_max_chiral):
    global _IK_MAP, _IK14_MAP, _STRICT_ALPHA, _STRICT_BETA, _STRICT_GAMMA, _USE_STRICT
    global _G_MAX_AMIDE, _G_MAX_ESTER, _G_MAX_CARBOXY
    global _G_MAX_MW, _G_MAX_HEAVY, _G_MAX_RINGS, _G_MAX_HETERO, _G_MAX_OXY, _G_MAX_CHIRAL
    _IK_MAP = ik_map; _IK14_MAP = ik14_map
    _STRICT_ALPHA = alpha_smarts; _STRICT_BETA = beta_smarts; _STRICT_GAMMA = gamma_smarts
    _USE_STRICT = use_strict
    _G_MAX_AMIDE   = g_max_amide
    _G_MAX_ESTER   = g_max_ester
    _G_MAX_CARBOXY = g_max_carboxy
    _G_MAX_MW      = g_max_mw
    _G_MAX_HEAVY   = g_max_heavy
    _G_MAX_RINGS   = g_max_rings
    _G_MAX_HETERO  = g_max_hetero
    _G_MAX_OXY     = g_max_oxy
    _G_MAX_CHIRAL  = g_max_chiral
    RDLogger.DisableLog("rdApp.*")

def classify_strict(m: Chem.Mol):
    mm = SaltRemover().StripMol(m, dontRemoveEverything=True)
    mm = Chem.RemoveHs(mm)
    hits = []
    if _STRICT_ALPHA and mm.HasSubstructMatch(_STRICT_ALPHA): hits.append("alpha")
    if _STRICT_BETA  and mm.HasSubstructMatch(_STRICT_BETA):  hits.append("beta")
    if _STRICT_GAMMA and mm.HasSubstructMatch(_STRICT_GAMMA): hits.append("gamma")
    return hits[0] if hits else None, hits   # priority: alpha > beta > gamma

def _process_chunk(lines, idx_id, idx_sm,
                   use_strict, free_only, remove_salts, standardize,
                   allow_conn14, keep_all_ccd):
    rows=[]; n_total=n_parsed=0
    for line in lines:
        n_total+=1
        parts=line.rstrip("\n").split("\t")
        if len(parts)<=max(idx_id, idx_sm): continue
        cid=parts[idx_id]; smi=parts[idx_sm]
        m=Chem.MolFromSmiles(smi)
        if not m: continue
        n_parsed+=1

        # classification
        if use_strict:
            clz, hits = classify_strict(m)
            class_ok = clz is not None
            multi_hit = (len(hits) > 1)
        else:
            clz = classify_path_based(m, free_only=True if free_only else False)
            class_ok = clz is not None
            multi_hit = False

        if not class_ok:
            continue  # audit only alpha/beta/gamma candidates

        # global hard filters (functional groups + size/complexity)
        ok, reasons, n_amide, n_ester, n_carb, heavy, mw, rings, hetero, oxy, chiral = passes_global_filters(m)
        if multi_hit:
            reasons.append("MULTI_CLASS(" + ",".join(hits) + ")")

        # InChIKey + CCD
        ik = mol_to_inchikey_robust(m, remove_salts=remove_salts, standardize=standardize)
        ccd_ids=[]; match_type=""
        if ik:
            if ik in _IK_MAP:
                ccd_ids=_IK_MAP[ik]; match_type="exact"
            elif allow_conn14 and ik[:14] in _IK14_MAP and ik[:14]:
                ccd_ids=_IK14_MAP[ik[:14]]; match_type="conn14"
        if not keep_all_ccd and ccd_ids:
            ccd_ids=[ccd_ids[0]]

        rows.append({
            "CHEMBL_ID": cid,
            "SMILES": smi,
            "Class": clz,                      # alpha / beta / gamma
            "OK": 1 if ok else 0,
            "FailReasons": ";".join(reasons) if reasons else "",
            "n_amide": n_amide,
            "n_ester": n_ester,
            "n_carboxyl": n_carb,
            "heavy": heavy,
            "MW": round(mw, 4),
            "rings": rings,
            "hetero": hetero,
            "oxy": oxy,
            "chiral_centers": chiral,
            "InChIKey": ik,
            "CCD": ";".join(ccd_ids) if ccd_ids else "",
            "CCD_match_type": match_type
        })
    return rows, n_total, n_parsed

# ----------------- main workflow -----------------
def main():
    global VERBOSE
    ap=argparse.ArgumentParser(description="Audit + filter alpha/beta/gamma amino-acid monomers from ChEMBL chemreps, map to CCD (parallel).")
    ap.add_argument("--chemreps", required=True, help="chembl_XX_chemreps.txt.gz")
    ap.add_argument("--ccd", required=True, help="components-pub.sdf.gz")
    ap.add_argument("--out", default="AA_monomers_audited.xlsx")
    # behavior controls
    ap.add_argument("--no-strict-smarts", action="store_true", help="Disable strict SMARTS (fallback to topology path classification; looser)")
    ap.add_argument("--acid-only", action="store_true", help="Carboxyl group must be -C(=O)OH only (no deprotonated form)")
    ap.add_argument("--free-only", action="store_true", help="Path-based only: require at least one H on amine nitrogen")
    ap.add_argument("--no-conn14", action="store_true", help="Disable CCD 14-char connectivity-layer fallback")
    ap.add_argument("--keep-all-ccd", action="store_true", help="Keep all CCD candidates (default: pick the first)")
    ap.add_argument("--no-remove-salts", action="store_true")
    ap.add_argument("--no-standardize", action="store_true")
    # functional-group thresholds
    ap.add_argument("--max-amide-bonds", type=int, default=1, help="Allowed number of C(=O)N in the whole molecule (default 1)")
    ap.add_argument("--max-esters", type=int, default=0, help="Allowed number of C(=O)O-C in the whole molecule (default 0)")
    ap.add_argument("--max-carboxyl", type=int, default=2, help="Allowed number of carboxylic acids/carboxylates (default 2)")
    # size/complexity thresholds (enabled when >0; <=0 disables that criterion)
    ap.add_argument("--max-mw", type=float, default=-1.0, help="Molecular weight upper bound, e.g. 350 (<=0 means no limit)")
    ap.add_argument("--max-heavy", type=int, default=-1, help="Heavy atom count upper bound (<=0 means no limit)")
    ap.add_argument("--max-rings", type=int, default=-1, help="Ring count upper bound (<=0 means no limit)")
    ap.add_argument("--max-hetero", type=int, default=-1, help="Hetero atom count upper bound (<=0 means no limit)")
    ap.add_argument("--max-oxy", type=int, default=-1, help="Oxygen atom count upper bound (<=0 means no limit)")
    ap.add_argument("--max-chiral", type=int, default=-1, help="Chiral center count upper bound (<=0 means no limit)")
    # parallel and logging
    ap.add_argument("--sample", type=int, default=0, help="Only process first N lines (debug)")
    ap.add_argument("--progress-every", type=int, default=200000)
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--jobs", type=int, default=max(1,(os.cpu_count() or 2)-1))
    ap.add_argument("--chunk-size", type=int, default=20000)
    args=ap.parse_args()
    VERBOSE=bool(args.verbose)

    print(">> start aa_from_chemreps_audit_strict.py", flush=True)
    vlog(f"[ARGS] {args}", force=True)

    chemreps=Path(args.chemreps); ccd=Path(args.ccd)
    if not chemreps.exists(): print(f"[ERROR] chemreps not found: {chemreps}", flush=True); sys.exit(1)
    if not ccd.exists(): print(f"[ERROR] CCD not found: {ccd}", flush=True); sys.exit(1)

    # CCD mapping
    print("Loading CCD mapping...", flush=True)
    t0=time.time()
    ik_map, ik14_map = load_ccd_maps(str(ccd), build_conn14=(not args.no_conn14))
    print(f"  CCD exact InChIKey count: {len(ik_map)}; conn14: {len(ik14_map) if ik14_map else 0}", flush=True)
    vlog(f"[CCD] Build time {time.time()-t0:.1f}s", force=True)

    # prepare SMARTS
    alpha_s, beta_s, gamma_s = build_strict_smarts(acid_only=args.acid_only)
    use_strict = not args.no_strict_smarts

    # read chemreps header
    with gzip.open(str(chemreps), "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    hmap={h.strip().lower():i for i,h in enumerate(header)}
    def find_col(cands):
        for c in cands:
            if c in hmap: return hmap[c]
        return None
    idx_id=find_col(["chembl_id","molecule_chembl_id","chemblid"])
    idx_sm=find_col(["canonical_smiles","smiles","standard_smiles"])
    if idx_id is None or idx_sm is None:
        print("[ERROR] chemreps header mismatch:", ", ".join(header), flush=True); sys.exit(1)

    # parallel params
    remove_salts = not args.no_remove_salts
    standardize = not args.no_standardize
    allow_conn14 = not args.no_conn14
    keep_all_ccd = bool(args.keep_all_ccd)

    rows_all=[]; total_read=total_parsed=0
    step=max(10000, int(args.progress_every))
    chunk_size=max(1000, int(args.chunk_size))
    jobs=max(1, int(args.jobs))

    print(f"Scanning chemreps (strict monomer α/β/γ audit). Parallel {jobs} processes, chunk {chunk_size} lines...", flush=True)

    futures=[]
    with gzip.open(str(chemreps), "rt", encoding="utf-8", errors="replace") as fh:
        _=fh.readline()
        buf=[]; submitted=0
        if jobs==1:
            _init_worker(ik_map, ik14_map, alpha_s, beta_s, gamma_s, use_strict,
                         args.max_amide_bonds, args.max_esters, args.max_carboxyl,
                         args.max_mw, args.max_heavy, args.max_rings, args.max_hetero, args.max_oxy, args.max_chiral)
            for line in fh:
                if args.sample and total_read>=args.sample: break
                buf.append(line); total_read+=1
                if len(buf)>=chunk_size:
                    rows, nt, np = _process_chunk(
                        buf, idx_id, idx_sm,
                        use_strict, args.free_only,
                        remove_salts, standardize,
                        allow_conn14, keep_all_ccd
                    )
                    rows_all.extend(rows); total_parsed+=np
                    buf=[]
                    if total_read%step==0: vlog(f"  Progress: {total_read} lines...", force=True)
            if buf:
                rows, nt, np = _process_chunk(
                    buf, idx_id, idx_sm,
                    use_strict, args.free_only,
                    remove_salts, standardize,
                    allow_conn14, keep_all_ccd
                )
                rows_all.extend(rows); total_parsed+=np
        else:
            with ProcessPoolExecutor(
                max_workers=jobs,
                initializer=_init_worker,
                initargs=(ik_map, ik14_map, alpha_s, beta_s, gamma_s, use_strict,
                          args.max_amide_bonds, args.max_esters, args.max_carboxyl,
                          args.max_mw, args.max_heavy, args.max_rings, args.max_hetero, args.max_oxy, args.max_chiral)
            ) as ex:
                for line in fh:
                    if args.sample and total_read>=args.sample: break
                    buf.append(line); total_read+=1
                    if len(buf)>=chunk_size:
                        fut = ex.submit(
                            _process_chunk, list(buf), idx_id, idx_sm,
                            use_strict, args.free_only, remove_salts, standardize,
                            allow_conn14, keep_all_ccd
                        )
                        futures.append(fut); submitted+=1; buf=[]
                        if total_read%step==0: vlog(f"  Submitted chunks: {submitted}; total lines: {total_read}", force=True)
                if buf:
                    fut = ex.submit(
                        _process_chunk, list(buf), idx_id, idx_sm,
                        use_strict, args.free_only, remove_salts, standardize,
                        allow_conn14, keep_all_ccd
                    )
                    futures.append(fut); submitted+=1
                done=0
                for fut in as_completed(futures):
                    rows, nt, np = fut.result()
                    rows_all.extend(rows); total_parsed+=np
                    done+=1
                    if VERBOSE or done%max(1, math.ceil(len(futures)/10))==0:
                        vlog(f"  Completed chunks: {done}/{len(futures)}; parsed: {total_parsed}", force=True)

    print(f"Audit finished: total read {total_read} lines; candidate records {len(rows_all)} (all are α/β/γ candidates).", flush=True)

    # --- write Excel ---
    df = pd.DataFrame(rows_all)
    out = Path(args.out)
    print("Writing output to Excel...", flush=True)
    with pd.ExcelWriter(out, engine="openpyxl") as w:
        if not df.empty:
            df_all = df.drop_duplicates(subset=["CHEMBL_ID","SMILES","Class"])
            df_ok  = df_all[df_all["OK"]==1].copy()
            df_bad = df_all[df_all["OK"]==0].copy()

            cols = ["CHEMBL_ID","SMILES","Class","OK","FailReasons",
                    "n_amide","n_ester","n_carboxyl","heavy","MW","rings","hetero","oxy","chiral_centers",
                    "InChIKey","CCD","CCD_match_type"]

            df_ok.sort_values(["Class","CHEMBL_ID"])[cols].to_excel(w, index=False, sheet_name="Monomers_OK")
            df_bad.sort_values(["Class","CHEMBL_ID"])[cols].to_excel(w, index=False, sheet_name="Rejected_with_Reasons")
            df_all.sort_values(["Class","CHEMBL_ID"])[cols].to_excel(w, index=False, sheet_name="All_Annotated")

            # summary
            by_cls = df_ok["Class"].value_counts().rename_axis("Class").reset_index(name="Count_OK")
            hits = df_ok.assign(_hit=df_ok["CCD"].ne("").astype(int)).groupby("Class")["_hit"].sum().rename("CCD_Hits").reset_index()
            summary_by_class = pd.merge(by_cls, hits, on="Class", how="left").fillna({"CCD_Hits":0}).astype({"CCD_Hits":int})
            summary_by_class.to_excel(w, index=False, sheet_name="Summary_by_Class")

            summary = pd.DataFrame({
                "Metric": [
                    "Total_lines_read",
                    "Candidates(alpha+beta+gamma)",
                    "OK_total",
                    "CCD_Matched_in_OK"
                ],
                "Count": [
                    int(total_read),
                    int(len(df_all)),
                    int(len(df_ok)),
                    int(df_ok["CCD"].ne("").sum())
                ]
            })
            summary.to_excel(w, index=False, sheet_name="Summary")
        else:
            pd.DataFrame(columns=["CHEMBL_ID","SMILES","Class","OK","FailReasons","InChIKey","CCD","CCD_match_type"]).to_excel(
                w, index=False, sheet_name="All_Annotated"
            )
            pd.DataFrame({"Metric":["Total_lines_read","Candidates(alpha+beta+gamma)","OK_total","CCD_Matched_in_OK"],
                          "Count":[int(total_read),0,0,0]}).to_excel(w, index=False, sheet_name="Summary")

    print(f"OK -> {out.resolve()}", flush=True)

if __name__ == "__main__":
    main()