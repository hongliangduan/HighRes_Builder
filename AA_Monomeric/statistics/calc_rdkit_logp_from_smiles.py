#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors


def calc_logp(smiles: str):
    """Compute cLogP with RDKit; return None if parsing fails."""
    if pd.isna(smiles):
        return None
    smiles = str(smiles).strip()
    if not smiles:
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    return Crippen.MolLogP(mol)


def calc_mw(smiles: str):
    """
    Compute molecular weights with RDKit:
    - RDKit_MolWt: average molecular weight (commonly referred to as molecular weight)
    - RDKit_ExactMolWt: exact mass
    """
    if pd.isna(smiles):
        return None, None
    smiles = str(smiles).strip()
    if not smiles:
        return None, None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None

    mw = Descriptors.MolWt(mol)          # average molecular weight
    exact_mw = Descriptors.ExactMolWt(mol)  # exact mass
    return mw, exact_mw


def main():
    # ===== Places you need to edit =====
    # Input table: must contain at least CHEMBL_ID and SMILES columns
    input_file = "chembl_smiles.xlsx"   # change to your filename

    # Output filename
    output_file = "chembl_with_rdkit_logp_mw.xlsx"
    # ========================

    # Read data
    if input_file.lower().endswith((".xls", ".xlsx")):
        df = pd.read_excel(input_file)
    else:
        df = pd.read_csv(input_file)

    # Unified column name lookup
    cols = {c.lower(): c for c in df.columns}
    if "smiles" not in cols:
        raise ValueError("SMILES column not found; please confirm the column name is 'SMILES'")
    if "chembl_id" not in cols:
        raise ValueError("CHEMBL_ID column not found; please confirm the column name is 'CHEMBL_ID'")

    smiles_col = cols["smiles"]
    chembl_col = cols["chembl_id"]

    # Compute logP
    df["RDKit_cLogP"] = df[smiles_col].apply(calc_logp)

    # Compute molecular weights
    mw_results = df[smiles_col].apply(calc_mw)
    df["RDKit_MolWt"] = mw_results.apply(lambda x: x[0])        # average molecular weight
    df["RDKit_ExactMolWt"] = mw_results.apply(lambda x: x[1])   # exact mass

    # Sort by CHEMBL_ID (optional)
    df = df.sort_values(chembl_col).reset_index(drop=True)

    # Write to Excel
    df.to_excel(output_file, index=False)
    print(f"RDKit cLogP and molecular weights written to: {output_file}")


if __name__ == "__main__":
    main()