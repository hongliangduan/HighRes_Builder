#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
from collections import deque
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors


def find_alpha_carbon(mol):
    """
    Find an amino-acid-like alpha carbon in the molecule:
    - a carbon atom that is connected to:
      * at least one N (single bond)
      * at least one carbonyl carbon C(=O)
    Returns: (ca_idx, n_idx, c_idx) or (None, None, None)
    """
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "C":
            continue
        ca_idx = atom.GetIdx()
        neighbors = atom.GetNeighbors()

        n_idx = None
        carbonyl_c_idx = None

        for nb in neighbors:
            sym = nb.GetSymbol()
            if sym == "N" and n_idx is None:
                n_idx = nb.GetIdx()
            elif sym == "C" and carbonyl_c_idx is None:
                # Is this carbon a C=O?
                for bond in nb.GetBonds():
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(nb)
                        if other.GetSymbol() == "O":
                            carbonyl_c_idx = nb.GetIdx()
                            break

        if n_idx is not None and carbonyl_c_idx is not None:
            return ca_idx, n_idx, carbonyl_c_idx

    return None, None, None


def analyze_one_smiles(smiles):
    """
    Use RDKit to compute for one SMILES:
      - RDKit_cLogP
      - RDKit_MolWt / RDKit_ExactMolWt
      - sidechain_count
      - has_aromatic_in_sidechain
      - sidechain_has_non_CNO
      - molecule_has_non_CNO
    If alpha carbon is not found or SMILES parsing fails, side-chain-related fields are None.
    """
    base_result = {
        "RDKit_cLogP": None,
        "RDKit_MolWt": None,
        "RDKit_ExactMolWt": None,
        "sidechain_count": None,
        "has_aromatic_in_sidechain": None,
        "sidechain_has_non_CNO": None,
        "molecule_has_non_CNO": None,
    }

    if pd.isna(smiles):
        return base_result

    smiles = str(smiles).strip()
    if not smiles:
        return base_result

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return base_result

    # First compute logP and molecular weights
    clogp = Crippen.MolLogP(mol)
    mw = Descriptors.MolWt(mol)
    exact_mw = Descriptors.ExactMolWt(mol)

    # Whether the whole molecule contains heavy atoms other than C/N/O
    has_non_cno_mol = any(
        (a.GetSymbol() not in ("C", "N", "O", "H")) and a.GetAtomicNum() > 1
        for a in mol.GetAtoms()
    )

    result = {
        "RDKit_cLogP": clogp,
        "RDKit_MolWt": mw,
        "RDKit_ExactMolWt": exact_mw,
        "sidechain_count": None,
        "has_aromatic_in_sidechain": None,
        "sidechain_has_non_CNO": None,
        "molecule_has_non_CNO": has_non_cno_mol,
    }

    # Find alpha carbon
    ca_idx, n_idx, c_idx = find_alpha_carbon(mol)
    if ca_idx is None:
        # Not an amino-acid-like backbone; skip side-chain analysis
        return result

    # Backbone atom set: N, CA, C(=O), and all O attached to that carbon
    backbone = {ca_idx, n_idx, c_idx}
    c_atom = mol.GetAtomWithIdx(c_idx)
    for nb in c_atom.GetNeighbors():
        if nb.GetSymbol() == "O":
            backbone.add(nb.GetIdx())

    # Side-chain roots from CA: exclude backbone and H
    ca_atom = mol.GetAtomWithIdx(ca_idx)
    side_roots = []
    for nb in ca_atom.GetNeighbors():
        j = nb.GetIdx()
        if j in backbone:
            continue
        if nb.GetSymbol() == "H":
            continue
        side_roots.append(j)
    side_roots = sorted(set(side_roots))

    # BFS over the whole side chain
    side_atoms = set()
    for root in side_roots:
        visited = {root}
        q = deque([root])
        while q:
            cur = q.popleft()
            side_atoms.add(cur)
            atom_cur = mol.GetAtomWithIdx(cur)
            for nb in atom_cur.GetNeighbors():
                j = nb.GetIdx()
                if j in backbone:
                    continue
                if nb.GetSymbol() == "H":
                    continue
                if j not in visited:
                    visited.add(j)
                    q.append(j)

    sidechain_count = len(side_roots)
    has_arom_side = any(mol.GetAtomWithIdx(i).GetIsAromatic() for i in side_atoms)
    has_non_cno_side = any(
        (mol.GetAtomWithIdx(i).GetSymbol() not in ("C", "N", "O", "H"))
        and mol.GetAtomWithIdx(i).GetAtomicNum() > 1
        for i in side_atoms
    )

    result.update(
        {
            "sidechain_count": sidechain_count,
            "has_aromatic_in_sidechain": has_arom_side,
            "sidechain_has_non_CNO": has_non_cno_side,
        }
    )
    return result


def main():
    # ===== Change this to your spreadsheet filename =====
    input_file = "chembl_smiles.xlsx"         
    output_file = "chembl_smiles_sidechain.xlsx"
    # ==================================

    # Read table (supports xlsx/csv)
    if input_file.lower().endswith((".xls", ".xlsx")):
        df = pd.read_excel(input_file)
    else:
        df = pd.read_csv(input_file)

    # Case-insensitive column name handling
    cols = {c.lower(): c for c in df.columns}
    if "smiles" not in cols:
        raise ValueError("SMILES column not found; please confirm the column name is 'SMILES'")
    smiles_col = cols["smiles"]

    # Compute row by row
    results = df[smiles_col].apply(analyze_one_smiles)
    res_df = pd.DataFrame(list(results))

    # Merge back to original table
    out_df = pd.concat([df, res_df], axis=1)

    # Save as Excel
    out_df.to_excel(output_file, index=False)
    print("Saved results to:", output_file)


if __name__ == "__main__":
    main()