from rdkit import Chem
import numpy as np
import pandas as pd
import osmordred as rd
import rdkit

def CalcOsmordred(mol: Chem.Mol) -> np.ndarray | None:
    try:
        # Calculate ExtendedTopochemicalAtom fingerprint
        result = rd.CalcExtendedTopochemicalAtom(mol)
        return np.atleast_1d(result)
    except Exception as e:
        print(f"Error processing molecule {mol}: {e}")
        return None

# Read input CSV
df = pd.read_csv("/osmordred/testosmordred.csv")
smiles_list = df["smiles"].tolist()

# Prepare results
results = []
for smiles in smiles_list:
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        fp = CalcOsmordred(mol)
        if fp is not None:
            # Save each fingerprint bit as a separate column
            result = {"smiles": smiles}
            for i, val in enumerate(fp):
                result[f"ExtendedTopochemicalAtom_{i+1}"] = val
            results.append(result)

# Create DataFrame and save
result_df = pd.DataFrame(results)
result_df.to_csv(f"testosmordred-{rdkit.__version__}.csv", index=False)
