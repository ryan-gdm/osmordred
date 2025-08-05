from rdkit import Chem
import osmordred as rd



if __name__ == "__main__":
    smiles = ['CCCO','CCCN','c1ccccc1']
    for s in smiles:
        mol = Chem.MolFromSmiles(s)
        results =  rd.CalcBEState(mol)
        print(list(results))
    
