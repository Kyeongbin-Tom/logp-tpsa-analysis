from rdkit import Chem
from rdkit.Chem import Descriptors

def smiles_to_mol(smiles):
    return Chem.MolFromSmiles(smiles)

def smiles_to_descriptors(smiles):
    mol = smiles_to_mol(smiles)
    if mol is None:
        return None

    return {
        "TPSA": Descriptors.TPSA(mol),
        "MolWt": Descriptors.MolWt(mol),
        "NumHAcceptors": Descriptors.NumHAcceptors(mol),
        "NumHDonors": Descriptors.NumHDonors(mol),
        "NumRotatableBonds": Descriptors.NumRotatableBonds(mol),
    }
