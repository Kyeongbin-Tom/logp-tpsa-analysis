from rdkit import Chem
from rdkit.Chem import Descriptors

def smiles_to_mol(smiles):
    """
    Convert a SMILES string to an RDKit Mol object.
    """
    return Chem.MolFromSmiles(smiles)

def calculate_tpsa(mol):
    """
    Calculate Topological Polar Surface Area (TPSA) from an RDKit Mol object.
    """
    return Descriptors.TPSA(mol)

def smiles_to_tpsa(smiles):
    """
    Convert SMILES to TPSA in one step.
    """
    mol = smiles_to_mol(smiles)
    if mol is None:
        return None
    return calculate_tpsa(mol)
