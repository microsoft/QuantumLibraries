from qdk_chemistry.widgets.jsmol_widget import JsmolWidget
from qdk_chemistry.geometry import Geometry

from rdkit.Chem import AllChem as Chem

class Molecule(Chem.Mol):
    """
    Molecule object for visualization and geometry generation
    """
    def __init__(self, mol: Chem.Mol):
        self.mol = mol

    @classmethod
    def from_smiles(cls, smiles: str):
        mol = Chem.MolFromSmiles(self.value.smiles)
        return cls(mol=mol)

    def geometry(self):
        return Geometry.from_mol(self.mol)

    def show(self):
        return JsmolWidget.from_mol(self.mol)
