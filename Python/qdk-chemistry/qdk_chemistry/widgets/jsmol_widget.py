from ..convert import mol_to_xyz

from jupyter_jsmol import JsmolView


class JsmolWidget(JsmolView):
    """Jupyter widget for JSMol molecular geometry visualization
    """
    def __init__(self, *args, **kwargs):
        self.default_info['color']='white'
        super().__init__(*args, **kwargs)

    @classmethod
    def from_mol(cls, mol: "Mol", num_confs: int = 10, *args, **kwargs):
        """Generate JsmolWidget object from RDKit molecule

        Args:
            mol (Mol): RDKit molecule to visualize
            num_confs (int, optional): Number of conformers to generate. Defaults to 10.

        Returns:
            JsmolWidet
        """
        xyz = mol_to_xyz(mol, num_confs=num_confs)
        return cls.from_str(xyz)
