# Copyright (c) Microsoft Corporation.
# Licensed under the MIT License.

"""Module for Jupyter widget that displays JSMol editor
"""
from ..geometry import mol_to_xyz

from jupyter_jsmol import JsmolView


class JsmolWidget(JsmolView):
    """Jupyter widget for JSMol molecular geometry visualization
    """
    def __init__(self, *args, **kwargs):
        self.default_info['color'] = 'white'
        super().__init__(*args, **kwargs)

    @classmethod
    def from_mol(cls, mol: "Mol", num_confs: int=10, *args, **kwargs):
        """Generate JsmolWidget object from RDKit molecule

        :param mol: RDKit molecule to visualize
        :type mol: Mol
        :param num_confs: Number of conformers to generate, defaults to 10
        :type num_confs: int, optional
        :return: Widget instance
        :rtype: JsmolWidget
        """
        xyz = mol_to_xyz(mol, num_confs=num_confs)
        return cls.from_str(xyz)
