# Copyright (c) Microsoft Corporation.
# Licensed under the MIT License.

"""Module for converting coordinates to XYZ-formatted data

The formatting of the .xyz file format is as follows:

    <number of atoms>
    comment line
    <element> <X> <Y> <Z>
    ...

Source: https://en.wikipedia.org/wiki/XYZ_file_format.
"""

from typing import Iterable, Tuple

def element_coords_to_xyz(name: str, x: float, y: float, z: float) -> str:
    """Convert element to XYZ formatted line

    :param name: Element name
    :type name: str
    :param x: X coordinate
    :type x: float
    :param y: Y coordinate
    :type y: float
    :param z: Z coordinate
    :type z: float
    :return: XYZ formatted string
    :rtype: str
    """
    return f"{name} {x} {y} {z}"


def coordinates_to_xyz(
        number_of_atoms: int, 
        charge: int, 
        coordinates: Iterable[Tuple[str, float, float, float]],
        title: str = "unnamed"
    ) -> str:
    """Convert coordinates to XYZ file formatted string.

    :param number_of_atoms: Number of atoms in the conformer
    :type number_of_atoms: int
    :param charge: Charge of the conformer
    :type charge: int
    :param coordinates: List of tuples with values element name, x, y and z coordinates
    :type coordinates: List[Tuple[str, float, float, float]]
    :param title: XYZ file title
    :type title: str
    :return: XYZ file format
    :rtype: str
    """
    result = [
        f"{number_of_atoms}",
        title
    ]

    for element in coordinates:
        # Convert to <element> <x> <y> <z> and add space to end
        result.append(f"{element_coords_to_xyz(*element)} ")

    if charge != 0:
        result.extend([
            "$set",
            f"chrg {charge}",
            "$end"
        ])

    return "\n".join(result)
