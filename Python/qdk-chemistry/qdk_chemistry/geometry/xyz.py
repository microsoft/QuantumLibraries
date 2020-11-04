from typing import List, Tuple

def coordinates_to_xyz(number_of_atoms: int, charge: int, coordinates: List[Tuple[str, float, float, float]]) -> str:
    """Convert coordinates to XYZ file formatted string.

    :param number_of_atoms: Number of atoms in the conformer
    :type number_of_atoms: int
    :param charge: Charge of the conformer
    :type charge: int
    :param coordinates: List of tuples with values element name, x, y and z coordinates
    :type coordinates: List[Tuple[str, float, float, float]]
    :return: XYZ file format
    :rtype: str
    """
    result = [
        f"{number_of_atoms}",
        "title"
    ]

    for element, x, y, z in coordinates:
        result.append(f"{element} {x} {y} {z} ")

    if charge != 0:
        result.extend([
            "$set",
            f"chrg {charge}",
            "$end"
        ])

    return "\n".join(result)
