from atom2seq.classes import Mol


def verify_mol(molecule: Mol):
    if len(molecule.get_atoms()) < 10:
        return False
    num_elements = {}
    req_atoms = ["H", "C", "O", "N"]
    valid_atoms = ["H", "C", "O", "N", "S"]
    for atom in molecule.get_atoms():
        if atom.symbol not in valid_atoms:
            return False
        if atom.symbol in num_elements.keys():
            num_elements[atom.symbol] += 1
        else:
            num_elements[atom.symbol] = 1
    for elt in req_atoms:
        if elt not in num_elements.keys():
            return False
    return True
