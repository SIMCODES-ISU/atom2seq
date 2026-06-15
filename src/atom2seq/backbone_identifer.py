from atom2seq.classes import Mol


def id_n_term(molecule: Mol) -> list[int]:
    out = []
    # loops over the atoms in the molecule.
    for i in len(molecule.get_atoms):
        atom = molecule.get_atoms[i]
        valid = False
        # only considers the atom if it is a nitrogen.
        if atom.symbol == "N":
            # makes list of the symbols of atoms bonded to nitrogen.
            bonds = molecule.get_bonded(i)
            bonded_syms = [molecule.get_atoms[bond].symbol for bond in bonds]
            # initializes counters and helper variables.
            found_invalid = False
            h_counter = 0
            c_counter = 0
            for sym in bonded_syms:
                # Each hydrogen adds one to h_counter, carbon to c_counter, and
                # any other atom triggers the found_invalid flag.
                if sym == "H":
                    h_counter += 1
                elif sym == "C":
                    c_counter += 1
                else:
                    found_invalid = True
            # if the found_invalid flag was not triggered, then we check if it
            # is in either a) Proline or b) any other AA.
            if not found_invalid:
                if (h_counter == 1) and (c_counter == 2):
                    valid = True
                elif c_counter == 1:
                    if h_counter == 2 or h_counter == 3:
                        valid = True
            if valid:
                out.append(i)
    return out


def find_atom(molecule: Mol, listy: list, sym: str) -> list[int]:
    idx = listy[-1]
    bonds = molecule.get_bonded(idx)
    out = []
    for i in bonds:
        if molecule.get_atoms()[i].symbol == sym:
            out.append(i)
    return out


def small_backbone_iter(
    molecule: Mol, listy: list, step: str
) -> list[list[int]]:  # noqa
    new_list = []
    for sublist in listy:
        next_atoms = find_atom(molecule, listy, step)
        for atom in next_atoms:
            new_list.append(sublist.append(atom))
    return new_list


def large_backbone_iter(
    molecule: Mol, listy: list
) -> tuple[list[list[int]], bool]:  # noqa
    plus_c = small_backbone_iter(molecule, listy, "C")
    if len(plus_c) != 0:
        plus_c = small_backbone_iter(molecule, plus_c, "C")
        if len(plus_c) != 0:
            plus_o = small_backbone_iter(molecule, plus_c, "O")
            if len(plus_o) != 0:
                plus_h = find_atom(molecule, plus_o, "H")
                if (len(plus_h) == 1) or (len(plus_h) == 0):
                    return plus_o, True
            plus_n = small_backbone_iter(molecule, plus_c, "N")
            if len(plus_n) != 0:
                return plus_n, False
            else:
                raise Exception(
                    "This molecule is not a protein as it has no backbone."
                )  # noqa
        else:
            raise Exception(
                "This molecule is not a protein as it has no backbone."
            )  # noqa
    else:
        raise Exception("This molecule is not a protein as it has no backbone.")  # noqa


def find_backbone(molecule: Mol) -> list[int]:
    n_termini = id_n_term(molecule)
    paths = [[n_terminus] for n_terminus in n_termini]
    found = False
    i = 0
    while not found:
        paths, found = large_backbone_iter(molecule, paths)
        i += 1
        if i > 2 * len(molecule.atoms):
            raise RecursionError("While loop maximum depth exceeded.")
    return paths[0]


def label_backbone(molecule: Mol) -> None:
    backbone_idcs = find_backbone(molecule)
    molecule.set_backbone(backbone_idcs)
    molecule.set_n_term(backbone_idcs[0])
