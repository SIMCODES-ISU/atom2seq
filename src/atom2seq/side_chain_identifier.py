from atom2seq.classes import Atom, Mol

aa_dict = {
    [1, 3, 0, 0, 0]: "A",
    [1, 3, 0, 0, 1]: "C",
    [2, 3, 0, 2, 0]: "D",
    [3, 5, 0, 2, 0]: "E",
    [7, 7, 0, 0, 0]: "F",
    [0, 1, 0, 0, 0]: "G",
    [4, 6, 2, 0, 0]: "H",
    [4, 9, 0, 0, 0]: "I/L",
    [4, 10, 1, 0, 0]: "K",
    [3, 7, 0, 0, 1]: "M",
    [2, 4, 1, 1, 0]: "N",
    [3, 5, 0, 0, 0]: "P",
    [3, 6, 1, 1, 0]: "Q",
    [4, 12, 3, 0, 0]: "R",
    [1, 3, 0, 1, 0]: "S",
    [2, 5, 0, 1, 0]: "T",
    [3, 7, 0, 0, 0]: "V",
    [9, 8, 1, 0, 0]: "W",
    [7, 7, 0, 1, 0]: "Y",
}


def id_side_chains(molecule: Mol) -> list[str]:
    old_idx = False
    current_idx = molecule.get_n_term()
    new_current_idx = False
    out = []
    while len(molecule.get_atoms()) > 1:
        # Finds the next step in the backbone (it will be the alpha carbon) and
        # deletes any hydrogens connected to this nitrogen.
        to_del = []
        for i in molecule.get_bonded(current_idx):
            if molecule.get_backbone()[i]:
                old_idx = current_idx
                new_current_idx = i
            if molecule.get_atoms()[i].symbol == "H":
                to_del.append(i)
        for i in to_del:
            molecule.del_atom(i)
        current_idx = new_current_idx
        # Deletes the atom we were just at.
        molecule.del_atom(old_idx)
        # Finds how many hydrogens are bonded to the alpha carbon. If it is 2,
        # this is a glycine and the R-group will be deleted, so we set the
        # R-group prematurely.
        numH = [
            molecule.get_atoms()[i].symbol
            for i in molecule.get_bonded(current_idx)  # noqa
        ].count("H")
        r_group = False
        r_group_idx = False
        if numH == 2:
            to_del = []
            for i in molecule.get_bonded(current_idx):
                if molecule.get_atoms()[i].symbol == "H":
                    to_del.append(i)
                elif molecule.get_backbone()[i]:
                    old_idx = current_idx
                    new_current_idx = i
            for i in to_del:
                molecule.del_atom(i)
            r_group = Mol(Atom("H", (0, 0, 0)), [])
        else:
            # Deletes all hydrogens bonded to the alpha carbon. Now, the only
            # two things bonded to the alpha carbon are a carbon with 3 bonds
            # (in the carboxylic acid) and the start of the R-group.
            to_del = []
            for i in molecule.get_bonded(current_idx):
                if molecule.get_atoms()[i].symbol == "H":
                    to_del.append(i)
                elif len(molecule.get_bonded(i)) != 3:
                    r_group_idx = i
                elif molecule.get_backbone()[i]:
                    old_idx = current_idx
                    new_current_idx = i
        for i in to_del:
            molecule.del_atom(i)
        current_idx = new_current_idx
        # Deletes the alpha carbon. Now, the r-group is seperated from the rest
        # of the protein, so we can extract it with split_submol().
        molecule.del_atom(old_idx)
        if r_group_idx:
            r_group = molecule.split_submol(r_group_idx)
        r_group_makeup = [0, 0, 0, 0, 0]
        for atom in r_group.get_atoms():
            if atom.symbol == "C":
                r_group_makeup[0] += 1
            elif atom.symbol == "H":
                r_group_makeup[1] += 1
            elif atom.symbol == "N":
                r_group_makeup[2] += 1
            elif atom.symbol == "O":
                r_group_makeup[3] += 1
            elif atom.symbol == "S":
                r_group_makeup[4] += 1
        symbol = aa_dict[r_group_makeup]
        real_symbol = symbol
        if symbol == "I/L":
            for i in range(len(r_group)):
                atom = r_group.get_atoms()[i]
                # Checks how many carbons each carbon is bonded to. If one of
                # them has 3, then it is isoleucine.
                if atom.symbol == "C":
                    bonded_carbons = [
                        j
                        for j in r_group.get_bonded(i)
                        if r_group.get_atoms()[j].symbol == "C"
                    ]
                    if len(bonded_carbons) == 3:
                        real_symbol = "I"
            if real_symbol == "I/L":
                real_symbol = "L"
        out.append(real_symbol)
        for i in molecule.get_bonded(current_idx):
            if molecule.get_backbone()[i]:
                old_idx = current_idx
                new_current_idx = i
        current_idx = new_current_idx
        for i in molecule.get_bonded(old_idx):
            if molecule.get_atoms()[i].symbol == "O":
                to_del = i
        molecule.del_atom(to_del)
        molecule.del_atom(old_idx)
    return out
