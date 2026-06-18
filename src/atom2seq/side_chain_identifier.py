from atom2seq.classes import Atom, Mol

aa_dict = {
    (1, 3, 0, 0, 0): "A",
    (1, 3, 0, 0, 1): "C",
    (2, 3, 0, 2, 0): "D",
    (3, 5, 0, 2, 0): "E",
    (7, 7, 0, 0, 0): "F",
    (0, 1, 0, 0, 0): "G",
    (4, 6, 2, 0, 0): "H",
    (4, 9, 0, 0, 0): "I/L",
    (4, 10, 1, 0, 0): "K",
    (3, 7, 0, 0, 1): "M",
    (2, 4, 1, 1, 0): "N",
    (3, 5, 0, 0, 0): "P",
    (3, 6, 1, 1, 0): "Q",
    (4, 12, 3, 0, 0): "R",
    (1, 3, 0, 1, 0): "S",
    (2, 5, 0, 1, 0): "T",
    (3, 7, 0, 0, 0): "V",
    (9, 8, 1, 0, 0): "W",
    (7, 7, 0, 1, 0): "Y",
}


def update_idx_del(molecule: Mol, current_idx: int, del_sym: str) -> int:
    to_del = []
    old_idx = False
    new_current_idx = False
    for i in molecule.get_bonded(current_idx):
        if molecule.get_backbone()[i]:
            old_idx = current_idx
            new_current_idx = i
        if molecule.get_atoms()[i].symbol == del_sym:
            to_del.append(i)
    for i in sorted(to_del, reverse=True):
        molecule.del_atom(i)
    molecule.del_atom(old_idx)
    if new_current_idx > old_idx:
        current_idx = new_current_idx - 1
    else:
        current_idx = new_current_idx
    return current_idx


def update_idx_delH_find_Rgroup(molecule: Mol, current_idx: int) -> tuple[int]:
    to_del = []
    old_idx = False
    new_current_idx = False
    r_group_idx = False
    for i in molecule.get_bonded(current_idx):
        if molecule.get_backbone()[i]:
            old_idx = current_idx
            new_current_idx = i
        elif len(molecule.get_bonded(i)) != 3:
            r_group_idx = i
        elif molecule.get_atoms()[i].symbol == "H":
            to_del.append(i)
    for i in sorted(to_del, reverse=True):
        molecule.del_atom(i)
    molecule.del_atom(old_idx)
    if new_current_idx > old_idx:
        current_idx = new_current_idx - 1
    else:
        current_idx = new_current_idx
    return (current_idx, r_group_idx)


def ile_or_leu(r_group):
    out = "I/L"
    for i in range(len(r_group)):
        atom = r_group.get_atoms()[i]
        # Checks how many carbons each carbon is bonded to. If one of
        # them has 3, then it is isoleucine.
        if atom.symbol == "C":
            bonded_carbons = [
                j
                for j in r_group.get_bonded(i)
                if r_group.get_atoms()[j].symbol == "C"  # noqa
            ]
            if len(bonded_carbons) == 3:
                out = "I"
    if out == "I/L":
        out = "L"
    return out


def id_side_chains(molecule: Mol) -> list[str]:
    print(f"{molecule=}")
    current_idx = molecule.get_n_term()
    out = []
    num_iter = 0
    while len(molecule.get_atoms()) > 1:
        num_iter += 1
        print(f"{num_iter=}, {len(molecule.get_atoms())=}")
        # Finds the next step in the backbone (it will be the alpha carbon) and
        # deletes any hydrogens connected to this nitrogen.
        current_idx = update_idx_del(molecule, current_idx, "H")
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
            current_idx = update_idx_del(molecule, current_idx, "H")
            r_group = Mol([Atom("H", (0, 0, 0))], [])
        else:
            # Deletes all hydrogens bonded to the alpha carbon. Now, the only
            # two things bonded to the alpha carbon are a carbon with 3 bonds
            # (in the carboxylic acid) and the start of the R-group.
            current_idx, r_group_idx = update_idx_delH_find_Rgroup(
                molecule, current_idx
            )
        if (not r_group) and (r_group_idx):
            print(f"{r_group_idx=}, {molecule=}")
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
        symbol = aa_dict[tuple(r_group_makeup)]
        real_symbol = symbol
        if symbol == "I/L":
            real_symbol = ile_or_leu(r_group)
        out.append(real_symbol)
        current_idx = update_idx_del(molecule, current_idx, "O")
    return out
