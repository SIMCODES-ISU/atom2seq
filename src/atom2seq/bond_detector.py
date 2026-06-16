import time

from atom2seq.classes import Mol

radii = {"H": 0.31, "O": 0.66, "N": 0.71, "C": 0.76, "S": 1.05}

max_bonds = {"H": 1, "O": 2, "N": 3, "C": 4, "S": 2}


def bond_mol(molecule: Mol) -> None:
    i = 0
    begin = time.time()
    bond_list = [0 for atom in molecule.get_atoms()]
    for idx1 in range(len(molecule.get_atoms())):
        for idx2 in range(len(molecule.get_atoms())):
            i += 1
            if (idx1 != idx2) and (not molecule.is_bond(idx1, idx2)):
                sym1 = molecule.get_atoms()[idx1].symbol
                sym2 = molecule.get_atoms()[idx2].symbol
                radii_sum = radii[sym1] + radii[sym2]
                coords1 = molecule.get_atoms()[idx1].coords
                coords2 = molecule.get_atoms()[idx2].coords
                if (bond_list[idx1] < max_bonds[sym1]) and (
                    bond_list[idx2] < max_bonds[sym2]
                ):
                    if (
                        (abs(coords1[0] - coords2[0]) <= 1.1 * radii_sum)
                        and (abs(coords1[1] - coords2[1]) <= 1.1 * radii_sum)
                        and (abs(coords1[2] - coords2[2]) <= 1.1 * radii_sum)
                    ):
                        if (
                            0.25
                            <= molecule.squared_dist(idx1, idx2)
                            <= 1.21 * radii_sum
                        ):
                            molecule.add_bond(idx1, idx2)
                            bond_list[idx1] += 1
                            bond_list[idx2] += 1
                elif bond_list[idx1] >= max_bonds[sym1]:
                    # print(f"Atom {idx1} already has {bond_list[idx1]} bonds")
                    i += len(molecule.get_atoms()) - idx2
                    break
            # if i % 1000000 == 0:
            #     num_mil = int(i / 1000000)
            #     end = time.time()
            #     print(
            #         f"{num_mil} million pairs checked in {end - begin}"
            #         f"seconds. At this rate, the estimated total time to completion is {360000 * (end - begin)}"  # noqa
            #     )
            #     begin = time.time()
        continue
