from atom2seq.classes import Atom, Mol


def parse_bean(filename):
    file = open(filename, "r")
    contents = file.read()
    file.close()
    contents = contents.strip()
    contents = contents.replace(" ", "")
    contents = contents.split("\\")
    contents = [line.split(",") for line in contents]
    contents = [
        [int(elt) if elt.isdigit() else elt for elt in listy]
        for listy in contents  # noqa
    ]
    atoms = [Atom(listy[0], tuple(listy[1:])) for listy in contents]
    return Mol(atoms, [])
