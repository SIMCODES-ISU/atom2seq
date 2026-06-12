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


def parse_xyz(filename):
    """Parses the xyz-coordinates of atoms in a molecule stored in XYZ format
    and returns a Mol object that has those atoms in it.

    Parameters:
        filename (str): The path to the file to parse.

    Returns:
        Mol: A molecule obect containing the atoms in the initial file.
    """
    # Opens the file and reads the lines of the file into a list, then closes
    # the file.
    file = open(filename, "r")
    contents = file.readlines()
    file.close()

    # Gets rid of newlines within the content.
    contents = [line.strip() for line in contents]

    # Checking if the first line is the number of atoms. If it is, remove
    # the first line and any blank lines that come after it.
    if contents[0][0].isdigit():
        contents.pop(0)
        # This method for removing blank strings was taken from
        # geeksforgeeks.org/python/python-remove-empty-strings-from-list-of-strings/
        contents = [line for line in contents if line]

    # Initialize a list of atoms.
    atoms = []
    for line in contents:
        # Split into [symbol, x, y, z], then pop the symbol.
        coords = line.split()
        symbol = coords.pop(0)
        # Converts the coordinates into floats.
        for i in range(len(coords)):
            coords[i] = float(coords[i])
        # Appends an atom containing the symbol and coordinates to the
        # list.
        atoms.append(Atom(symbol, tuple(coords)))

    return Mol(atoms, [])


def parse_pdb(filename):
    file = open(filename, "r")
    contents = file.readlines()
    file.close()
    contents = [line.strip() for line in contents if line[0:4] == "ATOM"]
    contents = [[line.split()[-1], *line.split()[-6:-3]] for line in contents]
    contents = [
        [int(elt) if elt.isdigit() else elt for elt in listy]
        for listy in contents  # noqa
    ]
    atoms = [Atom(line[0], tuple(line[1:])) for line in contents]
    return Mol(atoms, [])


def parse_cif(filename):
    file = open(filename, "r")
    contents = file.readlines()
    file.close()
    print(contents)
    contents = [line.strip() for line in contents if line[0:4] == "ATOM"]
    print(contents)
    contents = [[line.split()[2], *line.split()[8:11]] for line in contents]
    print(contents)
    contents = [
        [int(elt) if elt.isdigit() else elt for elt in listy]
        for listy in contents  # noqa
    ]
    print(contents)
    atoms = [Atom(line[0], tuple(line[1:])) for line in contents]
    return Mol(atoms, [])
