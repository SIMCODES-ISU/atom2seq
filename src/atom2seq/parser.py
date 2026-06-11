from atom2seq.classes import Atom, Mol


def parse(filename):
    """Parses the xyz-coordinates of atoms in a molecule stored in PDB format
    and returns a Mol object that has those atoms in it.

    Parameters:
        filename (str): The path to the file to parse.

    Returns:
        Mol: A molecule obect containing the atoms in the initial file.
    """
    # Opens the file and reads the lines of the file into a list, then closes
    # the file.
    file = open(filename, "r")
    raw_content = file.readlines()
    file.close()

    # Gets rid of newlines within the content.
    content = []
    for line in raw_content:
        content.append(line.strip())

    # Checking if the first line is the number of atoms. If it is, remove
    # the first line and any blank lines that come after it.
    if content[0][0].isdigit():
        content.pop(0)
        # This method for removing blank strings was taken from
        # geeksforgeeks.org/python/python-remove-empty-strings-from-list-of-strings/
        content = [line for line in content if line]

    # Initialize a list of atoms.
    atoms = []
    for line in content:
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
