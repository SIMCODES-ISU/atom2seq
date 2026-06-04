import math


class Atom:
    """A class representing an atom.

    Attributes:
        symbol (str): The atomic symbol of the atom.
        coords (tuple[float]): The XYZ coordinates of the atom in Å.

    Supported Special Methods:
        eq (==)
    """

    def __init__(self, symbol, coords):
        self.symbol = symbol
        self.coords = coords

    def __repr__(self):
        return f"Atom('{self.symbol}', {self.coords})"

    def __eq__(self, other):
        # Returns True if and only if both the symbol and the coordinates are
        # the same
        return (self.symbol == other.symbol) and (self.coords == other.coords)


class Mol:
    """A class representing a molecule.

    Attributes:
        atoms (list[Atom]): A list of atoms in the molecule.
        bonds (list[list]): A list of the bonds between atoms in the molecule.
            A bond is stored as [n, m, type] where n and m are the indices of
            the atoms the bond is between in atoms and type is 's' for single,
            'd' for double, or 'h' for hydrogen.

    Methods:
        is_bond(n: int, m: int) -> str: Returns 'n' if there is no bond between
            indices n and m, and returns the type of bond if there is.
        add_bond(n: int, m: int, type: str) -> None: Adds a bond of type type
            between indices n and m. Throws an AttributeError if there is
            already a bond there.
        del_bond(n: int, m: int) -> None: Deletes the bond between indices n
            and m. Does nothing if there is no bond.
        dist(n: int, m: int) -> float: Calculates the Euclidean distance
            between the atoms at indices n and m.

    Supported Special Methods:
        eq (==)
    """

    def __init__(self, atoms, bonds):
        self.atoms = atoms
        self._bonds = bonds

    def __repr__(self):
        return f"Mol({self.atoms}, {self._bonds})"

    def __eq__(self, other):
        if (len(self.atoms) != len(other.atoms)) or (
            len(self._bonds) != len(other.get_bonds())
        ):
            return False
        else:
            for i in range(len(self.atoms)):
                if self.atoms[i] != other.atoms[i]:
                    return False
            for j in range(len(self._bonds)):
                if self._bonds[j] != other.get_bonds()[j]:
                    return False
        return True

    def is_bond(self, n, m):
        """Detects whether the atoms at indices n and m are bonded.

        Parameters:
            n (int): One of the indices to check.
            m (int): The other index to check.

        Returns:
            str: One of 's', 'd', 'h', 'n'. If 'n', the atoms passed in are not
                bonded. Otherwise, they are bonded and the letter indicates the
                type of bond: 's' is single, 'd' is double, and 'h' is
                hydrogen.
        """
        # Sets up two tracking variables
        i = 0
        found = False
        # Sets the output to 'n', then enters a loop while it hasn't found a
        # bond between these atoms and it has't reached the last bond in the
        # list.
        out = "n"
        while (not found) and (i < len(self._bonds)):
            bond = self._bonds[i]
            if (bond[0] == n) and (bond[1] == m):
                out = bond[2]
                found = True
            elif (bond[1] == n) and (bond[0] == m):
                out = bond[2]
                found = True
            i += 1
        return out

    def del_bond(self, n, m):
        """Deletes a bond between the two given atoms.

        Parameters:
            n (int): One of the indices to check.
            m (int): The other index to check.
        """
        # Goes over each bond and checks if it is between the two given atoms.
        # Removes it if it is.
        for bond in self._bonds:
            if (bond[0] == n) and (bond[1] == m):
                self._bonds.remove(bond)
            elif (bond[1] == n) and (bond[0] == m):
                self._bonds.remove(bond)

    def add_bond(self, n, m, type):
        """Adds a bond of the given type between the two given atoms.

        Parameters:
            n (int): One of the indices to bond.
            m (int): The other index to bond.
            type (str): The type of bond that is between the two given atoms.
                Must be one of 's', 'd', and 'h'.
        """
        if type not in ["s", "d", "h"]:
            # Raises a ValueError if the bond type is invalid.
            raise ValueError(f"{type} is not a valid bond type.")
        elif self.is_bond(n, m) != "n":
            # Raises an AttributeError if there is already a bond between the
            # given atoms.
            raise AttributeError("You cannot add a bond between bonded atoms.")
        else:
            # If neither error is raised,
            self._bonds.append([min(n, m), max(n, m), type])

    def dist(self, n, m):
        """Calculates the Euclidean distance bewteen the given atoms.

        Parameters:
            n (int): One of the indices to check.
            m (int): The other index to check.
        """
        n_coords = self.atoms[n].coords
        m_coords = self.atoms[m].coords
        return math.sqrt(
            (n_coords[0] - m_coords[0]) ** 2
            + (n_coords[1] - m_coords[1]) ** 2
            + (n_coords[2] - m_coords[2]) ** 2
        )

    def get_bonds(self):
        return self._bonds
