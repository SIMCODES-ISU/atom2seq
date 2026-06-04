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
            Boolean: Whether the atoms are bonded.
        """
        # Sets the output to False, then enters a loop while it hasn't found a
        # bond between these atoms and it has't reached the last bond in the
        # list.
        i = 0
        found = False
        while (not found) and (i < len(self._bonds)):
            bond = self._bonds[i]
            if (bond[0] == n) and (bond[1] == m):
                found = True
            elif (bond[1] == n) and (bond[0] == m):
                found = True
            i += 1
        return found

    def del_bond(self, n, m):
        """Deletes a bond between the two given atoms.

        Parameters:
            n (int): One of the indices to check.
            m (int): The other index to check.

        Returns:
            None
        """
        # Goes over each bond and checks if it is between the two given atoms.
        # Removes it if it is.
        for bond in self._bonds:
            if (bond[0] == n) and (bond[1] == m):
                self._bonds.remove(bond)
            elif (bond[1] == n) and (bond[0] == m):
                self._bonds.remove(bond)

    def add_bond(self, n, m):
        """Adds a bond between the two given atoms.

        Parameters:
            n (int): One of the indices to bond.
            m (int): The other index to bond.

        Returns:
            None
        """
        if self.is_bond(n, m):
            # Raises an AttributeError if there is already a bond between the
            # given atoms.
            raise AttributeError("You cannot add a bond between bonded atoms.")
        else:
            self._bonds.append([min(n, m), max(n, m)])

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
