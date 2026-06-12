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
        self._atoms = atoms
        self._bonds = bonds
        self._backbone = [False for _ in self._atoms]
        self._aas = [False for _ in self._atoms]
        self._n_term = False

    def __repr__(self):
        return f"Mol({self._atoms}, {self._bonds})"

    def __eq__(self, other):
        if (len(self._atoms) != len(other.get_atoms())) or (
            len(self._bonds) != len(other.get_bonds())
        ):
            return False
        else:
            for i in range(len(self._atoms)):
                if self._atoms[i] != other.get_atoms()[i]:
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
        # Loops over every bond and checks whether that bond is between the two
        # given atoms
        i = 0
        while i < len(self._bonds):
            bond = self._bonds[i]
            if (bond[0] == n) and (bond[1] == m):
                return True
            elif (bond[1] == n) and (bond[0] == m):
                return True
            i += 1
        return False

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
        if n == m:
            raise ValueError("You cannot bond an atom with itself.")
        if self.is_bond(n, m):
            # Raises an AttributeError if there is already a bond between the
            # given atoms.
            raise AttributeError("You cannot bond two bonded atoms.")
        else:
            self._bonds.append([min(n, m), max(n, m)])

    def dist(self, n, m):
        """Calculates the Euclidean distance bewteen the given atoms.

        Parameters:
            n (int): One of the indices to check.
            m (int): The other index to check.
        """
        n_coords = self._atoms[n].coords
        m_coords = self._atoms[m].coords
        return math.sqrt(
            (n_coords[0] - m_coords[0]) ** 2
            + (n_coords[1] - m_coords[1]) ** 2
            + (n_coords[2] - m_coords[2]) ** 2
        )

    def get_bonds(self):
        return self._bonds

    def get_atoms(self):
        return self._atoms

    def del_atom(self, idx):
        self._atoms.pop(idx)
        self._aas.pop(idx)
        self._backbone.pop(idx)
        new_bonds = []
        for i in range(len(self._bonds)):
            bond = self._bonds[i]
            if (bond[0] != idx) and (bond[1] != idx):
                new_bonds.append(bond)
        for bond in new_bonds:
            if bond[0] > idx:
                bond[0] -= 1
            if bond[1] > idx:
                bond[1] -= 1
        self._bonds = new_bonds

    def get_bonded(self, idx):
        out = []
        for bond in self._bonds:
            if bond[0] == idx:
                out.append(bond[1])
            elif bond[1] == idx:
                out.append(bond[0])
        return out

    def find_submol(self, idx):
        current_idcs = {idx}
        out = {idx}
        done = False
        while not done:
            for i in current_idcs:
                for elt in self.get_bonded(i):
                    out.add(elt)
                if out == current_idcs:
                    done = True
            for elt in self.get_bonded(i):
                current_idcs.add(elt)
        return list(out)

    def del_submol(self, idx):
        submol_idcs = self.find_submol(idx)
        submol_idcs = sorted(submol_idcs, reverse=True)
        for idx in submol_idcs:
            self.del_atom(idx)

    def split_submol(self, idx):
        submol_idcs = self.find_submol(idx)
        atoms = []
        bonds = []
        for idx in submol_idcs:
            atoms.append(self._atoms[idx])
        for bond in self._bonds:
            if bond[0] in submol_idcs:
                new_bond = [
                    submol_idcs.index(bond[0]),
                    submol_idcs.index(bond[1]),
                ]  # noqa
                bonds.append(new_bond)
        return Mol(atoms, bonds)

    def set_n_term(self, idx):
        self._n_term = idx

    def get_n_term(self):
        return self._n_term

    def set_backbone(self, idx_list):
        for idx in idx_list:
            self._backbone[idx] = True

    def number_aas(self, idx_list, num_list):
        for i in range(len(idx_list)):
            idx, num = idx_list[i], num_list[i]
            self._aas[idx] = num

    def get_backbone(self):
        return self._backbone

    def get_aas(self):
        return self._aas
