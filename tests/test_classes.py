from atom2seq import classes

test_atom = classes.Atom("H", [0, 0, 0])

assert test_atom.symbol == "H"
assert test_atom.coords == [0, 0, 0]
assert repr(test_atom) == "Atom('H', [0, 0, 0])"

test_atoms = [
    classes.Atom("H", [0, 1, 0]),
    classes.Atom("H", [1, 0, 0]),
    classes.Atom("O", [0, 0, 0]),
]
test_bonds = [[0, 2, "s"], [1, 2, "s"]]
test_mol = classes.Mol(test_atoms, test_bonds)
test_mol2 = classes.Mol(test_atoms, test_bonds)

assert test_mol == test_mol2

assert test_mol._bonds == [[0, 2, "s"], [1, 2, "s"]]
assert test_mol.atoms == [
    classes.Atom("H", [0, 1, 0]),
    classes.Atom("H", [1, 0, 0]),
    classes.Atom("O", [0, 0, 0]),
]

assert (
    repr(test_mol) == "Mol([Atom('H', [0, 1, 0]), Atom('H', [1, 0, 0]), "
    "Atom('O', [0, 0, 0])], [[0, 2, 's'], [1, 2, 's']])"
)

assert test_mol.is_bond(0, 1) == "n"
assert test_mol.is_bond(2, 0) == "s"

test_mol.del_bond(0, 2)
assert test_mol.is_bond(0, 2) == "n"

test_mol.add_bond(0, 2, "d")
assert test_mol.is_bond(0, 2) == "d"

assert test_mol.dist(0, 2) == test_mol.dist(1, 2) == 1
assert test_mol.dist(0, 1) == 2**0.5
