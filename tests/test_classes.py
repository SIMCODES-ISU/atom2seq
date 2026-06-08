import pytest

from atom2seq import classes

test_atom = classes.Atom("H", [0, 0, 0])


def test_symbol():
    assert test_atom.symbol == "H"


def test_coords():
    assert test_atom.coords == [0, 0, 0]


def test_repr():
    assert repr(test_atom) == "Atom('H', [0, 0, 0])"


test_atoms = [
    classes.Atom("H", [0, 1, 0]),
    classes.Atom("H", [1, 0, 0]),
    classes.Atom("O", [0, 0, 0]),
]
test_bonds = [[0, 2], [1, 2]]
test_mol = classes.Mol(test_atoms, test_bonds)
test_mol2 = classes.Mol(test_atoms, test_bonds)


def test_eq():
    assert test_mol == test_mol2


def test_get_bonds():
    assert test_mol.get_bonds() == [[0, 2], [1, 2]]


def test_atoms():
    assert test_mol.atoms == [
        classes.Atom("H", [0, 1, 0]),
        classes.Atom("H", [1, 0, 0]),
        classes.Atom("O", [0, 0, 0]),
    ]


def test_mol_repr():
    assert (
        repr(test_mol) == "Mol([Atom('H', [0, 1, 0]), Atom('H', [1, 0, 0]), "
        "Atom('O', [0, 0, 0])], [[0, 2], [1, 2]])"
    )


def test_not_is_bond():
    assert not test_mol.is_bond(0, 1)


def test_is_bond():
    assert test_mol.is_bond(2, 0)


def test_del_bond():
    test_mol.del_bond(0, 2)
    assert not test_mol.is_bond(0, 2)


def test_add_bond():
    test_mol.add_bond(0, 2)
    assert test_mol.is_bond(0, 2)


def test_dist_1():
    assert test_mol.dist(0, 2) == test_mol.dist(1, 2) == 1


def test_dist_sqrt_2():
    assert test_mol.dist(0, 1) == pytest.approx(2**0.5)
