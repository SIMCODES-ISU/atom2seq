from atom2seq.classes import Atom, Mol
import pytest


@pytest.fixture
def hydrogen():
    return ("H", (0, 0, 0))


def test_symbol(hydrogen):
    test_atom = Atom(*hydrogen)
    assert test_atom.symbol == "H"


def test_coords(hydrogen):
    test_atom = Atom(*hydrogen)
    assert test_atom.coords == (0, 0, 0)


def test_atom_repr(hydrogen):
    test_atom = Atom(*hydrogen)
    assert repr(test_atom) == "Atom('H', (0, 0, 0))"


@pytest.fixture
def water():
    test_atoms = [
        Atom("H", [0, 1, 0]),
        Atom("H", [1, 0, 0]),
        Atom("O", [0, 0, 0]),
    ]
    test_bonds = [[0, 2], [1, 2]]
    return (test_atoms, test_bonds)


def test_mol_eq(water):
    test_mol = Mol(*water)
    test_mol2 = Mol(*water)
    assert test_mol == test_mol2


def test_get_bonds(water):
    test_mol = Mol(*water)
    assert test_mol.get_bonds() == [[0, 2], [1, 2]]


def test_get_atoms(water):
    test_mol = Mol(*water)
    assert test_mol.get_atoms() == [
        Atom("H", [0, 1, 0]),
        Atom("H", [1, 0, 0]),
        Atom("O", [0, 0, 0]),
    ]


def test_mol_repr(water):
    test_mol = Mol(*water)
    assert (
        repr(test_mol) == "Mol([Atom('H', [0, 1, 0]), Atom('H', [1, 0, 0]), "
        "Atom('O', [0, 0, 0])], [[0, 2], [1, 2]])"
    )


def test_not_is_bond(water):
    test_mol = Mol(*water)
    assert not test_mol.is_bond(0, 1)


def test_is_bond(water):
    test_mol = Mol(*water)
    assert test_mol.is_bond(2, 0)


def test_del_bond(water):
    test_mol = Mol(*water)
    test_mol.del_bond(0, 2)
    assert not test_mol.is_bond(0, 2)


def test_add_bond(water):
    test_mol = Mol(*water)
    test_mol.add_bond(0, 1)
    assert test_mol.is_bond(0, 1)


def test_dist_1(water):
    test_mol = Mol(*water)
    assert test_mol.dist(0, 2) == test_mol.dist(1, 2) == 1


def test_dist_sqrt2(water):
    test_mol = Mol(*water)
    assert test_mol.dist(0, 1) == pytest.approx(2**0.5)


def test_del_atom(water):
    test_mol = Mol(*water)
    test_mol.del_atom(0)
    assert (
        test_mol.get_atoms()
        == [
            Atom("H", [1, 0, 0]),
            Atom("O", [0, 0, 0]),
        ]
    ) and (test_mol.get_bonds() == [[0, 1]])


def test_get_bonded(water):
    test_mol = Mol(*water)
    assert test_mol.get_bonded(2) == [0, 1]


def test_get_n_term(water):
    test_mol = Mol(*water)
    assert not test_mol.get_n_term()


def test_set_n_term(water):
    test_mol = Mol(*water)
    test_mol.set_n_term(0)
    assert test_mol.get_n_term() == 0


def test_get_backbone(water):
    test_mol = Mol(*water)
    assert test_mol.get_backbone() == [False, False, False]


def test_set_backbone(water):
    test_mol = Mol(*water)
    test_mol.set_backbone([0, 2])
    assert test_mol.get_backbone() == [True, False, True]


def test_get_aas(water):
    test_mol = Mol(*water)
    assert test_mol.get_aas() == [False, False, False]


def test_number_aas(water):
    test_mol = Mol(*water)
    test_mol.number_aas([0, 2], [0, 1])
    assert test_mol.get_aas() == [0, False, 1]


def test_find_submol_1(water):
    test_mol = Mol(*water)
    assert test_mol.find_submol(0) == [0, 1, 2]


@pytest.fixture
def water_x2():
    test_atoms = [
        Atom("H", [0, 1, 0]),
        Atom("H", [1, 0, 0]),
        Atom("O", [0, 0, 0]),
        Atom("H", [0, 3, 0]),
        Atom("H", [1, 2, 0]),
        Atom("O", [0, 2, 0]),
    ]
    test_bonds = [[0, 2], [1, 2], [3, 5], [4, 5]]
    return (test_atoms, test_bonds)


def test_find_submol_2(water_x2):
    test_mol = Mol(*water_x2)
    assert test_mol.find_submol(0) == [0, 1, 2]


def test_del_submol(water, water_x2):
    test_mol = Mol(*water)
    test_mol_x2 = Mol(*water_x2)
    test_mol_x2.del_submol(3)
    assert test_mol_x2 == test_mol


def test_split_submol(water, water_x2):
    test_mol = Mol(*water)
    test_mol_x2 = Mol(*water_x2)
    new_mol = test_mol_x2.split_submol(0)
    assert new_mol == test_mol
