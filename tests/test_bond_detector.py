import pytest

from atom2seq.bond_detector import bond_mol, radii
from atom2seq.classes import Atom, Mol
from atom2seq.parsers import parse_xyz

prefix = __file__.removesuffix("test_bond_detector.py")


@pytest.fixture
def pairs():
    return [
        ["H", "O"],
        ["H", "N"],
        ["H", "C"],
        ["C", "N"],
        ["C", "O"],
        ["S", "C"],
        ["S", "H"],
    ]


def out_of_range(sym1, sym2):
    atom1 = Atom(sym1, (0, 0, 0))
    atom2 = Atom(sym2, (1.2 * (radii[sym1] + radii[sym2]), 0, 0))
    return Mol([atom1, atom2], [])


def at_range(sym1, sym2):
    atom1 = Atom(sym1, (0, 0, 0))
    atom2 = Atom(sym2, (1.1 * (radii[sym1] + radii[sym2]), 0, 0))
    return Mol([atom1, atom2], [])


def in_range(sym1, sym2):
    atom1 = Atom(sym1, (0, 0, 0))
    atom2 = Atom(sym2, (radii[sym1] + radii[sym2], 0, 0))
    return Mol([atom1, atom2], [])


def under_range(sym1, sym2):
    atom1 = Atom(sym1, (0, 0, 0))
    atom2 = Atom(sym2, (0.4, 0, 0))
    return Mol([atom1, atom2], [])


def test_out_of_range(pairs):
    for pair in pairs:
        mol = out_of_range(*pair)
        bond_mol(mol)
        assert mol.get_bonds() == []


def test_at_range(pairs):
    for pair in pairs:
        mol = at_range(*pair)
        bond_mol(mol)
        assert mol.get_bonds() == [[0, 1]]


def test_in_range(pairs):
    for pair in pairs:
        mol = in_range(*pair)
        bond_mol(mol)
        assert mol.get_bonds() == [[0, 1]]


def test_under_range(pairs):
    for pair in pairs:
        mol = under_range(*pair)
        bond_mol(mol)
        assert mol.get_bonds() == []


@pytest.fixture
def glycine():
    return parse_xyz(prefix + "assets/bond_detector_tests/glycine.xyz")


def test_glycine(glycine):
    bond_mol(glycine)
    glycine_bonds = [
        [5, 6],
        [4, 5],
        [4, 7],
        [0, 4],
        [0, 8],
        [0, 9],
        [1, 2],
        [1, 3],
        [0, 1],
    ]
    found_non_bond = False
    for bond in glycine.get_bonds():
        if bond not in glycine_bonds:
            found_non_bond = True
        else:
            glycine_bonds.remove(bond)
    print(glycine.get_bonds())
    assert (not found_non_bond) and (len(glycine_bonds) == 0)
