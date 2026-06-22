import pytest

from atom2seq import backbone_identifier as bid
from atom2seq.parsers import parse_xyz

prefix = __file__.removesuffix("test_backbone_identifier.py")


@pytest.fixture
def glycine():
    gly = parse_xyz(prefix + "assets/backbone_identifier_tests/glycine.xyz")
    gly.add_bond(0, 1)
    gly.add_bond(1, 2)
    gly.add_bond(2, 3)
    gly.add_bond(2, 5)
    gly.add_bond(4, 5)
    gly.add_bond(5, 6)
    gly.add_bond(5, 7)
    gly.add_bond(7, 8)
    gly.add_bond(7, 9)
    return gly


def test_id_n_term_G(glycine):
    assert bid.id_n_term(glycine) == [7]


def test_id_pro_n_term_G(glycine):
    assert bid.id_pro_n_term(glycine) == []


def test_find_atom_G(glycine):
    assert bid.find_atom(glycine, [7], "C") == [5]


def test_small_backbone_iter_G(glycine):
    assert bid.small_backbone_iter(glycine, [[7]], "C") == [[7, 5]]


def test_large_backbone_iter_G(glycine):
    assert bid.large_backbone_iter(glycine, [[7]]) == (
        [[7, 5, 2, 1], [7, 5, 2, 3]],
        True,
    )


def test_find_backbone_G(glycine):
    assert bid.find_backbone(glycine) == [7, 5, 2, 1]


def test_label_backbone_G(glycine):
    bid.label_backbone(glycine)
    assert glycine.get_backbone() == [
        False,
        True,
        True,
        False,
        False,
        True,
        False,
        True,
        False,
        False,
    ]


@pytest.fixture
def gas():
    gas = parse_xyz(prefix + "assets/backbone_identifier_tests/gas.xyz")
    gas.add_bond(0, 1)
    gas.add_bond(1, 2)
    gas.add_bond(1, 4)
    gas.add_bond(3, 4)
    gas.add_bond(4, 5)
    gas.add_bond(4, 7)
    gas.add_bond(6, 7)
    gas.add_bond(7, 9)
    gas.add_bond(8, 9)
    gas.add_bond(9, 12)
    gas.add_bond(10, 13)
    gas.add_bond(11, 12)
    gas.add_bond(12, 13)
    gas.add_bond(12, 16)
    gas.add_bond(13, 14)
    gas.add_bond(13, 17)
    gas.add_bond(15, 16)
    gas.add_bond(16, 19)
    gas.add_bond(18, 19)
    gas.add_bond(19, 22)
    gas.add_bond(20, 23)
    gas.add_bond(21, 22)
    gas.add_bond(22, 23)
    gas.add_bond(22, 27)
    gas.add_bond(23, 24)
    gas.add_bond(23, 28)
    gas.add_bond(24, 25)
    gas.add_bond(26, 27)
    gas.add_bond(27, 29)
    gas.add_bond(29, 30)
    return gas


def test_id_n_term_GAS(gas):
    print(gas)
    assert bid.id_n_term(gas) == [1]


def test_find_atom_GAS(gas):
    assert bid.find_atom(gas, [1], "C") == [4]


def test_small_backbone_iter_GAS(gas):
    assert bid.small_backbone_iter(gas, [[1]], "C") == [[1, 4]]


def test_large_backbone_iter_GAS(gas):
    assert bid.large_backbone_iter(gas, [[1]]) == (
        [[1, 4, 7, 9]],
        False,
    )


def test_find_backbone_GAS(gas):
    assert bid.find_backbone(gas) == [1, 4, 7, 9, 12, 16, 19, 22, 27, 26]


def test_label_backbone_GAS(gas):
    bid.label_backbone(gas)
    assert gas.get_backbone() == [
        False,
        True,
        False,
        False,
        True,
        False,
        False,
        True,
        False,
        True,
        False,
        False,
        True,
        False,
        False,
        False,
        True,
        False,
        False,
        True,
        False,
        False,
        True,
        False,
        False,
        False,
        True,
        True,
        False,
        False,
        False,
    ]
