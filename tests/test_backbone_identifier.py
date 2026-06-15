import pytest

from atom2seq import backbone_identifer as bid
from atom2seq.classes import Atom, Mol
from atom2seq.parsers import parse_xyz

prefix = __file__.removesuffix("test_backbone_identifier.py")


@pytest.fixture
def glycine():
    gly = parse_xyz(prefix + "assets/glycine.xyz")
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
def gay():
    glyAlaTyr = parse_xyz(prefix + "assets/gay.xyz")
    glyAlaTyr.add_bond(0, 1)
    glyAlaTyr.add_bond(1, 2)
    glyAlaTyr.add_bond(1, 4)
    glyAlaTyr.add_bond(3, 4)
    glyAlaTyr.add_bond(4, 5)
    glyAlaTyr.add_bond(4, 6)
    glyAlaTyr.add_bond(6, 7)
    glyAlaTyr.add_bond(6, 8)
    glyAlaTyr.add_bond(8, 11)
    glyAlaTyr.add_bond(9, 12)
    glyAlaTyr.add_bond(10, 11)
    glyAlaTyr.add_bond(11, 12)
    glyAlaTyr.add_bond(11, 15)
    glyAlaTyr.add_bond(12, 13)
    glyAlaTyr.add_bond(12, 14)
    glyAlaTyr.add_bond(15, 16)
    glyAlaTyr.add_bond(15, 19)
    glyAlaTyr.add_bond(17, 21)
    glyAlaTyr.add_bond(18, 22)
    glyAlaTyr.add_bond(19, 24)
    glyAlaTyr.add_bond(20, 25)
    glyAlaTyr.add_bond(21, 22)
    glyAlaTyr.add_bond(21, 26)
    glyAlaTyr.add_bond(22, 27)
    glyAlaTyr.add_bond(23, 24)
    glyAlaTyr.add_bond(24, 25)
    glyAlaTyr.add_bond(24, 33)
    glyAlaTyr.add_bond(25, 26)
    glyAlaTyr.add_bond(25, 30)
    glyAlaTyr.add_bond(26, 31)
    glyAlaTyr.add_bond(27, 28)
    glyAlaTyr.add_bond(27, 32)
    glyAlaTyr.add_bond(28, 29)
    glyAlaTyr.add_bond(31, 35)
    glyAlaTyr.add_bond(31, 32)
    glyAlaTyr.add_bond(32, 36)
    glyAlaTyr.add_bond(33, 34)
    glyAlaTyr.add_bond(33, 37)
    glyAlaTyr.add_bond(37, 38)
    return glyAlaTyr


def test_id_n_term_GAY(gay):
    assert bid.id_n_term(gay) == [1]


def test_find_atom_GAY(gay):
    assert bid.find_atom(gay, [1], "C") == [4]


def test_small_backbone_iter_GAY(gay):
    assert bid.small_backbone_iter(gay, [[1]], "C") == [[1, 4]]


def test_large_backbone_iter_GAY(gay):
    assert bid.large_backbone_iter(gay, [[1]]) == (
        [[1, 4, 6, 8]],
        False,
    )


def test_find_backbone_GAY(gay):
    assert bid.find_backbone(gay) == [1, 4, 6, 8, 11, 15, 19, 24, 33, 37]


def test_label_backbone_GAY(gay):
    bid.label_backbone(gay)
    assert gay.get_backbone() == [
        False,
        True,
        False,
        False,
        True,
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
        False,
        True,
        False,
        False,
        False,
        False,
        True,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        True,
        False,
        False,
        False,
        True,
        False,
    ]
