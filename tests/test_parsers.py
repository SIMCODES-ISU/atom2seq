import pytest

from atom2seq.classes import Atom, Mol
from atom2seq.parsers import parse_cif, parse_gjf, parse_pdb, parse_xyz

prefix = __file__.removesuffix("test_parsers.py")


@pytest.fixture
def water_mol():
    return Mol(
        [Atom("H", (1, 0, 0)), Atom("H", (0, 1, 0)), Atom("O", (0, 0, 0))], []
    )  # noqa


def test_gjf_parser(water_mol):
    assert parse_gjf(prefix + "assets/parser_tests/water.bean") == water_mol


def test_xyz_basic(water_mol):
    out_basic = parse_xyz(prefix + "assets/parser_tests/water.xyz")
    assert out_basic == water_mol


def test_xyz_no_num(water_mol):
    out_no_num = parse_xyz(
        prefix + "assets/parser_tests/water_no_number_of_atoms.xyz"
    )  # noqa
    assert out_no_num == water_mol


def test_xyz_extra_lines(water_mol):
    out_extra_lines = parse_xyz(
        prefix + "assets/parser_tests/water_extra_lines.xyz"
    )  # noqa
    assert out_extra_lines == water_mol


def test_xyz_no_lines(water_mol):
    out_no_lines = parse_xyz(prefix + "assets/parser_tests/water_no_lines.xyz")
    assert out_no_lines == water_mol


def test_cif_parser(water_mol):
    assert parse_cif(prefix + "assets/parser_tests/water.cif") == water_mol


def test_cif_extra_info(water_mol):
    assert (
        parse_cif(prefix + "assets/parser_tests/water_extra_info.cif")
        == water_mol  # noqa
    )  # noqa


def test_pdb_parser(water_mol):
    assert parse_pdb(prefix + "assets/parser_tests/water.pdb") == water_mol


def test_pdb_extra_info(water_mol):
    assert (
        parse_pdb(prefix + "assets/parser_tests/water_extra_info.pdb")
        == water_mol  # noqa
    )  # noqa
