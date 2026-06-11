from atom2seq.xyz_parser import parse_xyz
from atom2seq.classes import Atom, Mol

prefix = __file__.removesuffix("test_xyz_parser.py")


def test_basic():
    out_basic = parse_xyz(prefix + "assets/water.xyz")
    assert out_basic == Mol(
        [
            Atom("O", (0, 0, 0)),
            Atom("H", (0.758602, 0, 0.504284)),
            Atom("H", (0.758602, 0, -0.504284)),
        ],
        [],
    )


def test_no_num():
    out_no_num = parse_xyz(prefix + "assets/water_no_number_of_atoms.xyz")
    assert out_no_num == Mol(
        [
            Atom("O", (0, 0, 0)),
            Atom("H", (0.758602, 0, 0.504284)),
            Atom("H", (0.758602, 0, -0.504284)),
        ],
        [],
    )


def test_extra_lines():
    out_extra_lines = parse_xyz(prefix + "assets/water_extra_lines.xyz")
    assert out_extra_lines == Mol(
        [
            Atom("O", (0, 0, 0)),
            Atom("H", (0.758602, 0, 0.504284)),
            Atom("H", (0.758602, 0, -0.504284)),
        ],
        [],
    )


def test_no_lines():
    out_no_lines = parse_xyz(prefix + "assets/water_no_lines.xyz")
    assert out_no_lines == Mol(
        [
            Atom("O", (0, 0, 0)),
            Atom("H", (0.758602, 0, 0.504284)),
            Atom("H", (0.758602, 0, -0.504284)),
        ],
        [],
    )
