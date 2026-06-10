from atom2seq import parser
from atom2seq.classes import Atom, Mol

prefix = __file__.removesuffix("test_parser.py")


def test_basic():
    out_basic = parser.parse(prefix + "assets/water.xyz")
    assert out_basic == Mol(
        [
            Atom("O", (0, 0, 0)),
            Atom("H", (0.758602, 0, 0.504284)),
            Atom("H", (0.758602, 0, -0.504284)),
        ],
        [],
    )


def test_no_num():
    out_no_num = parser.parse(prefix + "assets/water_no_number_of_atoms.xyz")
    assert out_no_num == Mol(
        [
            Atom("O", (0, 0, 0)),
            Atom("H", (0.758602, 0, 0.504284)),
            Atom("H", (0.758602, 0, -0.504284)),
        ],
        [],
    )


def test_extra_lines():
    out_extra_lines = parser.parse(prefix + "assets/water_extra_lines.xyz")
    assert out_extra_lines == Mol(
        [
            Atom("O", (0, 0, 0)),
            Atom("H", (0.758602, 0, 0.504284)),
            Atom("H", (0.758602, 0, -0.504284)),
        ],
        [],
    )


def test_no_lines():
    out_no_lines = parser.parse(prefix + "assets/water_no_lines.xyz")
    assert out_no_lines == Mol(
        [
            Atom("O", (0, 0, 0)),
            Atom("H", (0.758602, 0, 0.504284)),
            Atom("H", (0.758602, 0, -0.504284)),
        ],
        [],
    )


# Adding a comment so I can commit this
