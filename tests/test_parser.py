from atom2seq import parser
from atom2seq.classes import Atom, Mol

prefix = __file__.removesuffix("test_parser.py")
print(prefix)

out_basic = parser.parse(prefix + "assets/water.xyz")


def test_basic():
    assert out_basic == Mol(
        [
            Atom("O", (0, 0, 0)),
            Atom("H", (0.758602, 0, 0.504284)),
            Atom("H", (0.758602, 0, -0.504284)),
        ],
        [],
    )


out_no_num = parser.parse(prefix + "assets/water_no_number_of_atoms.xyz")


def test_no_num():
    assert out_no_num == Mol(
        [
            Atom("O", (0, 0, 0)),
            Atom("H", (0.758602, 0, 0.504284)),
            Atom("H", (0.758602, 0, -0.504284)),
        ],
        [],
    )


out_extra_lines = parser.parse(prefix + "assets/water_extra_lines.xyz")


def test_extra_lines():

    assert out_extra_lines == Mol(
        [
            Atom("O", (0, 0, 0)),
            Atom("H", (0.758602, 0, 0.504284)),
            Atom("H", (0.758602, 0, -0.504284)),
        ],
        [],
    )


out_no_lines = parser.parse(prefix + "assets/water_no_lines.xyz")


def test_no_lines():
    assert out_no_lines == Mol(
        [
            Atom("O", (0, 0, 0)),
            Atom("H", (0.758602, 0, 0.504284)),
            Atom("H", (0.758602, 0, -0.504284)),
        ],
        [],
    )
