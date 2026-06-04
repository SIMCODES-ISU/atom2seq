from atom2seq import parser
from atom2seq.classes import Atom, Mol

out_basic = parser.parse("./assets/water.xyz")
assert out_basic == Mol(
    [
        Atom("O", (0, 0, 0)),
        Atom("H", (0.758602, 0, 0.504284)),
        Atom("H", (0.758602, 0, -0.504284)),
    ],
    [],
)

out_no_num = parser.parse("./assets/water_no_number_of_atoms.xyz")
assert out_no_num == Mol(
    [
        Atom("O", (0, 0, 0)),
        Atom("H", (0.758602, 0, 0.504284)),
        Atom("H", (0.758602, 0, -0.504284)),
    ],
    [],
)

out_extra_lines = parser.parse("./assets/water_extra_lines.xyz")
assert out_extra_lines == Mol(
    [
        Atom("O", (0, 0, 0)),
        Atom("H", (0.758602, 0, 0.504284)),
        Atom("H", (0.758602, 0, -0.504284)),
    ],
    [],
)

out_no_lines = parser.parse("./assets/water_no_lines.xyz")
assert out_no_lines == Mol(
    [
        Atom("O", (0, 0, 0)),
        Atom("H", (0.758602, 0, 0.504284)),
        Atom("H", (0.758602, 0, -0.504284)),
    ],
    [],
)
