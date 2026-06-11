from atom2seq.bean_parser import parse_bean
from atom2seq.classes import Atom, Mol

prefix = __file__.removesuffix("test_bean_parser.py")

water = prefix + "assets/water.bean"


def test_bean_parser():
    assert parse_bean(water) == Mol(
        [Atom("O", (0, 0, 0)), Atom("H", (0, 1, 0)), Atom("H", (1, 0, 0))], []
    )
