import pytest

from atom2seq.backbone_identifier import label_backbone
from atom2seq.parsers import parse_xyz
from atom2seq.side_chain_identifier import id_side_chains

# Will uncomment once merged into main with Bushman's testing_data branch
# from atom2seq.bond_detector import bond_mol
# from atom2seq.parsers import parse_nwc


prefix = __file__.removesuffix("test_side_chain_identifier.py")


aas = [
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
]


def test_aa_monomers():
    # Will uncomment when merged into main with Bushman's testing_data branch
    # for letter in aas:
    #     mol = parse_nwc(
    #         prefix + f"assets/side_chain_identifier_tests/mon{letter}.nwc"
    #     )  # noqa
    #     bond_mol(mol)
    #     label_backbone(mol)
    #     assert id_side_chains(mol) == [letter]
    pass


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
    label_backbone(gas)
    return gas


def test_gas(gas):
    assert id_side_chains(gas) == ["G", "A", "S"]
