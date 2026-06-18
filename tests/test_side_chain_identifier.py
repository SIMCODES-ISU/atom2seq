from atom2seq.backbone_identifier import label_backbone
from atom2seq.bond_detector import bond_mol
from atom2seq.parsers import parse_nwc
from atom2seq.side_chain_identifier import id_side_chains

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
    for letter in aas:
        mol = parse_nwc(
            prefix + f"assets/side_chain_identifier_tests/mon{letter}.nwc"
        )  # noqa
        bond_mol(mol)
        label_backbone(mol)
        assert id_side_chains(mol) == [letter]
