from atom2seq.detect_primary_structure import detect_primary_structure

prefix = __file__.removesuffix("glycine.py")

print(
    "This file contains information for the peptide ",
    detect_primary_structure(
        prefix + "tests/assets/bond_detector_tests/glycine.xyz", "xyz"
    ),
)
