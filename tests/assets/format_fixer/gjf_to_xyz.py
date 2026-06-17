for prefix in ["mono", "di", "tri"]:
    filebase = "{0}mer".format(prefix)
    with open(
        r"C:\Users\eb528353\OneDrive - Loras College\Desktop\atom2seq\tests\assets\{0}mers\{0}mers_gjf_files\{1}.gjf".format(  # noqa
            prefix, filebase
        )
    ) as f:
        content = f.read()
        # This is the path for any file I'm searching. Since I have monomers,
        # dimers, trimers, and anymer_gjf_files, what can I do to account for
        # all of them? E.g. can I use a wildcard or something similar?
        # WILL THIS WORK?
    cleaned_content = (
        content.replace(" ", "")
        .replace("\n", "")
        .replace("\\", "\n")
        .replace(",-", "    -")
        .replace(",", "     ")
    )  # noqa
"""
    This reads the file and pulls out, in this order, white spaces and
    replaces them with nothing; new lines and replaces them with nothing;
    backslashes and replaces them with a new line; commas immediately
    followed by a negative sign and replaces them with a tab and keeps the
    negative sign; and commas not immediately followed by a negative sign
    and replaces it with a tab and one space. This is all to format the gjf
    file into something more like an xyz file.
"""
#        print(cleaned_content)
"""
    This is to print everything to terminal.
"""

# Each coordinate has 11 significant figures.

# path = r"C:\Users\eb528353\OneDrive - Loras College\Desktop\atom2seq\tests\assets\{0}mers\{0}mers_xyz_files".format  # noqa
# filname = path + {} + "\\.xyz"

# with open("*.xyz", "w") as f:
#     print("Filename:", Filename, file=f)
# I want to print cleaned_content into a file with the same path as is in
# line 2, but with an .xyz extension.

# Also the goal of this script is to be able to pass a file through, e.g. in
# terminal and get
