with open(
    r"C:\Users\eb528353\OneDrive - Loras College\Desktop\atom2seq\tests\assets\trimers\trimers_gjf_files\trivaline.bean"  # noqa
) as f:
    content = f.read()
cleaned_content = (
    content.replace(" ", "")
    .replace("\n", "")
    .replace("\\", "\n")
    .replace(",-", "    -")
    .replace(",", "     ")
)  # noqa


print(cleaned_content)
