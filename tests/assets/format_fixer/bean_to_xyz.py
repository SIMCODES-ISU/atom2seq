with open(
    r"C:\Users\eb528353\OneDrive - Loras College\Desktop\atom2seq\tests\assets\dimers\dimers_gjf_files\dimethionine.bean"  # noqa
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
