def parse_pdb(file_name):
    file_name.open()
    lines = file_name.readlines()
    file_name.close()


[line.strip() for line in lines]
