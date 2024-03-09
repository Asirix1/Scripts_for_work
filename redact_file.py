file_name = "GSE237721_ENCFF498SJT_pairs_GRCh38.pairs"

with open(file_name, "r") as file:
    lines = file.readlines()

new_lines = []
for line in lines:
    if not line.startswith("#"):
        break
    new_lines.append(line)

with open(file_name, "w") as file:
    file.writelines(new_lines)