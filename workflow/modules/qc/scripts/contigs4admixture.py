import sys
import shutil


def generate_mapping(input_file, bim_file, output_file):
    conversion_dict = {}
    with open(input_file, "r") as f:
        for line in f:
            line = line.strip().split()
            conversion_dict[line[0]] = line[1]

    # Copy original bim file to a new file with ".orig" appended to its name
    orig_bim_file = bim_file + ".orig"
    shutil.copyfile(bim_file, orig_bim_file)

    # read bim file and replace the scaffold names with numbering 1:n (n = number of scaffolds)
    updated_lines = []
    with open(bim_file, "r") as f:
        for line in f:
            elements = line.strip().split("\t")
            scaffold = elements[0]
            if scaffold in conversion_dict:
                elements[0] = conversion_dict[scaffold]
            updated_lines.append("\t".join(elements))

    with open(output_file, "w") as f:
        for line in updated_lines:
            f.write(line + "\n")


input_file = snakemake.input.fai
bim_file = snakemake.input.bim
output_file = snakemake.output.bim
generate_mapping(input_file, bim_file, output_file)
