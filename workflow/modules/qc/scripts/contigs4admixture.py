import csv
import shutil


def load_admixture_map(map_file):
    conversion_dict = {}
    with open(map_file, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"plink_contig", "admixture_id"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError("contig_map.tsv must contain plink_contig and admixture_id columns")
        for row in reader:
            conversion_dict[row["plink_contig"]] = str(row["admixture_id"])
    return conversion_dict


def generate_mapping(map_file, bim_file, output_file):
    conversion_dict = load_admixture_map(map_file)

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


input_file = snakemake.input.contig_map
bim_file = snakemake.input.bim
output_file = snakemake.output.bim
generate_mapping(input_file, bim_file, output_file)
