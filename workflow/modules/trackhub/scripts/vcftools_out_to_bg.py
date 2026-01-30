from pathlib import Path


def chrom_dict(chrom_sizes_file):
    chroms = {}
    with open(chrom_sizes_file, "r") as f:
        for line in f:
            line = line.strip().split()
            chroms[line[0]] = int(line[1])
    return chroms


def parse_stat_file(stat_file, out_file, chrom_sizes):
    stat_file = Path(stat_file)
    file_type = stat_file.suffix
    window = int(stat_file.stem)

    with open(out_file, "w") as out:
        results = []
        with open(stat_file, "r") as inp:
            next(inp)
            for line in inp:
                line = line.strip().split()
                chrom = line[0]
                if chrom not in chrom_sizes:
                    continue
                else:
                    start = int(line[1])
                    end = start + (window - 1)
                    if end >= chrom_sizes[chrom]:
                        end = chrom_sizes[chrom] - 1

                    if file_type == ".Tajima":
                        value = line[3]
                    elif file_type == ".SNP-Density":
                        value = line[2]
                    elif file_type == ".Pi":
                        value = line[4]
                    else:
                        raise (ValueError(f"Unknown file type: {file_type}"))

                    results.append((chrom, start, end, value))

        sorted_results = sorted(results, key=lambda x: (x[0], x[1]))

        for chrom, start, end, value in sorted_results:
            print(f"{chrom}\t{start}\t{end}\t{value}\n", file=out)


def main():
    chrom_sizes = chrom_dict(snakemake.input["chrom_sizes"])
    parse_stat_file(
        stat_file=snakemake.input["stat_file"],
        out_file=snakemake.output[0],
        chrom_sizes=chrom_sizes,
    )


if __name__ == "__main__":
    main()
