import math

cov_thresh = {}
stdev_scale = 2
cov_thresh_stdev = 2

with open(snakemake.input["stats"]) as stats:
    for line in stats:
        if "mean" in line:
            continue

        fields = line.split()
        mean = float(fields[1])
        stdev = math.sqrt(mean)
        cov_thresh[fields[0]] = {
            "low": mean - (stdev * stdev_scale),
            "high": mean + (stdev * stdev_scale),
        }

with open(snakemake.output[0], "w") as output_file:
    output_file.write("chrom\tmin\tmax\n")
    for chrom, thresholds in cov_thresh.items():
        if chrom == "total":
            continue
        output_file.write(f"{chrom}\t{thresholds['low']}\t{thresholds['high']}\n")
