rule parse_bam_stats:
    input:
        coverage="results/qc_metrics/bam/{sample}_coverage.txt",
        flagstat="results/qc_metrics/bam/{sample}_flagstat.txt",
    output:
        json="results/qc_metrics/bam/{sample}.json",
    log:
        "logs/parse_bam_stats/{sample}.txt"
    run:
        import json
        
        # Parse coverage - weighted average across scaffolds
        num_sites = []
        depths = []
        covered_bases = 0
        
        with open(input.coverage) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                start, end = int(fields[1]), int(fields[2])
                num_sites.append(end - start + 1)
                depths.append(float(fields[6]))
                covered_bases += int(fields[4])
        
        total_sites = sum(num_sites)
        mean_depth = sum(d * n / total_sites for d, n in zip(depths, num_sites)) if total_sites > 0 else 0
        
        # Parse flagstat
        with open(input.flagstat) as f:
            lines = f.readlines()
        
        total_reads = int(lines[0].split()[0])
        num_dups = int(lines[4].split()[0])
        num_mapped = int(lines[6].split()[0])
        pct_mapped = float(lines[7].split()[0].strip("%")) if total_reads > 0 else 0
        pct_proper_paired = float(lines[14].split()[0].strip("%")) if total_reads > 0 else 0
        
        out = {
            "sample": wildcards.sample,
            "total_reads": total_reads,
            "num_mapped": num_mapped,
            "percent_mapped": pct_mapped,
            "num_duplicates": num_dups,
            "percent_duplicates": (num_dups / total_reads * 100) if total_reads > 0 else 0,
            "percent_properly_paired": pct_proper_paired,
            "mean_depth": mean_depth,
            "covered_bases": covered_bases,
        }
        
        with open(output.json, "w") as f:
            json.dump(out, f, indent=2)


rule parse_sentieon_stats:
    input:
        insert_file="results/qc_metrics/sentieon/{sample}_insert_metrics.txt",
    output:
        json="results/qc_metrics/sentieon/{sample}.json",
    log:
        "logs/parse_sentieon_stats/{sample}.txt"
    run:
        import json
        
        median_insert = None
        median_abs_dev = None
        
        with open(input.insert_file) as f:
            for i, line in enumerate(f):
                if i == 2:  # Data line
                    fields = line.strip().split()
                    if len(fields) >= 2:
                        median_insert = float(fields[0])
                        median_abs_dev = float(fields[1])
                    break
        
        out = {
            "sample": wildcards.sample,
            "median_insert_size": median_insert,
            "median_abs_dev_insert": median_abs_dev,
        }
        
        with open(output.json, "w") as f:
            json.dump(out, f, indent=2)


def combine_qc_input(wildcards):
    """Get all QC metric files based on tool choice."""
    inputs = {
        "fastp": expand("results/qc_metrics/fastp/{sample}.json", sample=SAMPLES_WITH_FASTQ),
        "bam": expand("results/qc_metrics/bam/{sample}.json", sample=SAMPLES_WITH_BAM),
    }
    if USE_SENTIEON:
        inputs["sentieon"] = expand("results/qc_metrics/sentieon/{sample}.json", sample=SAMPLES_WITH_BAM)
    return inputs


rule combine_qc_metrics:
    input:
        unpack(combine_qc_input),
    output:
        tsv="results/qc_metrics/qc_report.tsv",
    log:
        "logs/combine_qc_metrics.txt"
    run:
        import json
        
        # Load fastp stats
        fastp_stats = {}
        for fn in input.fastp:
            sample = os.path.basename(fn).replace(".json", "")
            with open(fn) as f:
                fastp_stats[sample] = json.load(f)
        
        # Load bam stats
        bam_stats = {}
        for fn in input.bam:
            sample = os.path.basename(fn).replace(".json", "")
            with open(fn) as f:
                bam_stats[sample] = json.load(f)
        
        # Load sentieon stats if available
        sentieon_stats = {}
        if hasattr(input, "sentieon"):
            for fn in input.sentieon:
                sample = os.path.basename(fn).replace(".json", "")
                with open(fn) as f:
                    sentieon_stats[sample] = json.load(f)
        
        # Build header
        header = [
            "sample",
            "reads_before_filtering",
            "reads_after_filtering",
            "fraction_passed",
            "total_reads",
            "percent_mapped",
            "num_duplicates",
            "percent_duplicates",
            "percent_properly_paired",
            "mean_depth",
            "covered_bases",
        ]
        if sentieon_stats:
            header.extend(["median_insert_size", "median_abs_dev_insert"])
        
        # Write report
        with open(output.tsv, "w") as f:
            f.write("\t".join(header) + "\n")
            
            for sample in SAMPLES_WITH_BAM:
                row = [sample]
                
                # Fastp stats (may not exist for bam input type)
                if sample in fastp_stats:
                    fs = fastp_stats[sample]
                    row.extend([
                        str(fs["summary"]["before_filtering"]["total_reads"]),
                        str(fs["summary"]["after_filtering"]["total_reads"]),
                        f"{fs['summary']['after_filtering']['total_reads'] / fs['summary']['before_filtering']['total_reads']:.4f}",
                    ])
                else:
                    row.extend(["NA", "NA", "NA"])
                
                # BAM stats
                bam = bam_stats[sample]
                row.extend([
                    str(bam["total_reads"]),
                    f"{bam['percent_mapped']:.2f}",
                    str(bam["num_duplicates"]),
                    f"{bam['percent_duplicates']:.2f}",
                    f"{bam['percent_properly_paired']:.2f}",
                    f"{bam['mean_depth']:.2f}",
                    str(bam["covered_bases"]),
                ])
                
                # Sentieon stats
                if sentieon_stats:
                    if sample in sentieon_stats:
                        sent = sentieon_stats[sample]
                        row.extend([
                            str(sent["median_insert_size"]) if sent["median_insert_size"] else "NA",
                            str(sent["median_abs_dev_insert"]) if sent["median_abs_dev_insert"] else "NA",
                        ])
                    else:
                        row.extend(["NA", "NA"])
                
                f.write("\t".join(row) + "\n")