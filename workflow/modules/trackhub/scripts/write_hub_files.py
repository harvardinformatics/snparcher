from os.path import basename
import shutil
import os

# https://genome.ucsc.edu/goldenPath/help/hgTracksHelp.html#UseOneFile
hub_text = """hub {genome}
shortLabel {genome} snpArcher Track Hub
longLabel {genome} snpArcher Track Hub
useOneFile on
descriptionUrl index.html
email {email}\n
genome {genome}\n"""

vcf_track_txt = """track VCF
bigDataUrl {vcf_file}
shortLabel VCF
longLabel VCF
visibility squish
html index.html
type vcfTabix\n"""

window_parent_txt = """track {track_type}
compositeTrack on
shortLabel {track_type}
longLabel {track_type}
color {color}
altColor 0,102,204
autoScale on
type bigWig
allButtonPair on
html index.html
visibility full\n"""

window_track_txt = """track {track_name}
parent {parent} on
bigDataUrl {data_url}
type bigWig
visibility {vis}
shortLabel {label}
longLabel {label}\n"""

allele_freq_txt = """track MinorAlleleFrequency
bigDataUrl {data_url}
type bigWig
color 88,85,120
altColor 0,102,204
autoScale on
visibility full
shortLabel Minor Allele Frequency
html index.html
longLabel Minor Allele Frequency\n"""

snp_depth_txt = """track SNPDepth
bigDataUrl {data_url}
type bigWig
color 120,172,145
altColor 0,102,204
autoScale on
visibility full
shortLabel SNP Depth
html index.html
longLabel SNP Depth\n"""

coverage_track_txt = """track NonCallableSites
bigDataUrl {cov_file}
shortLabel Non Callable Sites
type bigBed
longLabel Non Callable Sites
color 0,0,0
html index.html
visibility dense\n"""

COLORS = {
    "Tajima": "(70,130,180)",
    "SNP-Density": "(186,85,211)",
    "Pi": "(248,174,51)",
}

html = """<!DOCTYPE html>
<html>
<head>
  <title>snpArcher Track Hub Description</title>
  <style>
    body {
      font-family: Arial, sans-serif;
      margin: 20px;
    }
    h1 {
      font-size: 24px;
      font-weight: bold;
    }
    h2 {
      font-size: 20px;
      font-weight: bold;
      margin-top: 20px;
    }
    h3 {
      font-size: 16px;
      font-weight: bold;
      margin-top: 10px;
      margin-left: 20px;
    }
    p {
      font-size: 16px;
      margin-left: 40px;
    }
  </style>
</head>
<body>
  <h1>snpArcher Track Hub Description</h1>
  
  <h2>Introduction</h2>
  <p>To facilitate downstream data exploration and as an example of the module development components of this work, we
    developed a module to generate UCSC Genome Browser track files to explore population variation data (<a
        href="https://doi.org/10.1093/molbev/msad270">see paper for details</a>).</p>
  
  <h2>Track Descriptions</h2>
  
  <h3>Tajima’s D</h3>
  <p>This track provides windowed estimates of Tajima’s D, a population genetic statistic that measures the departure from neutral evolution in a DNA sequence.</p>
  
  <h3>SNP Density</h3>
  <p>This track displays the density of single nucleotide polymorphisms (SNPs) across the genome, showing regions with high or low levels of genetic variation.</p>
  
  <h3>Pi</h3>
  <p>The Pi track represents the average number of nucleotide differences per site between any two sequences in a population, providing an estimate of genetic diversity.</p>
  
  <h3>Minor Allele Frequency</h3>
  <p>This track shows the frequency of the less common allele at a SNP locus, providing insights into the genetic variation within a population.</p>
  
  <h3>SNP Depth</h3>
  <p>The SNP Depth track displays the number of reads or sequencing depth at each SNP position, indicating the coverage and quality of the variant calls.</p>
  
  <h3>Non Callable Sites</h3>
  <p>The Non Callable Sites track highlights regions in the genome that are considered non-callable, meaning that they have low sequencing coverage or other technical limitations that make it difficult to accurately determine genetic variation in those regions.</p>

</body>
</html>




"""


def human_format(num):
    num = float("{:.3g}".format(num))
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
    return "{}{}".format(
        "{:f}".format(num).rstrip("0").rstrip("."), ["", "K", "M", "B", "T"][magnitude]
    )


def main():
    file_types = snakemake.params["file_types"]  # noqa: F821
    email = snakemake.params["email"]  # noqa: F821
    trackhub_windows = snakemake.params["windows"]  # noqa: F821
    vcf_file = basename(snakemake.input["vcf"][0])  # noqa: F821
    cov_file = basename(snakemake.input["callable_sites"][0])  # noqa: F821
    freq_file = basename(snakemake.input["allele_freq"][0])  # noqa: F821
    depth_file = basename(snakemake.input["depth"][0])  # noqa: F821
    genome = snakemake.params["refGenome"]  # noqa: F821
    trackhub_file = snakemake.output["trackhub_file"]  # noqa: F821
    html_file = snakemake.output["html"]  # noqa: F821

    with open(html_file, "w") as f:
        f.write(html)

    with open(trackhub_file, "w") as out:
        print(hub_text.format(genome=genome, email=email), file=out)
        print(vcf_track_txt.format(vcf_file=vcf_file), file=out)
        print(coverage_track_txt.format(cov_file=cov_file), file=out)
        print(allele_freq_txt.format(data_url=freq_file), file=out)
        print(snp_depth_txt.format(data_url=depth_file), file=out)

        for file in file_types:
            print(window_parent_txt.format(track_type=file, color=COLORS[file]), file=out)
            for window in trackhub_windows:
                track_name = f"{file}_{human_format(window)}_bp_bins"
                label = f"{file}_{human_format(window)}_bp bins"
                url = f"{file}_{window}.bw"
                if window == 1000:
                    vis = "full"
                else:
                    vis = "hide"
                print(
                    window_track_txt.format(
                        track_name=track_name,
                        label=label,
                        parent=file,
                        data_url=url,
                        vis=vis,
                    ),
                    file=out,
                )


if __name__ == "__main__":
    main()
