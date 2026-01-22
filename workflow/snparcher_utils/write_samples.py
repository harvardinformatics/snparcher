import gzip
import shutil
import argparse
from pathlib import Path
from typing import TextIO
# User provides list of sample names, 1 per line
# User provides path to where fastq files are. Assume paired end and that file name has sample name in it
# User provides path to reference genome. Need to copy this to proper path


def read_sample_list(sample_fh: TextIO) -> list:
    return sample_fh.read().splitlines()


def find_sample_fastqs(samples: list, fastq_dir: Path) -> dict:
    """Searches fastq_dir for sample names and associates in a dict"""
    sample_fastq_paths = {}
    cant_find = []
    for samp_name in samples:
        fqs = sorted(
            list(fastq_dir.glob(f"*{samp_name}*"))
        )  # Hoping that sorting will make fq1 first.
        if len(fqs) != 2:
            cant_find.append(samp_name)
        else:
            sample_fastq_paths[samp_name] = fqs
    return sample_fastq_paths, cant_find


def copy_reference(ref: Path) -> str:
    exts = [".fna", ".fa", ".fasta"]
    for ext in exts:
        if ext in ref.name:
            ref_name = ref.name.split(ext)[0]
    if Path("..", "data", "genome", ref_name + ".fna").exists():
        return ref_name
    if not Path("../data/genome").exists():
        Path("../data/genome").mkdir(parents=True)
    if ref.suffix == ".gz":
        with gzip.open(ref, "rb") as f_in:
            with open(Path("..", "data", "genome", ref_name + ".fna"), "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
    else:
        shutil.copyfile(ref, Path("data", "genome", ref_name + ".fna"))
    return ref_name


def write_sample_sheet(sample_dict: dict, ref_name: str, ref_path: str, ncbi_ref: bool) -> None:
    """Writes the sample sheet"""
    with open(Path("../config", "samples.csv"), "w") as out:
        if ncbi_ref:
            out.write("BioSample,LibraryName,refGenome,Run,BioProject,fq1,fq2\n")
            for i, (k, v) in enumerate(sample_dict.items()):
                out.write(f"{k},lib_{k},{ref_name},{i},NaN,{v[0]},{v[1]}\n")
        else:
            out.write("BioSample,LibraryName,refGenome,refPath,Run,BioProject,fq1,fq2\n")
            for i, (k, v) in enumerate(sample_dict.items()):
                out.write(f"{k},lib_{k},{ref_name},{ref_path},{i},NaN,{v[0]},{v[1]}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Write sample files.")
    parser.add_argument(
        "-s", "--sample_list", dest="samp", required=True, help="Specify path to sample list"
    )
    parser.add_argument(
        "-f", "--fastq_dir", dest="fastq", required=True, help="Specify path to fastq dir"
    )
    parser.add_argument(
        "-c",
        "--copy",
        dest="copyref",
        required=False,
        default=False,
        help="Copy reference genome to data/genome dir and unzip.",
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-r",
        "--ref",
        dest="ref",
        help="Specify path to reference genome. Mutually exclusive with -a/--acc.",
    )
    group.add_argument(
        "-a",
        "--acc",
        dest="acc",
        help="Specify reference genome accession. Mutually exclusive with -r/--ref",
    )
    args = parser.parse_args()

    sample_list = args.samp
    fastq_dir = Path(args.fastq)

    with open(sample_list, "r") as f:
        samples = read_sample_list(f)

    sample_dict, cant_find = find_sample_fastqs(samples, fastq_dir)
    ncbi_ref = True

    if args.ref:
        ref = Path(args.ref)
        ncbi_ref = False
        if args.copyref:
            ref_name = copy_reference(ref)
            ref_path = "../data/genome/" + ref_name + ".fna"
        else:
            ref_name = ref.stem
            ref_path = args.ref
    else:
        ref_name = args.acc
        ref_path = ""

    write_sample_sheet(sample_dict, ref_name, ref_path, ncbi_ref)

    if cant_find:
        print("Couldnt' find fastqs for these files:")
        for name in cant_find:
            print(name)


if __name__ == "__main__":
    main()
