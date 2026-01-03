# Setting up snpArcher
## Environment Setup
First, you will need to have [Mamba](https://mamba.readthedocs.io/en/latest/mamba-installation.html#mamba-install) installed. Follow the link and use the "Fresh Install (recommended)" directions. 

Mamba is a faster version of conda. Conda is a package manager that makes it easy to create local environments with pre-configured versioning for your favorite packages. 

Once Mamba is installed, create a conda environment with snakemake. These are the only two dependencies you need for the pipeline to work, the workflow will create mamba environments for each rule, and there is no need to install each package separately. 

```
mamba create -c conda-forge -c bioconda -n snparcher "snakemake>=8" "python==3.11.4"
mamba activate snparcher
```
If you encounter issues, please see the [Snakemake docs](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for detailed installation instructions.

Next, clone the [snpArcher github repo](https://github.com/harvardinformatics/snpArcher) to your machine:
```
git clone https://github.com/harvardinformatics/snpArcher.git
cd snpArcher
```

## Creating a sample sheet
In order to determine what outputs to create, snpArcher requires sample sheet file. This comma-separated file contains the required sample metadata about your samples in order to run the workflow. At a minimum, the snpArcher pipeline requires that each sample have a unique sample name, a reference genome accession or a path to a fasta file, and an SRA accession, or path to two paired end fastq files.

Below are all of the accepted fields for a sample sheet:
| Field | Description |
| ---- | -------------|
| BioSample | The name of the sample. This will be the sample name in the final VCF |
| LibraryName | LibraryID for sample, this can be the same or different than BioSample |
| Run | The SRR for the sample, if applicable. If not, must be some **unique** value. It is often the lane number if samples are sequenced on multiple lanes. |
| refGenome | Reference genome accession, if applicable. *See note* |
| refPath | Path to local reference genome, if applicable. *See note* |
| BioProject | If applicable. Otherwise any value is acceptable. |
| fq1 | Optional if no SRR value in Run. Path to read 1 for sample |
| fq2 | Optional if no SRR value in Run. Path to read 2 for sample |
| SampleType | Optional. Triggers postproccesing module. Accepted values are 'include' or 'exclude' |
| lat | Optional. Decimal latitude for sample, required to generate map in QC dashboard. |
| long | Optional. Decimal longitude for sample, required to generate map in QC dashboard. |

```{note}
refGenome is always required. refPath specifying the path to a reference fasta file is optional, but when specified, a name for the assembly (in refGenome) must also be included. 

If you are using the same reference genome for all samples in your sample sheet, you can omit the refGenome and/or refPath column from the sample sheet and specify these fields in the config file. See [config setup below](#configuring-snparcher) for more details.
```

It is important to note that samples are proccessed together based on their `refGenome` metadata, so **all BioSamples that share a reference genome will ultimately end up in the same final vcf file.** If you are mapping multiple populations / species to a single reference genome, and want separate VCF files for each population / species, you will need to split your final vcf after the pipeline completes, or run multiple indpendent sample sheets in different results directories. 

If your reads (and, optionally, your local reference genome) are stored in somewhere seperate of the workflow (e.g.: a scratch disk) then you can specify the path to your reads using the `fq1` and `fq2` fields, and the location of your reference genome fasta in the `refPath` field. 

### Using data from NCBI SRA
If you'd like to reanalyze an existing NCBI SRA BioProject, please follow these instructions to quickly create a sample sheet.

1. Go to the BioProject overview web page on the SRA.
2. In the subheading `Project Data` there is a table with the columns `Resource Name` and `Number of Links`. Click the link in the `Number of Links` column in the `SRA Experiments` row. You will be redirected to a search results page.
3. Near the top of the search results page, click the link `Send results to Run Selector`
4. On the Run Selector page, you can select/deselect samples you'd like to include/exclude in your sample sheet by using the checkboxes.
5. Once you are done selecting samples, click the `Metadata` button in the `Download` column in the table near the middle of the page.
6. This will download a a comma separated file called `SraRunTable.txt`.
7. Open the file in the editor of your choice, and add a column named `refGenome`. In this column, enter the reference genome accession you want to use for every row in the sheet.
8. Save the sample sheet, it is now ready to use.

### Using local data
A python script `workflow/snparcher_utils/write_samples.py` is included to help write the sample sheet for you. In order to use this script, you must have organized all of your fastq files in to one directory. The script requies you provide a file with one sample per name that maps uniquely to a pair of fastq files in the afformentioned directory. The script also requires either a reference genome accession or path to reference fasta. 

```{note}
This script cannot currently handle multiple sequencing runs per sample. Please see below for how to handle this case.
```

Usage details: 

|Argument| Description|
| ------ | ---------- |
| `-s / --sample_list` | Path to a sample list. One sample per line |
| `-f / --fastq_dir` | Path to directory containing ALL fastq files. It is assumed that each fastq file will contain the sample name uniquely. |
| `-r / --ref` | Path to reference fasta. Mutually exclusive with -a|
| `-a / --acc` | NCBI accession of reference. Mutually exclusive with -r|

#### Handling samples with more than one pair of reads

In order to specify samples that were sequenced multiple times in your sample sheet, you must:
1. Create a duplicate row for each unit of sequencing
2. Ensure the `BioSample` value is the same across all rows for the sample.
3. Give each row a unique `Run` value. This allows snpArcher to collect all read pairs for a `BioSample`. All runs for a sample will be mapped separately to the genome and subsequently merged.
4. Give each row a unique `LibraryName` value, if applicable. Used for marking duplicates, `LibraryName` should be the same in cases where the same library was sequenced multiple times. If a sample had multiple libraries prepared for it, then `LibraryName` should be unique for each library.  See [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups) for more info. 

For example, consider we have 2 samples: `A` and `B`. `Sample A` was sequenced 3 times, 2 of which were derived from the same library prep, and the other was a unique library. `Sample B `was only sequenced once. Below is how the sample sheet would look in order to define these relationships. Note, only the relevant fields have been included.

| BioSample | LibraryName | Run |
| --------- | ----------- | --- |
| sample_A | lib_A_1 | 1 |
| sample_A | lib_A_1 | 2 |
| sample_A | lib_A_2 | 3 |
| sample_B | lib_B_1 | 4 |


## Configuring snpArcher

Workflow variables such as output file prefix, tool settings, and other options are set in `config/config.yaml`. Resource settings such as threads and memory are controlled per tool in the `profiles/default/config.yaml`.

### Core configuration
The following options in `config/config.yaml` must be set before running snpArcher:

| Option | Description | Type | Required | Default |
| ---- | -------------| ------ |------ | ------ |
| `samples` | Path to CSV sample sheet.| `str` | `True` | `None` |
| `final_prefix` | Prefix to name final output files with (e.g. VCF) | `str` | `True` | `None` |
| `intervals` |  Use SplitByN interval approach for GATK variant calling | `bool` | `True` | `True` |
| `sentieon` | Use Sentieon tools instead of GATK for variant calling | `bool` | `True` | `False` |
| `sentieon_lic` | If using Sentieon tools, provide license here | `str` | `True` if `sentieon==True`| `None` |
| `remote_reads` | Use remote storage provider for reads. | `bool` | `False` | `False`|
| `bigtmp` | Set a directory for TMP. Default is $TMPDIR env var | `str` | `False` | `None` |
| `cov_filter` | Use coverage thresholds for filtering callable sites | `bool`| `True` | `True` |
| `generate_trackhub` | Generate population genomics stats trackhub | `bool`| `True` | `True` |
| `trackhub_email` | Trackhubs require an email address | `str` | `True` if `generate_trackhub==True` | `None` |
| `refGenome` | Reference genome name or accession | `str` | `True` if not provided in sample sheet |  `None` |
| `refPath` | Path to reference genome if not using NCBI genome accession | `str` | `False` | `None` |
| `mark_duplicates` | Mark optical duplicates before variant calling. | `str` | `True` | `True` |
| `sort_reads` | Sort reads by read name before running adapter trimming. | `str` | `True` | `False` |

### Other options
The following options can be adjusted based on your needs and your dataset.

#### Variant Calling Options
| Option | Description | Type |
| ---- | -------------| ------ |
|`minNmer`| The minimum span of Ns to split reference genome at for interval generation | `int`|
|`num_gvcf_intervals` | The maximum number of GVCF intervals to create. Actual number of intervals may be less if reference genome is highly contiguous. | `int`|
|`db_scatter_factor` | Used to calculate number of DB intervals to create. `num_db_intervals = (scatter_factor * num_samples * num_gvcf_intervals)`. Recommend <1 | `float`|
| `minP` | Controls `--min-pruning` in GATK HaplotypeCaller. Recommend 1 for low coverage (<10x), 2 for high coverage (>10x) | `int` |
| `minD` | Controls `--min-dangling-branch-length` in GATK HaplotypeCaller. Recommend 1 for low coverage (<10x), 4 for high coverage (>10x) | `int` |
| `ploidy` | Ploidy for variant calling step. | `int` |

#### Callable Sites Options
| Option | Description | Type |
| ---- | -------------| ------ |
|`mappability_min`| Genomic regions with mappability score less than this will be removed from callable sites. | `int`|
|`mappability_k`| Kmer size to compute mappability. | `int`|
|`mappability_merge`| Merge passing mappability regions separated by this or fewer basepairs into a signle region | `int`|
|`cov_merge`| Merge passing coverage regions separated by this or fewer basepairs into a signle region | `int`|

#### Coverage Filtering Options
If `cov_filter` in the core options is set to `True`, then the following options can be adjusted to the user's needs. Coverage filtering can be handled 3 ways:

1. Hard upper and lower thresholds: regions with a mean coverage that falls within these thresholds are considered callable.

| Option | Description | Type |
| ---- | -------------| ------ |
|`cov_threshold_lower`| Lower coverage threshold| `int`|
|`cov_threshold_upper`| Upper coverage threshold| `int`|

2. Standard deviations: regions with a mean coverage that is within N standard deviations (assumes Poisson distribution) are considered callable.

| Option | Description | Type |
| ---- | -------------| ------ |
|`cov_threshold_stdev`| Number of standard deviations is considered callable | `int`|

3. Absolute scaling: Thresholds set by factor N. Lower bowund is (global mean coverage / N), upper bound (global mean coverage * N). A region is callable if its mean coverage is within these bounds.

| Option | Description | Type |
| ---- | -------------| ------ |
|`cov_threshold_rel`| Scaling factor for coverage threshold| `int`|

```{note}
In order to use one of the above coverage filtering approaches, you must set the options of the desired approach, and leave the others blank.
```

#### Module Options
Please refer to the [modules page](./modules.md) for each module's options.

### Resources
Compute resources (threads, memory, etc.) as well as Snakemake arguments are set by the workflow profile located in `workflow-profiles/default/config.yaml`. This profile is used by default when running Snakemake. To specify a different profile, use the `--workflow-profile` option in your Snakemake command.

In the profile you can set resources to be applied to all rules via the `default-resources` key. You can override this default per-rule by uncommenting that rule under the `set-threads` and/or `set-resources` key.

```{note}
We recommend you use `--workflow-profile` to set resources for the workflow run. To set Snakemake options such as executor, retries, etc, use `--profile`. Read [here](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) for more info.
```

#### Threads
The profile controls how many `threads` (or CPU cores) a rule can use via the `set-thread` key. We have provided reasonable defaults, though you may need to adjust depending on the resources available on your system/cluster.
```{note}
Many rules can only use 1 thread, and providing more threads **will not** decrease runtime or improve performance. Please refer to the `profiles/default/config.yaml` for details.
``` 

#### Memory and other resources
The profile controls how much memory and what other resources a rule can use via the `set-resources` key. When executing snpArcher on a cluster or the cloud, specifying memory is important as these environments will typically kill jobs that use more memory than they requested. 

Other resources, such as `slurm_partition`, `runtime`, etc. can also be set here if they are required by your cluster. We have provided a SLURM profile `profiles/slurm` that has the most common SLURM resources.

```{note}
Snakemake allows you to dynamically assign resources. We use the `attempt` keyword to specify memory. For example. `attempt * 2000` will provide 2GB on the first attempt of the rule, if the rule fails (out of memory) then on the second attempt it will be provided 4GB. This behavior requires the `-T/--retries` Snakemake option.
```









