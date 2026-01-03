# Examples
On this page you will find an example project scenario and how to setup and run it using snpArcher.

In this example, we have 10 resequenced individuals we would like to generate variant calls for. We will cover creating the sample sheet, selecting config options, and running the workflow.
## Directory structure
First, let's setup our directories as suggested in our [executing](./executing.md#optional-directory-setup) instructions. Let's assume we are working in a directory called `workdir/`, and the snpArcher repository has already been cloned there. We have also already created the `snparcher` conda env as instructed in the [setup docs](./setup.md#environment-setup).

1. Let's create a directory to organize this project and future ones, call it `projects`. Then, create a new directory for this project, we'll call it `secretarybird_reseq`. 
```
.
├── projects
│   └── secretarybird_reseq
└── snpArcher
    └── ...
```
```{note}
Not all files and directories are shown, only relevant ones. 
```
2. Copy the snpArcher config directory `snpArcher/config` to `projects/secretarybird_reseq`:
```
.
├── projects
│   └── secretarybird_reseq
│       └── config
│           └── config.yaml
└── snpArcher
    └── ...
```

3. Assume we already have all our sequence data and reference genome on our system, stored in a different location `/storage/data`. We do not need to move the raw data to our project directory. 
```{note}
We'll cover the cases using SRA data and refSeq genomes later on in this example.
```
## Sample sheet setup
Now we need to setup our sample sheet to inform snpArcher of our samples and their metadata. You can use any editor to create the sheet, as long as it is a CSV file. We will save the sample sheet in our project's config directory: `projects/secretarybird_reseq/samples.csv`. Below is the final sample sheet that we will use going forward, with explanations of each column following.

For a more comprehensive explanation of the sample sheet, please refer to [here](./setup.md#creating-a-sample-sheet) for more details.


### Final sample sheet
```
BioSample,LibraryName,Run,fq1,fq2,lat,long
bird_1,bird_1_lib,1,/storage/data/bird_1_R1.fq.gz,/storage/data/bird_1_R2.fq.gz,-8.758119,-36.280061
bird_2,bird_2_lib,2,/storage/data/bird_2_R1.fq.gz,/storage/data/bird_2_R2.fq.gz,-72.336165,35.751903
bird_3,bird_3_lib,3,/storage/data/bird_3_R1.fq.gz,/storage/data/bird_3_R2.fq.gz,-11.874137,-5.382251
bird_4,bird_4_lib,4,/storage/data/bird_4_R1.fq.gz,/storage/data/bird_4_R2.fq.gz,-73.235723,-145.261219
bird_5,bird_5_lib,5,/storage/data/bird_5_R1.fq.gz,/storage/data/bird_5_R2.fq.gz,88.08701,-52.658705
bird_6,bird_6_lib,6,/storage/data/bird_6_R1.fq.gz,/storage/data/bird_6_R2.fq.gz,69.640536,-12.971862
bird_7,bird_7_lib,7,/storage/data/bird_7_R1.fq.gz,/storage/data/bird_7_R2.fq.gz,18.608941,-100.485774
bird_8,bird_8_lib,8,/storage/data/bird_8_R1.fq.gz,/storage/data/bird_8_R2.fq.gz,-36.570632,-102.38721
bird_9,bird_9_lib,9,/storage/data/bird_9_R1.fq.gz,/storage/data/bird_9_R2.fq.gz,-88.592265,157.406505
bird_10,bird_10_lib,10,/storage/data/bird_10_R1.fq.gz,/storage/data/bird_10_R2.fq.gz,40.106437,-58.649016
```
### Description of Columns
1. **BioSample**: This is the name for the sample.
2. **LibraryName**: Identifier for the sample's sequencing library. This is especially important if you have samples that were sequenced multiple times across multiple lanes, which is not the case in this example. See [here](./setup.md#handling-samples-with-more-than-one-pair-of-reads) for more details.
3. **Run**: If we were using reads from the SRA, this is where the sample's SRR accession would go. However, since we have local data, this just has to be a unique value.
4. **fq1**: Path to the first read pair. Absolute paths are recommended. If we were using SRA data, this column should be omitted.
5. **fq1**: Path to the second read pair. Same note as fq1.
6. **lat**: Decimal latitude for the sample, used to generate map in QC module output.
6. **long**: Decimal longitude for the sample, used to generate map in QC module output. 

```{note}
If your project has multiple genomes, you can add the refPath and refGenome columns.
```

## Config file setup
Now that we've created our sample sheet, we need to edit the config file we copied earlier: `projects/secretarybird_reseq/config.yaml`. This file controls the main options for controlling snpArcher's outputs. Refer to the [setup section](./setup.md#configuring-snparcher) for more details. 

In our example we are using all of the default options. This will configure snpArcher to perform variant calling using GATK with the scatter-by-intervals approach. Also, we have set our reference genome name and path since we want to use the same genome for all samples in our sample sheet.

```
samples: "config/samples.csv" # path to the sample metadata CSV
final_prefix: "" # prefix for final output files
intervals: True #Set to True if you want to perform variant calling using interval approach.
sentieon: False #set to True if you want to use sentieon, False if you want GATK
sentieon_lic: "" #set to path of sentieon license
remote_reads: False # Set True if reads are in a location seperate from --default-remote-prefix.
bigtmp: "" #Set to a path with lots of free space to use for commands that require large amounts of temp space; defaults to system tmpdir if empty
cov_filter: True #set to True if you want to include coverage thresholds in the callable sites bed file (default uses mappability only)
generate_trackhub: True #Set to true if you want to generate a Genome Browser Trackhub. Dependent on postprocessing module.
trackhub_email: "hi@email.com"
##############################
# Variables you *might* need to change
##############################

# Set reference genome here if you would like to you use the same reference genome for all samples in sample sheet. See docs for more info.
refGenome: "bird_genome" # Name for reference genome
refPath: "/storage/data/bird.fa.gz"
```

## Profile setup
Snakemake uses profile YAML files to specify commonly used command line arguments, so you don't have to remember all of the arguments you need. Read more about profiles [here](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles). To specify a profile, you can use the `--workflow-profile` option when running Snakemake.

```
cp -r snpArcher/workflow-profiles projects/secretarybird_reseq
```

The profile also enables you to specify the compute resources any of snpArcher's rules can use. This is done via the YAML keys `default-resources`, `set-resources`, and `set-threads`. `default-resources` will apply to all rules, and `set-resources` can be applied to indiviudal rules, overriding what the default was set to. There is no way to set a default thread value. 

First, we will specify how many threads each rule can use. This is the same using the default or SLURM profile. Both profiles come with reasonable default thread values, but you may need to adjust based on your system or cluster. 

Let's say we wanted the alignment step (bwa mem) to use more threads:
```
# ...
set-threads:
  bwa_map: 16 # Changed from 8 to 16.
# ...
```
Next, we will specify memory and other resources. This step only applies if you are running on a SLURM cluster.

In our example cluster, we have two compute partitions, "short" and "long". So we want to put long running jobs on the "long" partition, and the rest on "short". Additionally, the "short" partition has a timelimit of 1 hour and "long" 10 hours, so we will specify that. 

First, lets specify the default resources:
```
default-resources:
  mem_mb: attempt * 2000
  mem_mb_reduced: (attempt * 2000) * 0.9 # Mem allocated to java for GATK rules (tries to prevent OOM errors)
  slurm_partition: "short" # This line was changed
  slurm_account: # Same as sbatch -A. Not all clusters use this.
  runtime: 60 # In minutes 
```
Then, lets modify the specific resources for the GATK HaplotypeCaller step:
```
set-resources:
# ... other rules
   bam2gvcf: # HaplotypeCaller <--- This line was uncommented
#     mem_mb: attempt * 2000
#     mem_mb_reduced: (attempt * 2000) * 0.9 # Mem allocated to java (tries to prevent OOM errors)
     slurm_partition: "long" # This line was changed
     runtime: 600 # This line was changed
```

## Running the workflow
We are now ready to run the workflow! From our working directory we can run the command:
```
snakemake -s snpArcher/workflow/Snakefile -d projects/secretarybird_reseq --workflow-profile projects/secretarybird_reseq/workflow-profiles/default
```
This instructs Snakemake to use snpArcher's workflow file, and to run in the project directory we setup using the config and sample sheet we setup there.

If we were on a SLURM cluster, we would add `--executor slurm` to our command:
```
snakemake --executor slurm -s snpArcher/workflow/Snakefile -d projects/secretarybird_reseq --workflow-profile projects/secretarybird_reseq/workflow-profiles/default
```