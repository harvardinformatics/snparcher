# Running snpArcher
## Setup
Please refer to our [setup instructions](./setup.md) to prepare the snpArcher environment and requisite files. 
## Test datasets
To test that your environment is properly setup, you can run a quick test with the following command:
```
snakemake -d .test/ecoli --cores 1 --use-conda --workflow-profile workflow-profiles/default
```
If this runs without errors, you are ready to go!
## Using the Dry-run option
Snakemake offers the `--dry-run (-n)` CLI option to perform a dry-run of the workflow to show what jobs would be run. We recommend doing this before executing snpArcher to ensure that the sample sheet was setup correctly, and Snakemake has correctly built the workflow DAG.
## Local Execution
Once you have setup the requisite configuration files and sample sheet, executing snpArcher on your local machine is as simple as running the Snakemake command with the number of cores you would like to use. For example, to use 8 cores you would run:
```
snakemake --cores 8 --use-conda --workflow-profile workflow-profiles/default
```

### Optional directory setup 
To maintain organization across many different projects, you may consider creating a new directory for each project you run using snpArcher. This way, each of your project directories will contain the configuration files used for that run. Below is an example directory structure:

```
.
├── snpArcher
├── project_1/
│   ├── config/
│   │   ├── config.yaml
│   │   └── samples.csv
│   ├── data
│   └── results
└── project_2/
    ├── config/
    │   ├── config.yaml
    │   └── samples.csv
    └── data
```

When creating a new directory for an analysis, ensure that you copy the `config` directory from the snpArcher directory to your new directory.

Then, to run snpArcher on `project_2` from our example, we would execute the command:
```
snakemake -s ./snpArcher/workflow/Snakefile -d ./project_2 <other CLI options>
```

## Cluster Execution
Snakemake [supports most cluster schedulers](https://snakemake.github.io/snakemake-plugin-catalog/) via executor plugins. Here, we provide documentation for SLURM, however please refer to Snakemake's documentation for further details on using other plugins.
### SLURM
#### Install plugin
To execute snpArcher on a SLURM cluster, you will need to install the [SLURM executor plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html) into the snpArcher environment.
```shell
conda activate snparcher
pip install snakemake-executor-plugin-slurm
```
#### Profile Setup
To specify resources for the workflow to SLURM, you must use a workflow profile. We have provided a SLURM profile template (`workflow-profiles/slurm/config.yaml`) which you can modify to specify SLURM partitions, memory allocation, etc. Please refer to the [profiles setup section](./setup.md#resources) for more details. Also see the [Snakemake documentation on profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) for more info.

An example SLURM profile specifies required and recommended Snakemake options:
```yaml
executor: slurm
use-conda: True
jobs: 100 # Have up to N jobs submitted at any given time
latency-wait: 20 # Wait N seconds for output files due to latency
retries: 3 # Retry jobs N times.
```

#### Running the workflow
Once you have modified the SLURM profile appropriately, you can run snpArcher with the following command:
```shell
snakemake --workflow-profile <path/to/profile-directory> <other options>
```
Depending on your cluster, you can run this command on the head node and Snakemake will submit jobs to the SLURM queue. You can also submit this command via `srun` or `sbatch`.

## Cloud Execution
```{warning}
Google Lifesciences execution is not supported by Snakemake versions >8. Please use Snakemake v7.32.4 if you would like to use this execution mode.
```
Like cluster execution, Snakemake [supports a number of cloud providers](https://snakemake.readthedocs.io/en/stable/executing/cloud.html). Here we provide documentation for executing using Snakemake's Google Lifesciences integration. Please refer to Snakemake's documentation for details on using other cloud providers.
### Google Lifesciences
Snakemake's integration with the Google Lifesciences (GLS) API allows you to easily run snpArcher on the Google Cloud Platform (GCP). Using this execution mode allows you to take advantage of hundreds or thousands of GCP virtual machine instances. Snakemake manages deploying instances, running jobs, and deleting instances with finished jobs. Furthermore, you can use preemptible instances which are offered at a large cost discount, but can only run for a maximum of 24 hours. 

We include profiles for GLS using preemptible instances so that you can get up and running quickly.
#### Google Credential Setup
In order to use the Google Lifesciences execution option, you must first setup your Google Cloud credentials. Please refer [here](https://snakemake.readthedocs.io/en/stable/executor_tutorial/google_lifesciences.html#credentials) for full details.
#### Data setup
To use this execution mode, you must have a Google Storage bucket with your raw data files. This can be achieved by using Google's web interface, or at the command line using [`gsutil`](https://cloud.google.com/storage/docs/gsutil). For example, if we have some data locally like so:
```
.
└── data/
    ├── raw_reads/
    │   ├── samp_1_R1.fq.gz
    │   ├── samp_1_R2.fq.gz
    │   ├── samp_2_R1.fq.gz
    │   ├── samp_2_R2.fq.gz
    │   ├── samp_3_R1.fq.gz
    │   ├── samp_3_R2.fq.gz
    │   ├── samp_4_R1.fq.gz
    │   ├── samp_4_R2.fq.gz
    │   └── ...
    └── ref_genome/
        └── genome.fa
```

We can copy this data to our bucket like so:
```
gsutil cp -r ./data gs://<bucket-name>
```

```{note}
When using cloud execution, do not include the bucket name in any path fields of the sample sheet, such as fq1, fq2, or refPath
```
```{note}
If you are using data hosted on NCBI, you do not need to upload those data to your bucket, snpArcher will handle this for you. However, you still need to create storage bucket to be used for the workflow.
```

Some users may want to store their raw reads in a separate bucket from where the workflow will store files. To do so, you can specify the remote prefix in `config/config.yaml`.

#### Running the workflow
Once your credentials and data are setup, you can run snpArcher using the included profile `profiles/google_lifesciences/config.yaml`. You should modify this profile to set resources for the workflow such as instance type, threads, preemptible rules, etc.  

To run the workflow, execute the following command:
```
snakemake --workflow-profile <GLS profile> --default-remote-prefix <bucket name>
```

As the workflow runs, Snakemake will print out logging information to the terminal. Please refer [here](https://snakemake.readthedocs.io/en/stable/executor_tutorial/google_lifesciences.html#step-5-debugging) for further details.
