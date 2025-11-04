# SAILOR
General use of SAILOR

## Navigating to your supercomputer space

Access IU's supercomputer, Quartz, using the remote desktop server, ThinLinc, on your personal computer.

Once your on the remote desktop, open the terminal (shown below):



Once at the terminal:

```console
module load conda
module load apptainer
```
To check that miniconda is loaded, run:

```console
module list
```

## Installation and getting set up

You only need to run these commands once:

```console
conda install -c bioconda sra-tools
```
Now we'll create a conda environment and load it with the programs/packages we'll need:

```console
conda create -n rnaseq -c conda-forge -c bioconda entrez-direct sra-tools==3.1.1 fastqc seqtk star samtools subread
conda activate rnaseq
```
Note: Use conda list to check your version of sra-tools and change if needed in above line.

## Transferring RNA-seq files to your working space

We'll use the ```scp``` command to move the FASTQ files from my workspace to your workspace.
In my terminal, I'll run the command
