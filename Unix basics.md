# SAILOR
General use of SAILOR

## Navigating to your supercomputer space

Access IU's supercomputer, Quartz, using the remote desktop server, ThinLinc, on your personal computer.

Once your on the remote desktop, open the terminal (shown below):

<img width="946" height="480" alt="remote-desktop-terminal" src="https://github.com/user-attachments/assets/84790ec6-286a-45c2-b946-fbfb6819f0e3" />

Once at the terminal, run the following commands to load preset modules:

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
In my terminal, I'll run the command:

```
scp For_Josh.tar.gz jwt10@quartz.uits.iu.edu:/N/slate/jwt10/
```
And you'll have to type in your password on my laptop for the transfer to go through. Then you'll be able to find the directory as ```/N/slate/jwt10```.

## Working with FASTQ files

To decompress the directory containing the FASTQ files, run:

```
tar -xzvf For_Josh.tar.gz
```

Now if you move into the directory and list out the contents, you should be able to see the 6 FASTQ files:

```
cd For_Josh
ls
```
