## Navigating to your supercomputer space

Access IU's supercomputer, Quartz, using the remote desktop server, ThinLinc, on your personal computer.

Once your on the remote desktop, open the terminal (shown below):

<br><br> <img width="946" height="480" alt="remote-desktop-terminal" src="https://github.com/user-attachments/assets/84790ec6-286a-45c2-b946-fbfb6819f0e3" /> <br><br>

Once at the terminal, you should see something like this:

<br><br> ![Terminal](https://github.com/user-attachments/assets/70abd7d8-e748-4f0b-bc0c-61423df3a14e) <br><br>

Once there, run the following commands to load preset modules:

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
Note: Use ```conda list``` to check your version of sra-tools and change if needed in above line.

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

Note: I recommend doing all of your work in one parent directory. For example, you can keep all of the FASTQ files in the ```For_Josh``` directory (or can also rename it), but keep that directory within a parent directory that you'll do all of the downstream steps in. You could name it something like ```rnaseq_analysis```.

## Sanity checks with the FASTQ files

Each FASTQ file is a replicate and there are three replicates each for wildtype (strain is called N2) and the _adr-2(-)_ mutant.

You can see how many reads are in each FASTQ file with the following command (replace the file name with the file you want to count the reads in):

```
wc -l yourfile.fastq
```

## Align reads

Next, you'll want to align the reads in the FASTQ files to the _C. elegans_ genome. But, first you'll need the reference genome for _C. elegans_, which you can download from WormBase.

Copy and paste this link into a browser to download the reference genome:
https://downloads.wormbase.org/releases/WS275/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS275.genomic.fa.gz

Move the file to your Documents folder on your laptop (this is for Windows, if you have a Mac, we'll have to talk to someone else in the lab about how to do this step).

Now, on your remote desktop, go to the ThinDrives folder (shown below) and then into the MyDocuments folder. This should show all of the files/folders in your laptop's Document folder.

<br><br> <img width="935" height="478" alt="remote-desktop-thiclincdrive" src="https://github.com/user-attachments/assets/3193f7e1-db27-4422-897c-c67173bf0f29" /> <br><br>

Move the reference genome file to your remote desktop Desktop (drag and drop).

Then, in the terminal, run this command to transfer the reference genome from the Desktop to your Slate space:

```
scp ~/Desktop/c_elegans.PRJNA13758.WS275.genomic.fa.gz /N/slate/jwt10/rnaseq_analysis
```
Repeat this process for the _C. elegans_ genome annotation file, found here:
https://downloads.wormbase.org/releases/WS275/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS275.canonical_geneset.gtf.gz

Make sure to unzip the reference files:

```
gunzip c_elegans.PRJNA13758.WS275.genomic.fa.gz
gunzip c_elegans.PRJNA13758.WS275.canonical_geneset.gtf.gz
```

Once your have both the reference genome and annotation files on your Slate workspace, you can move on to mapping the reads from the FASTQ files.

## Mapping reads using STAR

First, you'll have to create a genomic index. This only uses the reference genome and annotation files (not the RNA-seq data).

Create genomic index:

```console
STAR \
  --runThreadN 4 \  # Number of CPU threads used in parallel
  --runMode genomeGenerate \  # Buuld a genome index (as opposed to aligning reads)
  --genomeDir genome_index \  # Directory where the genome index will be stored (create this before performing this command)
  --genomeFastaFiles c_elegans.PRJNA13758.WS275.genomic.fa \  # Reference genome FASTA file that you provide
  --sjdbGTFfile c_elegans.PRJNA13758.WS275.canonical_geneset.gtf \  # Annotation file that you proivde; helps STAR handle exon-exon junctions
  --genomeSAindexNbases 12  # Length of the suffix array index; affects memory usage and speed
```
This will create various index files that STAR will use during the actual mapping process. The files will be stored in the ```genome_index``` directory that you specificed in the above code.

Now you can align your reads (in the FASTQ files) to the genome.

Align reads to genome:

```console

FASTQ=$GSF4254*

for FASTQ in ${FASTQ[@]}; do
  PREFIX=results/aligned/$(basename $FASTQ .fastq)_
  STAR \
    --runThreadN 8 \  # Uses 8 threads in parallel
    --outFilterMultimapNmax 1 \  # Set the maximum number of loci that a read is allowed to map to
    --outFilterScoreMinOverLread .66 \  # Set the minimum alignment score, scaled by read length
    --outFilterMismatchNmax 10 \  # Set the maximum number of mismatches allowed per read
    --runMode alignReads \
    --genomeDir genome_index \
    --readFilesIn $FASTQ \
    --outFileNamePrefix $PREFIX \
    --outSAMattributes All \
    --outSAMtype BAM SortedByCoordinate
done
```
You should determine how many reads are in each FASTQ file because this information is important for downstream steps like mapping (you'll want to know how many reads mapped out of the total number of reads you started with, etc.).
