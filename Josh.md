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
## Sanity checks with the FASTQ files

Each FASTQ file is a replicate and there are three replicates each for wildtype (strain is called N2) and the _adr-2(-)_ mutant.

You can see how many reads are in each FASTQ file with the following command (replace the file name with the file you want to count the reads in):

```
wc -l yourfile.fastq
```

## Align reads

Next, you'll want to align the reads in the FASTQ files to the _C. elegans_ genome. But, first you'll need the reference genome for _C. elegans_, which you can download from WormBase.

Set a link to the reference genome at WormBase:

```
ASSEMBLY='https://downloads.wormbase.org/releases/WS275/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS275.genomic.fa.gz'
```
Now use the curl command to download the file, followed by moving the file to 

Create genomic index:

```console
STAR \
  --runThreadN 4 \  # Number of CPU threads used in parallel
  --runMode genomeGenerate \  # Buuld a genome index (as opposed to aligning reads)
  --genomeDir genome/index \  # Directory where the genome index will be stored (create this before performing this command)
  --genomeFastaFiles genome/assembly.fasta \  # Reference genome FASTA file that you provide
  --sjdbGTFfile genome/annotation.gtf \  # Annotation file that you proivde; helps STAR handle exon-exon junctions
  --genomeSAindexNbases 12  # Length of the suffix array index; affects memory usage and speed
```

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
    #--outFilterMismatchNoverLmax .3 \  # Set the maximum number of mismatches allowed per read, scaled by read length. I'm omitting it because Nmax supersedes it.
    --runMode alignReads \
    --genomeDir genome/index \
    --readFilesIn $FASTQ \
    --outFileNamePrefix $PREFIX \
    --outSAMattributes All \
    --outSAMtype BAM SortedByCoordinate
done
```
You should determine how many reads are in each FASTQ file because this information is important for downstream steps like mapping (you'll want to know how many reads mapped out of the total number of reads you started with, etc.).
