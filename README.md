# SAILOR
General use of SAILOR

## Navigating to your supercomputer space

Access IU's supercomputer, Quartz, but using the remote desktop server, ThinLinc, on your personal computer.

[details]

Once at the terminal:

```console
module load miniconda
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

Note: #Use conda list to check your version of sra-tools and change if needed in above line.

Move into your Slate directory and make a directory for the dataset that you will be analyzing

cd /N/slate/lhkelley
mkdir GSF2848

## Analyzing GSF4254 (*glh-1* and *glh-2* short read RNA-seq)

```console
cd /N/slate/lhkelley/GSF4254/
conda activate rnaseq
```
### Create an index with STAR

Create a directory to store the index:

```console
mkdir -p genome/index
```
Create the STAR genome index:

Notes: 
--genomeSAindexNbases 12 was recommended by software.

```console
STAR \
  --runThreadN 4 \  # Number of CPU threads used in parallel
  --runMode genomeGenerate \  # Buuld a genome index (as opposed to aligning reads)
  --genomeDir genome/index \  # Directory where the genome index will be stored (create this before performing this command)
  --genomeFastaFiles genome/assembly.fasta \  # Reference genome FASTA file that you provide
  --sjdbGTFfile genome/annotation.gtf \  # Annotation file that you proivde; helps STAR handle exon-exon junctions
  --genomeSAindexNbases 12  # Length of the suffix array index; affects memory usage and speed
```
Additional details on ```--genomeSAindexNbases```:

The suffix array is like a dictionary, so that STAR can quickly look up reads. However, instead of searching the entire genome for the entire legnth of a read, the SA index allows STAR to look at the first *k* bases of the read.

If *k* is too small = the k-mers won't be unique enough and STAR will have to scan more of the entries in the SA (slow mapping)

If *k* is too large = the index becomes large, requiring more memory

This is the calculation: ```genomeSAindexNbases (k) ≈ log2​(genomeLength)/2 − 1```

The *C. elegans* genome is roughly 100 Mb, so using 12 is appropriate.

Outputs will be in genome/index.

### Align the reads

Still in :

```console
cd /N/slate/lhkelley/GSF4254/
conda activate rnaseq
```

Create a variable so that all files that start with GSF4254 (should only be FASTQ files) will be used as input in the STAR mapping command:

```console
FASTQ=$GSF2848*
```



```console
for FASTQ in ${FASTQ[@]}; do
  PREFIX=results/aligned/$(basename $FASTQ .fastq)_
  STAR \
    --runThreadN 8 \  # Uses 8 threads in parallel
    --outFilterMultimapNmax 1 \  # Set the maximum number of loci that a read is allowed to map to (how STAR should handle multi-mapping reads)
    --outFilterScoreMinOverLread .66 \
    --outFilterMismatchNmax 10 \
    --outFilterMismatchNoverLmax .3 \
    --runMode alignReads \
    --genomeDir genome/index \
    --readFilesIn $FASTQ \
    --outFileNamePrefix $PREFIX \
    --outSAMattributes All \
    --outSAMtype BAM SortedByCoordinate
done
```

```--outFilterMultimapNmax```:
This parameter sets the maximum number of loci that a read is allowed to map to. You can input a range as well.

If a read maps to =< N loci --> STAR keeps the read and reports all alignments (up to N)

If a read maps to > N loci --> STAR discards the read entirely (not reported in BAM)

```--outFilterMultimapNmax 1``` would only keep uniquely mapped reads and discards read that map to more than 1 locus. Typically, in RNA-seq experiments where you might want to look at differential expression, you would want to only keep unique reads.

```--outFilterMultimapNmax 10``` will keep reads that map to up to 10 loci. They will all be reported in the BAM file.
