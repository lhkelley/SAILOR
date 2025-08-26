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
    --outFilterMultimapNmax 1 \  # Set the maximum number of loci that a read is allowed to map to
    --outFilterScoreMinOverLread .66 \  # Set the minimum alignment score, scaled by read length
    --outFilterMismatchNmax 10 \  # Set the maximum number of mismatches allowed per read
    --outFilterMismatchNoverLmax .3 \
    --runMode alignReads \
    --genomeDir genome/index \
    --readFilesIn $FASTQ \
    --outFileNamePrefix $PREFIX \
    --outSAMattributes All \
    --outSAMtype BAM SortedByCoordinate
done
```

####```--outFilterMultimapNmax```:

This parameter sets the maximum number of loci that a read is allowed to map to. You can input a range as well.

If a read maps to =< N loci --> STAR keeps the read and reports all alignments (up to N)

If a read maps to > N loci --> STAR discards the read entirely (not reported in BAM)

```--outFilterMultimapNmax 1``` would only keep uniquely mapped reads and discards read that map to more than 1 locus. Typically, in RNA-seq experiments where you might want to look at differential expression, you would want to only keep unique reads.

```--outFilterMultimapNmax 10``` will keep reads that map to up to 10 loci. They will all be reported in the BAM file.

When **identifying edit sites** you do not want a read that maps to multiple locations, so set this parameter to 1.

####```--outFilterScoreMinOverLread```

```minScore = readLength × outFilterScoreMinOverLread```

```readLength``` = length of the read in bases
```outFilterScoreMinOverLread``` = fraction of read length required to be aligned well

The score that STAR assigns to each read is based on matches, mismatches, and gap. Matches are positive and mismatches/gaps are negative. A perfect alignment would be a score that equals the read's length. The length of the reads depends on your RNA-seq experimental design (AKA what you told the sequencing facility). Common read lengths for single end RNA-seq experiments are 75 and 100 bp.

If ```--outFilterScoreMinOverLread``` is set t0 0.66, that means if you have 100 bp reads, then ```100 x 0.66 = 66```, so a read must have a score of at least 66 to be kept, meaning two-thirds of the read align.

####```--outFilterMismatchNmax```

This sets the maximum number of mismatches allowed per read. Setting this to 10 means that you're allowing 10% of the read have mismatches.

For identifying and quantifying editing sites, it's important to include mismatches because **edits are mismatches to the genome.** Not sure how to empirically determine this parameter yet.
