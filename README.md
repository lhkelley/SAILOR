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
    #--outFilterMismatchNoverLmax .3 \  # Set the maximum number of mismatches allowed per read, scaled by read length. I'm omitting it because Nmax supersedes it.
    --runMode alignReads \
    --genomeDir genome/index \
    --readFilesIn $FASTQ \
    --outFileNamePrefix $PREFIX \
    --outSAMattributes All \
    --outSAMtype BAM SortedByCoordinate
done
```

#### ```--outFilterMultimapNmax```:

This parameter sets the maximum number of loci that a read is allowed to map to. You can input a range as well.

If a read maps to =< N loci --> STAR keeps the read and reports all alignments (up to N)

If a read maps to > N loci --> STAR discards the read entirely (not reported in BAM)

```--outFilterMultimapNmax 1``` would only keep uniquely mapped reads and discards read that map to more than 1 locus. Typically, in RNA-seq experiments where you might want to look at differential expression, you would want to only keep unique reads.

```--outFilterMultimapNmax 10``` will keep reads that map to up to 10 loci. They will all be reported in the BAM file.

When **identifying edit sites** you do not want a read that maps to multiple locations, so set this parameter to 1.

#### ```--outFilterScoreMinOverLread```

```minScore = readLength × outFilterScoreMinOverLread```

```readLength``` = length of the read in bases

```outFilterScoreMinOverLread``` = fraction of read length required to be aligned well

The score that STAR assigns to each read is based on matches, mismatches, and gap. Matches are positive and mismatches/gaps are negative. A perfect alignment would be a score that equals the read's length. The length of the reads depends on your RNA-seq experimental design (AKA what you told the sequencing facility). Common read lengths for single end RNA-seq experiments are 75 and 100 bp.

If ```--outFilterScoreMinOverLread``` is set t0 0.66, that means if you have 100 bp reads, then ```100 x 0.66 = 66```, so a read must have a score of at least 66 to be kept, meaning two-thirds of the read align.

#### ```--outFilterMismatchNmax``` and ```---outFilterMismatchNoverLmax```

```Nmax``` sets the maximum number of mismatches allowed per read. This is a strict cap. Setting this to 10 means that you're allowing 10% of the read have mismatches.

For identifying and quantifying editing sites, it's important to include mismatches because **edits are mismatches to the genome.** Not sure how to empirically determine this parameter yet.

```NoverLmax``` sets the maximum number of mismatches by a relative cap, scaling with read length. This is a functional cap. If you have 100 bp reads and set this parameter to 0.3, then ```100 x 0.3 = 30``` so reads with 30 mismatches would be allowed through.

However, for short read lengths (i.e. 100 bp), if the ```Nmax``` parameter is smaller than  ```NoverLmax```, ```Nmax`` supersedes ```NoverLmax```.

### Checking the alignment

After running the STAR align, the BAMs and output files should be in ```results/align```.

```console
cd results/align
```
Check the number of uniquely aligned reads for each sample. Here is a command that will pull this information from each of the STAR output log files and will combine it into one table:

```console
awk 'BEGIN{OFS="\t"; print "Sample","UniquelyMappedReads","UniquelyMappedReadsPercent"}
FNR==1{sample=FILENAME; sub(/.*\//,"",sample); sub(/\.Log\.final\.out$/,"",sample)}
/Uniquely mapped reads number/ {num=$NF}
/Uniquely mapped reads %/ {pct=$NF; print sample,num,pct}' *_Log.final.out > GSF4254_STAR_uniquely_mapped_summary.txt
```
### SAILOR (through FLARE)

Make a directory for SAILOR:

```console
mkdir sailor
```
This should contain:
* Your ```bam``` files (merged or unmerged, see below)
* A ```bed``` file of known SNPs
* ```workflow_sailor``` directory containing scripts and Snakefile

Now you need to create a .json file that will be used to specify the parameters for the workflow (input files, reference files, etc.). Example below:

```console
{
  "samples_path":"/N/slate/lhkelley/GSF4254/sailor/",
  "samples": [
    "GSF4254-N2.merged.bam",  # List bam files here
  ],
  "reverse_stranded":true,  # This depends on your library prep and sequencing machine
  "reference_fasta": "/N/slate/lhkelley/GSF4254/genome/assembly.fasta",
  "known_snps": "/N/slate/lhkelley/GSF4254/sailor/c.elegans.WS275.snps.nostrand.sorted.bed",
  "edit_type": "AG",  # The type of edits you're looking for (A to I edits would be AG here)
  "output_dir": "/N/slate/lhkelley/GSF4254/results/fromSAILOR_N2results"
}
```
Then you run the command:

```console
snakemake \
  --snakefile /N/slate/lhkelley/GSF2848/sailor/workflow_sailor/Snakefile \
  --configfile /N/slate/lhkelley/GSF2848/sailor/run-sailor-GSF2848-adr2.json \
  --use-singularity \
  --singularity-args "\ 
    --bind /N/slate/lhkelley/GSF2848 \
    --bind /N/slate/lhkelley/GSF2848/genome \
    --bind /N/slate/lhkelley/GSF2848/sailor/workflow_sailor/scripts \
    --bind /N/slate/lhkelley/GSF2848/sailor/workflow_sailor \
    --bind /N/slate/lhkelley/GSF2848/ranSailor_adr2" \
  -j1
```
The ```singularity-args``` allows you to list extra arguments to pass to the Singularity command. The ```--bind``` mounts a host directory inside the Singularity container, which is necessary because the containers have their own filtesystem, so any input files, scripts, etc. must be accessible from inside the container.

The ```-j1``` parameter specifies to run one job at a time.

#### Brief overview of what the Snakefile is doing

* Indexes the input BAM files using ```samtools```
* Split the strands in the BAM files into forward and reverse
* Remove duplicate reads with ```samtools rmdup``` [need to change to ```samtools markdup```]
* Indexes the reads
* Generates bigWif coverage files
* Filters the reads based on edit type, junction overhang, etc.
* Creates pileups using ```samtools mpileup```
* Calls variants (SNVs) with ```bcftools call_snvs```
* Process and filters variant calls
* Ranks the RNA editing sites (??)
* Combines forward and reverse strands into same file
* Produces a bedgraph of editing fractions

The output will be in the directory that you specificed above (last of the ```-bind``` arguments) and will contain many things:

```console
1_split_strands
3_index_reads
4_filter_reads
5_pileup_reads
6_vcfs
7_scored_outputs
8_bw_and_bam
9_edit_fraction_bedgraphs
subsampled.bam.combined.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed
```
But the final out are the ```ranked.bed``` files in the ```7_scored_outputs``` directory. There will be separate files for the forward and reverse strands.

### SAILOR annotation (through FLARE)

Now that we know the 

