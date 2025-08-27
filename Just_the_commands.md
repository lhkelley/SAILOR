# Identifying edit sites via SAILOR (through FLARE)

## Aliging reads to the genome (STAR)

Create index:

```console
STAR \
  --runThreadN 4 \  # Number of CPU threads used in parallel
  --runMode genomeGenerate \  # Buuld a genome index (as opposed to aligning reads)
  --genomeDir genome/index \  # Directory where the genome index will be stored (create this before performing this command)
  --genomeFastaFiles genome/assembly.fasta \  # Reference genome FASTA file that you provide
  --sjdbGTFfile genome/annotation.gtf \  # Annotation file that you proivde; helps STAR handle exon-exon junctions
  --genomeSAindexNbases 12  # Length of the suffix array index; affects memory usage and speed
```

Align the reads:

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

## Identify variant sites with SAILOR

```.json``` file:

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

Snakemake command:

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
