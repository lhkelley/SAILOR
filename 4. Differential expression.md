# Differential expression

Count features:
```
BAMS=$(find ./results/aligned -name "*\.bam")

featureCounts \
  -a genome/annotation.gtf \
  -o results/counts/counts_16Sept25.tsv \
  -t gene \
  -g gene_id \
  --largestOverlap \
  --readExtension3 150 \
  --primary \
  -s 2 \
  -T 8 \
  ${BAMS}
  
done
```
