# Comparison of ```samtools``` packages ```rmdup``` and ```markdup ```

```rmdup``` identifies PCR duplicates from pair-end reads. It identifies "pairs of reads where multiple reads align to the same exact start position in the genome, and the reverse read on the 3′ end maps at the exact same location (i.e. external mapping coordinates are identical)" (Ebbert et al. 2016). This program takes the read pairs with the highest mapping quality score and removes the other read pairs (thus, removing data completely).

For single-end reads, it identifies duplicates based on identical external 5' coordinates.

```markdup``` identifies "read pairs with the same orientation that have the exact same 5′ start position in the mapping. It takes into account clipping on the 5′ end of the read and makes calculations based on where the 5′ start position would be if the entire read had mapped to the reference" (Ebbert et al. 2016). This program does not remove reads, but flags the duplicates with the SAM flag 1024.

For single-end reads, it identifies duplicates based on the 5' mapping position and oreintation and uses soft-clipping to determine the true 5' position. So, if the following are the same for a group of reads, they are considered duplicates: the 5' end of the reads are the same, they have the same orientation (forward vs reverse), and they map to the same chromosome.

