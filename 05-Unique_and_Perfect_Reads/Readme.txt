find_unique_perfect.py is a script used gather statistics about unique and perfectly mapped reads.

Perfectly mapped reads are the one that exactly match the reference (no mismatches; no indels). Uniquely mapped reads are the ones that have only one primary alignment (no secondary alignment at all) or the reads do have secondary alignment but the primary alignment score is significantly greater than the secondary alignment score

RUN  : samtools view <INPUT_1.bam> | python find_unique_perfect.py

Here <INPUT_1.bam> is the bam file for which you want to calculate the statistics.

NOTE : This bam file should be genarated by bowtie2 and must represent reads mapped in a single end (unpairred) mode.
