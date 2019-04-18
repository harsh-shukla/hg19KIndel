
In order to build annotation for hg19Kindel reGene.txt.gz for hg19 was downloaded from UCSC table browser.

1)

This files was given as input to convert_cordinates.py which lifted over coordinates of transcripts to hg19Kindel reference.This code is very similar to /variants_liftover/liftover_vcf_cordinates.py except for few minor tweaks.

 RUN    : python convert_cordinates.py <hg19_refGene.txt> <INDELs_Replaced_1000Genome.vcf.gz> <hg19K.chr_length.txt> <hg19Kindel.chr_length.txt>

 OUTPUT : hg19Kindel_genes.refSeq - intermediary files representing  transcripts in hg19Kindel
          UnMapped.refseq         - list of transcripts that could,nt be converted

2)  

As mentioned in the paper an additional step of error correction had to be done in order make the frame information consistent with the changes we made to the reference. This is done by the script error_correct.py

 RUN    : python error_correct.py <hg19Kindel_genes.refSeq>

 OUTPUT : hg19Kindel_refGene.txt   - best possible refGene file representing the transcripts in hg19Kindel 

