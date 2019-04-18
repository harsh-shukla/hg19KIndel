
1) extract_major_Allele.py

This python script takes a vcf file as input and extracts the variants having alternate allele(AA) frequency greater than 0.5.In case of multi-allelic sites the alternate allele having the frequency > 0.5 is selectively extracted and ALT field of the vcf file is updated accordingly. The resulting vcf file generated has only one ALT allele(the major allele) for a given REF allele

2) replace_variants.py

This python script takes a fasta file and a vcf file as input and replaces for all the records (variants), the REF with the ALT field to give a modified/updated fasta file.

