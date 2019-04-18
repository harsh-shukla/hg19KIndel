import re
import sys

############################################################## USAGE ###########################################################
#
# Run   : python generate_chain_file_hg19Kindeltohg19.py <INPUT_VCF_FILE.vcf> <INDELs_Replaced_1000Genome.vcf.gz> <hg19K.chr_length.txt> <hg19Kindel.chr_length.txt> 
#   
# INPUT :
#  <INPUT_VCF_FILE.vcf>                - variants called using hg19Kindel as reference  
#  <INDELs_Replaced_1000Genome.vcf.gz> - the vcf containing all the indels we replaced in hg19 to create hg19Kindel
#  <hg19K.chr_length.txt>              - tab seperated file of chromosomes and their correspoding lengths in hg19K/hg19 reference
#  <hg19Kindel.chr_length.txt>         - tab seperated file of chromosomes and their correspoding lengths in hg19Kindel reference
#
# Output :  1) INPUT_VCF_FILE.vcf.CONVERTED.vcf  - all the variants with positions lifted to hg19 coordinate frame
#           2) INPUT_VCF_FILE.vcf.UNMAPPED.vcf   - the variants that could not be mapped. (That particualr position is deleted in hg19)
#
#  DESCRIPTION :
#       This python script liftsover coordiantes in an input vcf called using hg19Kindel as reference to hg19 coordinate frame  

#  NOTE        :  All the above input files (barring the input vcf) are present in Supplementary_Files folders in the current repository
##################################################################################################################################



## Input ##


f_vcf = open(sys.argv[1],"r")

hg19K_chrom_size_file      =  open(sys.argv[3],"r")
hg19Kindel_chrom_size_file =  open(sys.argv[4],"r")


##  Output  ###

coll = re.split("/",sys.argv[1])

fname = coll[-1]


filename = fname + ".CONVERTED.vcf"
fout = open(filename,"w")

filename2 = fname + ".UNMAPPED.vcf"
funmapped = open(filename2,"w")

hg19K_chromSize = {}
hg19Kindel_chromSize={}


for line in hg19K_chrom_size_file:

    cols = re.split("\t|\n",line)

    chr_curr = cols[0]
    chr_size = int(cols[1])  
  
    hg19K_chromSize[chr_curr] = chr_size


for line in hg19Kindel_chrom_size_file:

    cols = re.split("\t|\n",line)

    chr_curr = cols[0]
    chr_size = int(cols[1])  
  
    hg19Kindel_chromSize[chr_curr] = chr_size


hg19K_chrom_size_file.close()
hg19Kindel_chrom_size_file.close()

include_chr=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"]

#include_chr=["chr1"]

# cannot gurantee the converetd cordinates when looked at the sequence will have the excat sequence in hg19 especially in case of indels also in case of SNPS

current_chr=""
map_current={}


for line in f_vcf:

    if line[0:3]!="chr":
        fout.write(line)
        continue


    cols = re.split("\t|\n",line)
    chrm_refSeq = cols[0]

    if chrm_refSeq not in include_chr:        
        fout.write(line)
        continue  

    #print cols

    if chrm_refSeq!=current_chr:

        del map_current
        map_current={}

        ## generate the mapping dictionary for every position in 

        f_variant_vcf = open(sys.argv[2],"r")

        map_k=1
        map_kindel=1
        displaced = 0
        chromosome_size_hg19       = hg19K_chromSize[chrm_refSeq]
        chromosome_size_hg19Kindel = hg19Kindel_chromSize[chrm_refSeq]

        #print chromosome_size_hg19,chromosome_size_hg19Kindel

        for line in f_variant_vcf:

            if line[0:3]!="chr":   # process the header records write them as they are
                continue                 

            col1 = re.split("\t|\n",line)
            chr_vcf    = col1[0]
       
            if chr_vcf!=chrm_refSeq:
                continue


            old_pos = int(col1[1]) # No -1 since all is vcf 1 based cordinate
            ref     = col1[3]
            alt     = col1[4]

            len_ref = len(ref)
            len_alt = len(alt)


            while (map_k <= old_pos):

                map_current[map_kindel] = map_k

                map_k       = map_k + 1
                map_kindel  = map_kindel + 1
                                         

            if len_ref > len_alt:
                
                counter=0
                while(counter < (len_alt-1)):
                    map_current[map_kindel] = map_k 
                    map_k       = map_k + 1
                    map_kindel  = map_kindel + 1
                    counter = counter + 1
  
                map_k = map_k + (len_ref - len_alt)

            if len_alt > len_ref:

                counter=0
                while(counter < (len_ref-1)):
                    map_current[map_kindel] = map_k 
                    map_k       = map_k + 1
                    map_kindel  = map_kindel + 1
                    counter = counter + 1

                map_kindel = map_kindel + (len_alt - len_ref)

            #if map_k > 10000:
            #    break
            


        while(map_kindel <= chromosome_size_hg19Kindel):

            map_current[map_kindel] = map_k

            map_k       = map_k + 1
            map_kindel  = map_kindel + 1
   
            #print map_k
            #break

        #print  str((map_kindel - 1)),"\t",map_current[(map_kindel - 1)]
        #print  chromosome_size_hg19Kindel,"\t",map_current[chromosome_size_hg19Kindel]

        if map_current[(map_kindel - 1)] == chromosome_size_hg19:

            print " Sanity checked passed for ",chrm_refSeq


        f_variant_vcf.close()

        #for key in map_current:      
        #    print key,"\t",map_current[key]              



    pos = int(cols[1])

    try:

        mapped_pos = map_current[pos]

        write_line = chrm_refSeq + "\t" + str(mapped_pos) + "\t" + cols[2] + "\t" + cols[3] + "\t" + cols[4] + "\t" + cols[5] + "\t" + cols[6] + "\t" + cols[7] + "\t" + cols[8] + "\t" + cols[9] + "\n"

        fout.write(write_line)

    except KeyError,e:
        funmapped.write(line) 
 
  

    current_chr=chrm_refSeq

    #break

    





