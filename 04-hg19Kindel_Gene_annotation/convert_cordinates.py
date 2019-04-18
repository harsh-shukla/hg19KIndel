import re
import sys


## Input ##


f_refSeq = open(sys.argv[1],"r")



hg19K_chrom_size_file      =  open(sys.argv[3],"r")
hg19Kindel_chrom_size_file =  open(sys.argv[4],"r")


##  Output  ###
filename = "hg19Kindel_genes.refSeq"
fout = open(filename,"w")

funmapped = open("UnMapped.refSeq","w")


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


current_chr=""
map_current={}

for line in f_refSeq:  # Write the header line
    fout.write(line)

    write_line = "Error_Base" + "\t" + line
    funmapped.write(write_line)

    break
 

for line in f_refSeq:

    cols = re.split("\t|\n",line)
    chrm_refSeq = cols[2]
    #print chrm

    if chrm_refSeq not in include_chr:        
        fout.write(line)
        continue  

    if chrm_refSeq!=current_chr:

        del map_current
        map_current={}

        ## generate the mapping dictionary for every position in 

        f_vcf = open(sys.argv[2],"r")

        map_k=0
        map_kindel=0
        displaced = 0
        chromosome_size_hg19       = hg19K_chromSize[chrm_refSeq]
        chromosome_size_hg19Kindel = hg19Kindel_chromSize[chrm_refSeq]

        #print chromosome_size_hg19,chromosome_size_hg19Kindel

        for line in f_vcf:

            if line[0:3]!="chr":   # process the header records write them as they are
                continue                 

            col1 = re.split("\t|\n",line)
            chr_vcf    = col1[0]
       
            if chr_vcf!=chrm_refSeq:
                continue


            old_pos = int(col1[1]) - 1
            ref     = col1[3]
            alt     = col1[4]

            len_ref = len(ref)
            len_alt = len(alt)

            while (map_k <= old_pos):

                map_current[map_k] = map_kindel

                map_k       = map_k + 1
                map_kindel  = map_kindel + 1
                                         

            if len_ref > len_alt:

                counter=0
                while(counter < (len_alt-1)):
                    map_current[map_k] =map_kindel 
                    map_k       = map_k + 1
                    map_kindel  = map_kindel + 1
                    counter = counter + 1

                map_k_old = map_k
                map_k = map_k + (len_ref - len_alt)
                                                      
                for cmap in xrange(map_k_old,map_k):
                    map_current[cmap]=map_kindel


            if len_alt > len_ref:

                counter=0
                while(counter < (len_ref-1)):
                    map_current[map_k] = map_kindel
                    map_k       = map_k + 1
                    map_kindel  = map_kindel + 1
                    counter = counter + 1

                map_kindel = map_kindel + (len_alt - len_ref)                    


            

    
            #if map_k > 100000:
            #    break
            #break
 
      
        while(map_k < chromosome_size_hg19):

            map_current[map_k] = map_kindel

            map_k       = map_k + 1
            map_kindel  = map_kindel + 1
   
            #print map_k
            #break


        if map_current[(map_k - 1)] == (chromosome_size_hg19Kindel - 1):

            print " Sanity checked passed for ",chrm_refSeq
                
        #print     


        #for key in map_current:      
        #    print key,"\t",map_current[key]




    #bin    name    chrom   strand  txStart txEnd   cdsStart        cdsEnd  exonCount       exonStarts      exonEnds        score  name2   cdsStartStat    cdsEndStat      exonFrames    

    bin_rec       = cols[0]
    name          = cols[1]
#   chrm_refSeq   = cols[2]
    strand        = cols[3]
 
    txStart       = int(cols[4])    
    txEnd         = int(cols[5])
    cdsStart      = int(cols[6])
    cdsEnd        = int(cols[7])
    exonCount     = int(cols[8])

    exonStarts    = cols[9]
    exonEnds      = cols[10]
    score         = cols[11]
    name2         = cols[12]
    cdsStartStat  = cols[13]
    cdsEndStat    = cols[14]
    exonFrames    = cols[15]


    # change the cordinates here

    try:

        txStart_mapped        = map_current[txStart]
        txEnd_mapped          = map_current[txEnd]    
        cdsStart_mapped       = map_current[cdsStart] 
        cdsEnd_mapped         = map_current[cdsEnd]

        cols_exonStarts = re.split(",",exonStarts)
        cols_exonEnds   = re.split(",",exonEnds)

        if ( (len(cols_exonStarts) - 1) != exonCount ) or ( (len(cols_exonEnds) - 1) != exonCount ):
            print name

        exonStarts_mapped = ""
        exonEnds_mapped   = "" 

        for j in xrange(0,exonCount):

            current_starts = int(cols_exonStarts[j])
            current_ends   = int(cols_exonEnds[j])

            current_starts_mapped = map_current[current_starts]
            current_ends_mapped   = map_current[current_ends]

            exonStarts_mapped = exonStarts_mapped + str(current_starts_mapped) + ","
            exonEnds_mapped   = exonEnds_mapped + str(current_ends_mapped) + ","


        write_line = bin_rec + "\t" + name + "\t" + chrm_refSeq + "\t" + strand + "\t" + str(txStart_mapped) + "\t" + str(txEnd_mapped) + "\t" + str(cdsStart_mapped) + "\t" + str(cdsEnd_mapped) + "\t" + str(exonCount) + "\t" + exonStarts_mapped + "\t" + exonEnds_mapped + "\t" + score + "\t" + name2 + "\t" +  cdsStartStat + "\t" + cdsEndStat + "\t" + exonFrames + "\n"


        fout.write(write_line)        


    except KeyError,e:
        
        print str(e)
        write_line = str(e) + "\t" + line
        funmapped.write(write_line) 
 
 
         
 

    current_chr=chrm_refSeq


    #break 



fout.close()
funmapped.close()







