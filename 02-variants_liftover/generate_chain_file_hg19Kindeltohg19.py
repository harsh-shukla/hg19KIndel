import re
import sys

############################################################## USAGE ###########################################################
#
# Run   : python generate_chain_file_hg19Kindeltohg19.py <INDELs_Replaced_1000Genome.vcf.gz> <hg19K.chr_length.txt> <hg19Kindel.chr_length.txt> 
#   
# INPUT :  
#  <INDELs_Replaced_1000Genome.vcf.gz> - the vcf containing all the indels we replaced in hg19 to create hg19Kindel
#  <hg19K.chr_length.txt>              - tab seperated file of chromosmes and their correspoding lengths in hg19K/hg19 reference
#  <hg19Kindel.chr_length.txt>         - tab seperated file of chromosmes and their correspoding lengths in hg19Kindel reference

#  OUTPUT :
#  <OUTPUT_1.fasta> - hg19Kindetohg19.over.chain file 
#
#  DESCRIPTION :
#       This python script does quasi-global-alignment between hg19 and hg19Kindel for every chromosome ; finds corresponding syntenic alignment blocks and outputs hg19Kindeltohg19.over.chain file 

#  NOTE        :  All the above input files are present in Supplementary_Files folders in the repository
##################################################################################################################################




vcf_f = open(sys.argv[1],"r")

hg19K_chrom_size_file      =  open(sys.argv[2],"r")
hg19Kindel_chrom_size_file =  open(sys.argv[3],"r")

fout = open("hg19Kindel_to_hg19.over.chain","w")

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


current_chr = ""
chain_id=0 
flag=0


total_ungapped_alignment = 0
total_dt=0
total_dq=0

for line in vcf_f:


    if line[0:1]=="#":
        continue
 
    
    cols = re.split("\t|\n",line)
    chrm = cols[0]

    if chrm != current_chr:   # A new chain alignment has started

        ## process the end of the older chromosome
 
        if flag==0:
            flag=1
       
        else:				                	 # process here older alignment
            pos = tSize - 1
            diff = pos - prev_pos
            size = diff - correction

            str_print = str(size)
            fout.write(str_print)
            fout.write("\n\n")

            total_ungapped_alignment = total_ungapped_alignment +  size
            calculated_size_target = total_ungapped_alignment + total_dt
            calculated_size_query  = total_ungapped_alignment + total_dq

            print "\n\n"
            print " Actual size of ",current_chr," of hg19K          :    ",tSize
            print " Actual size of ",current_chr," of hg19Kindel     :    ",qSize

            print " Calculated size of ",current_chr," of hg19K      :    ",calculated_size_target
            print " Calculated size of ",current_chr," of hg19Kindel :    ",calculated_size_query
            print "\n\n"

            total_ungapped_alignment = 0
            total_dt=0
            total_dq=0
            
        ## start the processing of the new chain block

        score   = 10000000                                       # chain score
        tName   = chrm                         			 # chromosome (reference sequence)
        tSize   = hg19K_chromSize[chrm]         		 # chromosome size (reference sequence)
        tStrand = "+"                           		 # strand (reference sequence)
        tStart  = 0                      		         # alignment start position (reference sequence)
        tEnd    = hg19K_chromSize[chrm]            	         # alignment end position (reference sequence)
        qName   = chrm                       			 # chromosome (query sequence)
        qSize   = hg19Kindel_chromSize[chrm] 			 # chromosome size (query sequence)
        qStrand = "+"						 # strand (query sequence)
        qStart  = 0						 # alignment start position (query sequence)
        qEnd    = hg19Kindel_chromSize[chrm]			 # alignment end position (query sequence)
        chain_id= chain_id + 1					 # chain ID

        #str_print = "chain" + " " + str(score) + " " + tName + " " + str(tSize) + " " + tStrand + " " + str(tStart) + " " + str(tEnd) + " " + qName + " " +  str(qSize) + " " + qStrand + " " + str(qStart) + " " + str(qEnd) + " " + str(chain_id)
        str_print  = "chain" + " " + str(score) + " " + qName + " " + str(qSize) + " " + qStrand + " " + str(qStart) + " " + str(qEnd) + " " + tName + " " +  str(tSize) + " " + tStrand + " " + str(tStart) + " " + str(tEnd) + " " + str(chain_id)  
        fout.write(str_print)
        fout.write("\n")

        pos  = int(cols[1]) - 1  # since vcf file is 1 based cordinate and we will be working with 0 based

        len_ref = len(cols[3])
        len_alt = len(cols[4]) 

        size = pos + 1
        dt   = len_ref - 1
        dq   = len_alt - 1

        #str_print = str(size) + "\t" +  str(dt) + "\t" + str(dq)
        str_print = str(size) + "\t" +  str(dq) + "\t" + str(dt)
        fout.write(str_print)
        fout.write("\n")

        prev_pos = pos
        correction = dt
        current_chr = chrm
        
        total_ungapped_alignment = total_ungapped_alignment +  size
        total_dt = total_dt + dt
        total_dq = total_dq + dq

        continue


    pos  = int(cols[1]) - 1  # since vcf file is 1 based cordinate and we will be working with 0 based

    len_ref = len(cols[3])
    len_alt = len(cols[4]) 

    if len_ref==1 and len_alt==1:
        continue  

    diff = pos - prev_pos
    size = diff - correction

    dt   = len_ref - 1
    dq   = len_alt - 1

    if dt > 0 and dq > 0:

        if dt > dq:
            size = size + dq
            w_dt = dt - dq
            w_dq = 0

        if dq > dt:
            size = size + dt
            w_dq = dq - dt          
            w_dt = 0 

    #str_print = str(size) + "\t" +  str(dt) + "\t" + str(dq)
    str_print = str(size) + "\t" +  str(w_dq) + "\t" + str(w_dt)
    fout.write(str_print)
    fout.write("\n")

    prev_pos = pos
    current_chr = chrm
    correction = dt

    total_ungapped_alignment = total_ungapped_alignment +  size
    total_dt = total_dt + dt
    total_dq = total_dq + dq



pos = tSize - 1
diff = pos - prev_pos
size = diff - correction

str_print = str(size)
fout.write(str_print)
fout.write("\n\n")

total_ungapped_alignment = total_ungapped_alignment +  size
calculated_size_target = total_ungapped_alignment + total_dt
calculated_size_query  = total_ungapped_alignment + total_dq

print "\n\n"
print " Actual size of ",current_chr," of hg19K          :    ",tSize
print " Actual size of ",current_chr," of hg19Kindel     :    ",qSize

print " Calculated size of ",current_chr," of hg19K      :    ",calculated_size_target
print " Calculated size of ",current_chr," of hg19Kindel :    ",calculated_size_query
print "\n\n"






