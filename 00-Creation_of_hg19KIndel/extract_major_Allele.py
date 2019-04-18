
############################################################## USAGE ###########################################################
#
#  Run   : python extract_major_Allele.py <INPUT.vcf> <OUTPUT_filename.vcf>
#   
#  INPUT :  
#  <INPUT.vcf> can be files filtered for a particular variant type like SNP,Indels etc
#
#  OUTPUT :
#  <OUTPUT_filename.vcf> - the filename of the output vcf ; will be created in same folder
#
#  DESCRIPTION :
#  This python script takes a vcf file as input and extracts the variants having alternate allele(AA) frequency greater than 0.5.In case of multi-allelic sites the alternate allele having the frequency > 0.5 is selectively extracted and ALT field of the vcf file is updated accordingly. The resulting vcf file generated has only one ALT allele(the major allele) for a given REF allele

##################################################################################################################################
      
import re
import sys


f_vcf = open(sys.argv[1],"r")

filename = sys.argv[2]

fout = open(filename,"w")


for line in f_vcf:

    if line[0:3]!="chr":   # process the header records write them as they are

        #print line
        fout.write(line)
        continue


    #print line

    col1 = re.split("\t|\n",line)
 
    col2 = re.split(";",col1[7])
  
    total_info_fields=len(col2)  

    #variant_index = total_info_fields - 1
    #print col2[variant_index]

    keys=[]
    values=[]
    

    for i in xrange(0,total_info_fields):

        col_temp = re.split("=",col2[i])
 
        try:

            values.append(col_temp[1])
            keys.append(col_temp[0])
            
        
        except IndexError:

            #print line
            #sys.exit()   
            continue 


    #print keys
    #print values
    flag_AF=-1 

    infomap = dict(zip(keys,values))
 
    try:
        #print infomap["VT"]
   
        current_AF = float(infomap["AF"])

        if current_AF > 0.5:
            fout.write(line)

    except ValueError:    # process the multiple alleles to find the major allele
        #print line
        #print keys
        #print values 
        #sys.exit()

 

        alt_freqs = infomap["AF"]
        #print alt_freqs

        AF = re.split(",",alt_freqs)

        number_of_AFs = len(AF)
        max_AF=0.0
     
        for j in xrange(0,number_of_AFs):

            current_AF_j = float(AF[j])

            if current_AF_j > max_AF:

                max_AF = current_AF_j
                flag_AF= j

        #print max_AF
        
        if max_AF > 0.5:

            AA = re.split(",",col1[4])
            major_allele = AA[flag_AF]

            #print flag_AF,col1[4],major_allele
       
            #col1[4] = str(major_allele)

            str_print = col1[0] + "\t" + col1[1] + "\t" + col1[2] + "\t" + col1[3] + "\t" + str(major_allele) + "\t" + col1[5] + "\t" + col1[6] +"\t" + col1[7] + "\n"       

            #print flag_AF,col1[4],major_allele,str_print 

            fout.write(str_print)
  
            #print str_print
            #break


    del keys
    del values
    del infomap

    #break

fout.close()
