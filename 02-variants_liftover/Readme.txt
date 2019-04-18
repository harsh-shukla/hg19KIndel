-- Generating and using the chain file hg19Kindeltohg19

  1) generate_chain_file_hg19Kindeltohg19.py

     This python script does quasi-global-alignment between hg19 and hg19Kindel for every chromosome ; finds corresponding syntenic alignment blocks and outputs hg19Kindeltohg19.over.chain file 

  2) liftover_vcf_coordinates_using_chainfile.py

   Description : Python code to liftover coordinates of variants called (vcf file) on hg19Kindel to hg19. 

   Usage  : python convert_coordinates_vcf.py hg19Kindeltohg19.over.chain <INPUT_VCF_FILE.vcf>

   Input  : <INPUT_VCF_FILE.vcf> variants called using hg19Kindel as reference

   Output : 1) INPUT_VCF_FILE.vcf.CONVERTED.vcf- all the variants with positions lifted to hg19 coordinate frame
            2) INPUT_VCF_FILE.vcf.UNMAPPED.vcf - the variants that could not be mapped. (That particualr position is deleted in hg19)

   Dependency : PyLiftover 
                PyLiftover is a library for quick and easy conversion of genomic (point) coordinates between different assemblies.It 
                uses the same logic and coordinate conversion mappings as the UCSC liftOver tool.
                Ref:https://pypi.org/project/pyliftover/

   Note : The above python code converts genomic point coordinates only!!!. It in no way checks whether the reference allele in the
          target genome is consistent after liftover and neither does it update REF or ALT allele to reflect the change in genome. 
          hg19Kindeltohg19.over.chain file can be used with existing liftover tools such as UCSC liftOver,LiftoverVcf from 
          Picard,CrossMap etc. The given chain file has been verified by using it with LiftoverVcf from Picard Tools


-- Using the final <INDELs_Replaced_1000Genome.vcf>/<combined.vcf> to liftover cordiantes

   1) liftover_vcf_cordinates.py

      This python script liftsover coordiantes in an input vcf called using hg19Kindel as reference to hg19 coordinate frame using the final vcf containing all the indels we replaced in hg19 to make hg19Kindel.
