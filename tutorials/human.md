## Tutorial
This tutorial provides a step-by-step walkthrough for calling genetic variants and performing related analyses. It is divided into three main sections: installing the required tools, downloading the necessary data, and running the variant-calling pipeline followed by downstream analysis. The workflow is demonstrated using human represented by a single sample. Variant calling for other species follows the same procedure, differing only in species names.<br>    
The tutorial assumes you are working in a Unix/Linux environment. Several directories will be generated throughout the process, and their structure is shown in Fig. 1. To follow along, you will need to install Python and Java-specifically Python 3.6 or later, and JDK 1.8 (required for GATK). Although the tutorial references tool versions used at the time it was written, you may use the latest available versions.<br>    
In the examples, the ‘$’ symbol indicates the command prompt, and lines beginning with ‘#’ represent comments that should be omitted when running commands. All executable commands are shown in italics.

![](https://user-images.githubusercontent.com/63629577/209596455-a7696db8-98a9-483e-a39c-ebd71579813e.png)   
*Fig. 1 : The overall structure of the directories.*

<br>

## Part I: Install tools
1.	Create a directory  
Create the directory “tools” in your home directory.   
___$mkdir tools___  

2.	Go to the directory “tools”, and download and install the following tools.

    *	__BWA__: https://sourceforge.net/projects/bio-bwa/files/
    *	__Samtools__: https://github.com/samtools/samtools/releases/
    *	__Picard__: https://github.com/broadinstitute/picard/releases/
    *	__GATK__: https://github.com/broadinstitute/gatk/releases/    
    
    (note) BWA and Samtools may require libraries (e.g., bzip2-devel, ncurses-devel, xz-devel, zlib-devel, and curl-devel, etc) installed depending on the open source linux operating system.    
3.	Download and install BWA(Burrows-Wheeler Aligner) using the following commands.  
___$wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.12.tar.bz2___       &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# download  
___$bunzip2 bwa-0.7.12.tar.bz2___  	        &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# unzip and untar file  
___$tar xvf bwa-0.7.12.tar___  
___$mv bwa-0.7.12  bwa___   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# change directory name  
___$cd bwa___    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# go to directory bwa and install BWA  
___$make___  
___$make install___  

    (note) The up-to-date versions of bwa and bwa2 are bwa-0.7.17 (Nov 7, 2017, https://sourceforge.net/projects/bio-bwa/files/) and bwa-mem2-2.3 (Jun 30, 2025, https://github.com/bwa-mem2/bwa-mem2/releases/), respectively. 

4.	Download and install Samtools using the following commands.   
___$wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2___   
___$bunzip2 samtools-1.16.1.tar.bz2___		&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# unzip and untar file   
___$tar xvf samtools-1.16.1.tar___   
___$mv samtools-1.16.1 	 samtools___		&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# change directory name   
___$cd samtools___  				&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# go to samtools and install samtools   
___$make___    
___$make install___   
       (note) The most recent releases is samtools-1.22.1 (Jul 15, 2025) available from its respective project page: https://github.com/samtools/samtools/releases/ 
     
5.	Download picard using the following commands.   
___$mkdir picard___			&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# create a directory under directory tools   
___$cd picard___  			&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# go to directory picard   
___$wget https://github.com/broadinstitute/picard/releases/download/2.26.0/picard.jar___   
  
    (note) Make sure JDK version 1.8 has been installed.   
    (note) Note: The most recent release is picard.3.4.0 (Apr 14, 2025) available from its  respective project page: https://github.com/broadinstitute/picard/releases/       


6.	Download and install GATK using the following command.  
___$wget https://github.com/broadinstitute/gatk/releases/download/4.6.2.0/gatk-4.6.2.0.zip___      	 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# download GATK   
___$unzip gatk-4.6.2.0.zip___      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# unzip and untar file   
___$mv gatk-4.6.2.0 gatk___    
___$cd gatk___      
___$mv gatk-package-4.6.0-local.jar  GenomeAnalysisTk.jar___         &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# change file name   
    (note) The most recent release is gatk-4.6.2.0 (Apr 15, 2025) available from its  respective project page: https://github.com/broadinstitute/gatk/releases/    
	
<br>    

## Part II: Data download
1. Create the required directories  
    a. If you are working with human data, create a directory named "human" in your home directory make a directory.    
___$mkdir human___      
    b. In the "human" directory, create two subdirectories, "data" and "module" (see Fig. 1).    
   ___$cd human___    
   ___$mkdir data module___      
    c. Navigate to the "data" directory and create the three sub-directories: "fastq", "ref", and "db (see Fig. 1).    
   ___$cd data___    
   ___$mkdir fastq ref db___    

2. Download the following datasets into the directory “data” directory.  
    *	FASTQ
    *	reference sequence
    *	databases of variants

3.	Navigate to the “fastq” directory and download FASTQ file
    (note) You can download more than one sample for variant calling, but in this tutorial we will focus on using just a single sample.      <br>       
    a. Visit the website: https://www.internationalgenome.org/data-portal/sample.       <br>
  	
    b. Search for the sample “HG00096” (see Fig. 2).    <br>
	
  	![image](https://user-images.githubusercontent.com/63629577/209597435-7c156350-bb4a-4d1d-9b73-220ea83d35ff.png)     
    *Fig. 2: https://www.internationalgenome.org/data-portal/sample.*        <br>
	
    c. Select “HG00096” under the 1 matching sample results.
  	
    d. Under Data Types, choose “sequence,” and under Technologies, select “Low coverage WGS.” For sample HG00096, you will see six available FASTQ files (see Fig. 3).    
  	![image](https://user-images.githubusercontent.com/63629577/209597483-24b1a42b-becb-40e6-af57-b8bf25a463e8.png)   
    *Fig. 3: Result of searching HG00096.*
  	
    e. Navigate to the directory “fastq” and download the matching data (FASTQ) files.    
   ___$cd fastq___    
   ___$wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR062/SRR062634/SRR062634_1.fastq.gz___    
   ___$wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR062/SRR062634/SRR062634_2.fastq.gz___    
   ___$wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR062/SRR062635/SRR062635_1.fastq.gz___    
   ___$wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR062/SRR062635/SRR062635_2.fastq.gz___    	   
   ___$wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR062/SRR062641/SRR062641_1.fastq.gz___	       
   ___$wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR062/SRR062641/SRR062641_2.fastq.gz___    	   

    f. Combine the FASTQ files and rename the resulting merged files.    
    __$zcat SRR062634_1.fastq.gz SRR062635_1.fastq.gz SRR062641_1.fastq.gz | gzip -c > HG00096_1.fastq.gz__    
    __$zcat SRR062634_2.fastq.gz SRR062635_2.fastq.gz SRR062641_2.fastq.gz | gzip -c > HG00096_2.fastq.gz__    

4.	Navigate to the “ref” directory and download the reference sequence
   ___$cd ref___    	    
  ___$wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa___    	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# download human reference sequence   
  
5.	Navigate to the “db” directory and download two variant databases: dbSNP and pseudoDB.  <br>
    a.	Download dbSNP of human and assign a new name to it.   
      ___$wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz___  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# download   
      ___$mv 00-All.vcf.gz      dbSNP_dbSNP.vcf.gz___        &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# change DB name    

    b.	Download the pseudoDB of human at
      https://zenodo.org/record/7488070/files/human_pseudoDB.vcf.gz?download=1   
      
<br><br><br>
## Part III: Variant calling with analysis
1.	Download “gatk.py” module from the github repository into directory “tools”.   
	___$curl -L -O https://github.com/infoLab204/pseudo_DB/raw/main/gatk.py___  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# download “gatk.py” module   
<br>

2.	Go to the directory “tools” and import the module as follows.   

     ___>>>import  gatk___        &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# import the “gatk.py” module   
  
    (note) The “gatk.py” module contains the following functions:   
    *	___set_wd( )___: set working directory   
    *	___pre_align( )___: create files from reference sequence for alignment   
    *   ___align_fastq( )___: align FASTQ to reference sequence     
    *	___pseudo_db( )___: construct pseudo-database    
    *	___qs_recal( )___: recalibrate base quality score   
    *	___variant_call( )___: call genetic variants   
    *	___error_rate()____ : estimate error rate of sample   
    *	___qs_model( )___: estimate model-adjusted base quality score   

    (note) execute the above functions at directory “tools”.
<br>

3.	Create subdirectories under directory “module”. 
  
    ```
      Format: gatk.set_wd(“species_name”)   
    ```
    
     ___>>>gatk.set_wd(“human”)___          &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# create subdirectories   

    The list of subdirectories created under directory “module”:   
    *	___align___ : results of aligning FASTQ to reference   
    *	___error___ : result of estimating sample error rate   
    *	___machine___ : result of recalibrating machine-provided base quality score    
    *	___model___ : result of estimating model-adjusted base quality score   
    *	___variants___ : result of genetic variant calling   
    
<br>

4.	Create file names for the alignment under directory “ref”.    

    ```
    Format: gatk.pre_align(“species_name”, “reference_file”)   
    ```
  
	 ___>>>gatk.pre_align(“human”, “GRCh38_full_analysis_set_plus_decoy_hla.fa”)___   
    
    The following files are created in the directory “ref”:
    *	GRCh38_full_analysis_set_plus_decoy_hla.fa.amb
    *	GRCh38_full_analysis_set_plus_decoy_hla.fa.ann
    *	GRCh38_full_analysis_set_plus_decoy_hla.fa.bwt
    *	GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
    *	GRCh38_full_analysis_set_plus_decoy_hla.fa.pac
    *	GRCh38_full_analysis_set_plus_decoy_hla.fa.sa
    *	GRCh38_full_analysis_set_plus_decoy_hla.dict 
<br>

5.	Align FASTQ file of single sample (Case 1) or all samples (Case 2) to the reference.    

     ___Case 1 : align single sample___     
     
    ```
    Format: gatk.align_fastq(“species_name”, “reference”, “sample_name”)   
    ```
    ___>>>gatk.align_fastq(“human”, “GRCh38_full_analysis_set_plus_decoy_hla.fa”,”HG00096”)___   
    
    (note) Files HG00096_aligned.bam and HG00096_aligned.bai are created in the directory “align” with a sample HG00096.    <br>        

    <br>

     ___Case 2 : align all samples___           
    
    ```
    Format: gatk.align_fastq(“species_name”, “reference”)   
    ```
    
    ___>>>gatk.align_fastq(“human”, “GRCh38_full_analysis_set_plus_decoy_hla.fa”)___   
    
    (note) Files human_aligned.bam and human_aligned.bai are created in the directory “align” with all human samples.    
<br>


6.	Construct a pseudo database.   
    ```
    Format: gatk.pseudo_db(“species_name”, “reference”)   
    ```
    (note) Constructed pseudoDB used all samples in the “align” directory.        
    <br>
    
    ___>>>gatk.pseudo_db(“human”,“GRCh38_full_analysis_set_plus_decoy_hla.fa”)___       
    
    (note) File “human_pseudoDB.vcf” and “human_pseudoDB.vcf.idx” are created in the directory “db”.    
    
<br>

7.	Recalibrate machine-provided base quality score from single sample (Case 1) or all samples (Case 2).    

     ___Case 1 : recalibrate single sample___   

    ```
	  Format: gatk.qs_recal(“species_name”, “reference”, “name of database”, “db_type”, “sample_name”)   
    ```
    
    (note) The argument “db_type” can be either “dbSNP” or “pseudoDB”   <br><br>
    
     ___>>>gatk.qs_recal(“human”, “GRCh38_full_analysis_set_plus_decoy_hla.fa”, “dbSNP_b151.vcf”, “dbSNP”, ”HG00096”)___     
     
     (note) Files HG00096_dbSNP_recalibrated.bam and HG00096_dbSNP_recalibrated.bai are created in the directory “machine”.  <br><br> 

     ___>>>gatk.qs_recal(“human”,“GRCh38_full_analysis_set_plus_decoy_hla.fa”, “human_pseudoDB.vcf”, “pseudoDB”, ”HG00096”)___      
     
     (note) Files HG00096_pseudoDB_recalibrated.bam and HG00096_pseudoDB_recalibrated.bai are created in the directory “machine”.  <br><br> 
     
     
     ___Case 2 : recalibrate all samples___   
 
    ```
	  Format: gatk.qs_recal(“species_name”, “reference”, “name of database”, “db_type”)   
    ```
      
     ___>>>gatk.qs_recal(“human”, “GRCh38_full_analysis_set_plus_decoy_hla.fa”, “dbSNP_b151.vcf”, “dbSNP”)___     
     
     (note) Files human_dbSNP_recalibrated.bam and human_dbSNP_recalibrated.bai are created in the directory “machine”.  <br><br> 

     ___>>>gatk.qs_recal(“human”,“GRCh38_full_analysis_set_plus_decoy_hla.fa”, “human_pseudoDB.vcf”, “pseudoDB”)___    
     
     (note) Files human_pseudoDB_recalibrated.bam and human_pseudoDB_recalibrated.bai are created in the directory “machine”.   

    
<br>

8.	Call genetic variants.   

    ```
	  Format: gatk.variant_call(“species_name”, “reference”, “db_type”)   
    ```
  
     ___>>>gatk.variant_call(“human”,“GRCh38_full_analysis_set_plus_decoy_hla.fa”, “dbSNP”)___  
  
     (note) Files “human_dbSNP_variant_calling.vcf” and “human_dbSNP_variant_calling.vcf.idx” are created in the directory “variants”.   <br><br>
  
     ___>>>gatk.variant_call(“human”,“GRCh38_full_analysis_set_plus_decoy_hla.fa”, “pseudoDB”)___   
  
    (note) FIles “human_pseudoDB_variant_calling.vcf” and “human_pseudoDB_variant_calling.vcf.idx” are created in the directory “variants”.   
    
<br>

 9.	Estimate sample error rate   
  
    ```
	  Format: gatk.error_rate(“species_name”, “sample_name”, “reference”, “name of database”, “db_type”)   
    ```
    ___>>>gatk.error_rate(“human”,“HG00096”, “GRCh38_full_analysis_set_plus_decoy_hla.fa”, “dbSNP_b151.vcf”, “dbSNP”)___   
    
    (note) File “HG00096_dbSNP_erate” is created in the directory “error”.   <br><br>
    
    ___>>>gatk.error_rate(“human”,“HG00096”, “GRCh38_full_analysis_set_plus_decoy_hla.fa”, “human_pseudoDB.vcf”, “pseudoDB”)___   
    
    (note) File “HG00096_pseudoDB_erate” is created in the directory “error”.
  
  <br>

10.	Estimate model-adjusted base quality score.   

    ```
    Format: gatk.qs_model(“species_name”, “sample_name”, “db_type”)   
    ```
  
      ___>>>gatk.qs_model(“human”,“HG00096”, “dbSNP”)___   
  
      (note) File “HG00096_dbSNP_qs” is created in the directory “model”   <br><br>
  
      ___>>>gatk.qs_model(“human”,“HG00096”, “pseudoDB”)___  
  
      (note) File “HG00096_pseudoDB_qs” is created in directory “model”   
  
<br><br>
####  End of tutorial  

