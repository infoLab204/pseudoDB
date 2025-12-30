## Tutorial
This tutorial provides a step-by-step walkthrough for calling genetic variants and performing related analyses. It is divided into three main sections: installing the required tools, downloading the necessary data, and running the variant-calling pipeline followed by downstream analysis. The workflow is demonstrated using sheep represented by a single sample.<br>   

The tutorial assumes you are working in a Unix/Linux environment. Several directories will be generated throughout the process, and their structure is shown in Fig. 1. To follow along, you will need to install Python and Java-specifically Python 3.6 or later, and JDK 1.8 (required for GATK). Although the tutorial references tool versions used at the time it was written, you may use the latest available versions.<br>    

In the examples, the ‘$’ symbol indicates the command prompt, and lines beginning with ‘#’ represent comments that should be omitted when running commands. All executable commands are shown in italics.

![](https://github.com/user-attachments/assets/2763063c-48c1-4f1d-83a6-f3c43a7db296)   
*Fig. 1 : The overall structure of the directories.*

<br>

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
    a. If you are working with sheep data, create a directory named "sheep" in your home directory.    
___$mkdir sheep___      <br>

    b. In the "sheep" directory, create two sub-directories, "data" and "module" (see Fig. 1).    
   ___$cd sheep___    
   ___$mkdir data module___      <br>
   
    c. Navigate to the "data" directory and create the three sub-directories: "fastq", "ref", and "db" (see Fig. 1).    
   ___$cd data___    
   ___$mkdir fastq ref db___    

2. Download the following datasets into the directory “data” directory.  
    *	FASTQ
    *	reference sequence
    *	databases of variants

3.	Navigate to the “fastq” directory and download FASTQ file
    (note) You can download more than one sample for variant calling, but in this tutorial we will focus on using just a single sample.      <br>       
    a. Visit the website: https://www.ebi.ac.uk/ena/browser/view/PRJNA160933.       <br>
  	
    b. Search for the sample "SRR501898".    <br>
	
    c. Navigate to the directory “fastq” and download the matching data (FASTQ) files.    
   ___$cd fastq___    
   ___$wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR501/SRR501898/SRR501898_1.fastq.gz___    
   ___$wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR501/SRR501898/SRR501898_2.fastq.gz___    


4.	Navigate to the “ref” directory and download the reference sequence    
   ___$cd ref___    	    
  ___$wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Ovis_aries/latest_assembly_versions/GCA_000298735.2_Oar_v4.0/GCA_000298735.2_Oar_v4.0_genomic.fna.gz___    	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# download sheep reference sequence   
  __$gzip -d GCA_000298735.2_Oar_v4.0_genomic.fna.gz__  # unzip
  __$ mv GCA_000298735.2_Oar_v4.0_genomic.fna oar_ref_Oar_v4.0_full.fa__  # change reference name

6.	Navigate to the “db” directory and download two variant databases: dbSNP and pseudoDB.  <br>
    a.	Download dbSNP of human and assign a new name to it.   
      ___$wget https://ftp.ncbi.nih.gov/snp/organisms/archive/sheep_9940/VCF/00-All.vcf.gz___  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# download   
      ___$mv 00-All.vcf.gz      sheep_dbSNP.vcf.gz___        &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# change DB name    

    b.	Download the pseudoDB of human    
      https://zenodo.org/record/7488425/files/sheep_pseudoDB.vcf.gz?download=1       
      
<br>

## Part III: Variant calling with analysis
1.	Download "gatk.py" module from the github repository into directory "tools".   
    ```
    $curl -L -O https://github.com/infoLab204/pseudoDB/pipeline/raw/main/gatk.py # download "gatk.py" module   
    ```

2.	Go to the directory "tools" and import the module as follows.   
    ```
    import  gatk       # import the "gatk.py" module   
    ```  
    (note) The "gatk.py" module contains the following functions:   
    *	___set_wd( )___: set working directory   
    *	___pre_align( )___: create files from reference sequence for alignment   
    *   ___align_fastq( )___: align FASTQ to reference sequence     
    *	___pseudo_db( )___: construct pseudo-database    
    *	___qs_recal( )___: recalibrate base quality score   
    *	___variant_call( )___: call genetic variants   
    *	___error_rate()____ : estimate error rate of sample   
    *	___qs_model( )___: estimate model-adjusted base quality score   

    (note) execute the above functions at directory "tools".

3.	Create subdirectories under directory "module".    <br>

       <b>Format: gatk.set_wd("species_name")</b>       
	
    ```
    gatk.set_wd("sheep") 
    ```

    The list of subdirectories created under directory "module":   
    *	__align__ : results of aligning FASTQ to reference   
    *	___error___ : result of estimating sample error rate   
    *	___machine___ : result of recalibrating machine-provided base quality score    
    *	___model___ : result of estimating model-adjusted base quality score   
    *	___variants___ : result of genetic variant calling   
    
<br>

4.	Create file names for the alignment under directory "ref".  <br>  


    <b>Format: gatk.pre_align("species_name", "reference_file")   </b>

    ```
	 gatk.pre_align("sheep", "oar_ref_Oar_v4.0_full.fa")   
    ```    
    The following files are created in the directory "ref":
    *	oar_ref_Oar_v4.0_full.fa.amb
    *	oar_ref_Oar_v4.0_full.fa.ann
    *	oar_ref_Oar_v4.0_full.fa.bwt
    *	oar_ref_Oar_v4.0_full.fa.fai
    *	oar_ref_Oar_v4.0_full.fa.pac
    *	oar_ref_Oar_v4.0_full.fa.sa
    *	oar_ref_Oar_v4.0_full.dict 
<br>

5.	Align FASTQ file of single samples to the reference.    <br>

    <b> Format: gatk.align_fastq("species_name", "reference", "sample_name") </b>       
    ```
    gatk.align_fastq("sheep", "oar_ref_Oar_v4.0_full.fa","SRR501898")       
    ``` 
    (note) Files SRR501898_aligned.bam and HG00096_aligned.bai are created in the directory “align” with a sample SRR501898.    <br>        

    <br>

6.	Construct a pseudo database.   

    <b>Format: gatk.pseudo_db("species_name", "reference") </b>     
    (note) Constructed pseudoDB used all samples in the “align” directory.        
    <br>
    ``` 
    gatk.pseudo_db("sheep","oar_ref_Oar_v4.0_full.fa")       
    ``` 
    (note) File "human_pseudoDB.vcf.gz" and "human_pseudoDB.vcf.gz.tbi" are created in the directory "db".    
    
<br>

7.	Recalibrate base quality score from samples    

	  <b>Format: gatk.qs_recal("species_name", "reference", "db_type", "sample_name")   </b>

    
    (note) The argument "db_type" can be either "dbSNP" or "pseudoDB"   <br><br>
     ```
     gatk.qs_recal("sheep", "oar_ref_Oar_v4.0_full.fa", "dbSNP", "SRR501898")     
     ```
     (note) Files SRR501898_dbSNP_recalibrated.bam and SRR501898_dbSNP_recalibrated.bai are created in the directory "machine".  <br><br> 
     ```
     gatk.qs_recal("sheep","oar_ref_Oar_v4.0_full.fa", "pseudoDB", "SRR501898")          
     ```     
     (note) Files SRR501898_pseudoDB_recalibrated.bam and SRR501898_pseudoDB_recalibrated.bai are created in the directory "machine".  <br><br> 
     
  <br>

8.	Call genetic variants.   

    a. Call variants per-sample
	  <b>Format: gatk.variant_call("species_name", "reference", "db_type","sample_name")  </b> 

    ```
    gatk.variant_call("sheep","oar_ref_Oar_v4.0_full.fa", "dbSNP","SRR501898")
    ```  
     (note) Files "SRR501898_dbSNP.g.vcf.gz" and "SRR501898_dbSNP.g.vcf.gz.tbi" are created in the directory "variants".   <br><br>
    ```
    gatk.variant_call("sheep","oar_ref_Oar_v4.0_full.fa", "pseudoDB","SRR501898") 
    ```
    (note) FIles "SRR501898_pseudoDB.g.vcf.gz" and "SRR501898_pseudoDB.g.vcf.gz.tbi" are created in the directory "variants".
  	
    b. Joint-Call Cohort
  	  <b>Format: gatk.variant_joint_call("species_name", "reference", "db_type")  </b>

    ```
    gatk.variant_joint_call("sheep","oar_ref_Oar_v4.0_full.fa", "dbSNP")
    ```  
     (note) Files "sheep_dbSNP_variants.vcf.gz" and "sheeo_dbSNP_variants.vcf.gz.tbi" are created in the directory "variants".   <br><br>
    ```
    gatk.variant_call("sheep","oar_ref_Oar_v4.0_full.fa", "pseudoDB") 
    ```
    (note) FIles "sheep_pseudoDB_variants.vcf.gz" and "sheep_pseudoDB_variants.vcf.gz.tbi" are created in the directory "variants".
<br>

 9.	Estimate sample error rate   
  
	  <b> Format: gatk.error_rate("species_name", "sample_name", "reference", "name of database", "db_type")   </b>
    ```
    gatk.error_rate("sheep","SRR501898", "oar_ref_Oar_v4.0_full.fa", "sheep_dbSNP.vcf.gz", "dbSNP")   
    ```    
    (note) File "SRR501898_dbSNP_erate" is created in the directory "error".   <br><br>
    ```    
    gatk.error_rate("sheep","SRR501898", "oar_ref_Oar_v4.0_full.fa", "sheep_pseudoDB.vcf.gz", "pseudoDB")   
    ```    
    (note) File "SRR501898_pseudoDB_erate" is created in the directory "error".
  
  <br>

10.	Estimate model-adjusted base quality score.   


    <b> Format: gatk.qs_model("species_name", "sample_name", "db_type")</b>   

    ``` 
    gatk.qs_model("sheep","SRR501898", "dbSNP")   
    ```  
      (note) File "SRR501898_dbSNP_qs" is created in the directory "model"   <br><br>
    ```  
    gatk.qs_model("sheep","SRR501898", "pseudoDB") 
    ```  
      (note) File "SRR501898_pseudoDB_qs" is created in directory "model"   
  
<br><br>
####  End of tutorial  
