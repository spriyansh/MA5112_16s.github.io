#!/usr/bin/env nextflow

/*
Nextflow script for 16S metagenomic analysis of amplicon sequence variants (ASVs) using DADA2

Group 4 AUTHORS:
Niall Garvey
Ben Nolan
Breena McCormack
Priyansh Srivastava

View on GitHub: https://spriyansh.github.io/MA5112_16s.github.io/

*/

/*About the dataset
Link to metatdata: https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP315998&o=acc_s%3Aa
Instrument: Illumina MiSeq
Isolation_source: gut
LibraryLayout: PAIRED
LibrarySelection: PCR
LibrarySource: METAGENOMIC
Organism: insect gut metagenome
Platform: ILLUMINA
ReleaseDate: 2021-04-23

Link to study: https://link.springer.com/article/10.1007/s00248-021-01756-1
*/


/*

Contributions

Ben Nolan:
1. Created Docker container with R, DADA2, multiqc, fastqc: https://hub.docker.com/repository/docker/bennolan/dada2_mf 
2. Created singularity image: /data/containers/dada2_mf.img

Breena McCormack:

Niall Garvey:

Priyansh Srivastava:
1. Integrated the nextflow-script
2. Selection of Dataset
3. Process designed-> filterAndTrim and Dereplication
*/

/*Saving path of the compressed fastq files*/
params.reads = "/home/vega/Semester-2/genomicsDataAnalysis-II/dataForAnalysis/fastqFiles/*_R[1,2].fastq.gz"

/*Reading channel; if the directory is empty, an error message is printed; read files are sent to two different process*/
Channel
        .fromFilePairs(params.reads)
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
        .into {reads_ch; dada2ReadPairs}
        

/*Author: Ben Nolan*/
process FastQC {

	/*Declaration of directory for publishing results*/
	publishDir "output/Process_1_fastqcBeforeTrim", mode: "link"
	
	/*Name for the process*/
	label 'FastQC'
        
        /*Input from Reading channel*/
        input:
        set val(key), file(reads) from reads_ch
        
        /*Output is sent through fastqc_ch*/
        output:
        file "*_fastqc.{zip,html}" into fastqc_ch
        
        /*Shell command for running fastQC*/
        script:
        """
        fastqc -q $reads
        """  
}

/*Author: Ben Nolan*/
process MultiQC {

	/*Declaration of directory for publishing results*/
	publishDir "output/Process_2_multiqcBeforeTrim", mode: "link"
	
	/*Name for the process*/	
	label 'MultiQC'
        
        /*Dynamic input from FastQC-process*/
        input:
        file("FastQC/*") from fastqc_ch.collect()
        
        /*Output is sent through multiqc_report; report is explicitly saved*/
        output:
        file "multiqc_report.html" into multiqc_report
        file "multiqc_data"
        
        /*Shell command for running multiQC*/
        script:
        """
        multiqc .
        """
}

/*Author: Priyansh Srivastava*/
process filterAndTrim {

	/*Declaration of directory for publishing results*/	
	publishDir "output/Process_3_FilterAndTrim", mode: "link"
	
	/*Name for the process*/	
	label 'FilterTrim'
	
	/*Input from Reading channel*/
	input:
	set pairId, file(reads) from dada2ReadPairs
	
	/*Saving and sending to multiple process all at once*/
	output:
	set val(pairId), "*_R1f.fastq.gz", "*_R2f.fastq.gz" into filteredReads /*to Derep*/
	file "*_R1f.fastq.gz" into forReads /*to Fastqc2*/
	file "*_R2f.fastq.gz" into revReads /*to Fastqc2*/
	file "*.trimmed.txt" into trimTracking
	
	/*Running R commands by calling R's interpreter*/
	script:
	"""
	#!/usr/bin/env Rscript
	library(dada2)
	
	# Saving as "out" object
	out <- filterAndTrim(fwd = "${reads[0]}", # Forward Reads
			filt = paste0("${pairId}", "_R1f.fastq.gz"), # Roots for forward reads
    			rev = "${reads[1]}",  # Reverse Reads
    			filt.rev = paste0("${pairId}", "_R2f.fastq.gz"), # Roots for reverse reads
    			maxN = 0, # Removing Ns
    			maxEE = c(2,2), # Maximum expected errors, using EE = sum(10^(-Q/10))
    			truncQ = 0, # Trim start
    			truncLen = 245, # Trim stop
    			rm.phix=FALSE, # Not discarding seqs that match against phix db
    			multithread=TRUE, # Parallel compute
    			compress = T) # Compressing
    	
    	# Writting Text file
    	write.csv(out, paste0("${pairId}", ".trimmed.txt"))
    	"""
}

/*Author: */
process FastQC_2 {

	/*Declaration of directory for publishing results*/
	publishDir "output/Process_4_fastqcAfterTrim", mode: "link"

	/*Name for the process*/
        label 'FastQC_2'
        
        /*Input from trimmed reads*/
        input:
        file(readsforF) from forReads.collect()
        file(readsrevF) from revReads.collect()
        
        /*Output is sent through fastqc_ch2*/
        output:
        file "*_fastqc.{zip,html}" into fastqc_ch2
        
        /*Shell command for running fastQC*/
        script:
        """
        fastqc -q $readsforF $readsrevF
        """  
}

/*Author: */
process MultiQC_2 {

	/*Declaration of directory for publishing results*/
	publishDir "output/Process_5_multiqcAfterTrim", mode: "link"
       
	/*Name for the process*/
	label 'MultiQC_2'
        
        /*Dynamic input from FastQC_2 process*/
	input:
	file("FastQC_F/*") from fastqc_ch2.collect()
        
        /*Output is sent through multiqc_report2; report is explicitly saved*/
	output:
	file "multiqc_report.html" into multiqc_report2
	file "multiqc_data"
        
        /*Shell command for running multiQC*/
	script:
	"""
	multiqc .
	"""
}

/*Author: Priyansh Srivastava*/
process Dereplication {

	/*Declaration of directory for publishing results*/
	publishDir "output/Process_6_dereplication", mode: "link"
	
	/*Name for the process*/
	label "Dereplication"
	
	/*Reading the trimmed files collected from filterAndTrim process*/
	input:
	set val(pairId), file(filtFor), file(filtRev) from filteredReads
	
	/* Saving and sending output to two process simultaneously*/
	output:
	file "derepf.RDS" into derepFor
	file "derepr.RDS" into derepRev
	file "derepf.RDS" into derepFordada
	file "derepr.RDS" into derepRevdada
	
	/*Running R commands by calling R's interpreter*/
	script:
	"""
	#!/usr/bin/env Rscript
	# Calling dada library
	library(dada2)
	
	# Dereplication of trimmed sequences
	derepF <- derepFastq("${filtFor}")
	derepR <- derepFastq("${filtRev}")
	
	# Saving as R data structure
    	saveRDS(derepF, "derepf.RDS")
    	saveRDS(derepR, "derepr.RDS")
    	"""
}

/*Author: */
process ErrModel{

	/*Declaration of directory for publishing results*/
	publishDir "output/Process_7_errorModel", mode: "link"
	
	/*Name for the process*/
	label "ErrorEstimation"
	
	/*Input from Dereplication process*/
	input:
	file (dFs) from derepFor
	file (dRs) from derepRev
	
	/*Saving to R's data structure; passing into two channels*/
	output:
	file "errorsF.RDS" into errorsFor
	file "errorsR.RDS" into errorsRev
	
	/*Running R commands by calling R's interpreter*/
	script:
	"""
	#!/usr/bin/env Rscript
	# Calling dada library
	library(dada2)	
	
	# Setting seed to make it reproducible
	set.seed(100)
	
	# Reading the dereplicated sequences
	derepF <- readRDS("${dFs}")
	derepR <- readRDS("${dRs}")

	# Calculating errors
	errMF <- dada(derepF, err=NULL, selfConsist = T)
	errMR <- dada(derepR, err=NULL, selfConsist = T)
	
	# Saving as R data structure
	saveRDS(errMF, "errorsF.RDS")
	saveRDS(errMR, "errorsR.RDS")
    	"""
}

/*Author: */
process dada2{

	/*Declaration of directory for publishing results*/	
	publishDir "output/Process_8_dada", mode: "link"
	
	/*Name for the process*/
	label "dadaAlgo"
	
	/*Input from ErrModel and Dereplication process*/
	input:
    	file (errMF) from errorsFor
    	file (errMR) from errorsRev
    	file (dFs) from derepFordada
	file (dRs) from derepRevdada

    	output:
    	file "merged.RDS" into mergedReads

	/*Running R commands by calling R's interpreter*/
	script:
	"""
	#!/usr/bin/env Rscript
	# Calling dada library
	library(dada2)
	
	# Reading the dereplicated sequences; for forward and reverse reads
	dF <- readRDS("${dFs}")
	dR <- readRDS("${dRs}")
	
	# Reading the error models; for forward and reverse reads
	eMF <- readRDS("${errMF}")
	eMR <- readRDS("${errMR}")
	
	# Running the DADA by passing the error rates and dereplicated sequences   
    	ddF <- dada(dF, err=eMF)
	ddR <- dada(dR, err=eMR)
	
	# Merging forward and reverse reads
	merger <- mergePairs(ddF, dF, ddR, dR, trimOverhang = TRUE)
	
	# Saving as R data structure
	saveRDS(merger, "merged.RDS")
    	"""
}

/*Author: */
process ASVnoChimera{
	
	/*Declaration of directory for publishing results*/
	publishDir "output/Process_9_ASV", mode: "link"
	
	/*Name for the process*/
	label "ASVnoChimera"
	
	/*Input from dada2 process*/
	input:
	file (merged) from mergedReads
	
	/*Saving to CSV and RDS, no further process*/
	output:
	file "ASV-Table.csv"
	file "ASV.RDS"
	
	/*Running R commands by calling R's interpreter*/
	script:
	"""
	#!/usr/bin/env Rscript
	# Calling dada library
	library(dada2)
	
	# Reading the R data-structure of mereged amplicons
	mergers <- readRDS("${merged}")
	
	# Extarcting Amplicon Sequence Variants
	seqtab.all <- makeSequenceTable(mergers)

	# Removing chimera
	seqtab <- removeBimeraDenovo(seqtab.all)
	
	# Saving as a CSV file
	write.csv(t(seqtab),"ASV-Table.csv")
	
	# Saving as R data structure
	saveRDS(seqtab, "ASV.RDS")	
	"""
}
