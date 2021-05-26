#!/usr/bin/env nextflow

params.reads = "/home/vega/Semester-2/genomicsDataAnalysis-II/dataForAnalysis/fastqFiles/*_R[1,2].fastq.gz"

Channel
        .fromFilePairs(params.reads)
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
        .into {reads_ch; dada2ReadPairs}
        
process FastQC {

	publishDir "output/Process_1_fastqcBeforeTrim", mode: "link"
		
	label 'FastQC'
        
        input:
        set val(key), file(reads) from reads_ch
        
        output:
        file "*_fastqc.{zip,html}" into fastqc_ch
        
        script:
        """
        fastqc -q $reads
        """  
}

process MultiQC {

	publishDir "output/Process_2_multiqcBeforeTrim", mode: "link"
	
	label 'MultiQC'
        
        input:
        file("FastQC/*") from fastqc_ch.collect()
        
        output:
        file "multiqc_report.html" into multiqc_report
        file "multiqc_data"
        
        script:
        """
        multiqc .
        """
}

process filterAndTrim {
	
	publishDir "output/Process_3_FilterAndTrim", mode: "link"
	
	input:
	set pairId, file(reads) from dada2ReadPairs
	
	output:
	set val(pairId), "*_R1f.fastq.gz", "*_R2f.fastq.gz" into filteredReads /*Derep*/
	file "*_R1f.fastq.gz" into forReads /*Fastqc2*/
	file "*_R2f.fastq.gz" into revReads /*Fastqc2*/
	file "*.trimmed.txt" into trimTracking

    script:
    """
    #!/usr/bin/env Rscript
    library(dada2)
    
    
    out <- filterAndTrim(fwd = "${reads[0]}", # Roots for forward reads
    			filt = paste0("${pairId}", "_R1f.fastq.gz"),
    			
    			rev = "${reads[1]}", # Roots for reverse reads
    			filt.rev = paste0("${pairId}", "_R2f.fastq.gz"),
    	maxN = 0,
    	maxEE = c(2,2),
    	truncQ = 0, truncLen = 245,
    	rm.phix=FALSE, multithread=TRUE,
    	compress = T)
    	
    	write.csv(out, paste0("${pairId}", ".trimmed.txt"))
    """
}

process FastQC_2 {

	publishDir "output/Process_4_fastqcAfterTrim", mode: "link"

        label 'FastQC_2'
        
        input:
        file(readsforF) from forReads.collect()
        file(readsrevF) from revReads.collect()
        
        output:
        file "*_fastqc.{zip,html}" into fastqc_ch2
        
        script:
        """
        fastqc -q $readsforF $readsrevF
        """  
}

process MultiQC_2 {

	publishDir "output/Process_5_multiqcAfterTrim", mode: "link"
       
	label 'MultiQC_F'
        
        input:
        file("FastQC_F/*") from fastqc_ch2.collect()
        
        output:
        file "multiqc_report.html" into multiqc_report2
        file "multiqc_data"
        
        script:
        """
        multiqc .
        """
}

process Dereplication {

	publishDir "output/Process_6_dereplication", mode: "link"
	
	label "Dereplication"
	
	input:
	set val(pairId), file(filtFor), file(filtRev) from filteredReads
	
	output:
	file "derepf.RDS" into derepFor
	file "derepr.RDS" into derepRev
	file "derepf.RDS" into derepFordada
	file "derepr.RDS" into derepRevdada
	
	script:
	"""
	#!/usr/bin/env Rscript
	library(dada2)
	
	
	derepF <- derepFastq("${filtFor}")
	derepR <- derepFastq("${filtRev}")
	
    
    	saveRDS(derepF, "derepf.RDS")
    	saveRDS(derepR, "derepr.RDS")
    	"""
}

process ErrModel{

	publishDir "output/Process_7_errorModel", mode: "link"
	
	label "ErrorEstimation"
	
	input:
	file (dFs) from derepFor
	file (dRs) from derepRev
	
	output:
	file "errorsF.RDS" into errorsFor
	file "errorsR.RDS" into errorsRev
	
	script:
	"""
	#!/usr/bin/env Rscript
	library(dada2)	
	
	set.seed(100)
	
	derepF <- readRDS("${dFs}")
	derepR <- readRDS("${dRs}")

	errMF <- dada(derepF, err=NULL, selfConsist = T)
	errMR <- dada(derepR, err=NULL, selfConsist = T)
	
	saveRDS(errMF, "errorsF.RDS")
	saveRDS(errMR, "errorsR.RDS")
    	"""
}

process dada2{
	
	publishDir "output/Process_8_dada", mode: "link"
	
	label "dadaAlgo"
	
	input:
    	file (errMF) from errorsFor
    	file (errMR) from errorsRev
    	file (dFs) from derepFordada
	file (dRs) from derepRevdada

    	output:
    	file "merged.RDS" into mergedReads

	script:
	"""
	#!/usr/bin/env Rscript
	library(dada2)
	
	dF <- readRDS("${dFs}")
	dR <- readRDS("${dRs}")
	
	eMF <- readRDS("${errMF}")
	eMR <- readRDS("${errMR}")
    
	ddF <- dada(dF, err=eMF)
	ddR <- dada(dR, err=eMR)
	
	merger <- mergePairs(ddF, dF, ddR, dR, trimOverhang = TRUE)
	saveRDS(merger, "merged.RDS")
    	"""
}

process ASVnoChimera{

	publishDir "output/Process_9_ASV", mode: "link"
	
	label "ASVnoChimera"
	
	input:
	file (merged) from mergedReads
	
	output:
	file "ASV-Table.csv"
	file "ASV.RDS"
	
	script:
	"""
	#!/usr/bin/env Rscript
	library(dada2)
	
	mergers <- readRDS("${merged}")
	
	seqtab.all <- makeSequenceTable(mergers)

	seqtab <- removeBimeraDenovo(seqtab.all)
	
	write.csv(t(seqtab),"ASV-Table.csv")
	saveRDS(seqtab, "ASV.RDS")	
	"""
}
