process FastQC {
    tag "FastQC on $sample.name using $task.cpus CPUs and $task.memory memory"
    container "staphb/fastqc:0.12.1"
   	// container 'registry.gitlab.ics.muni.cz:443/450402/btk_k8s:16'	
	label "s_cpu"
	label "m_mem"
	publishDir "${params.outDirectory}/QC/", mode:'copy'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path ("*")

    script:
    """
    fastqc $reads -o ./
    rm -r ?
    """
}

process BamQC {
	tag "BamQC on $sample.name using $task.cpus CPUs and $task.memory memory"
    // container 'registry.gitlab.ics.muni.cz:443/450402/btk_k8s:16'

	container 'quay.io/biocontainers/mulled-v2-b0664646864bfdb46c5343b1b2b93fc05adb4b77:39a005770a3e30fb6aa3bf424b57ddf52bae7ece-0'	
	
	input:
	tuple val(sample), path(bam), path(bai)

	output:
    tuple val(sample), path ("*")

	script:
	"""
	samtools flagstat $bam > ${sample.name}.flagstat
	samtools stats $bam > ${sample.name}.samstats
	picard BedToIntervalList I=${sample.bedFile} O=${sample.name}.interval_list SD=${params.ref}.dict
	picard CollectHsMetrics I=$bam BAIT_INTERVALS=${sample.name}.interval_list TARGET_INTERVALS=${sample.name}.interval_list R=${params.ref}.fasta O=${sample.name}.aln_metrics
	"""
}


process MultiQC {
    tag "MultiQC on $runName using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/QC/", mode:'copy'
    container 'multiqc/multiqc:v1.22.1'
	label "s_cpu"
	label "s_mem"

    input:
	tuple val(runName), path(qcResults)

    output:
    path "${runName}.QCreport.html"

    script:
    //     Rscript --vanilla ${params.coverstat} ${launchDir}/coverage $run
    // Rscript --vanilla ${params.covcompare} $run ${params.seqstats}
    // cat PerExon.html | perl  -pe 's/[^[:ascii:]]//g;' > PerExon_mqc.html
    // echo "log_filesize_limit: 30000000" > multiqc_config.yaml
    """
    multiqc . -n ${runName}.QCreport.html
    """
}