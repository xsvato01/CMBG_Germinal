// process BamQC {
// 	tag "first QC on $sample.name using $task.cpus CPUs and $task.memory memory"
// 	container 'registry.gitlab.ics.muni.cz:443/450402/btk_k8s:16'	
	
// 	input:
// 	tuple val(sample), path(bam)

// 	output:
// 	path "*"

// 	script:
// 	"""
// 	samtools flagstat $bam > ${sample.name}.flagstat
// 	samtools stats $bam > ${sample.name}.samstats
// 	picard BedToIntervalList I=${params.covbed} O=${sample.name}.interval_list SD=${params.ref}.dict
// 	picard CollectHsMetrics I=$bam BAIT_INTERVALS=${sample.name}.interval_list TARGET_INTERVALS=${sample.name}.interval_list R=${params.ref}.fasta O=${sample.name}.aln_metrics
// 	"""
// }


// process MULTIQC {
// 	tag "MultiQC using $task.cpus CPUs and $task.memory memory"
// 	publishDir "${params.outDirectory}/multiqc_reports/", mode:'copy'
// 	label "smallest_process"
// 	container "staphb/multiqc:1.19"
// 	// container 'registry.gitlab.ics.muni.cz:443/450402/tp53_nf:5'

// 	input:
// 	path '*'

// 	output:
// 	path '*.html'

// 	script:
// 	"""
// 	export LC_ALL=C.UTF-8
// 	export LANG=C.UTF-8
// 	multiqc . -n MultiQC-"`date +"%d-%m-%Y"`".html
// 	"""
// }