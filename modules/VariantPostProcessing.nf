
process GatherVCFPanels {
	// collapse files from individual intervals into one file per panel
	tag "GatherVCFPanels on $panelName using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/mergedVCFs/", mode:'copy'
	container 'broadinstitute/gatk:4.2.3.0'
	label "s_cpu"
	label "m_mem"
	
	input:
	tuple val(panelName), val(intervalNames), path(vcfsPanelIntervals)

	output:
	tuple val(panelName), path("${panelName}.vcf.gz")

	script:
	input = vcfsPanelIntervals.collect{"-I ${it}"}.join(' ')
	"""
	gatk --java-options "-Xmx${task.memory.toMega()}m" \
	GatherVcfs \
		$input \
		--REORDER_INPUT_BY_FIRST_VARIANT \
		-O ${panelName}.vcf.gz
	"""
}

process FilterVCFPanels {
	// filter variants for bedfile for respective panel
    // applicable for virtual panels mostly
	tag "FilterVCFPanels on $panelName using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/mergedVCFs/", mode:'copy'
	container 'registry.gitlab.ics.muni.cz:443/450402/btk_k8s:16'
	label "s_cpu"
	label "xxs_mem"
	
	input:
	tuple val(panel), path(mergedVcf)

	output:
	tuple val(panel), path("${panelName}.filtered.vcf.gz")

	script:
	"""
	bedtools intersect -header -a $mergedVcf -b ${panel.bedFile} | bgzip > ${panelName}.filtered.vcf.gz
	"""
}