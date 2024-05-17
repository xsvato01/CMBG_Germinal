
process AnnotateVariantsVEP {
	tag "AnnotateVariantsVEP on $panelName: $intervalName using $task.cpus CPUs $task.memory"
	container "ensemblorg/ensembl-vep:release_108.0"
	publishDir "${params.outDirectory}/Annotate/", mode:'copy'
	label "s_cpu"
	label "l_mem"

	input:
	tuple val(panelName), val(intervalName), path(vcfRegion)
	
	output:
	tuple val(panelName), val(intervalName), path("${panelName}_${intervalName}.vep.vcf")

	script:
	"""
	vep \
		-i $vcfRegion \
		--cache \
		--cache_version 108 \
		--format vcf \
		--dir_cache $params.vep \
		--fasta ${params.ref}.fasta \
		--merged \
		--mane_select \
		--offline \
		--vcf \
		--everything \
		-o ${panelName}_${intervalName}.vep.vcf
	"""	
}

process FilterVEPTranscripts {
	tag "FilterVEPTranscripts on $panelName: $intervalName using $task.cpus CPUs $task.memory"
	container "ensemblorg/ensembl-vep:release_108.0"
	publishDir "${params.outDirectory}/Annotate/", mode:'copy'
	label "s_cpu"
	label "l_mem"

	input:
	tuple val(panelName), val(intervalName), path(vcfRegion)
	
	output:
	tuple val(panelName), val(intervalName), path("${panelName}_${intervalName}.MANE.vep.vcf")

	script:
	"""
	filter_vep \
		--input_file $vcfRegion \
		--output_file ${panelName}_${intervalName}.MANE.vep.vcf \
		--format vcf \
		--only_matched \
		--filter "(CANONICAL is YES) or (MANE_SELECT is YES) or (MANE_PLUS_CLINICAL is YES)"
	"""	
}