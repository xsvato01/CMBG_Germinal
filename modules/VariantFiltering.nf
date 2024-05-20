
process FilterPanelSNPs {
	// this is applicable for smaller panels
	tag "FilterPanelSNPs on $panelName: $intervalName using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/FilteredSNPs/", mode:'copy'
	label "gatk"
	label "s_cpu"
	label "s_mem"
	
	input:
	tuple val(panelName), val(intervalName), path(genotypedVcfInterval), path(tbi)

	output:
	tuple val(panelName), val(intervalName), path("${intervalName}_${panelName}.filtered.SNPs.vcf.gz")

	script:
	"""
	gatk --java-options -Xmx${task.memory.toMega()}m \
	SelectVariants \
		-V ${genotypedVcfInterval} \
		-select-type SNP \
		-O snps.vcf.gz

	gatk --java-options -Xmx${task.memory.toMega()}m \
	VariantFiltration \
		-V snps.vcf.gz \
		-filter "QD < 2.0" --filter-name "QD2" \
		-filter "QUAL < 30.0" --filter-name "QUAL30" \
		-filter "SOR > 3.0" --filter-name "SOR3" \
		-filter "FS > 60.0" --filter-name "FS60" \
		-filter "MQ < 40.0" --filter-name "MQ40" \
		-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
		-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
		-O ${intervalName}_${panelName}.filtered.SNPs.vcf.gz
	"""
}

process FilterPanelIndels {
	// this is applicable for smaller panels
	tag "FilterPanelIndels on $panelName: $intervalName using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/FilteredINDELs/", mode:'copy'
	label "gatk"
	label "s_cpu"
	label "s_mem"
	
	input:
	tuple val(panelName), val(intervalName), path(genotypedVcfInterval), path(tbi)

	output:
	tuple val(panelName),val(intervalName), path("${intervalName}_${panelName}.filtered.INDELs.vcf.gz")

	script:
	"""
	gatk --java-options -Xmx${task.memory.toMega()}m \
	SelectVariants \
		-V ${genotypedVcfInterval} \
		-select-type INDEL \
		-O indels.vcf.gz

	gatk --java-options -Xmx${task.memory.toMega()}m \
	VariantFiltration \
		-V indels.vcf.gz \
		-filter "QD < 2.0" --filter-name "QD2" \
		-filter "QUAL < 30.0" --filter-name "QUAL30" \
		-filter "FS > 200.0" --filter-name "FS200" \
		-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
		-O ${intervalName}_${panelName}.filtered.INDELs.vcf.gz
	"""
}


process MergeFilteredVariants {
	// this is applicable for smaller panels
	tag "MergeFilteredVariants on $panelName: $intervalName using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/FilteredMerged/", mode:'copy'
	label "gatk"
	label "s_cpu"
	label "s_mem"
	
	input:
	tuple val(panelName), val(intervalName), path(filteredSNPvcf), path(filteredINDELSvcf)

	output:
	tuple val(panelName), val(intervalName), path("${intervalName}_${panelName}.filtered.merged.vcf.gz")

	script:
	"""
	gatk --java-options -Xmx${task.memory.toMega()}m \
	MergeVcfs \
		--INPUT ${filteredSNPvcf} \
		--INPUT ${filteredINDELSvcf} \
		--OUTPUT filtered.vcf.gz

	gatk --java-options -Xmx${task.memory.toMega()}m \
	SelectVariants \
		-R ${params.ref}.fasta \
		-V filtered.vcf.gz \
		--exclude-filtered \
		-O ${intervalName}_${panelName}.filtered.merged.vcf.gz
	"""
}