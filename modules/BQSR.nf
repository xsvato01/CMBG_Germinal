
process BQSR {
	tag "BaseRecalibration on $sample.name-${intervalBed.getSimpleName()} using $task.cpus CPUs and $task.memory memory"
	label "gatk"
	label "s_cpu"
	label "s_mem"

	input:
	tuple val(panelName),val(sample), path(bam), path(bai), path(intervalBed)

	output:
	tuple val(sample), path("*.recal.table")

	script:
	"""
	gatk --java-options -Xmx${task.memory.toMega()}m \
	BaseRecalibrator \
		-I ${bam} \
		-O ${intervalBed.getSimpleName()}_${sample.name}.recal.table \
		-L ${intervalBed} \
		--tmp-dir . \
		-R ${params.ref}.fasta \
		--known-sites ${params.gatkKnownSites}/hapmap_3.3.hg38.vcf.gz \
		--known-sites ${params.gatkKnownSites}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
		--known-sites ${params.gatkKnownSites}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
		--known-sites ${params.gatkKnownSites}/dbsnp_146.hg38.vcf.gz \
		--verbosity INFO
	"""
}


process GatherRecalibrationTables {
	tag "GatherRecalibrationTables on $sample.name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/reCals/", mode:'copy'
	label "gatk"
	label "s_cpu"
	label "m_mem"

	input:
	tuple val(sample), path(intervals)

	output:
	tuple val(sample), path("${sample.name}.recal.table")

	script:
	input = intervals.collect{"-I ${it}"}.join(' ')
	"""
	gatk --java-options -Xmx${task.memory.toMega()}m  \
	GatherBQSRReports \
		$input \
		-O ${sample.name}.recal.table \
	"""
}

process ApplyBQSR {
	tag "ApplyBQSR on $sample.name-${intervalBed.getSimpleName()} using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/reCals/", mode:'copy'
	label "gatk"
	label "s_cpu"
	label "s_mem"

	input:
	tuple val(panelName), val(sample), path(recal_table), path(bam), path(bai), path(intervalBed)

	output:
	tuple val(sample), path("${intervalBed.getSimpleName()}_${sample.name}.recal.bam")

	script:
	"""
	gatk --java-options -Xmx${task.memory.toMega()}m \
	ApplyBQSR \
		-R ${params.ref}.fasta \
		--input ${bam} \
		--output ${intervalBed.getSimpleName()}_${sample.name}.recal.bam \
		-L ${intervalBed} \
		--bqsr-recal-file ${recal_table}
	"""
}


process MergeRecalibratedBams {
	tag "MergeRecalibratedBams on $sample.name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/bams/", mode:'copy'
	label "gatk"
	label "m_cpu"
	label "xxs_mem"

	input:
	tuple val(sample), path(bams)

	output:
	tuple val(sample), path("${sample.name}.recal.bam"), path("${sample.name}.recal.bam.bai")

	script:
	"""
	samtools merge --threads ${task.cpus} ${sample.name}.recal.bam ${bams}
	samtools index ${sample.name}.recal.bam
	"""
}


process GatherGenotypedVCFPanels {
	// collapse files from individual intervals into one file per panel
	tag "GatherGenotypedVCFPanels on $panelName using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/mergedVCFs/", mode:'copy'
	label "gatk"
	label "s_cpu"
	label "m_mem"
	
	input:
	tuple val(panelName), val(intervalNames), path(vcfsPanelIntervals), path(tbisPanelIntervals)

	output:
	tuple val(panelName), path("${panelName}.vcf.gz"), path("${panelName}.vcf.gz.tbi")

	script:
	input = vcfsPanelIntervals.collect{"-I ${it}"}.join(' ')
	"""
	gatk --java-options "-Xmx${task.memory.toMega()}m" \
	GatherVcfs \
		$input \
		--REORDER_INPUT_BY_FIRST_VARIANT \
		-O ${panelName}.vcf.gz

	gatk --java-options -Xmx${task.memory.toMega()}m \
	IndexFeatureFile \
		-I ${panelName}.vcf.gz
	"""
}

