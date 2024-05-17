
process BaseRecalibration {
	tag "BaseRecalibration on $sample.name-${intervalBed.getSimpleName()} using $task.cpus CPUs and $task.memory memory"
	container 'broadinstitute/gatk:4.2.3.0'
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
	container 'broadinstitute/gatk:4.2.3.0'
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

process ApplyRecalibration {
	tag "ApplyRecalibration on $sample.name-${intervalBed.getSimpleName()} using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/reCals/", mode:'copy'
	container 'broadinstitute/gatk:4.2.3.0'
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
	container 'broadinstitute/gatk:4.2.3.0'
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

//#####

// process VAR_RECALL {
// 	// NOT APLICABLE FOR LOW NUMBER OF SAMPLES
// 	// https://gatk.broadinstitute.org/hc/en-us/articles/360036510892-VariantRecalibrator
// 	tag "VAR_RECALL on $panelName: $intervalName using $task.cpus CPUs and $task.memory memory"
// 	container 'broadinstitute/gatk:4.2.3.0'
// 	label "s_cpu"
// 	label "m_mem"

// 	input:
// 	tuple val(panelName), val(intervalName), path(genotypedVcfInterval), path(tbi)

// 	output:
// 	tuple val(panelName), val(intervalName), path("${intervalName}_${panelName}*")

// 	script:
// 	"""
// 	gatk --java-options "-Xmx${task.memory.toMega()}m -XX:ParallelGCThreads=${task.cpus}" \
// 	VariantRecalibrator \
// 		-R ${params.ref}.fasta \
// 		-V ${genotypedVcfInterval} \
// 		--resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${params.gatkKnownSites}/hapmap_3.3.hg38.vcf.gz \
// 		--resource:omni,known=false,training=true,truth=false,prior=12.0 ${params.gatkKnownSites}/1000G_omni2.5.hg38.vcf.gz \
// 		--resource:1000G,known=false,training=true,truth=false,prior=10.0 ${params.gatkKnownSites}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
// 		--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.gatkKnownSites}/dbsnp_146.hg38.vcf.gz  \
// 		-an QD -an ReadPosRankSum -an FS -an SOR -an DP \
// 		-mode SNP \
// 		--max-gaussians 1 \
// 		-O ${intervalName}_${panelName}.recal \
// 		--tranches-file ${intervalName}_${panelName}.tranches \
// 		--rscript-file ${intervalName}_${panelName}.plots.R
// 	"""
// }


// process APPLY_VQSR {
// 	// NOT APLICABLE FOR LOW NUMBER OF SAMPLES
// 	// https://gatk.broadinstitute.org/hc/en-us/articles/360037056912-ApplyVQSR
// 	tag "APPLY_VQSR on $panelName: $intervalName using $task.cpus CPUs and $task.memory memory"
// 	container 'broadinstitute/gatk:4.2.3.0'
// 	label "s_cpu"
// 	label "m_mem"

// 	input:
// 	tuple val(panelName), val(intervalName), path(var_recall_files)

// 	output:
// 	tuple val(panelName), val(intervalName), path("${intervalName}_${panelName}*")

// 	script:
// 	"""
// 	gatk --java-options "-Xmx${task.memory.toMega()}m -XX:ParallelGCThreads=${task.cpus}" \
// 	ApplyVQSR \
// 		-R ${params.ref}.fasta \
// 		-V ${intervalName}_${panelName}.vcf.gz \
// 		-O ${intervalName}_${panelName}.recalled.vcf.gz  \
// 		--ts_filter_level 99.0 \
// 		--tranches-file ${intervalName}_${panelName}.tranches \
// 		--recal-file ${intervalName}_${panelName}.recal \
// 		-mode SNP
// 	"""
// }