
process MakeSitesOnlyVcfs {
	// should speed up VQSR
	// https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
	tag "MakeSitesOnlyVcfs on $panelName using $task.cpus CPUs and $task.memory memory"
	label "gatk"
	label "s_cpu"
	label "s_mem"

	input:
	tuple val(panelName), path(vcfPanel)

	output:
	tuple val(panelName), path("${panelName}_sitesOnly.vcf.gz"), path("${panelName}_sitesOnly.vcf.gz.tbi")

	script:
	"""
	gatk --java-options "-Xmx${task.memory.toMega()}m -XX:ParallelGCThreads=${task.cpus}" \
	MakeSitesOnlyVcf \
        -I $vcfPanel \
        -O ${panelName}_sitesOnly.vcf.gz
	"""
}

process VQSR_SNPs {
	// NOT APLICABLE FOR LOW NUMBER OF SAMPLES
	// minimum 30 wes / 1 wgs
	// https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
	// https://gatk.broadinstitute.org/hc/en-us/articles/360036510892-VariantRecalibrator
	tag "VQSR_SNPs on $panelName using $task.cpus CPUs and $task.memory memory"
	label "gatk"
	label "s_cpu"
	label "s_mem"

	input:
	tuple val(panelName), path(vcf), path(tbi)

	output:
	tuple val(panelName), path("${panelName}_SNPs.recal"), path("${panelName}_SNPs.recal.idx"), path("${panelName}_SNPs.tranches")

	script:
	"""
	gatk --java-options "-Xmx${task.memory.toMega()}m -XX:ParallelGCThreads=${task.cpus}" \
	VariantRecalibrator \
		-R ${params.ref}.fasta \
		-V $vcf \
		--trust-all-polymorphic \
		-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 \
		-tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 \
		-tranche 97.0 -tranche 90.0 \
		-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
		-mode SNP \
		--max-gaussians 6 \
		--resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${params.gatkKnownSites}/hapmap_3.3.hg38.vcf.gz \
		--resource:omni,known=false,training=true,truth=false,prior=12.0 ${params.gatkKnownSites}/1000G_omni2.5.hg38.vcf.gz \
		--resource:1000G,known=false,training=true,truth=false,prior=10.0 ${params.gatkKnownSites}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
		--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.gatkKnownSites}/dbsnp_146.hg38.vcf.gz  \
		-O ${panelName}_SNPs.recal \
		--tranches-file ${panelName}_SNPs.tranches
	"""
}

process VQSR_INDELs {
	// NOT APLICABLE FOR LOW NUMBER OF SAMPLES
	// minimum 30 wes / 1 wgs
	// https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
	// https://gatk.broadinstitute.org/hc/en-us/articles/360036510892-VariantRecalibrator
	tag "VQSR_INDELs on $panelName using $task.cpus CPUs and $task.memory memory"
	// container 'broadinstitute/gatk:4.5.0.0'
	label "gatk"
	label "s_cpu"
	label "s_mem"

	input:
	tuple val(panelName), path(vcf), path(tbi)

	output:
	tuple val(panelName), path("${panelName}_INDELs.recal"), path("${panelName}_INDELs.recal.idx"), path("${panelName}_INDELs.tranches")

	script:
	"""
	gatk --java-options "-Xmx${task.memory.toMega()}m -XX:ParallelGCThreads=${task.cpus}" \
	VariantRecalibrator \
		-R ${params.ref}.fasta \
		-V $vcf \
		-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 \
		-tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 \
		-tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
		-an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
		-mode INDEL \
		--max-gaussians 4 \
		--resource:mills,known=false,training=true,truth=true,prior=12 ${params.gatkKnownSites}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
		--resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${params.gatkKnownSites}/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
		--resource:dbsnp,known=true,training=false,truth=false,prior=2 ${params.gatkKnownSites}/dbsnp_146.hg38.vcf.gz \
		-O ${panelName}_INDELs.recal \
		--tranches-file ${panelName}_INDELs.tranches
	"""
}


process ApplySNPsVQSR {
	// NOT APLICABLE FOR LOW NUMBER OF SAMPLES
	// https://gatk.broadinstitute.org/hc/en-us/articles/360037056912-ApplyVQSR
	tag "ApplySNPsVQSR on $panelName using $task.cpus CPUs and $task.memory memory"
	label "gatk"
	label "s_cpu"
	label "s_mem"

	input:
	tuple val(panelName), path(vcf), path(tbi), path(snpsRecalibration), path(snpsRecalibrationIdx), path(snpsTranches)

	output:
	tuple val(panelName), path("${panelName}.SNP.recalled.vcf.gz")
	//, path("${panelName}.SNP.recalled.vcf.gz.tbi")

	script:
	"""
	gatk --java-options "-Xmx${task.memory.toMega()}m -XX:ParallelGCThreads=${task.cpus}" \
	ApplyVQSR \
		-R ${params.ref}.fasta \
		-V $vcf \
		--recal-file ${snpsRecalibration} \
		--tranches-file ${snpsTranches} \
		--truth-sensitivity-filter-level 99.7 \
		--create-output-variant-index true \
		-mode SNP \
		-O ${panelName}.SNP.recalled.vcf.gz
	"""
}

process ApplyINDELsVQSR {
	// NOT APLICABLE FOR LOW NUMBER OF SAMPLES
	// https://gatk.broadinstitute.org/hc/en-us/articles/360037056912-ApplyVQSR
	tag "ApplyINDELsVQSR on $panelName using $task.cpus CPUs and $task.memory memory"
	label "gatk"
	label "s_cpu"
	label "s_mem"

	input:
	tuple val(panelName), path(vcf), path(tbi), path(indelsRecalibration), path(snpsRecalibrationIdx), path(indelsTranches)

	output:
	tuple val(panelName), path("${panelName}.INDEL.recalled.vcf.gz")
	//, path("${panelName}.INDEL.recalled.vcf.gz.tbi")

	script:
	"""
	gatk --java-options "-Xmx${task.memory.toMega()}m -XX:ParallelGCThreads=${task.cpus}" \
	ApplyVQSR \
		-R ${params.ref}.fasta \
		-V $vcf \
		--recal-file ${indelsRecalibration} \
		--tranches-file ${indelsTranches}\
		--truth-sensitivity-filter-level 99.7 \
		--create-output-variant-index true \
		-mode INDEL \
		-O ${panelName}.INDEL.recalled.vcf.gz
	"""
}


process MergeINDELsSNPsVQSR {
	// NOT APLICABLE FOR LOW NUMBER OF SAMPLES
	// https://gatk.broadinstitute.org/hc/en-us/articles/360037056912-ApplyVQSR
	tag "MergeINDELsSNPsVQSR on $panelName using $task.cpus CPUs and $task.memory memory"
	label "gatk"
	label "s_cpu"
	label "m_mem"

	input:
	tuple val(panelName), path(snpVcf), path(indelVcf)

	output:
	tuple val(panelName), path("${panelName}.mergedVQSR.vcf.gz"), path("${panelName}.mergedVQSR.vcf.gz.tbi")

	script:
	"""
	gatk --java-options -Xmx${task.memory.toMega()}m \
	MergeVcfs \
		--INPUT ${snpVcf} \
		--INPUT ${indelVcf} \
		--OUTPUT ${panelName}.mergedVQSR.vcf.gz

	gatk --java-options -Xmx${task.memory.toMega()}m \
	IndexFeatureFile \
		-I ${panelName}.mergedVQSR.vcf.gz
	"""
}

process VcfToIntervals {
	tag "VcfToIntervals on $panelName: ${intervalList.getSimpleName()} using $task.cpus CPUs and $task.memory memory"
	label "gatk"
	label "s_cpu"
	label "s_mem"

	input:
	tuple val(panelName), path(vcf), path(tbi), path(intervalList)

	output:
	tuple val(panelName), val("${intervalList.getSimpleName()}"), path("${intervalList.getSimpleName()}.${panelName}.VQSRed.vcf.gz")

	script:
	"""
	gatk --java-options -Xmx${task.memory.toMega()}m \
	SelectVariants \
		-R ${params.ref}.fasta \
		-V $vcf \
		-L $intervalList \
		-O ${intervalList.getSimpleName()}.${panelName}.VQSRed.vcf.gz
	"""
}
