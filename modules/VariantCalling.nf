
process HaplotypeCallerGVCF {
	tag "HaplotypeCallerGVCF on $sample.name-${intervalBed.getSimpleName()} using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/VCFs/", mode:'copy'
	container 'broadinstitute/gatk:4.2.3.0'
	label "m_cpu"
	label "m_mem"

	input:
	tuple val(panelName), val(sample), path(bam), path(bai), path(intervalBed)

	output:
	tuple val(panelName), val("${intervalBed.getSimpleName()}"), path("$sample.name-${intervalBed.getSimpleName()}.g.vcf.gz"), path("$sample.name-${intervalBed.getSimpleName()}.g.vcf.gz.tbi")

	script:
	"""
	gatk --java-options "-Xmx${task.memory.toMega()}m -XX:ParallelGCThreads=${task.cpus}" \
	HaplotypeCaller  \
		-R ${params.ref}.fasta \
		-I $bam \
		-L ${intervalBed} \
		-O $sample.name-${intervalBed.getSimpleName()}.g.vcf.gz \
		--create-output-variant-index \
		-ERC GVCF
	"""	
}

process BuildDB {
	// build DB per panel type (its respective samples) and interval
	tag "BuildDB on $panelName:${intervalBed.getSimpleName()} using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/DBs/", mode:'copy'
	container 'broadinstitute/gatk:4.2.3.0'
	label "s_cpu"
	label "s_mem"

	input:
	tuple val(panelName), val(intervalName), path(intervalBed), path(vcfsSamples), path(tbiSamples)

	output:
	tuple val(panelName), val(intervalName), path("${panelName}-${intervalBed.getSimpleName()}.db")

	script:
	input = vcfsSamples.collect{"-V ${it}"}.join(' ')
	"""
	gatk --java-options "-Xmx${task.memory.toMega()}m -XX:ParallelGCThreads=${task.cpus}" \
	GenomicsDBImport \
		$input \
		--genomicsdb-workspace-path ${panelName}-${intervalBed.getSimpleName()}.db \
		-L $intervalBed
	"""
}

process GenotypeGVCFsIntervals {
	tag "GenotypeGVCFsIntervals on $panelName: $intervalName using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/genotyped/", mode:'copy'
	container 'broadinstitute/gatk:4.2.3.0'
	label "s_cpu"
	label "m_mem"
	
	input:
	tuple val(panelName), val(intervalName), path(db)

	output:
	tuple val(panelName), val(intervalName), path("${panelName}_${intervalName}.vcf.gz"), path("${panelName}_${intervalName}.vcf.gz.tbi")


	script:
	"""
	gatk --java-options "-Xmx${task.memory.toMega()}m -XX:ParallelGCThreads=${task.cpus}" \
	GenotypeGVCFs \
		-R ${params.ref}.fasta \
		-V gendb://${db} \
		-O ${panelName}_${intervalName}.vcf.gz
	"""
}


// process BEAGLE_IMPUTATION {
// 	tag "BEAGLE_IMPUTATION on $panelName:$intervalName using $task.cpus CPUs and $task.memory memory"
// 	publishDir "${params.outDirectory}/phased/", mode:'copy'
// 	container 'quay.io/biocontainers/beagle:5.4_22Jul22.46e--hdfd78af_0'
// 	label "s_cpu"
// 	label "m_mem"
	
// 	input:
// 	tuple val(panelName), val(intervalName), path(vcfInterval)

// 	output:
// 	tuple val(panelName), path("*")

// 	script:
// 	splited = intervalName.split("_")
// 	chrom = splited[0]
// 	"""
// 	echo $chrom
// 	beagle \
// 		gt=$vcfInterval \
// 		out=${panelName}_${intervalName}_phased
// 		map=${params.plinkMap}/plink.${chrom}.GRCh38.map \
// 		chrom=$chrom:${splited[1]}-${splited[2]} \
// 		nthreads=$task.cpus
// 	"""
// }