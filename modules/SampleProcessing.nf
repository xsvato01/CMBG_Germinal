
process VCFsPerSample {
	tag "VCFsPerSample on $sample.name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/VCFsPerSample/", mode:'copy'
	label "gatk"
	label "s_cpu"
	label "m_mem"
	
	input:
	tuple val(sample), path(mergedVcf)

	output:
	tuple val(sample), path("${sample.name}.vcf.gz")

	script:
	"""
	gatk --java-options -Xmx${task.memory.toMega()}m \
	IndexFeatureFile \
		-I $mergedVcf

	gatk --java-options -Xmx${task.memory.toMega()}m \
	SelectVariants \
		-R ${params.ref}.fasta \
		-V ${mergedVcf} \
		-sn $sample.name \
		-O ${sample.name}.vcf.gz
	"""
}


process FilterSampleVCFs {
	tag "FilterSampleVCFs_INDELs on $sample.name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/VCFsPerSample/", mode:'copy'
	container 'staphb/bcftools:1.20'
	label "s_cpu"
	label "xxs_mem"
	
	input:
	tuple val(sample), path(vcf)

	output:
	tuple val(sample), path("${sample.name}.filtered.vcf.gz")

	script:
	"""
	bcftools view -i 'AC>0 && GQ>=20 && FILTER="PASS"' $vcf -o ${sample.name}.filtered.vcf.gz
	"""
}


process CreateTextFile {
	tag "CreateTextFile on $sample.name using $task.cpus CPUs and $task.memory memory"
	container 'registry.gitlab.ics.muni.cz:443/450402/btk_k8s:16'
	publishDir "${params.outDirectory}/txt/", mode:'copy'
	label "s_cpu"
	label "xs_mem"
	
	input:
	tuple val(sample), path(vcfgz)

	output:
	tuple val(sample), path("${sample.name}.csv")
	
	script:
	"""
	bcftools view $vcfgz > temp.vcf
	python ${params.vcf2txt} simple --build GRCh38 -i temp.vcf -t ${sample.name} -o ${sample.name}.csv  
	"""
}