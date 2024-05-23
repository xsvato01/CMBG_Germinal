process CollectBasecalled {
	tag "CollectBasecalled on $sample.name using $task.cpus CPUs and $task.memory memory"
	container "ubuntu:22.04"
	label "s_cpu"
	label "xxs_mem"

	input:
	val(sample)

	output:
	tuple val(sample), path("*.fastq.gz")

	script:
	"""
	cp  /mnt/share/710000-CEITEC/713000-cmm/713003-pospisilova/base/sequencing_results/primary_data/*${sample.run}/raw_fastq/${sample.name}* ./
	"""
} 

process CreateIntervals {
	tag "CreateIntervals on $panelName using $task.cpus CPUs and $task.memory memory"
	container 'pegi3s/bedtools:2.31.0'
	label "s_cpu"
	label "xxs_mem"

	input:
	tuple val(panelName), val(panel)

	output:
	tuple val(panelName), path("*.intervals")

	script:

	"""
	bedtools merge -d 20000000 -i $panel.bedFile | awk -v FS="\t" -v OFS="\t" '{print \$1 ":" \$2+1 "-" \$3}' > intervals.temp
	while IFS= read -r line; do
		echo "\$line" > "\$(echo "\$line" | sed 's/[:\\-]/_/g').intervals"
	done < intervals.temp
	"""
}

process AlignBams {
	tag "AlignBams on $sample.name using $task.cpus CPUs and $task.memory memory"
	container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:a34558545ae1413d94bde4578787ebef08027945-0"
	publishDir "${params.outDirectory}/mapped/", mode:'copy'
	label "l_cpu"
	label "l_mem"

	input:
	tuple val(sample), path(reads)

	output:
	tuple val(sample), path("${sample.name}.sorted.bam"), path("${sample.name}.sorted.bai")

	script:
	rg = "\"@RG\\tID:${sample.name}\\tSM:${sample.name}\\tLB:${sample.name}\\tPL:ILLUMINA\""
	"""
	bwa mem -R ${rg} -t ${task.cpus} ${params.refindex} $reads \
	| samtools view -Sb -o - - | samtools sort -o ${sample.name}.sorted.bam
	samtools index ${sample.name}.sorted.bam ${sample.name}.sorted.bai	
	"""
}


process MarkDuplicates {
	tag "Mark duplicates on $sample.name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/markDups/", mode:'copy'
	label "gatk"
	label "m_cpu"
	label "l_mem"

	input:
	tuple val(sample), path(bam), path(bai)

	output:
	tuple val(sample), path("${sample.name}.md.bam"), path("${sample.name}.md.bai")
	tuple val(sample), path("${sample.name}.bam.metrics")


	script:
	"""
	gatk --java-options "-Xmx${task.memory.toMega()}m -XX:ParallelGCThreads=${task.cpus}" \
	MarkDuplicates \
		--INPUT $bam \
		--METRICS_FILE ${sample.name}.bam.metrics \
		--ASSUME_SORT_ORDER coordinate \
		--CREATE_INDEX true \
		--OUTPUT ${sample.name}.md.bam
	"""
}