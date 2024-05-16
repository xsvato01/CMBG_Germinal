process REFORMAT_PARAMS {
	tag "Reformating $references.panel using $task.cpus CPUs $task.memory"
	container "ubuntu:22.04"
	label "s_cpu"
	label "xxs_mem"

	input:
	val references

	output:
	tuple val(references), val(references.name)

	
	""" """ //this is not an error
}



process REFORMAT_SAMPLE {
	tag "Reformating $sample.name using $task.cpus CPUs $task.memory"
	container "ubuntu:22.04"
	label "s_cpu"
	label "xxs_mem"

	input:
	val sample

	output:
		tuple val(sample.name),val(sample)

	""" """ //this is not an error


}


process COLLECT_BASECALLED {
	tag "COLLECT_BASECALLED on $sample.name using $task.cpus CPUs and $task.memory memory"
	container "ubuntu:22.04"
	label "s_cpu"
	label "xxs_mem"

	input:
	val(sample)

	output:
	tuple val(sample), path("*.fastq.gz")

	script:
	"""
	echo COLLECT_BASECALLED $sample.name
	cp  /mnt/share/710000-CEITEC/713000-cmm/713003-pospisilova/base/sequencing_results/primary_data/*${sample.run}/raw_fastq/${sample.name}* ./
	"""
} 


process TRIMMING {
	tag "trimming on $sample.name using $task.cpus CPUs and $task.memory memory"
	
	input:
	tuple val(sample), path(reads)

	output:
	tuple val(sample), path("*.fastq.gz")

	script:
	""" 
	cutadapt -m 50 -o ${sample.name}.trimmed.R1.fastq.gz -p ${sample.name}.trimmed.R2.fastq.gz $reads
	"""
}

process ALIGN {
	tag "ALIGN on $sample.name using $task.cpus CPUs and $task.memory memory"
	container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:a34558545ae1413d94bde4578787ebef08027945-0"
	publishDir "${params.outDirectory}/mapped/", mode:'copy'
	label "m_cpu"
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

// process FIRST_QC {
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

process MARK_DUPLICATES {
	tag "Mark duplicates on $sample.name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/markDups/", mode:'copy'
	container 'broadinstitute/gatk:4.2.3.0'
	label "s_cpu"
	label "l_mem"

	input:
	tuple val(sample), path(bam), path(bai)

	output:
	tuple val(sample), path("${sample.name}.md.bam"), path("${sample.name}.md.bai")
	path "${sample.name}.bam.metrics"

	script:
	"""
	gatk --java-options -Xmx${task.memory.toGiga()}g \
	MarkDuplicates \
		--INPUT $bam \
		--METRICS_FILE ${sample.name}.bam.metrics \
		--ASSUME_SORT_ORDER coordinate \
		--CREATE_INDEX true \
		--OUTPUT ${sample.name}.md.bam
	"""
}


process CREATE_INTERVALS {
	tag "CREATE_INTERVALS on $panelName using $task.cpus CPUs and $task.memory memory"
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

process BASE_RECALIBRATOR {
	tag "BASE_RECALIBRATOR on $sample.name-${intervalBed.getSimpleName()} using $task.cpus CPUs and $task.memory memory"
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


process GATHER_BSQR {
	tag "GATHER_BSQR on $sample.name using $task.cpus CPUs and $task.memory memory"
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

process APPLY_BSQR {
	tag "APPLY_BSQR on $sample.name-${intervalBed.getSimpleName()} using $task.cpus CPUs and $task.memory memory"
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


process MERGE_BAMRECALLS {
	tag "MERGE_BAMRECALLS on $sample.name using $task.cpus CPUs and $task.memory memory"
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

process HAPLOTYPECALLER_GVCF {
	tag "HAPLOTYPECALLER_GVCF on $sample.name-${intervalBed.getSimpleName()} using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/VCFs/", mode:'copy'
	container 'broadinstitute/gatk:4.2.3.0'
	label "s_cpu"
	label "m_mem"

	input:
	tuple val(panelName), val(sample), path(bam), path(bai), path(intervalBed)

	output:
	tuple val(panelName), val("${intervalBed.getSimpleName()}"), path("$sample.name-${intervalBed.getSimpleName()}.g.vcf.gz"), path("$sample.name-${intervalBed.getSimpleName()}.g.vcf.gz.tbi")
	path "*"

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

process DB_IMPORT {
	// build DB per panel type (its respective samples) and interval
	tag "DB_IMPORT on $panelName:${intervalBed.getSimpleName()} using $task.cpus CPUs and $task.memory memory"
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

process GENOTYPE_GVCFs_INTERVALS {
	tag "GENOTYPE_GVCFs_INTERVALS on $panelName: $intervalName using $task.cpus CPUs and $task.memory memory"
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



process FILTER_PANEL_SNPs {
	// this is applicable for smaller panels
	tag "FILTER_PANEL_SNPs on $panelName: $intervalName using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/FilteredSNPs/", mode:'copy'
	container 'broadinstitute/gatk:4.2.3.0'
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

process FILTER_PANEL_INDELs {
	// this is applicable for smaller panels
	tag "FILTER_PANEL_INDELs on $panelName: $intervalName using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/FilteredINDELs/", mode:'copy'
	container 'broadinstitute/gatk:4.2.3.0'
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


process MERGE_FILTERED {
	// this is applicable for smaller panels
	tag "MERGE_FILTERED on $panelName: $intervalName using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/FilteredMerged/", mode:'copy'
	container 'broadinstitute/gatk:4.2.3.0'
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

process VAR_RECALL {
	// NOT APLICABLE FOR LOW NUMBER OF SAMPLES
	// https://gatk.broadinstitute.org/hc/en-us/articles/360036510892-VariantRecalibrator
	tag "VAR_RECALL on $panelName: $intervalName using $task.cpus CPUs and $task.memory memory"
	container 'broadinstitute/gatk:4.2.3.0'
	label "s_cpu"
	label "m_mem"

	input:
	tuple val(panelName), val(intervalName), path(genotypedVcfInterval), path(tbi)

	output:
	tuple val(panelName), val(intervalName), path("${intervalName}_${panelName}*")

	script:
	"""
	gatk --java-options "-Xmx${task.memory.toMega()}m -XX:ParallelGCThreads=${task.cpus}" \
	VariantRecalibrator \
		-R ${params.ref}.fasta \
		-V ${genotypedVcfInterval} \
		--resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${params.gatkKnownSites}/hapmap_3.3.hg38.vcf.gz \
		--resource:omni,known=false,training=true,truth=false,prior=12.0 ${params.gatkKnownSites}/1000G_omni2.5.hg38.vcf.gz \
		--resource:1000G,known=false,training=true,truth=false,prior=10.0 ${params.gatkKnownSites}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
		--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.gatkKnownSites}/dbsnp_146.hg38.vcf.gz  \
		-an QD -an ReadPosRankSum -an FS -an SOR -an DP \
		-mode SNP \
		--max-gaussians 1 \
		-O ${intervalName}_${panelName}.recal \
		--tranches-file ${intervalName}_${panelName}.tranches \
		--rscript-file ${intervalName}_${panelName}.plots.R
	"""
}


process APPLY_VQSR {
	// NOT APLICABLE FOR LOW NUMBER OF SAMPLES
	// https://gatk.broadinstitute.org/hc/en-us/articles/360037056912-ApplyVQSR
	tag "APPLY_VQSR on $panelName: $intervalName using $task.cpus CPUs and $task.memory memory"
	container 'broadinstitute/gatk:4.2.3.0'
	label "s_cpu"
	label "m_mem"

	input:
	tuple val(panelName), val(intervalName), path(var_recall_files)

	output:
	tuple val(panelName), val(intervalName), path("${intervalName}_${panelName}*")

	script:
	"""
	gatk --java-options "-Xmx${task.memory.toMega()}m -XX:ParallelGCThreads=${task.cpus}" \
	ApplyVQSR \
		-R ${params.ref}.fasta \
		-V ${intervalName}_${panelName}.vcf.gz \
		-O ${intervalName}_${panelName}.recalled.vcf.gz  \
		--ts_filter_level 99.0 \
		--tranches-file ${intervalName}_${panelName}.tranches \
		--recal-file ${intervalName}_${panelName}.recal \
		-mode SNP
	"""
}

process ANNOTATE_VEP {
	tag "ANNOTATE_VEP on $panelName: $intervalName using $task.cpus CPUs $task.memory"
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

process FILTER_VEP_TRANSCRIPTS {
	tag "FILTER_VEP_TRANSCRIPTS on $panelName: $intervalName using $task.cpus CPUs $task.memory"
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

process GATHER_VCF_PANELS {
	// collapse files from individual intervals into one file per panel
	tag "GATHER_VCF_PANELS on $panelName using $task.cpus CPUs and $task.memory memory"
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

process FILTER_VCF_PANELS {
	// filter variants for bedfile for respective panel
	tag "FILTER_VCF_PANELS on $panel.panel using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/mergedVCFs/", mode:'copy'
	container 'registry.gitlab.ics.muni.cz:443/450402/btk_k8s:16'	

	label "s_cpu"
	label "m_mem"
	
	input:
	tuple val(panel), path(mergedVcf)

	output:
	tuple val(panel), path("${panel.panel}.filtered.vcf.gz")

	script:
	"""
	bedtools intersect -header -a $mergedVcf -b ${panel.bedFile} | bgzip > ${panel.panel}.filtered.vcf.gz
	"""
}


process VCFs_PER_SAMPLE {
	tag "VCFs_PER_SAMPLE on $sample.name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/VCFsPerSample/", mode:'copy'
	container 'broadinstitute/gatk:4.2.3.0'
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


process FILTER_SAMPLE {
	tag "FILTER_SAMPLE_INDELs on $sample.name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/VCFsPerSample/", mode:'copy'
	container 'staphb/bcftools:1.20'
	label "s_cpu"
	label "s_mem"
	
	input:
	tuple val(sample), path(vcf)

	output:
	tuple val(sample), path("${sample.name}.filtered.vcf.gz")

	script:
	"""
	bcftools view -i 'AC>0 && GQ>=20' $vcf -o ${sample.name}.filtered.vcf.gz
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

process CREATE_TXT {
	tag "CREATE_TXT on $sample.name using $task.cpus CPUs and $task.memory memory"
	container 'registry.gitlab.ics.muni.cz:443/450402/btk_k8s:16'
	publishDir "${params.outDirectory}/txt/", mode:'copy'
	label "s_cpu"
	label "m_mem"
	
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

workflow {

Samples = channel.fromList(params.samples)
Panels = channel.fromList(params.panels).map({[it.panel,it]})

SamplesPanels = Samples.map({sample -> [ getPanel(sample.name), sample ]})

SamplesWithPanels = SamplesPanels.combine(Panels, by: 0) //match parameters to individual samples
		.map({ sample -> sample[1]+(sample[2]) }) //concatenate sample info with panel info


panelsUnique = SamplesWithPanels.map{[it.panel,it]}.groupTuple().map({[it[0],it[1][0]]}) //return just one panel info per tuple
intervals = CREATE_INTERVALS(panelsUnique)
// intervals.map({it -> it[0].map({itNested -> return[it[1], itnNested]})}).view()
NamedIntervals = intervals.flatMap { tuple ->
	tuple[1].collect { val -> [tuple[0], val] } //explode individual intervals to get [panelName, interval.list]
}

rawfastq = COLLECT_BASECALLED(SamplesWithPanels)
sortedBams	= ALIGN(rawfastq)
(dupBams, markTable) = MARK_DUPLICATES(sortedBams)

SamplesIntervals = dupBams
	.map({[it[0].panel, it[0], it[1], it[2]]})
	.combine(NamedIntervals, by: 0)

RecalTables = BASE_RECALIBRATOR(SamplesIntervals)
	.groupTuple() //[sample, [inverval1, ...invervalX]]
GatheredTables = GATHER_BSQR(RecalTables)
BSQRTablesBams = GatheredTables.join(dupBams)
	.map({[it[0].panel, it[0], it[1], it[2], it[3]]})
	.combine(NamedIntervals, by: 0) //join with markDup bams and back to intervals
RecalBamsToMerge = APPLY_BSQR(BSQRTablesBams).groupTuple()
MergedRecalBams = MERGE_BAMRECALLS(RecalBamsToMerge)
	.map({[it[0].panel, it[0], it[1], it[2]]})
	.combine(NamedIntervals, by: 0)

(RawVcfs, _) = HAPLOTYPECALLER_GVCF(MergedRecalBams) //haplotype per sample and interval
PerPanelRegionVcfs = RawVcfs.groupTuple(by:[0,1]) //group samples for same panel and interval
// DBsPanelRegion = DB_IMPORT(PerPanelRegionVcfs)

PerPanelRegionBedVcfs = NamedIntervals
	.map({[it[0], it[1].getSimpleName(), it[1]]})
	.join(PerPanelRegionVcfs, by: [0,1]) //match panel name and interval name with interval.list
DBsPanelRegion = DB_IMPORT(PerPanelRegionBedVcfs)

GenotypedPanelRegion = GENOTYPE_GVCFs_INTERVALS(DBsPanelRegion)
FilteredSNPs = FILTER_PANEL_SNPs(GenotypedPanelRegion)
FilteredINDELs = FILTER_PANEL_INDELs(GenotypedPanelRegion)

FilteredJoined = FilteredSNPs.join(FilteredINDELs, by:[0,1])
MergedFiltered = MERGE_FILTERED(FilteredJoined)

IntervalVCFsAnnotated = ANNOTATE_VEP(MergedFiltered)
ManeTranscriptsRegion = FILTER_VEP_TRANSCRIPTS(IntervalVCFsAnnotated)

GenotypedPerPanel = GATHER_VCF_PANELS(ManeTranscriptsRegion.groupTuple())
	.join(Panels)
	.map{[it[2], it[1]]}

// // 	// BEAGLE_IMPUTATION(GenotypedPanelRegion.first())
FilteredVcfPerPanel = FILTER_VCF_PANELS(GenotypedPerPanel)

SamplesMergedVCF = SamplesWithPanels //get [sample, panel.vcf.gz, panel.vcf.gz.tbi]
	.map{[it.panel,it]}
	.combine(FilteredVcfPerPanel)
	.map{[it[1], it[3]]}

VcfPerSample = VCFs_PER_SAMPLE(SamplesMergedVCF)
FilteredSampleVCFs = FILTER_SAMPLE(VcfPerSample)
CREATE_TXT(FilteredSampleVCFs)
}


def getPanel(name) {
	if (name[0] == 'd') {
		return "derma"
	} else if (name[0] == 'a') {
		return "atero"
	} else {
		return "Default"
	}
}

def removeSame(nm,lst) {
	def list=[]
	for (int i = 0; i < lst.size(); i++) {
		if (lst[i] != nm) {
			list.add(lst[i])
			}
		}
	return(list)
}