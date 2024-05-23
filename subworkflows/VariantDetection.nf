include {
    BQSR;
    GatherRecalibrationTables;
    ApplyBQSR;
    MergeRecalibratedBams;
	GatherGenotypedVCFPanels } from "${params.projectDirectory}/modules/BQSR"

include {
	MakeSitesOnlyVcfs;
	VQSR_SNPs;
	VQSR_INDELs;
	ApplySNPsVQSR;
	ApplyINDELsVQSR;
	MergeINDELsSNPsVQSR;
	VcfToIntervals } from "${params.projectDirectory}/modules/VQSR"

include {
    HaplotypeCallerGVCF;
    BuildDB;
    GenotypeGVCFsIntervals } from "${params.projectDirectory}/modules/VariantCalling"

include {
    FilterPanelSNPs;
    FilterPanelIndels;
    MergeFilteredVariants } from "${params.projectDirectory}/modules/VariantFiltering"

include {
    AnnotateVariantsVEP;
    FilterVEPTranscripts } from "${params.projectDirectory}/modules/Annotation"

include {
    GatherVCFPanels;
    FilterVCFPanels } from "${params.projectDirectory}/modules/VariantPostProcessing"

include {
    VCFsPerSample;
    FilterSampleVCFs;
    CreateTextFile } from "${params.projectDirectory}/modules/SampleProcessing"
    
workflow VariantDetection {
    take:
    dupBams
    panels
    samplesWithPanels
    namedIntervals

	main:
    

	if (params.skipBQSR) {
		finalBams = dupBams
			.map({ [it[0].panelName, it[0], it[1], it[2]] })
			.combine(namedIntervals, by: 0)
		}
		else {
		samplesIntervals = dupBams
			.map({ [it[0].panelName, it[0], it[1], it[2]] })
			.combine(namedIntervals, by: 0)

		recalTables = BQSR(samplesIntervals).groupTuple()
		gatheredTables = GatherRecalibrationTables(recalTables)
		bsqrTablesBams = gatheredTables
			.join(dupBams)
			.map({ [it[0].panelName, it[0], it[1], it[2], it[3]] })
			.combine(namedIntervals, by: 0)

		recalBamsToMerge = ApplyBQSR(bsqrTablesBams).groupTuple()
		finalBams = MergeRecalibratedBams(recalBamsToMerge)
			.map({ [it[0].panelName, it[0], it[1], it[2]] })
			.combine(namedIntervals, by: 0)
		}

    rawVcfs = HaplotypeCallerGVCF(finalBams)
	perPanelRegionVcfs = rawVcfs.groupTuple(by: [0, 1])

	perPanelRegionBedVcfs = namedIntervals
        .map({ [it[0], it[1].getSimpleName(), it[1]] })
        .join(perPanelRegionVcfs, by: [0, 1])

    dbsPanelRegion = BuildDB(perPanelRegionBedVcfs)
	genotypedPanelRegion = GenotypeGVCFsIntervals(dbsPanelRegion)

	if (params.skipVQSR) {
        // use hard filtering
		filteredSNPs = FilterPanelSNPs(genotypedPanelRegion)
		filteredINDELs = FilterPanelIndels(genotypedPanelRegion)
		filteredJoined = filteredSNPs.join(filteredINDELs, by: [0, 1])
		vcfToAnnotate = MergeFilteredVariants(filteredJoined)

	} else {	
        gatheredPanelVcf = GatherGenotypedVCFPanels(genotypedPanelRegion.groupTuple())
        // sitesOnlyPanelVcf = MakeSitesOnlyVcfs(gatheredPanelVcf)
        vqsrSNPs = gatheredPanelVcf.join(VQSR_SNPs(gatheredPanelVcf))
        vqsrINDELs = gatheredPanelVcf.join(VQSR_INDELs(gatheredPanelVcf))

        recalibratedSNPs = ApplySNPsVQSR(vqsrSNPs)//.view()
        recalibratedINDELs = ApplyINDELsVQSR(vqsrINDELs)

        recalibratedVariants = MergeINDELsSNPsVQSR(recalibratedSNPs.join(recalibratedINDELs))

        vcfToAnnotate = VcfToIntervals(recalibratedVariants
            .combine(namedIntervals, by: 0))
	}
	

    intervalVCFsAnnotated = AnnotateVariantsVEP(vcfToAnnotate)
	maneTranscriptsRegion = FilterVEPTranscripts(intervalVCFsAnnotated)

    genotypedPerPanel = GatherVCFPanels(maneTranscriptsRegion.groupTuple())
        .join(panels)
        .map{ [it[2], it[1]] }

	filteredVcfPerPanel = FilterVCFPanels(genotypedPerPanel)

	samplesMergedVCF = samplesWithPanels
        .map { [it.panel, it] }
        .combine(filteredVcfPerPanel)
        .map { [it[1], it[3]] }

	vcfPerSample = VCFsPerSample(samplesMergedVCF)
	filteredSampleVCFs = FilterSampleVCFs(vcfPerSample)
	CreateTextFile(filteredSampleVCFs)

    emit:
	vcfPerSample
}