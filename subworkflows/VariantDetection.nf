include {
    CollectBasecalled;
    AlignBams;
    MarkDuplicates } from "${params.projectDirectory}/modules/Preprocessing.nf"

include {
    BaseRecalibration;
    GatherRecalibrationTables;
    ApplyRecalibration;
    MergeRecalibratedBams } from "${params.projectDirectory}/modules/Recalibration.nf"

include {
    HaplotypeCallerGVCF;
    BuildDB;
    GenotypeGVCFsIntervals } from "${params.projectDirectory}/modules/VariantCalling.nf"

include {
    FilterPanelSNPs;
    FilterPanelIndels;
    MergeFilteredVariants } from "${params.projectDirectory}/modules/VariantFiltering.nf"

include {
    AnnotateVariantsVEP;
    FilterVEPTranscripts } from "${params.projectDirectory}/modules/Annotation.nf"

include {
    GatherVCFPanels;
    FilterVCFPanels } from "${params.projectDirectory}/modules/VariantPostProcessing.nf"

include {
    VCFsPerSample;
    FilterSampleVCFs;
    CreateTextFile } from "${params.projectDirectory}/modules/SampleProcessing.nf"
    
workflow VariantDetection {
    take:
    panels
    samplesWithPanels
    namedIntervals

	main:
    rawFastq = CollectBasecalled(samplesWithPanels)
    sortedBams = AlignBams(rawFastq)
    (dupBams, markTable) = MarkDuplicates(sortedBams)
    samplesIntervals = dupBams
        .map({ [it[0].panel, it[0], it[1], it[2]] })
        .combine(namedIntervals, by: 0)

    recalTables = BaseRecalibration(samplesIntervals).groupTuple()
	gatheredTables = GatherRecalibrationTables(recalTables)
	bsqrTablesBams = gatheredTables
        .join(dupBams)
        .map({ [it[0].panel, it[0], it[1], it[2], it[3]] })
        .combine(namedIntervals, by: 0)

    recalBamsToMerge = ApplyRecalibration(bsqrTablesBams).groupTuple()
	mergedRecalBams = MergeRecalibratedBams(recalBamsToMerge)
        .map({ [it[0].panel, it[0], it[1], it[2]] })
        .combine(namedIntervals, by: 0)

    rawVcfs = HaplotypeCallerGVCF(mergedRecalBams)
	perPanelRegionVcfs = rawVcfs.groupTuple(by: [0, 1])

	perPanelRegionBedVcfs = namedIntervals
        .map({ [it[0], it[1].getSimpleName(), it[1]] })
        .join(perPanelRegionVcfs, by: [0, 1])

    dbsPanelRegion = BuildDB(perPanelRegionBedVcfs)
	genotypedPanelRegion = GenotypeGVCFsIntervals(dbsPanelRegion)

	filteredSNPs = FilterPanelSNPs(genotypedPanelRegion)
	filteredINDELs = FilterPanelIndels(genotypedPanelRegion)
    filteredJoined = filteredSNPs.join(filteredINDELs, by: [0, 1])
	mergedFiltered = MergeFilteredVariants(filteredJoined)

    intervalVCFsAnnotated = AnnotateVariantsVEP(mergedFiltered)
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

}