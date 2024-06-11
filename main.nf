include { CreateIntervals } from "${params.projectDirectory}/modules/Preprocessing"
include {
    CollectBasecalled;
    AlignBams;
    MarkDuplicates } from "${params.projectDirectory}/modules/Preprocessing"
include { VariantDetection } from "${params.projectDirectory}/subworkflows/VariantDetection"
include { QualityControl } from "${params.projectDirectory}/subworkflows/QualityControl"


workflow {
	samples = channel.fromList(params.samples)
	panels = channel.fromList(params.panels).map({ [it.panelName, it] })
	samplesPanels = samples.map({ sample -> [getPanel(sample.name), sample] })
	samplesWithPanels = samplesPanels
		.combine(panels, by: 0) //match parameters to individual samples
		.map({ sample -> sample[1] + (sample[2]) }) //concatenate sample info with panel info
	panelsUnique = samplesWithPanels
		.map { [it.panelName, it] }
		.groupTuple()
		.map({ [it[0], it[1][0]] })

	intervals = CreateIntervals(panelsUnique)
	namedIntervals = intervals
		.flatMap { tuple -> tuple[1].collect { val -> [tuple[0], val] } }
	
	namedIntervals
	.map { it-> it[1] }
	.view()
	.count()
    .view{"____________Intervals count_____________: $it"}

	rawFastq = CollectBasecalled(samplesWithPanels)
    sortedBams = AlignBams(rawFastq)
    (dupBams, dupMetrics) = MarkDuplicates(sortedBams)

	// VariantDetection(dupBams, panels, samplesWithPanels, namedIntervals)

	QualityControl(rawFastq, dupBams, dupMetrics)

	VariantDetection(dupBams, panels, samplesWithPanels, namedIntervals)
}

def getPanel(name) {
	if (name[0] == 'd') {
		return "derma"
	} else if (name[0] == 'a') {
		return "atero"
	} else if (name[0] == 'n') {
		return "neuro"
	} else if (name[0] == 'x') {
		return "exom"
	} else {
		return "wgs"
	}
}
