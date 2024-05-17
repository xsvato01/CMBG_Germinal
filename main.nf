include { CreateIntervals } from "${params.projectDirectory}/modules/Preprocessing.nf"
include { VariantDetection } from "${params.projectDirectory}/subworkflows/VariantDetection"

workflow {
	samples = channel.fromList(params.samples).view()
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

	VariantDetection(panels, samplesWithPanels, namedIntervals)
}

def getPanel(name) {
  if (name[0] == 'd') {
    return "derma"
  } else if (name[0] == 'a') {
    return "atero"
  } else if (name[0] == 'n') {
    return "neuro"
  } else {
    return "Default"
  }
}
