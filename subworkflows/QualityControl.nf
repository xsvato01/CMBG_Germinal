include {
    FastQC;
    BamQC;
    MultiQC } from "${params.projectDirectory}/modules/QC"

    
workflow QualityControl {
    take:
    rawFastq
    dupBams
    dupMetrics


	main:
    fastQCs = FastQC(rawFastq)
    bamQCs = BamQC(dupBams)
    joinedQCs = bamQCs
    .join(fastQCs)
    .join(dupMetrics)
    // .map{[it[0].run, it[1].collect(), it[2].collect(), it[3]]}
        .map{[it[0].run, [it[1], it[2], it[3]]]}
        // .view()
    .groupTuple()
    .map{[it[0], it[1].flatten()]}

    // .collect()

    MultiQC(joinedQCs)


}