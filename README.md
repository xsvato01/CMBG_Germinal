# Nextflow GATK4 Germline pipeline for CMBG @2024

**Make sure to get the correct config parameters**. This is a private repository, you will need to create secrets to pull this repository by Nextflow:

append to `~/.nextflow/scm`:
`providers {
    github {
        user = 'UCO'
        password = 'YOUPERSONALTOKEN'
        }
    }`

## If your docker image is private and hosted on Gitlab Muni:

Run docker login registry.gitlab.ics.muni.cz:443 -u <username> -p <token> on an ubuntu machine. Find the ~/.docker/config.json file. Run base64 config.json and copy the output into dockerconfig-secret.yaml (in git repo). Finally run kubectl --namespace [your new namespace] apply -f dockerconfig-secret.yaml to add the secret to the kuba cluster.
This secret has to be added to nextflow.config:
`process {
    pod = [[imagePullSecret:'gitlab-svaton-secret'], ...]
    }`

## Samplesheet

The pipeline takes .json samplesheet in this format:
`{
  "samples": [
    {
      "run": "Next_qJ_MGI",
      "name": "dHL2308qJ"
    }
  ]
}
`

## Optional parameters

The pipeline takes .json samplesheet in this format:
`--skipBQSR` to skip recalibration of .bam files (save considerate amount of resources and time)
`--skipVQSR` to skip recalibrating of .vcf files

@Jan Svaton 2024
