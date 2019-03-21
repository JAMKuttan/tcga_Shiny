#!/usr/bin/env nextflow

// Path to an input file, or a pattern for multiple inputs
// Note - $baseDir is the location of this workflow file main.nf

// Define Input variables
params.run = "TRUE"

process copy {

  publishDir "$baseDir/output", mode: 'copy'

  input:

  params.run

  output:

  file("**/Clinical/**") into outPathClinical
  file("**/Expression/**") into outPathExpression

  script:

  """
  cp -r /project/shared/bicf_workflow_ref/tcga_shiny_app_data/* .
  """
}
