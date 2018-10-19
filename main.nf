#!/usr/bin/env nextflow

/*
 * SET UP CONFIGURATION VARIABLES
 */
bam = Channel
		.fromPath(params.bam)
		.ifEmpty { exit 1, "${params.bam} not found.\nPlease specify --bam option (--bed bamfile)"}

// if(params.bai) {
// 	bai = Channel
// 		.fromPath(params.bai)
// 		.ifEmpty { exit 1, "${params.bai} not found.\nYou can remove your --bai option (--bai baifile) and it will be produced for you automatically"}
// }

Channel
		.fromPath(params.genome)
		.ifEmpty { exit 1, "${params.genome} not found.\nPlease specify --genome option (--genome fastafile)"}
		.into { fastaToFai; fastaToHipSTR }

if(params.fai) {
	fai = Channel
			.fromPath(params.fai)
			.ifEmpty { exit 1, "${params.fai} not found.\nYou can remove your --fai option (--fai faifile) and it will be produced for you automatically"}
}

bed = Channel
    .fromPath(params.bed)
    .ifEmpty { exit 1, "${params.bed} not found.\nPlease specify --bed option (--bed bedfile)"}

		// Header log info
log.info """=======================================================
		HipSTR
======================================================="""
def summary = [:]
summary['Pipeline Name']    = 'HipSTR'
summary['Bam file']         = params.bam
summary['Bed file']         = params.bed
summary['Reference genome'] = params.genome
if(params.fai) summary['Fasta Index'] = params.fai
summary['Output dir']       = params.outdir
summary['Working dir']      = workflow.workDir
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


process preprocess_bam{

  tag "${bam}"
	container 'lifebitai/samtools'

  input:
  file bam from bam

  output:
  set file("ready/${bam}"), file("ready/${bam}.bai") into completeChannel

  script:
  """
  mkdir ready
  [[ `samtools view -H ${bam} | grep '@RG' | wc -l`   > 0 ]] && { mv $bam ready;}|| { java -jar /picard.jar AddOrReplaceReadGroups \
  I=${bam} \
  O=ready/${bam} \
  RGID=${params.rgid} \
  RGLB=${params.rglb} \
  RGPL=${params.rgpl} \
  RGPU=${params.rgpu} \
  RGSM=${params.rgsm};}
  cd ready ;samtools index ${bam};
  """
}

if(!params.fai) {
  process preprocess_fai {

      input:
      file fasta from fastaToFai

      output:
      file("${fasta}.fai") into fai

      script:
      """
      samtools faidx $fasta
      """
  }
}

process hipstr {
	container 'lifebitai/hipstr'

	publishDir "${params.outdir}", mode: 'copy'

	input:
	set file(bam), file(bai) from completeChannel
	file fasta from fastaToHipSTR
	file fai from fai
	file bed from bed

	output:
	file('output.*') into results

	script:
	"""
	HipSTR \
	--bams ${bam} \
	--fasta ${fasta} \
	--regions ${bed} \
	--str-vcf output.vcf.gz \
	--log output.log \
	--viz-out output.viz.gz \
	"""
}

workflow.onComplete {
	println ( workflow.success ? "\nHipSTR is done!" : "Oops .. something went wrong" )
}
