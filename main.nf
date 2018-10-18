params.genome = "s3://repeat-expansion/Reference/hs37d5.fa"
params.bam = "s3://repeat-expansion/Bams/HG00457.mapped.ILLUMINA.bwa.CHS.exome.20121211.bam"
params.bed= "s3://repeat-expansion/HipSTR/GRCh37.hipstr_reference.bed"

genome_file = file(params.genome)
genome_index = file(params.genome+".fai")
// bam = file(params.bam)
bai_file = file(params.bam+".bai")
bed_file = file(params.bed)

// Params for the Read Group Line to be added just in case its needed.
params.rgid=4
params.rglb="lib1"
params.rgpl="illumina"
params.rgpu="unit1"
params.rgsm=20

Channel
		.fromPath(params.bam)
		.ifEmpty { exit 1, "${params.bam} not found"}
		.set{bamChannel}

process preprocess_bam{

  tag "${bam}"
	container 'lifebitai/samtools'

  input:
  file(bam) from bamChannel

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

process hipstr {
	container 'lifebitai/hipstr'

	publishDir 'results'

	input:
	set file(bam), file(bai) from completeChannel
	file('genome.fa') from genome_file
	file('genome.fa.fai') from genome_index
	file('hipstr_reference.bed') from bed_file

	output:
	file('output.*') into results

	script:
	"""
	HipSTR \
	--bams ${bam} \
	--fasta genome.fa \
	--regions hipstr_reference.bed \
	--str-vcf output.vcf.gz \
	--log output.log \
	--viz-out output.viz.gz \
	"""
}

workflow.onComplete {
	println ( workflow.success ? "\nHipSTR is done!" : "Oops .. something went wrong" )
}
