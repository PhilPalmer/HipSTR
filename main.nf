params.genome = "s3://repeat-expansion/Reference/hs37d5.fa"
params.bam = "s3://repeat-expansion/Bams/HG00457.mapped.ILLUMINA.bwa.CHS.exome.20121211.bam"
params.bed= "s3://repeat-expansion/HipSTR/GRCh37.hipstr_reference.bed"
params.depth = "25"

genome_file = file(params.genome)
genome_index = file(params.genome+".fai")
bam_file = file(params.bam)
bai_file = file(params.bam+".bai")
bed_file = file(params.bed)
min = params.minreads

process expansionhunter {
	publishDir 'results'

	input:
	file('aln.bam') from bam_file
	file('aln.bam.bai') from bai_file
	file('genome.fa') from genome_file
	file('genome.fa.fai') from genome_index
	file('hipstr_reference.bed') from bed
	val(min) from min

	output:
	file('output.*') into results
	
	script:
	"""
	HipSTR
	--bams aln.bam \
	--fasta genome.fa \
	--regions hipstr_reference.bed \
	--str-vcf output.vcf.gz \
	--log output.log \
	--viz-out output.viz.gz \
	--min-reads $min \
	--def-stutter-model
	"""
}

workflow.onComplete {
	println ( workflow.success ? "\nHipSTR is done!" : "Oops .. something went wrong" )
}