params.index_dir="/mnt/d/Studies/Internship/nextflow_pipeline/ref"
params.ref="chr10.fa"
params.fastq="/mnt/d/Studies/Internship/nextflow_pipeline/fastq/*_{R1,R2}*"

params.bam_dir="/mnt/d/Studies/Internship/nextflow_pipeline/BAM"

process mapping {

publishDir("${params.bam_dir}", mode: 'copy')

input:
 path index_dir
 val ref
 tuple val(sample_id), path(fastq)

output:
 path"*"

script:
"""
bwa mem -t 16 ${index_dir}/${ref} ${fastq} | samtools view -h -b -o ${sample_id}.bam -
"""
}

workflow {

index_ch=Channel.fromPath(params.index_dir)
ref_ch=Channel.of(params.ref)

fastq_ch=Channel.fromFilePairs(params.fastq)

mapping(index_ch,ref_ch,fastq_ch)
mapping.out.view()

}
