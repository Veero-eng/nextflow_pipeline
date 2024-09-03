params.ref="/mnt/d/Studies/Internship/nextflow_pipeline/ref/chr10.fa"
params.index_dir="/mnt/d/Studies/Internship/nextflow_pipeline/ref"

process index {

publishDir("${params.index_dir}", mode: 'copy')

input:
 path genome

output:
path "*"

script:
"""
bwa index $genome

"""

}

workflow {

ref_ch=Channel.fromPath(params.ref)

index(ref_ch)


}
