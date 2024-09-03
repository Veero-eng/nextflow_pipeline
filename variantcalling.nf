params.ref = "/mnt/d/Studies/Internship/nextflow_pipeline/ref/chr10.fa"
params.sorted_bam_dir="/mnt/d/Studies/Internship/nextflow_pipeline/BAM/sorted"
params.vcf_dir="/mnt/d/Studies/Internship/nextflow_pipeline/VCF"

process variant_calling {

    publishDir("${params.vcf_dir}", mode: 'copy')

    input:
        path sorted_bam

    output:
        path "*"

    script:
        """
        freebayes -f ${params.ref} -i ${sorted_bam} > ${sorted_bam.baseName}.vcf
        """
}

workflow {

    sorted_bam_ch = Channel.fromPath(params.sorted_bam_dir + "/*.sorted.bam")

    variant_calling(sorted_bam_ch)
    variant_calling.out.view()
}
