params.filtered_vcf_dir="/mnt/d/Studies/Internship/nextflow_pipeline/VCF/filtered"
params.annotated_vcf_dir="/mnt/d/Studies/Internship/nextflow_pipeline/VCF/annotated"
params.snpEff_db = 'GRCh37.75'

process variant_annotation {

    publishDir("${params.annotated_vcf_dir}", mode: 'copy')

    input:
        path filtered_vcf

    output:
        path "*"

    script:
        """
        snpEff ${params.snpEff_db} ${filtered_vcf} > ${filtered_vcf.baseName}.annotated.vcf 
        """
}

workflow {

    filtered_vcf_ch = Channel.fromPath(params.filtered_vcf_dir + "/*.filtered.vcf")

    variant_annotation(filtered_vcf_ch)
    variant_annotation.out.view()
}
