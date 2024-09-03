params.vcf_dir="/mnt/d/Studies/Internship/nextflow_pipeline/VCF"
params.filtered_vcf_dir="/mnt/d/Studies/Internship/nextflow_pipeline/VCF/filtered"

process variant_filtering {

    publishDir("${params.filtered_vcf_dir}", mode: 'copy')

    input:
        path vcf

    output:
        path "*"

    script:
        """
        bcftools filter -i 'QUAL>20 && INFO/DP>10' -o ${vcf.baseName}.filtered.vcf ${vcf}
        """
}

workflow {

    vcf_ch = Channel.fromPath(params.vcf_dir + "/*.vcf")

    variant_filtering(vcf_ch)
    variant_filtering.out.view()
}
