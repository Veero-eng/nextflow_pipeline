params.ref = "/mnt/d/Studies/Internship/nextflow_pipeline/ref/chr10.fa"
params.normal_bam = "/mnt/d/Studies/Internship/nextflow_pipeline/somatic_vc/mutect2/normal_aln_sort_rg.bam"
params.tumor_bam = "/mnt/d/Studies/Internship/nextflow_pipeline/somatic_vc/mutect2/tumor_aln_sort_rg.bam"
params.strelka = "/mnt/d/Studies/Internship/nextflow_pipeline/somatic_vc/strelka/install/bin/configureStrelkaSomaticWorkflow.py"
params.output_dir = "/mnt/d/Studies/Internship/nextflow_pipeline/somatic_vc"

process strelka {
    publishDir("${params.output_dir}", mode: 'copy')
    input:
    path strelka
    path normal_bam
    path tumor_bam
    path ref
    output:
    path "*"
    script:
    """
    ${params.strelka} --normalBam ${params.normal_bam} --tumorBam ${params.tumor_bam} --referenceFasta ${params.ref} --runDir strelka_analysis
    """
}

workflow {
    strelka_ch = Channel.fromPath(params.strelka)
    normal_bam_ch = Channel.fromPath(params.normal_bam)
    tumor_bam_ch = Channel.fromPath(params.tumor_bam)
    ref_ch = Channel.fromPath(params.ref)

    strelka(strelka_ch, normal_bam_ch, tumor_bam_ch, ref_ch)
    strelka.out.view()
}