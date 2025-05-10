params.strelka_somatic_dir = "/mnt/d/Studies/Internship/nextflow_pipeline/somatic_vc/strelka_analysis"

process run_dir {
    input:
    path strelka_somatic_dir
    output:
    path strelka_somatic_dir
    script:
    """
    /mnt/d/Studies/Internship/nextflow_pipeline/somatic_vc/strelka_analysis/runWorkflow.py -m local -j 8
    """
}

workflow {
run_dir(params.strelka_somatic_dir)
}
