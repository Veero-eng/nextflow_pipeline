params.bam_dir="/mnt/d/Studies/Internship/nextflow_pipeline/BAM"
params.sorted_bam_dir="/mnt/d/Studies/Internship/nextflow_pipeline/BAM/sorted"

process sorting {

    publishDir("${params.sorted_bam_dir}", mode: 'copy')

    input:
        path bam

    output:
        path "*"

    script:
        """
        samtools sort -o ${bam.baseName}.sorted.bam ${bam}
        samtools index ${bam.baseName}.sorted.bam
        """
}

workflow {

    bam_ch = Channel.fromPath(params.bam_dir + "/*.bam")

    sorting(bam_ch)
    sorting.out.view()
}
