params.ref = "/mnt/d/Studies/Internship/nextflow_pipeline/ref/chr10.fa"
params.index_dir = "/mnt/d/Studies/Internship/nextflow_pipeline/ref"
params.fastq = "/mnt/d/Studies/Internship/nextflow_pipeline/fastq/*{R1,R2}*"
params.qc_report = "/mnt/d/Studies/Internship/nextflow_pipeline/fastqc_report"
params.bam_dir = "/mnt/d/Studies/Internship/nextflow_pipeline/BAM"
params.sorted_bam_dir = "/mnt/d/Studies/Internship/nextflow_pipeline/BAM/sorted"
params.vcf_dir = "/mnt/d/Studies/Internship/nextflow_pipeline/VCF"
params.filtered_vcf_dir = "/mnt/d/Studies/Internship/nextflow_pipeline/VCF/filtered"
params.annotated_vcf_dir = "/mnt/d/Studies/Internship/nextflow_pipeline/VCF/annotated"
params.snpEff_db = 'GRCh37.75'

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

process QC {
    publishDir("${params.qc_report}", mode: 'copy')
    input:
    path fastq
    output:
    path "*"
    script:
    """
    fastqc $fastq
    """
}

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
    bwa mem -t 16 ${ref} ${fastq} | samtools view -h -b -o > ${sample_id}.bam - 
    """
}

process sorting {
    publishDir("${params.sorted_bam_dir}", mode: 'copy')
    input:
    path bam
    output:
    path "*"
    script:
    """
    samtools sort -m 4G -@ 2 -o ${bam.baseName}.sorted.bam ${bam}
    samtools index ${bam.baseName}.sorted.bam
    """
}

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
    ref_ch = Channel.fromPath(params.ref)
    index(ref_ch)
    index.out.view()

    fastq_ch = Channel.fromPath(params.fastq)
    QC(fastq_ch)
    QC.out.view()

    index_ch=Channel.fromPath(params.index_dir)
    ref_ch=Channel.of(params.ref)
    fastq_ch=Channel.fromFilePairs(params.fastq)
    mapping(index_ch,ref_ch,fastq_ch)
    mapping.out.view()

    bam_ch = Channel.fromPath(params.bam_dir + "/*.bam")
    sorting(bam_ch)
    sorting.out.view()

    sorted_bam_ch = Channel.fromPath(params.sorted_bam_dir + "/*.sorted.bam")
    variant_calling(sorted_bam_ch)
    variant_calling.out.view()

    vcf_ch = Channel.fromPath(params.vcf_dir + "/*.vcf")
    variant_filtering(vcf_ch)
    variant_filtering.out.view()

    filtered_vcf_ch = Channel.fromPath(params.filtered_vcf_dir + "/*.filtered.vcf")
    variant_annotation(filtered_vcf_ch)
    variant_annotation.out.view()
}