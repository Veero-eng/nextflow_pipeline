params.read_groups_dir = "/mnt/d/Studies/Internship/nextflow_pipeline/Somatic_variant_calling"
params.indexed_bam_dir = "/mnt/d/Studies/Internship/nextflow_pipeline/Somatic_variant_calling"
params.mutect2_svc_dir = "/mnt/d/Studies/Internship/nextflow_pipeline/Somatic_variant_calling"
params.output_dir = "/mnt/d/Studies/Internship/nextflow_pipeline/Somatic_variant_calling"
params.ref = "/mnt/d/Studies/Internship/nextflow_pipeline/ref/chr10.fa"
params.normal_bam = "/mnt/d/Studies/Internship/nextflow_pipeline/somatic_vc/mutect2/normal_aln_sort.bam"
params.tumor_bam = "/mnt/d/Studies/Internship/nextflow_pipeline/somatic_vc/mutect2/tumor_aln_sort.bam"

process add_read_groups_normal {
  publishDir("${params.read_groups_dir}", mode: 'copy')
  input:
    path normal_bam
  output:
    path "*"

  script:
  """
  java -jar /mnt/d/Studies/Internship/nextflow_pipeline/tools/gatk-4.6.0.0/gatk-package-4.6.0.0-local.jar AddOrReplaceReadGroups -I ${params.normal_bam} -O normal_aln_sort_rg.bam -RGID 1 -RGLB Lipoma-N15 -RGPL Illumina -RGPU unit1 -RGSM normal
  """
}

process add_read_groups_tumor {
  publishDir("${params.read_groups_dir}", mode: 'copy')
  input:
    path tumor_bam
  output:
    path "*"

  script:
  """
  java -jar /mnt/d/Studies/Internship/nextflow_pipeline/tools/gatk-4.6.0.0/gatk-package-4.6.0.0-local.jar AddOrReplaceReadGroups -I ${params.tumor_bam} -O tumor_aln_sort_rg.bam -RGID 2 -RGLB Lipoma-T19 -RGPL Illumina -RGPU unit2 -RGSM tumor
  """
}

process index_normal_bam {
  publishDir("${params.indexed_bam_dir}", mode: 'copy')
  input:
    path normal_bam
  output:
    path "*"

  script:
  """
  samtools index ${normal_bam}
  """
}

process index_tumor_bam {
  publishDir("${params.indexed_bam_dir}", mode: 'copy')
  input:
    path tumor_bam
  output:
    path "*"

  script:
  """
  samtools index ${tumor_bam}
  """
}

process mutect2 {
  publishDir("${params.mutect2_svc_dir}", mode: 'copy')
  input:
    path tumor_bam
    path normal_bam
    path ref
  output:
    path "*"

  script:
  """
  java -jar /mnt/d/Studies/Internship/nextflow_pipeline/tools/gatk-4.6.0.0/gatk-package-4.6.0.0-local.jar Mutect2 -R ${params.ref} -I /mnt/d/Studies/Internship/nextflow_pipeline/Somatic_variant_calling/tumor_aln_sort_rg.bam -I /mnt/d/Studies/Internship/nextflow_pipeline/Somatic_variant_calling/normal_aln_sort_rg.bam -O somatic_variants_mutect2.vcf
  """
}

workflow {
  normal_bam_ch = Channel.fromPath(params.normal_bam)
  tumor_bam_ch = Channel.fromPath(params.tumor_bam)
  
  add_read_groups_normal_ch = Channel.fromPath(params.read_groups_dir + "/*_rg.bam")
  add_read_groups_normal(normal_bam_ch)

  add_read_groups_tumor_ch = Channel.fromPath(params.read_groups_dir + "/*_rg.bam")
  add_read_groups_tumor(tumor_bam_ch)

  indexed_bam_dir_ch = Channel.fromPath(params.indexed_bam_dir + "/*_rg.bam.bai")
  index_normal_bam(add_read_groups_normal.out)
  index_tumor_bam(add_read_groups_tumor.out)

  mutect2_svc_dir_ch = Channel.fromPath(params.mutect2_svc_dir + "/*.vcf")
  mutect2(add_read_groups_tumor.out, add_read_groups_normal.out, Channel.fromPath(params.ref))
  mutect2.out.view()
}