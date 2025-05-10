params.somatic_sniper_dir = "/mnt/d/Studies/Internship/nextflow_pipeline/tools/VarScan/VarScan_files/somatic_sniper_output"
params.varscan_dir = "/mnt/d/Studies/Internship/nextflow_pipeline/tools/VarScan/VarScan_files/varscan_output"
params.normal_bam = "/mnt/d/Studies/Internship/nextflow_pipeline/tools/VarScan/VarScan_files/normal_aln_sort_rg.bam"
params.tumor_bam = "/mnt/d/Studies/Internship/nextflow_pipeline/tools/VarScan/VarScan_files/tumor_aln_sort_rg.bam"
params.ref = "/mnt/d/Studies/Internship/nextflow_pipeline/tools/VarScan/VarScan_files/chr10.fa"
params.index_dir = "/mnt/d/Studies/Internship/nextflow_pipeline/tools/VarScan/VarScan_files"

process somatic_sniper {
  publishDir("${params.somatic_sniper_dir}", mode: 'copy')
  input:
    path normal_bam
    path tumor_bam
    path ref
  output:
    path "*"

    script:
    """
    bam-somaticsniper -f ${ref} ${tumor_bam} ${normal_bam} somatic_sniper_output.txt
    """
}

process varscan {
  publishDir("${params.varscan_dir}", mode: 'copy')
  input:
    path normal_bam
    path tumor_bam
    path ref
  output:
    path "*"

    script:
    """
    samtools mpileup -f ${params.ref} ${normal_bam} > mynormalData.mpileup
    samtools mpileup -f ${params.ref} ${tumor_bam} > mytumorData.mpileup
    java -jar ${params.index_dir}/VarScan.v2.3.9.jar somatic mynormalData.mpileup mytumorData.mpileup varscan_output
     """
}

workflow {
  normal_bam_ch = Channel.fromPath(params.normal_bam)
  tumor_bam_ch = Channel.fromPath(params.tumor_bam)
  ref_ch = Channel.fromPath(params.ref)
  somatic_sniper_ch = somatic_sniper(normal_bam_ch, tumor_bam_ch, ref_ch)
  varscan_ch = varscan(normal_bam_ch, tumor_bam_ch, ref_ch)
}
