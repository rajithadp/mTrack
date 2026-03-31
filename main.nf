#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input_dir = "results" // Folder containing the pipeline outputs
params.sample    = "sample_alias"
params.region    = "chr11:67641087-71654474"
params.outdir    = "plots"

process PLOT_VARIATION {
    tag "${sample_id}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(sv), path(cnv), path(meth), path(cov)
    val region

    output:
    path "${sample_id}_visualization.pdf"

    script:
    """
    plot_gviz_tracks.R \\
        --sample ${sample_id} \\
        --region "${region}" \\
        --sv_vcf ${sv} \\
        --cnv_vcf ${cnv} \\
        --meth_bw ${meth} \\
        --cov_bg ${cov} \\
        --out ${sample_id}_visualization.pdf
    """
}

workflow {
    // Collect the specific files based on the pipeline's naming convention
    def sv_vcf  = file("${params.input_dir}/${params.sample}.wf_sv.vcf.gz")
    def cnv_vcf = file("${params.input_dir}/${params.sample}.wf_cnv.vcf.gz")
    def meth_bw = file("${params.input_dir}/${params.sample}.wf_mods.5mC.bw")
    def cov_bg  = file("${params.input_dir}/${params.sample}.per-base.bedgraph.gz")

    data_ch = Channel.of([params.sample, sv_vcf, cnv_vcf, meth_bw, cov_bg])

    PLOT_VARIATION(data_ch, params.region)
}