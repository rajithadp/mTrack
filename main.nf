#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define input parameters (default to null, user should provide them via CLI)
params.sample_id   = "sample_1"
params.region      = "chr1:1000000-2000000" // Default region to plot
params.sv_vcf      = null
params.cnv_vcf     = null
params.meth_bed    = null
params.cov_bw      = null
params.outdir      = "results"

process VISUALIZE_TRACKS {
    publishDir "${params.outdir}/plots", mode: 'copy'
    
    input:
    val sample_id
    val region
    path sv_vcf
    path cnv_vcf
    path meth_bed
    path cov_bw

    output:
    path "${sample_id}_tracks.pdf"

    script:
    """
    # Run the R script located in the bin/ directory
    plot_gviz_tracks.R \\
        --sample ${sample_id} \\
        --region ${region} \\
        --sv ${sv_vcf} \\
        --cnv ${cnv_vcf} \\
        --meth ${meth_bed} \\
        --cov ${cov_bw} \\
        --out ${sample_id}_tracks.pdf
    """
}

workflow {
    // Check if required inputs are provided
    if (!params.sv_vcf || !params.cnv_vcf || !params.meth_bed || !params.cov_bw) {
        error "Please provide all required input tracks (--sv_vcf, --cnv_vcf, --meth_bed, --cov_bw)"
    }

    sv_ch   = file(params.sv_vcf)
    cnv_ch  = file(params.cnv_vcf)
    meth_ch = file(params.meth_bed)
    cov_ch  = file(params.cov_bw)

    VISUALIZE_TRACKS(params.sample_id, params.region, sv_ch, cnv_ch, meth_ch, cov_ch)
}