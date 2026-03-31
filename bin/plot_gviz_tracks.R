#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Gviz)
  library(GenomicRanges)
  library(VariantAnnotation)
})

option_list = list(
  make_option(c("--sample"), type="character", help="Sample Alias"),
  make_option(c("--region"), type="character", help="Region (chr:start-end)"),
  make_option(c("--sv_vcf"), type="character", help="wf_sv.vcf.gz"),
  make_option(c("--cnv_vcf"), type="character", help="wf_cnv.vcf.gz"),
  make_option(c("--meth_bw"), type="character", help="wf_mods.5mC.bw"),
  make_option(c("--cov_bg"), type="character", help="per-base.bedgraph.gz"),
  make_option(c("--out"), type="character", help="Output PDF name")
)
opt = parse_args(OptionParser(option_list=option_list))

# Parse coordinates
reg <- strsplit(opt$region, "[:-]")[]
chrom <- reg; s <- as.numeric(reg); e <- as.numeric(reg)
gen <- "hg38" # Standard for wf-human-variation

# 1. Axis & Ideogram
axTrack <- GenomeAxisTrack()
idxTrack <- IdeogramTrack(genome = gen, chromosome = chrom)

# 2. SV Track (Arrows)
# We parse the VCF specifically for the requested region to keep it light
sv_vcf <- readVcf(opt$sv_vcf, gen)
sv_gr <- rowRanges(sv_vcf)
svTrack <- AnnotationTrack(sv_gr, name = "SVs", chromosome = chrom, 
                           shape = "arrow", fill = "red", col = "red")

# 3. Coverage Track (Step-line from BedGraph)
covTrack <- DataTrack(range = opt$cov_bg, genome = gen, chromosome = chrom,
                      type = "s", name = "Coverage", fill = "gray", col = "black")

# 4. Methylation Track (Histogram from BigWig)
methTrack <- DataTrack(range = opt$meth_bw, genome = gen, chromosome = chrom,
                       type = "h", name = "5mC %", col = "black", fill = "black")

# 5. CNV Track
cnv_vcf <- readVcf(opt$cnv_vcf, gen)
cnv_gr <- rowRanges(cnv_vcf)
cnvTrack <- AnnotationTrack(cnv_gr, name = "CNVs", chromosome = chrom, 
                            fill = "blue", col = "darkblue")

# Plotting
pdf(opt$out, width = 12, height = 10)
plotTracks(list(idxTrack, axTrack, svTrack, covTrack, methTrack, cnvTrack),
           from = s, to = e, chromosome = chrom,
           background.title = "white", col.title = "black", 
           cex.title = 0.8, justification.title = "left")
dev.off()