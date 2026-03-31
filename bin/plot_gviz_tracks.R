#!/usr/bin/env Rscript

library(optparse)
library(Gviz)
library(GenomicRanges)
library(rtracklayer)

option_list = list(
  make_option(c("--sample"), type="character", default=NULL, help="Sample ID"),
  make_option(c("--region"), type="character", default=NULL, help="Genomic region (e.g., chr1:1000-2000)"),
  make_option(c("--sv"), type="character", default=NULL, help="Path to SV VCF"),
  make_option(c("--cnv"), type="character", default=NULL, help="Path to CNV VCF/BED"),
  make_option(c("--meth"), type="character", default=NULL, help="Path to Methylation BedGraph"),
  make_option(c("--cov"), type="character", default=NULL, help="Path to Coverage BigWig"),
  make_option(c("--out"), type="character", default="output.pdf", help="Output PDF file")
)
opt = parse_args(OptionParser(option_list=option_list))

# Parse Region
chrom <- strsplit(opt$region, ":")[]
coords <- strsplit(strsplit(opt$region, ":")[], "-")[]
start_pos <- as.numeric(coords)
end_pos <- as.numeric(coords)

# 1. Base Tracks
itrack <- IdeogramTrack(genome = "hg38", chromosome = chrom)
gtrack <- GenomeAxisTrack()

# 2. Coverage Track (BigWig)
cov_track <- DataTrack(range = opt$cov, genome = "hg38", type = "l", 
                       chromosome = chrom, name = "Coverage", fill.mountain = c("blue", "blue"))

# 3. Methylation Track (BedGraph)
meth_track <- DataTrack(range = opt$meth, genome = "hg38", type = "p", 
                        chromosome = chrom, name = "Methylation %", col = "red")

# 4. SV / CNV Tracks (Assuming BED format for simplicity in this template)
# Note: If these are raw VCFs, you will need to parse them with VariantAnnotation first.
sv_track <- AnnotationTrack(range = opt$sv, genome = "hg38", name = "SVs", chromosome = chrom)
cnv_track <- AnnotationTrack(range = opt$cnv, genome = "hg38", name = "CNVs", chromosome = chrom)

# Plot to PDF
pdf(opt$out, width = 10, height = 8)
plotTracks(list(itrack, gtrack, cov_track, meth_track, cnv_track, sv_track), 
           from = start_pos, to = end_pos, main = paste("Multi-Omics:", opt$sample))
dev.off()