#!/usr/bin/env Rscript
# annotate vcf file with titan seg file

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("rtracklayer"));
suppressPackageStartupMessages(library("VariantAnnotation"));
suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library("GenomicRanges"));

if (!interactive()) {
  options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList <- list(
                make_option("--facetsFile", default = NULL, type = "character", action = "store", help ="Facets CNCF file"),
                make_option("--outFile", default = NULL, type = "character", action = "store", help ="targeted interval bed"))

parser <- OptionParser(usage = "%prog [options] [vcf file]", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
  cat("Need vcf file\n\n")
  print_help(parser);
  stop();
} else if (is.null(opt$facetsFile)) {
  cat("Need facets CNCF file\n\n")
  print_help(parser);
  stop();
} else if (is.null(opt$outFile)) {
  cat("Need output file\n\n")
  print_help(parser);
  stop();
}

# read file file
facetsSeg <- read.table(opt$facetsFile, sep = '\t', header = T)
facetsGr <- with(facetsSeg, GRanges(seqnames = chrom,
                                    ranges = IRanges(start = loc.start, end = loc.end),
                                    CF = cf, TCN = tcn, LCN = lcn,
                                    CF_EM = cf.em, TCN_EM = tcn.em, LCN_EM = lcn.em))

# output file
outfn <- opt$outFile
null <- suppressWarnings(file.remove(outfn))
out <- file(outfn, open = 'a')

#input file
vcfFile <- arguments$args[1]
temp <- tempfile()
zipped <- bgzip(vcfFile, temp)
idx <- indexTabix(temp, "vcf")
tab <- TabixFile(zipped, idx, yieldSize = 8000)

open(tab)
while(nrow(vcf <- readVcf(tab, 'hg19'))) {
  info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "Float", Description = "Facets cellular fraction", row.names = "facetsCF"))
  info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "Integer", Description = "Facets total copy number", row.names = "facetsTCN"))
  info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "Integer", Description = "Facets lesser copy number", row.names = "facetsLCN"))
  info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "Float", Description = "Facets cellular fraction (EM)", row.names = "facetsCF_EM"))
  info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "Integer", Description = "Facets total copy number (EM)", row.names = "facetsTCN_EM"))
  info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "Integer", Description = "Facets lesser copy number (EM)", row.names = "facetsLCN_EM"))

  seqlevels(vcf) <- sub('chr', '', seqlevels(vcf))
  ol <- findOverlaps(rowRanges(vcf), facetsGr, select = 'first')
  if (sum(!is.na(ol)) > 0) {
    info(vcf)$facetsCF[!is.na(ol)] <- facetsGr$CF[ol[!is.na(ol)]]
    info(vcf)$facetsTCN[!is.na(ol)] <- facetsGr$TCN[ol[!is.na(ol)]]
    info(vcf)$facetsLCN[!is.na(ol)] <- facetsGr$LCN[ol[!is.na(ol)]]
    info(vcf)$facetsCF_EM[!is.na(ol)] <- facetsGr$CF_EM[ol[!is.na(ol)]]
    info(vcf)$facetsTCN_EM[!is.na(ol)] <- facetsGr$TCN_EM[ol[!is.na(ol)]]
    info(vcf)$facetsLCN_EM[!is.na(ol)] <- facetsGr$LCN_EM[ol[!is.na(ol)]]
  }
  writeVcf(vcf, out)
}
close(tab)
close(out)
