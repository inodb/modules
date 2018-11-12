#!/usr/bin/env Rscript
#### turn segmented copy number data to gene-based copy number with findOverlaps
## define HomDel as TCN=0, loss as TCN<ploidy, gain as TCN>ploidy, amp as TCN>=ploidy+4
## where ploidy= mode of TCN
### some variant of the below, also need one for the breast panel, IMPACT310 and exome

#---------------
# initialization
#---------------

# load base libraries
for (lib in c("optparse","RColorBrewer","GenomicRanges","plyr","dplyr", 
              "stringr","tidyr","magrittr","foreach", "rtracklayer","grid","rlist", "RMySQL")) {
    suppressPackageStartupMessages(library(lib, character.only=TRUE))
}

#--------------
# parse options
#--------------

optList <- list(
				make_option("--outFile", default = NULL, help = "output file"),
				make_option("--mysqlHost", default = '10.0.200.48', help = "MySQL server hostname"),
				make_option("--mysqlPort", default = 38493, help = "MySQL server port"),
				make_option("--mysqlUser", default = 'embl', help = "MySQL server username"),
				make_option("--mysqlPassword", default = NULL, help = "MySQL server password"),
				make_option("--mysqlDb", default = 'homo_sapiens_core_75_37', help = "MySQL server database"),
				make_option("--genesFile", default = NULL, help = "list of genes to include (hgnc symbols)"),
				make_option("--annotFile", default = "~/share/reference/annotation_gene_lists/geneCN.txt", help = "file with annotations to replace MySQL sever query"))
parser <- OptionParser(usage = "%prog [options] [facets files]", option_list = optList)

arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

if (length(arguments$args) < 1) {
	cat("Need facets output files\n")
	print_help(parser)
	stop()
} else if (is.null(opt$outFile)) {
	cat("Need output prefix\n")
	print_help(parser)
	stop()
} else {
	facetsFiles <- arguments$args
}

if (is.null(opt$annotFile)) {
	connect <- function() dbConnect(MySQL(), host = opt$mysqlHost, port = opt$mysqlPort, user = opt$mysqlUser, password = opt$mysqlPassword, dbname = opt$mysqlDb)
	cat('Connecting to ensembl ... ')
	mydb <- connect()
	on.exit(dbDisconnect(mydb))

	query <- "select r.name as chrom,
	g.seq_region_start as start,
	g.seq_region_end as end,
	x.display_label as hgnc,
	k.band as band
	from gene as g
	join seq_region as r on g.seq_region_id = r.seq_region_id
	join xref as x on g.display_xref_id = x.xref_id
	left join karyotype k on g.seq_region_id = k.seq_region_id
	and ((g.seq_region_start >= k.seq_region_start and g.seq_region_start <= k.seq_region_end)
	or (g.seq_region_end >= k.seq_region_start and g.seq_region_end <= k.seq_region_end))
	where x.external_db_id = 1100;"
	repeat {
		rs <- try(dbSendQuery(mydb, query), silent = T)
		if (is(rs, "try-error")) {
			cat("Lost connection to mysql db ... ")
			mydb <- connect()
			cat("reconnected\n")
		} else {
			break
		}
	}
	genes <- dbFetch(rs, -1)
} else {
	genes = read.csv(file=opt$annotFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
}
cat(paste("Found", nrow(genes), "records\n"))

genes %<>% filter(chrom %in% as.character(c(1:22, "X", "Y"))) %>%
		   filter(!duplicated(hgnc)) %>%
		   arrange(as.integer(chrom), start, end)

if (!is.null(opt$genesFile)) {
	g <- scan(opt$genesFile, what = 'character')
	genes %<>% filter(hgnc %in% g)
	absentGenes <- g[!g %in% genes$hgnc]
	if (length(absentGenes) > 0) {
		print("Unable to find", length(absentGenes), "in database\n");
		cat(absentGenes, sep = '\n');
	}
}

cat(paste("Filtering to", nrow(genes), "records\n"))

genesGR <- genes %$% GRanges(seqnames = chrom, ranges = IRanges(start, end), band = band, hgnc = hgnc)
			
mm <- lapply(facetsFiles, function(f) {
    #load(f)
    #tab <- fit$cncf
	#tab$chrom[which(tab$chrom==23)] <- "X"

	#tabGR <- tab %$% GRanges(seqnames = chrom, ranges = IRanges(start, end))
	#mcols(tabGR) <- tab %>% select(cnlr.median:lcn.em)

	#fo <- findOverlaps(tabGR, genesGR)

	#df <- as.data.frame(cbind(mcols(genesGR)[subjectHits(fo),], mcols(tabGR)[queryHits(fo),]))
	#df %<>% group_by(hgnc) %>% top_n(1, abs(cnlr.median))

	#ploidy <- table(df$tcn.em)
	#ploidy <- as.numeric(names(ploidy)[which.max(ploidy)])

	#df$GL <- 0
	#df$GL[df$tcn.em < ploidy] <- -1
	#df$GL[df$tcn.em == 0] <- -2
	#df$GL[df$tcn.em > ploidy] <- 1
	#df$GL[df$tcn.em >= ploidy + 4] <- 2

	#df %>% select(hgnc, GL) %>% ungroup

    print(paste("reading table", f, sep=" "));
    cn.table <- read.table(f,sep="\t",header=T)
	cn.table$chrom[which(cn.table$chrom==23)] <- "X"
	cn.table$chrom[which(cn.table$chrom==24)] <- "Y"

	tabGR <- cn.table %$% GRanges(seqnames = chrom, ranges = IRanges(start, end))
	mcols(tabGR) <- cn.table %>% select(cnTotal:date)

	fo <- findOverlaps(tabGR, genesGR)

	df <- as.data.frame(cbind(mcols(genesGR)[subjectHits(fo),], mcols(tabGR)[queryHits(fo),]))
	df %<>% group_by(hgnc) %>% top_n(1, abs(nBraw))

    # use cnTotal for copy number
    # and add ploidy
    df$GL <- df$cnTotal
	df$GL2 <- df$Ploidy

    # copy number change
    ploidy <- unique(df$Ploidy) # should be single value
	df$GL3 <- 0
	df$GL3[df$cnTotal < ploidy] <- -1
	df$GL3[df$cnTotal == 0] <- -2
	df$GL3[df$cnTotal > ploidy] <- 1
	df$GL3[df$cnTotal >= ploidy + 4] <- 2

	df %>% select(hgnc, GL, GL2, GL3) %>% ungroup
})
names(mm) <- facetsFiles
for (f in facetsFiles) {
	n <- sub('\\..*', '', sub('.*/', '', f))
	colnames(mm[[f]])[2:4] <- paste(n, c("cnTotal", "Ploidy", "cnChange"), sep="_")
}

mm <- left_join(genes, join_all(mm, type = 'full', by="hgnc")) %>% arrange(as.integer(chrom), start, end)
write.table(mm, file=opt$outFile, sep="\t", row.names=F, na="", quote=F)

