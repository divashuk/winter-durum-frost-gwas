# Этап 5.3 (R-часть): фильтрация выравниваний кДНК Svevo на интервалы QTL Langdon и сбор аннотации.
suppressPackageStartupMessages(library(data.table))
OUT <- "/work/project/gwas_frost"
G   <- "/work/scratch/genes"
p <- fread(file.path(G, "hits.paf"), header = FALSE, fill = TRUE, select = 1:12,
           col.names = c("tx", "qlen", "qs", "qe", "strand", "target", "tlen", "ts", "te", "nmatch", "alnlen", "mapq"))
p[, `:=`(identity = nmatch / alnlen, qcov = (qe - qs) / qlen)]
p <- p[identity >= 0.98 & qcov >= 0.90]
p[, gene := sub("\\.[0-9]+$", "", tx)]
p <- p[order(-nmatch)][, .SD[1], by = .(gene, target)]
p[, `:=`(QTL = sub("\\|.*", "", target), region = sub(".*\\|", "", target))]
p[, `:=`(chr = sub(":.*", "", region), int_start = as.numeric(sub(".*:([0-9]+)-.*", "\\1", region)))]
p[, `:=`(langdon_start = int_start + ts - 1, langdon_end = int_start + te - 1)]
h <- data.table(line = readLines(file.path(G, "cdna_headers.txt")))
h[, tx := sub(" .*", "", line)]
h[, description := ifelse(grepl("description:", line), sub(".*description:", "", line), "")]
h[, svevo_pos := sub(".*chromosome:Svevo.v1:([^ ]+).*", "\\1", line)]
h[, biotype := sub(".*gene_biotype:([^ ]+).*", "\\1", line)]
gff <- fread(file.path(G, "genes.gff3"), header = FALSE, sep = "\t", quote = "")
gff[, gene := sub(".*ID=gene:([^;]+).*", "\\1", V9)]
gff[, gff_description := ifelse(grepl("description=", V9), utils::URLdecode(sub(".*description=([^;]+).*", "\\1", V9)), "")]
p <- merge(p, h[, .(tx, svevo_pos, biotype, description)], by = "tx", all.x = TRUE)
p <- merge(p, unique(gff[, .(gene, gff_description)]), by = "gene", all.x = TRUE)
p[, annotation := fifelse(nzchar(gff_description), gff_description, fifelse(nzchar(description), description, "нет описания в аннотации Svevo"))]
res <- p[order(QTL, langdon_start), .(QTL, chr, langdon_start, langdon_end, svevo_gene = gene, svevo_transcript = tx, svevo_pos, biotype,
                                       identity = round(identity, 4), coverage = round(qcov, 3), annotation)]
fwrite(res, file.path(OUT, "tables/05_genes_in_qtl.csv"))
print(res[, .(n_genes = .N, n_with_description = sum(annotation != "нет описания в аннотации Svevo")), by = QTL])
