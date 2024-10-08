library(matrixStats)    # 1.3.0
library(dplyr)    # 1.1.4


# functions to calculate the fold of IQR
cal.fold <- function(x, Q1, Q3, IQR) {
    if (x > Q3) {
        return((x - Q3) / IQR)
    } else if (x >= Q1) {
        return(0)
    } else {
        return((x - Q1) / IQR)
    }
}

cal.folds.rowwise <- function(a.row) {
    mapply(cal.fold, a.row[5:length(a.row)], a.row[1], a.row[3], a.row[4])
}


# function to get the P-Value of Fisher's exact test of independence
get.fisher.p <- function(a.row, coldata, a.cofactor) {    
    coldata.tb <- as_tibble(coldata)
    coldata.tb$Outlier <- a.row

    count.tb <- coldata.tb %>% group_by(.data[[a.cofactor]], Outlier) %>% summarise(Count = n(), .groups = 'drop')
    count.m <- reshape(as.data.frame(count.tb), v.names = 'Count',
                       timevar = 'Outlier', idvar = a.cofactor, direction = 'wide')
    count.m[is.na(count.m)] <- 0
    
    fisher.res <- fisher.test(count.m[, -1])

    return(fisher.res$p.value)
}


# function to get the counts of OO-OO, OO-NO, NO-OO, and NO-NO
get.counts <- function(a.row, outlier.m) {
    a.row.df <- as.data.frame(t(outlier.m[a.row, ]))
    colnames(a.row.df) <- c('Gene1', 'Gene2')
    a.row.df <- as.data.frame(a.row.df %>% group_by(Gene1, Gene2) %>% summarise(Count = n(), .groups = 'drop'))
    
    a.row.m <- matrix(rep(0, 4), nr = 2, dimnames = list(c('OO', 'NO'), c('OO', 'NO')))
    for (i in 1:nrow(a.row.df)) {
        a.row.m[a.row.df[i, 'Gene2'], a.row.df[i, 'Gene1']] <- a.row.df[i, 'Count']
    }
    
    return(as.vector(a.row.m))
}


# function to get Kendall's tau and the P-Value of tau after removing outlier pairs
get.tau.P <- function(a.row, outlier.m, filtered.m) {
    a.row.outlier <- outlier.m[a.row, ]
    a.row.NO <- a.row.outlier == 'NO'
    a.row.NO.samp.IDs <- colnames(a.row.NO[, colSums(a.row.NO) == 2])
    
    a.row.NO.filtered <- filtered.m[a.row, a.row.NO.samp.IDs]
    a.row.NO.tau.P <- cor.test(a.row.NO.filtered[a.row[1], ], a.row.NO.filtered[a.row[2], ], method = 'kendall')
    
    return(c(a.row.NO.tau.P$estimate, a.row.NO.tau.P$p.value))
}


# function to analyze outlier
analyze.outlier <- function(TPM.df, coldata.df, high.similarity.genes, chromosomes, gene.metadata.cols, cofactors,
                            TPM.cutoff, fold.cutoff, p.cutoff, OO.pair.cutoff, out.prefix) {
    TPM.m <- as.matrix(TPM.df[, rownames(coldata.df)])

## filtration: TPM > TPM.cutoff in at least one sample, on autosome or X chromosome,
##     protein-coding, and without high sequence similarity
    filtered.index <- rowMaxs(TPM.m, useNames = T) > TPM.cutoff
    filtered <- TPM.m[filtered.index, ]
    filtered.gene.metadata <- TPM.df[filtered.index, gene.metadata.cols]

    filtered.index <- (filtered.gene.metadata$Chr %in% chromosomes) &
        (filtered.gene.metadata$Genebiotype == 'protein_coding') &
        (!(rownames(filtered.gene.metadata) %in% high.similarity.genes))
    filtered <- filtered[filtered.index, ]
    filtered.gene.metadata <- filtered.gene.metadata[filtered.index, ]

## get Q1, median, Q3, and IQR
    filtered.stats <- as.data.frame(rowQuantiles(filtered, useNames = T)[, c(2, 3, 4)])
    colnames(filtered.stats) <- c('Q1', 'Median', 'Q3')
    filtered.stats$IQR <- filtered.stats$Q3 - filtered.stats$Q1

## get the fold of IQR
    filtered.fold <- t(apply(cbind(filtered.stats, filtered), 1, cal.folds.rowwise))

## get over outlier (OO), non-outlier(NO), and under outlier (UO)
    outlier <- matrix('NO', nr = nrow(filtered), nc = ncol(filtered),
                      dimnames = list(rownames(filtered), colnames(filtered)))
    outlier[filtered.fold > fold.cutoff & filtered > TPM.cutoff] <- 'OO'
    outlier[filtered.fold < -fold.cutoff] <- 'UO'

    outlier.index <- rowSums(outlier != 'NO') != 0
    outlier <- outlier[outlier.index, ]
    outlier.gene.metadata <- filtered.gene.metadata[outlier.index, ]

    if (nrow(outlier) > 0) {
        cofactor.num <- length(cofactors)
        if (cofactor.num > 0) {
## get the P-Value of Fisher's exact test of independence
            outlier.ps <- matrix(nr = nrow(outlier), nc = cofactor.num,
                                 dimnames = list(rownames(outlier), cofactors))
            for (cofactor in cofactors) {
                outlier.ps[, cofactor] <- apply(outlier, 1, get.fisher.p, coldata.df, cofactor)
            }
            filtered.stats[rownames(outlier.ps), cofactors] <- outlier.ps
            
            outlier.indep.index <- rowAlls(outlier.ps > p.cutoff)
            outlier.indep <- outlier[outlier.indep.index, ]
            outlier.indep.gene.metadata <- outlier.gene.metadata[outlier.indep.index, ]
        } else {
            outlier.indep <- outlier
            outlier.indep.gene.metadata <- outlier.gene.metadata
        }
        
        if (nrow(outlier.indep) > 0) {
## outlier output
            filtered.round <- round(filtered[rownames(outlier.indep), ])
            outlier.indep.out <- data.frame(row.names = rownames(outlier.indep))
            for (samp.ID in rownames(coldata.df)) {
                outlier.indep.out[, samp.ID] <- gsub('NO-', '', paste(outlier.indep[, samp.ID],
                                                                      filtered.round[, samp.ID], sep = '-'))
            }
            GeneID.df <- data.frame(GeneID = rownames(outlier.indep))
            write.table(cbind(GeneID.df, outlier.indep.gene.metadata, outlier.indep.out),
                        paste0(out.prefix, 'outlier.tsv'), quote = F, sep = '\t', row.names = F)

## correlated OOs
### choose outlier: without UO, and with at least three OOs
            outlier.chosen <- outlier.indep[rowSums(outlier.indep == 'UO') == 0, ]
            outlier.chosen <- outlier.chosen[rowSums(outlier.chosen == 'OO') >= OO.pair.cutoff, ]
            
            if (nrow(outlier.chosen) > 0) {
### choose outlier pair: at least three OO-OO pairs, and NO-NO for all others
                genes.chosen <- rownames(outlier.chosen)
                outlier.pair.chosen <- data.frame(GeneID1 = rep(genes.chosen, each = length(genes.chosen)),
                                                  GeneID2 = rep(genes.chosen, length(genes.chosen))) %>%
                    filter(GeneID1 < GeneID2)
                outlier.pair.chosen[, 3:6] <- t(apply(outlier.pair.chosen, 1, get.counts, outlier.chosen))
                colnames(outlier.pair.chosen)[3:6] <- c('OO_OO', 'OO_NO', 'NO_OO', 'NO_NO')
                outlier.pair.chosen <- outlier.pair.chosen %>%
                    filter(OO_OO >= OO.pair.cutoff, OO_NO == 0, NO_OO == 0)

                if (nrow(outlier.pair.chosen) > 0) {
### get Kendall's tau and the P-Value of tau after removing outlier pairs
                    filtered.chosen <- filtered[genes.chosen, ]
                    outlier.pair.chosen[, 7:8] <- t(apply(outlier.pair.chosen[, 1:2], 1, get.tau.P,
                                                          outlier.chosen, filtered.chosen))
                    colnames(outlier.pair.chosen)[7:8] <- c('Tau', 'P')

                    outlier.pair.chosen$P_adj <- p.adjust(outlier.pair.chosen$P, method = 'bonferroni')

### correlated OOs output
                    outlier.indep.gene.metadata.chosen <- cbind(GeneID.df,
                        outlier.indep.gene.metadata)[genes.chosen, c('GeneID', 'Genename', 'Chr', 'Start')]
                    outlier.pair.chosen <- inner_join(outlier.pair.chosen, outlier.indep.gene.metadata.chosen,
                                                      by = c('GeneID1' = 'GeneID'))
                    outlier.pair.chosen <- inner_join(outlier.pair.chosen, outlier.indep.gene.metadata.chosen,
                                                      by = c('GeneID2' = 'GeneID'))
                    colnames(outlier.pair.chosen)[10:15] <- c('Genename1', 'Chr1', 'Start1',
                                                              'Genename2', 'Chr2', 'Start2')
                    outlier.pair.chosen <- outlier.pair.chosen %>% select(GeneID1, Genename1, Chr1, Start1,
                        GeneID2, Genename2, Chr2, Start2, OO_OO, OO_NO, NO_OO, NO_NO, Tau, P, P_adj)
                    write.table(outlier.pair.chosen, paste0(out.prefix, 'outlier_pair.tsv'),
                                quote = F, sep = '\t', row.names = F)
                }
            }
        }
    }
}


# input files
## TPM file with gene metadata
TPM.df <- read.table('data/GTEx_brain_TPM.tsv', header = T, sep = '\t', row.names = 'GeneID', check.names = F)

## sample annotation file
coldata.df <- read.table('data/GTEx_brain_coldata.tsv', header = T, sep = '\t', row.names = 'SAMPID')

## genes with high similarity paralogs
high.similarity.genes <- read.table('data/human_genes_wParalogs.list')[, 1]


# input parameters
chromosomes <- paste0(rep('chr', 23), c(as.character(1:22), 'X'))    ## autosomes and X chromosome
gene.metadata.cols <- 1:6    ## gene metadata columns
cofactors <- c('Sex', 'AGE', 'DTHHRDY')    ## cofactors

TPM.cutoff <- 5    ## minimum expression level for the highest sample or over outliers (OOs)
fold.cutoff <- 5    ## cutoff for the fold of IQR
p.cutoff <- 0.1    ## marginally significant
OO.pair.cutoff <- 3    ## cutoff for the number of OO-OO pairs

out.prefix <- 'data/GTEx_brain_'    ## prefix of output file names


analyze.outlier(TPM.df, coldata.df, high.similarity.genes, chromosomes, gene.metadata.cols, cofactors,
                TPM.cutoff, fold.cutoff, p.cutoff, OO.pair.cutoff, out.prefix)

