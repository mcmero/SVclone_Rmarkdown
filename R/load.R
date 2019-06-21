load_hiclus_snvs <- function(basename, basedir) {
    mix <- basename
    ccubeResultRDataFile = paste0(basedir, basename, '_out/',
                                  '/ccube_out/snvs/', basename, '_ccube_snv_results.RData')
    load(ccubeResultRDataFile)
    
    snvs <- data.frame(data='snvs', snvRes$ssm, sv=snvRes$ssm$mutation_id, mix=mix)
    snvs <- attach_001_ground_truth_snvs(wd_truth, snvs)
    
    snvin <- read.delim(paste0(basedir, basename, '_out/', basename, '_filtered_snvs.tsv'))
    snvin$mutation_id <- paste(snvin$chrom, snvin$pos, sep='_')
    snvs <- inner_join(snvs, snvin[,c('chrom', 'mutation_id', 'gtype')], by='mutation_id')
    
    snvs$true_ccf <- 1
    if(mix == '4clus') {
        snvs[snvs$sample%in%'bM_only' & snvs$chrom %% 2 == 1, 'true_ccf'] <- 0.2
        snvs[snvs$sample%in%'bM_only' & snvs$chrom %% 2 == 0, 'true_ccf'] <- 0.6
        snvs[snvs$sample%in%'gM_only' & snvs$chrom %% 2 == 1, 'true_ccf'] <- 0.4
    } else {
        snvs[snvs$sample%in%'bM_only' & snvs$chrom %% 2 == 1, 'true_ccf'] <- 0.8
        snvs[snvs$sample%in%'bM_only' & snvs$chrom %% 2 == 0, 'true_ccf'] <- 0.6
        snvs[snvs$sample%in%'gM_only' & snvs$chrom %% 2 == 1, 'true_ccf'] <- 0.2
        snvs[snvs$sample%in%'gM_only' & snvs$chrom %% 2 == 0, 'true_ccf'] <- 0.4
    }
    
    snvs$adjusted_support <- snvs$var_counts
    snvs$adjusted_depth <- snvs$var_counts + snvs$ref_counts
    snvs$norm_cn <- as.numeric(snvs$normal_cn)
    
    x <- t(data.frame(apply(snvs, 1, get_true_sc_cn)))
    snvs$true_cn <- x[,1]
    snvs$wtotal_cn <- x[,2]
    snvs$pv <- x[,3]
    snvs$best_ccf <- x[,4]
    snvs$best_ccf[snvs$best_ccf > 2] <- 2
    
    return(snvs)
}

load_hiclus_svs <- function(basename, basedir) {
    mix <- basename
    ccubeResultRDataFile = paste0(paste0(basedir, basename, '_out/'), '/ccube_out/', 
                                  basename, '_ccube_sv_results.RData')
    load(ccubeResultRDataFile)
    svs <- data.frame(data='svs', doubleBreakPtsRes$ssm, mix=mix)
    
    svin <- read.delim(paste0(basedir, basename, '_out/', '/', basename, '_filtered_svs.tsv'))
    svin$sv <- paste(svin$chr1, svin$pos1, svin$dir1, svin$chr2, svin$pos2, svin$dir2, sep=':')
    svin$mutation_id <- paste(paste(svin$chr1, svin$pos1, svin$dir1, sep=':'), 
                              paste(svin$chr2, svin$pos2, svin$dir2, sep=':'), sep='_')
    svs <- inner_join(svs, svin, by='mutation_id')
    svs <- attach_001_ground_truth_svs(wd_truth, svs)
    
    pp <- read.delim(paste0(basedir, basename, '_out/purity_ploidy.txt'))
    pur <- pp$purity; pl <- pp$ploidy
    af <- 1 - (pur / pl)
    gains <- svs$classification%in%c('DUP','INTDUP')
    svs$norm1_adjusted <- svs$norm1; svs$norm2_adjusted <- svs$norm2;
    svs$norm1_adjusted[gains] <- svs$norm1_adjusted[gains] * af
    svs$norm2_adjusted[gains] <- svs$norm2_adjusted[gains] * af
    svs$norm_cn <- 2; svs$norm_cn[svs$chr1%in%c('X','Y','MT')] <- 1
    
    svs$true_ccf <- 1
    if(mix == '4clus') {
        svs[svs$sample%in%'bM_only' & svs$chr1 %% 2 == 1, 'true_ccf'] <- 0.2
        svs[svs$sample%in%'bM_only' & svs$chr1 %% 2 == 0, 'true_ccf'] <- 0.6
        svs[svs$sample%in%'gM_only' & svs$chr1 %% 2 == 1, 'true_ccf'] <- 0.4
    } else {
        svs[svs$sample%in%'bM_only' & svs$chr1 %% 2 == 1, 'true_ccf'] <- 0.8
        svs[svs$sample%in%'bM_only' & svs$chr1 %% 2 == 0, 'true_ccf'] <- 0.6
        svs[svs$sample%in%'gM_only' & svs$chr1 %% 2 == 1, 'true_ccf'] <- 0.2
        svs[svs$sample%in%'gM_only' & svs$chr1 %% 2 == 0, 'true_ccf'] <- 0.4
    }
    
    rownames(svs) <- 1:nrow(svs)
    svs$mean_ccf <- apply(svs[,c('ccube_ccf1','ccube_ccf2')], 1, mean)
    
    x <- t(data.frame(apply(svs, 1, get_true_sc_cn_by_side, pur=pur, side=1)))
    svs$true_cn1 <- x[,1]
    svs$wtotal_cn1 <- x[,2]
    svs$pv1 <- x[,3]
    svs$best_ccf1 <- x[,4]
    
    x <- t(data.frame(apply(svs, 1, get_true_sc_cn_by_side, pur=pur, side=2)))
    svs$true_cn2 <- x[,1]
    svs$wtotal_cn2 <- x[,2]
    svs$pv2 <- x[,3]
    svs$best_ccf2 <- x[,4]
    return(svs)
}

get_all_sv_meta_3clus <- function(ssm, p1, resultFolder, wd_truth) {
    basename <- paste('001bM_p', p1, '_001gM_p', 1-p1, sep='')
    sv_filt_file <- paste0(resultFolder, '/', basename, '/', basename, '_filtered_svs.tsv')
    pp_file <- paste0(resultFolder, '/', basename, '/purity_ploidy.txt')
    
    svs <- read.delim(sv_filt_file, sep='\t')
    svs$mutation_id <- paste(paste(svs$chr1, svs$pos1, svs$dir1, sep=':'), 
                             paste(svs$chr2, svs$pos2, svs$dir2, sep=':'), sep='_')
    svs <- inner_join(ssm, svs, by='mutation_id')
    svs$sv <- paste(svs$chr1, svs$pos1, svs$dir1, svs$chr2, svs$pos2, svs$dir2, sep=':')
    svs <- attach_001_ground_truth_svs(wd_truth, svs)
    
    svs$true_ccf <- 1
    svs[svs$sample%in%'bM_only','true_ccf'] <- p1
    svs[svs$sample%in%'gM_only','true_ccf'] <- 1-p1
    
    pp <- read.delim(pp_file); pur <- pp$purity; pl <- pp$ploidy
    af <- 1 - (pur / pl)
    gains <- svs$classification%in%c('DUP','INTDUP')
    svs$norm1_adjusted <- svs$norm1; svs$norm2_adjusted <- svs$norm2;
    svs$norm1_adjusted[gains] <- svs$norm1_adjusted[gains] * af
    svs$norm2_adjusted[gains] <- svs$norm2_adjusted[gains] * af
    svs$norm_cn <- 2; svs$norm_cn[svs$chr1%in%c('X','Y','MT')] <- 1
    rownames(svs) <- 1:nrow(svs)
    
    x <- t(data.frame(apply(svs, 1, get_true_sc_cn_by_side, pur=pur, side=1)))
    svs$true_cn1 <- x[,1]
    svs$wtotal_cn1 <- x[,2]
    svs$pv1 <- x[,3]
    svs$best_ccf1 <- x[,4]
    
    x <- t(data.frame(apply(svs, 1, get_true_sc_cn_by_side, pur=pur, side=2)))
    svs$true_cn2 <- x[,1]
    svs$wtotal_cn2 <- x[,2]
    svs$pv2 <- x[,3]
    svs$best_ccf2 <- x[,4]
    
    return(svs)
}

get_all_snv_meta_3clus <- function(ssm, p1, resultFolder, wd_truth) {
    snvs <- data.frame(data='snvs', ssm, sv=ssm$mutation_id, mix=mix, p1=p1, p2=1-p1)
    snvs$total_cn <- snvs$major_cn + snvs$minor_cn
    
    snvs <- attach_001_ground_truth_snvs(wd_truth, snvs)
    
    snvin <- read.delim(paste0(resultFolder, '/', basename, '/', basename, '_filtered_snvs.tsv'))
    snvin$mutation_id <- paste(snvin$chrom, snvin$pos, sep='_')
    snvs <- inner_join(snvs, snvin[,c('mutation_id', 'gtype')], by='mutation_id')
    
    snvs$true_ccf <- 1
    snvs[snvs$sample%in%'bM_only','true_ccf'] <- p1
    snvs[snvs$sample%in%'gM_only','true_ccf'] <- 1-p1
    
    snvs$adjusted_support <- snvs$var_counts
    snvs$adjusted_depth <- snvs$var_counts + snvs$ref_counts
    
    snvs$norm_cn <- snvs$normal_cn
    x <- t(data.frame(apply(snvs, 1, get_true_sc_cn)))
    snvs$true_cn <- x[,1]
    snvs$wtotal_cn <- x[,2]
    snvs$pv <- x[,3]
    snvs$best_ccf <- x[,4]
    
    return(snvs)
}

proc_pyclone_results <- function(pycloneFolder, purity, burnIn=1000) {
    # load trace and mpear label
    traceFile <- dir(paste0(pycloneFolder, "/trace"),
                     pattern = "cellular_prevalence", full.names = T)
    paramsTrace <- read.delim(traceFile, stringsAsFactors = F, header = F)
    idx <- as.character(paramsTrace[1,])
    ssm <- data.frame(mutation_id = idx , stringsAsFactors = F)
    
    id <- do.call(rbind, strsplit(as.character(ssm$mutation_id), "_", fixed = T))
    
    paramsTrace <- as.matrix( paramsTrace[-1:-(burnIn+1), ])
    class(paramsTrace) <- "numeric"
    
    mpearFile <- paste0(pycloneFolder, "/pyclone_mpear.tsv")
    
    if (file.exists(mpearFile)) {
        mpear <- read.delim(paste0(pycloneFolder, "/pyclone_mpear.tsv"), stringsAsFactors = F)
    } else {
        labelTrace <- read.delim(gzfile(paste0(pycloneFolder, 
                                               "/trace/labels.tsv.bz2")), stringsAsFactors = F)
        mpearLabels <- compute_mpear_label(labelTrace[-1:-burnIn,])
        mpear <- data.frame(mutation_id = colnames(labelTrace), cluster_id = mpearLabels)
        write.table(mpear, file = paste0(pycloneFolder, '/pyclone_mpear.tsv'),
                    row.names = F, sep = "\t", quote = F)
        rm(labelTrace)
    }
    
    allData <- GetCcfFromLabel(t(paramsTrace), idx, mpear$cluster_id, func='mean')
    write.table(allData, file = paste0(pycloneFolder, '/pyclone_results.tsv'),
                row.names = F, sep = "\t", quote = F)
}