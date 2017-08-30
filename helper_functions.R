library(ggplot2)
library(gridExtra)
library(data.table)
library(plyr)
library(GenomicRanges)
library(gtools)
library(RColorBrewer)

subclonal_cutoff <- 0.9

adjust_vafs <- function(dat, pur) {
    mut_prop <- dat$prop_chrs_bearing_mutation
    vaf <- dat$adjusted_vaf
    frac <- dat$cn_frac
  
    cn_v <- dat$most_likely_variant_copynumber * frac
    cn_r <- dat$most_likely_ref_copynumber * (1 - frac)
    cn_r <- sapply(cn_r,function(x){if(x==1){0}else{x}})
    
    pn <- (1 - pur) * 2
    pr <- pur * cn_r
    pv <- pur * cn_v
    
    norm_const <- pn + pr + pv
    pv <- pv / norm_const

    prob <- pv * mut_prop
    adj_vaf <- (1 / prob) * vaf
    adj_vaf[adj_vaf > 2] <- 2
    
    return(adj_vaf)
}

get_pv <- function(pur, phi, cn_v, cn_r, mut_prop, frac) {
    cn_v <- cn_v * frac
    cn_r <- cn_r * (1 - frac)
    
    pn <- (1 - pur) * 2
    pr <- pur * cn_r
    pv <- pur * cn_v
    
    norm_const <- pn + pr + pv
    pv <- pv / norm_const
    pv <- phi * pv * mut_prop

    return(pv)
}

add_combos <- function(combos, vmaj, vmin, ref_cn, frac) {
    vtotal = vmaj + vmin
    if (vtotal == 0) {
        mu_v <- 0
    } else {
        for (i in 1:vmaj) {
            mu_v <- 1*i / vtotal
            combos <- rbind(combos, data.frame(cn_r=ref_cn, cn_v=vtotal, mu_v=mu_v, frac=frac))
        }
    }
    return(combos)
}

get_combos <- function(sc, pur, target=1) {
    sc <- strsplit(as.character(sc), '\\|')[[1]]
    if (length(sc) > 1) {
        sc1 <- as.numeric(strsplit(sc[1], ',')[[1]])
        sc2 <- as.numeric(strsplit(sc[2], ',')[[1]])
        m1 <- sc1[1]; n1 <- sc1[2]; t1 <- sc1[1]+sc1[2]; f1 <- sc1[3]
        m2 <- sc2[1]; n2 <- sc2[2]; t2 <- sc2[1]+sc2[2]; f2 <- sc2[3]
        combos <- add_combos(NULL, m1, n1, t2, f1)
        combos <- add_combos(combos, m2, n2, t1, f2)
        combos <- unique(combos)
    } else {
        sc1 <- as.numeric(strsplit(sc, ',')[[1]])
        m <- sc1[1]; n <- sc1[2]; f <- sc1[3]; t <- m+n
        combos <- add_combos(NULL, m, n, t, f)
    }
    combos <- unique(combos)
    for (j in 1:nrow(combos)) {
        combos$pv[j] <- get_pv(pur, target, combos$cn_v[j], combos$cn_r[j], combos$mu_v[j], combos$frac[j])
    }
    return(combos)
}

adjust_vafs_pyclone <- function(dat, pur) {
    cn <- dat$tumour_copynumber
    mut_prop <- dat$multiplicity
    prob <- (cn * mut_prop * pur) / (2 * (1 - pur) + cn * pur)
    vaf <- dat$vaf
    adj_vaf <- (1 / prob) * vaf
    adj_vaf[adj_vaf>2] <- 2
    return(adj_vaf)
}

get_frac <- function(x, snvs) {
    variant_cn <- x['most_likely_variant_copynumber']
    side <- x['preferred_side']
    gtype <- x['gtype']
    if (!snvs) {
        if(side==0) {
            gtype <- x['gtype1']
        } else {
            gtype <- x['gtype2']
        }
    }
    sc <- strsplit(gtype, '\\|')[[1]]
    if (length(sc) == 0){return(1)}
    sc <- strsplit(sc,',')
    sc1 <- as.numeric(sc[[1]])
    sc1 <- c(sum(sc1[1:2]),sc1[3])
    if(length(sc) > 1) {
        sc2 <- as.numeric(sc[[2]])
        sc2 <- c(sum(sc2[1:2]),sc2[3])
        if (sc1[1] == variant_cn) {
            return(sc1[2])
        } else {
            return(sc2[2])
        }
    }
    return(sc1[2])
}

load_pyclone_sc <- function(pyc_wd, pur_wd, p1, p2) {
    base <- paste('001bM_p', p1/100, '_001gM_p', p2/100, sep='')
    mix <- paste(p1, p2, sep='-')
    filename_sc <- paste(pyc_wd, base,'/',base,'_subclonal_structure.txt.gz',sep='')
    pur_file <- paste(pur_wd, base,'_cellularity_ploidy.txt',sep='')
    sc <- read.delim(file=gzfile(filename_sc))
    pp <- read.delim(file=pur_file)
    sc$proportion <- sc$proportion/pp$cellularity
    sc$mix <- paste(mix, 'pyc', sep='_')
    colnames(sc)[3] <- 'CCF'
    sc <- sc[, c('mix','cluster','n_ssms','CCF')]
    sc$variant_proportion <- sc$n_ssms/sum(sc$n_ssms)
    return(sc)
}

load_pyclone_points <- function(pyc_wd, p1, p2) {
    base <- paste('001bM_p', p1/100, '_001gM_p', p2/100, sep='')
    mix <- paste(p1, p2, sep='-')
    
    res_file <- paste(pyc_wd, base, '_pyclone_results_table.csv', sep='')
    dat <- read.delim(res_file, stringsAsFactors = F, sep=',')
    
    ccf <- data.frame(sv=paste(dat$chr, dat$pos, sep='_'), mix=paste(mix, 'pyc', sep='_'),
                      adjusted_support=dat$var_counts, adjusted_depth=dat$var_counts+dat$ref_counts,
                      CCF=dat$ccf, VAF=dat$vaf, cluster=dat$cluster_id, sample=NA)
    return(ccf)
}

load_coclus_sc <- function(snv_wd, p1, p2) {
    base <- paste('001bM_p', p1/100, '_001gM_p', p2/100, sep='')
    mix <- paste(p1, p2, sep='-')
    
    sc_file <- paste(snv_wd, base, '/best_run_coclus_post_assign/', base, '_subclonal_structure.txt', sep='')
    sc <- read.delim(sc_file, stringsAsFactors = F, sep='\t')
    
    sc <- data.frame(mix=mix, cluster=sc$cluster, 
                     n_ssms=sc$n_ssms, CCF=sc$CCF, variant_proportion=sc$proportion)
    return(sc)
}

get_data_svc2 <- function(mixes, wd, snvs=FALSE) {
    ccf_all <- NULL
    sc_all <- NULL
    for(i in 1:nrow(mixes)) {
        base_name <- paste('001bM_p', mixes[i,1]/100, '_001gM_p', mixes[i,2]/100, sep = '')
        mix <- paste(mixes[i,1], mixes[i,2], sep='-')
        # min_bic_run <- 'best_run_svs'        
        min_bic_run <- get_best_run(wd, base_name, 'svc_IC')
        min_bic_run <- paste('/', min_bic_run, '/', sep='')
        
        scs_file <-  paste(wd, base_name, '/', min_bic_run, '/', base_name, '_subclonal_structure.txt', sep = '')
        scs <- read.table(scs_file, sep = '\t', header = T)
        scs <- scs[scs$n_ssms>1, ]
        scs <- scs[order(scs$CCF, decreasing = T), ]
        scs$new_cluster <- 1:nrow(scs)

        sv_df <- read.table(paste(wd, base_name, '/', base_name, '_filtered_svs.tsv', sep=''), 
                            header=T, sep='\t', stringsAsFactors=F)
        pur <- read.table(paste(wd, base_name, '/purity_ploidy.txt', sep=''), header=T, sep='\t', stringsAsFactors=F)$purity

        cc_file <-  paste(wd, base_name, '/', min_bic_run, '/', base_name, '_cluster_certainty.txt', sep = '')
        cc <- read.table(cc_file, sep = '\t', stringsAsFactors = F, header = T)
        
        mlcn_file <- paste(wd, base_name, '/', min_bic_run, '/', base_name, '_most_likely_copynumbers.txt', sep='')
        mlcn <- read.table(mlcn_file, header = T, sep = '\t', stringsAsFactors = F)

        merge_cols <- c('chr1', 'pos1', 'dir1', 'chr2', 'pos2', 'dir2')
        dat <- merge(sv_df, mlcn, by.x=merge_cols, by.y=merge_cols)
        dat <- merge(dat, cc, by.x=merge_cols, by.y=merge_cols)
        dat$cn_frac <- apply(dat, 1, function(x){get_frac(x,snvs)})
        dat <- cbind(dat, CCF=adjust_vafs(dat, pur))
        
        dat$cluster <- NA
        for (j in 1:nrow(scs)) {
            dat[dat$most_likely_assignment==scs$cluster[j], 'cluster'] <- scs$new_cluster[j]
        }
        
        sv_ids <- paste(dat$chr1, dat$pos1, dat$dir1, 
                        dat$chr2, dat$pos2, dat$dir2, sep=':')
        ccfs <- data.frame(sv = sv_ids, mix = paste(mix, 'svs', sep='_'),
                           adjusted_support = dat$adjusted_support, adjusted_depth = dat$adjusted_depth, 
                           CCF = dat$CCF, VAF = dat$adjusted_vaf, cluster = dat$cluster)
        ccf_all <- rbind(ccf_all, ccfs)
    
        scs$cluster <- scs$new_cluster
        scs <- scs[,1:4]
        scs <- data.frame(mix = paste(mix, 'svs', sep='_'), scs[-3], variant_proportion=scs$n_ssms / sum(scs$n_ssms))
        sc_all <- rbind(sc_all, scs)
    }
    ccf_all <- ccf_all[!duplicated(ccf_all),]
    return(list(sc_all, ccf_all))
}

get_runs <- function(wd) {
    runs <- c()
    for (dir in list.dirs(wd)) {
        cur_dir <- strsplit(dir,'/')[[1]]
        if (length(cur_dir) <= 1)
            next
        cur_dir <- cur_dir[length(cur_dir)]
        if (substring(cur_dir,1,3) == 'run') {
            runs <- c(runs,cur_dir)
        }
    }
    return(unique(runs))
}

get_ic_table <- function(wd, base_name, runs, allowed_ics = c('BIC', 'AIC', 'AICc', 'svc_IC'), clus_penalty = 4, snvs = FALSE) {
    snv_pref <- ''
    samp_dir <- paste(wd, base_name, sep='/')

    if (snvs) {snv_pref <- 'snvs/'}
    ic_table <- NULL
    for (run in runs) {
        ic <- read.table(paste(samp_dir, '/', run, '/', snv_pref, base_name, '_fit.txt', sep=''), sep='\t', header=F, stringsAsFactors = F)        
        sc <- read.table(paste(samp_dir, '/', run, '/', snv_pref, base_name, '_subclonal_structure.txt', sep=''), sep='\t', header=T, stringsAsFactors = F)
        
        accept_solution <- !(nrow(sc) == 1 & sc[1,'CCF'] < 0.9)
        ic$V3 <- accept_solution

        ic <- ic[ic$V1%in%allowed_ics,]
        ic <- cbind(run=run, ic)
        ic_table <- rbind(ic_table, ic)
    }
    ic_table$run <- factor(ic_table$run, levels=mixedsort(unique(as.character(ic_table$run))))
    return(ic_table)
}

get_best_run <- function(wd, base_name, ic_metric, snvs = FALSE, clus_penalty = 4) {
    runs <- get_runs(wd)
    ic <- get_ic_table(wd, base_name, runs, clus_penalty = clus_penalty, snvs = snvs)
    
    min_ic <- min(ic[ic$V1==ic_metric,'V2'])
    best_run <- as.character(ic[ic$V2==min_ic & ic$V1==ic_metric,'run'])[1]
    return(best_run)
}

get_min_dic_run <- function(wd, id, runs) {
    ic_table <- get_ic_table(wd, id, runs)
    ic_table <- ic_table[ic_table$DIC>0, ]
    min_ic_run <- as.character(ic_table[min(ic_table$DIC)==ic_table$DIC,'run'])
    return(min_ic_run)
}

get_data_svc3 <- function(mixes, wd) {
    ccf_all <- NULL
    sc_all <- NULL
    for(i in 1:nrow(mixes)) {
        base_name <- paste('001bM_p', mixes[i,1]/100, '_001gM_p', mixes[i,2]/100, sep = '')
        mix <- paste(mixes[i,1], mixes[i,2], sep='-')
        
        samp_dir <- paste(wd, base_name, sep='/')
        runs <- get_runs(samp_dir)
        min_dic_run <- get_min_dic_run(wd, base_name, runs)
        
        print(paste('picking', min_dic_run, 'for', base_name))
        
        scs_file <-  paste(samp_dir, '/', min_dic_run, '/', base_name, '_subclonal_structure.txt', sep = '')
        scs <- read.table(scs_file, sep = '\t', header = T)
        scs <- scs[scs$n_ssms>1, ]
        scs <- scs[order(scs$CCF, decreasing = T), ]
        scs$new_cluster <- 1:nrow(scs)
        
        sv_df <- read.table(paste(wd, base_name, '/', base_name, '_filtered_svs.tsv', sep=''), 
                            header=T, sep='\t', stringsAsFactors=F)
        pur <- read.table(paste(wd, base_name, '/purity_ploidy.txt', sep=''), header=T, sep='\t', stringsAsFactors=F)$purity
        
        cc_file <-  paste(wd, base_name, '/', min_dic_run, '/', base_name, '_cluster_certainty.txt', sep = '')
        cc <- read.table(cc_file, sep = '\t', stringsAsFactors = F, header = T)
        
        mlcn_file <- paste(wd, base_name, '/', min_dic_run, '/', base_name, '_most_likely_copynumbers.txt', sep='')
        mlcn <- read.table(mlcn_file, header = T, sep = '\t', stringsAsFactors = F)
        
        merge_cols <- c('chr1', 'pos1', 'dir1', 'chr2', 'pos2', 'dir2')
        dat <- merge(sv_df, mlcn, by.x=merge_cols, by.y=merge_cols)
        dat <- merge(dat, cc, by.x=merge_cols, by.y=merge_cols)
        dat$cn_frac <- apply(dat, 1, function(x){get_frac(x)})
        dat <- cbind(dat, CCF=adjust_vafs(dat, pur))
        
        dat$cluster <- NA
        for (j in 1:nrow(scs)) {
            dat[dat$most_likely_assignment==scs$cluster[j], 'cluster'] <- scs$new_cluster[j]
        }
        
        ccfs <- data.frame(sv = apply(dat[,c('chr1','pos1','chr2','pos2')], 1, paste, collapse=':'), 
                           mix = paste(mix, 'svs', sep='_'), CCF = dat$CCF, cluster = dat$cluster)
        ccf_all <- rbind(ccf_all, ccfs)
        
        scs$cluster <- scs$new_cluster
        scs <- scs[,1:4]
        scs <- data.frame(mix = paste(mix, 'svs', sep='_'), scs[-3], variant_proportion=scs$n_ssms / sum(scs$n_ssms))
        sc_all <- rbind(sc_all, scs)
    }
    return(list(sc_all, ccf_all))
}

attach_001_ground_truth_svs <- function(wd_truth, ccf_all, ff=15) {
    bm_svs <- read.table(paste(wd_truth, "001bM_merged_list_svinfo.txt", sep = ''), sep='\t', stringsAsFactors = F, header = T)
    gm_svs <- read.table(paste(wd_truth, "001gM_merged_list_svinfo.txt", sep = ''), sep='\t', stringsAsFactors = F, header = T)
    
    # check SVs that don't match exactly
    bm_svs$id <- paste(bm_svs$chr1, bm_svs$pos1, bm_svs$dir1, 
                       bm_svs$chr2, bm_svs$pos2, bm_svs$dir2, sep=':')
    gm_svs$id <- paste(gm_svs$chr1, gm_svs$pos1, gm_svs$dir1,
                       gm_svs$chr2, gm_svs$pos2, gm_svs$dir2, sep=':')
    
    both <- c()
    for (i in 1:nrow(bm_svs)) {
        bm_chr1 <- bm_svs[i,'chr1']
        bm_chr2 <- bm_svs[i,'chr2']
        bm_pos1 <- bm_svs[i,'pos1']
        bm_pos2 <- bm_svs[i,'pos2']
        bm_dir1 <- bm_svs[i,'dir1']
        bm_dir2 <- bm_svs[i,'dir2']
        if (bm_svs$support[i] == 0) {
            next
        }
        for (j in 1:nrow(gm_svs)) {
            if (gm_svs$support[j] == 0) {
                next
            }            
            gm_chr2 <- gm_svs[j,'chr2']
            gm_chr1 <- gm_svs[j,'chr1']
            gm_pos1 <- gm_svs[j,'pos1']
            gm_pos2 <- gm_svs[j,'pos2']
            gm_dir1 <- gm_svs[j,'dir1']
            gm_dir2 <- gm_svs[j,'dir2']
            if (gm_chr1 == bm_chr1 & gm_chr2 == bm_chr2) {                
                if (bm_dir1 == gm_dir1 & bm_dir2 == gm_dir2) {
                    if (abs(gm_pos1 - bm_pos1) < ff & abs(gm_pos2 - bm_pos2) < ff) {
                        both <- c(both, bm_svs$id[i], gm_svs$id[j])
                    }
                }
            }
        }
    }
    both <- unique(both)
    bm_uniq <- bm_svs$id[(!bm_svs$id%in%both) & bm_svs$support>0]
    gm_uniq <- gm_svs$id[(!gm_svs$id%in%both) & gm_svs$support>0]
    
    ccf_all$sample <- NA
    ccf_all[ccf_all$sv%in%bm_uniq, 'sample'] <- 'bM_only'
    ccf_all[ccf_all$sv%in%gm_uniq, 'sample'] <- 'gM_only'
    ccf_all[ccf_all$sv%in%both, 'sample'] <- 'shared'
    ccf_all <- ccf_all[!is.na(ccf_all$sample), ]
    
    return(ccf_all)
}

attach_001_ground_truth_snvs <- function(snv_wd, ccfs) {
    bm_snvs <- suppressWarnings(fread(paste(snv_wd, "001bM_SS_mut_merged.call.stats.txt", sep = ''), 
                     sep='\t', stringsAsFactors = F, header = T))
    gm_snvs <- suppressWarnings(fread(paste(snv_wd, "001gM_SS_mut_merged.call.stats.txt", sep = ''), 
                     sep='\t', stringsAsFactors = F, header = T))
    bm <- bm_snvs$Row.names
    gm <- gm_snvs$Row.names
    
    both <- intersect(bm, gm)
    gm_uniq <- gm[!gm%in%bm]
    bm_uniq <- bm[!bm%in%gm]
    
    ccfs[ccfs$sv%in%bm_uniq, 'sample'] <- 'bM_only'
    ccfs[ccfs$sv%in%gm_uniq, 'sample'] <- 'gM_only'
    ccfs[ccfs$sv%in%both, 'sample'] <- 'shared'
    ccfs[is.na(ccfs$sample), 'sample'] <- 'unmatched'
    unmatched <- ccfs[ccfs$sample == 'unmatched', ]
    ccfs <- ccfs[ccfs$sample!='unmatched', ]
    
    return(ccfs)
}

get_truth <- function(sc_all, ccf_all) {
    truth <- data.frame(table(ccf_all[,c('mix','sample')]))
    truth$mix <- sapply(truth$mix, as.character)
    truth$sample <- sapply(truth$sample, as.character)
    
    truth_sv <- ccf_all
    truth_sv$mix <- sapply(truth_sv$mix, as.character)
    truth_sv$sample <- sapply(truth_sv$sample, as.character)
    
    truth$CCF <- NA
    truth$cluster <- NA
    for (i in 1:nrow(truth)) {
        this_mix <- truth$mix[i]
        this_sample <- truth$sample[i]
        sam_split <- strsplit(this_mix,'-')[[1]]
        sam_split[2] <- strsplit(sam_split[2],'_')[[1]][1]
        sam_split <- as.numeric(sapply(sam_split,as.numeric))
        if (this_sample == 'bM_only') {
            real_ccf <- sam_split[1]/100
            truth$CCF[i] <- real_ccf
            truth_sv$CCF[truth_sv$mix == this_mix & truth_sv$sample == this_sample] <- real_ccf
            if(sam_split[1] > sam_split[2]){
                truth$cluster[i] <- 2
                truth_sv$cluster[truth_sv$mix == this_mix & truth_sv$sample == this_sample] <- 2
            } else {
                truth$cluster[i] <- 3
                truth_sv$cluster[truth_sv$mix == this_mix & truth_sv$sample == this_sample] <- 3
            }
        } else if (this_sample == 'gM_only') {
            real_ccf <- sam_split[2]/100
            truth$CCF[i] <- real_ccf
            truth_sv$CCF[truth_sv$mix == this_mix & truth_sv$sample == this_sample] <- real_ccf
            if(sam_split[2] > sam_split[1]){
                truth$cluster[i] <- 2
                truth_sv$cluster[truth_sv$mix == this_mix & truth_sv$sample == this_sample] <- 2
            } else {
                truth$cluster[i] <- 3
                truth_sv$cluster[truth_sv$mix == this_mix & truth_sv$sample == this_sample] <- 3
            }
        } else {
            truth$CCF[i] <- 1.0
            truth_sv$CCF[truth_sv$mix == this_mix & truth_sv$sample == this_sample] <- 1.0
            truth$cluster[i] <- 1
            truth_sv$cluster[truth_sv$mix == this_mix & truth_sv$sample == this_sample] <- 1
        }
    }
    truth$mix <- paste(truth$mix, 'truth', sep='_')
    truth_sv$mix <- paste(truth_sv$mix, 'truth', sep='_')
    truth <- truth[, c('mix', 'cluster', 'Freq', 'CCF')]
    colnames(truth)[3] <- 'n_ssms'

    #get ssm proportion
    truth$variant_proportion <- apply(truth, 1, function(x){as.numeric(x['n_ssms']) / sum(truth[truth$mix == x['mix'], 'n_ssms'])})
    return(list(truth, truth_sv))
}

get_purity <- function(wd, id) {
    pur <- read.table(paste(wd, '/purity_ploidy.txt', sep=''),
                      header=T, sep='\t', stringsAsFactors=F)$purity
    return(pur)
}

get_purity_mixes <- function(wd, p1, p2, postfix='') {
    base_name <- paste('001bM_p', p1/100, '_001gM_p', p2/100, sep = '')
    pur <- read.table(paste(wd, base_name, postfix, '/purity_ploidy.txt', sep=''), 
               header=T, sep='\t', stringsAsFactors=F)$purity
    return(pur)
}

get_ploidy_mixes <- function(wd, p1, p2, postfix='') {
    base_name <- paste('001bM_p', p1/100, '_001gM_p', p2/100, sep = '')
    pl <- read.table(paste(wd, base_name, postfix, '/purity_ploidy.txt', sep=''), 
                      header=T, sep='\t', stringsAsFactors=F)$ploidy
    return(pl)
}

get_gtype <- function(x) {
    if(as.numeric(x['frac1_A'])==1) {
        return(paste(x['nMaj1_A'], x['nMin1_A'], x['frac1_A'],sep=','))
    } else {
        return(paste(paste(x['nMaj1_A'], x['nMin1_A'], x['frac1_A'],sep=','),
                     paste(as.character(x['nMaj2_A']), as.character(x['nMin2_A']), as.character(x['frac2_A']),sep=','),sep='|'))
    }
}

get_matched_cnv_locs <- function(wd, p1, p2, truth_sv, side=0) {
    mix <- paste(p1, p2, sep='-')
    base_name <- paste('001bM_p', p1/100, '_001gM_p', p2/100, sep = '')
    
    pur <- get_purity_mixes(wd, p1, p2)
    pl <- get_ploidy_mixes(wd, p1, p2)
    
    sv_file <- paste(wd, base_name, '/', base_name, '_filtered_svs_post_assign.tsv', sep='')
    x <- read.delim(sv_file, header=T, stringsAsFactors = F)
    x$sv <- paste(x$chr1, x$pos1, x$dir1, x$chr2, x$pos2, x$dir2, sep=':')
    x <- suppressWarnings(merge(x, truth_sv, by='sv', all.x=T))
    
    off <- 100000
    x_sv <- NULL
    for (i in 1:nrow(x)) {
        class <- x[i, 'classification']
        if (side==0) {
            # if (x[i,'preferred_side']==0) {
            chr   <- x[i, 'chr1']
            pos   <- x[i, 'pos1']
            gtype <- x[i, 'gtype1']
            if (class == 'INTRX') {
                if (x[i, 'dir1'] == '-') {pos <- pos + off}
                else {pos <- pos - off}
            } else {
                pos <- pos - off
            }
            x_sv <- rbind(x_sv, data.frame(chr=chr, pos=pos, gtype=gtype, truth=x$sample[i]))
        } else {
            chr   <- x[i, 'chr2']
            pos   <- x[i, 'pos2']
            gtype <- x[i, 'gtype2']
            if (class == 'INTRX') {
                if (x[i, 'dir1'] == '-') {pos <- pos + off}
                else {pos <- pos - off}
            } else {
                pos <- pos + off
            }
            x_sv <- rbind(x_sv, data.frame(chr=chr, pos=pos, gtype=gtype, truth=x$sample[i]))
        }
    }
    af <- as.numeric(sapply(x$classification, function(x){if(x%in%c('DUP','INTDUP')){return(1-(pur/pl))}else{return(1)}}))
    x_sv$adjusted_norm <- (x$norm1*af)
    if(side==1){x_sv$adjusted_norm <- x$norm2*af}
    x_sv$adjusted_support <- x$adjusted_support
    x_sv$adjusted_depth <- x_sv$adjusted_norm+x_sv$adjusted_support
    x_sv$adjusted_vaf <- x_sv$adjusted_support / x_sv$adjusted_depth
    x_sv$id <- x$sv
    return(x_sv)
}

get_snv_wcnv <- function(wd, ccf_all, base_name, mix) {
    x_snv <- ccf_all[ccf_all$mix==mix,]
    x_snv$chr <- sapply(x_snv$sv, function(x){strsplit(as.character(x),'_')[[1]][1]})
    x_snv$pos <- as.numeric(sapply(x_snv$sv, function(x){strsplit(as.character(x),'_')[[1]][2]}))
    
    sc <- read.delim(paste(wd, base_name, '_subclones.txt', sep=''), stringsAsFactors = F)
    sc_gx <- GRanges(seqnames = sc$chr, ranges = IRanges(start = sc$startpos, end = sc$endpos))
    snv_gx <- GRanges(seqnames = x_snv$chr, ranges = IRanges(start = as.numeric(x_snv$pos)-1, end = as.numeric(x_snv$pos)))
    
    bb_hits <- sc[subjectHits(findOverlaps(snv_gx, sc_gx)),]
    gtypes <- NULL;for(i in 1:nrow(bb_hits)){gtypes <- c(gtypes, get_gtype(bb_hits[i,]))}
    x_snv <- x_snv[queryHits(findOverlaps(snv_gx, sc_gx)),]
    x_snv$gtype <- gtypes
    colnames(x_snv)[c(6,8)] <- c('adjusted_vaf', 'truth')
    
    return(x_snv)
}

compare_cnv_truth <- function(x_sv, wd_truth) {
    # compare "truth" vs. assigned copy-numbers
    bm_sc <- read.delim(paste(wd_truth, '001bM_subclones.txt', sep=''), header=T, stringsAsFactors = F)[,1:15]
    gm_sc <- read.delim(paste(wd_truth, '001gM_subclones.txt', sep=''), header=T, stringsAsFactors = F)[,1:15]
    
    bm_gx <- GRanges(seqnames = bm_sc$chr, ranges = IRanges(start = bm_sc$startpos, end = bm_sc$endpos))
    gm_gx <- GRanges(seqnames = gm_sc$chr, ranges = IRanges(start = gm_sc$startpos, end = gm_sc$endpos))
    
    x_gx <- GRanges(seqnames = x_sv$chr, ranges = IRanges(start = x_sv$pos, end = x_sv$pos+1))
    x_sv[queryHits(findOverlaps(x_gx, bm_gx)),'bM_hit'] <- subjectHits(findOverlaps(x_gx, bm_gx))
    x_sv[queryHits(findOverlaps(x_gx, gm_gx)),'gM_hit'] <- subjectHits(findOverlaps(x_gx, gm_gx))
    
    x_sv$cnv_true <- NA
    for (i in 1:nrow(x_sv)) {
        sc <- as.character(x_sv[i,'gtype'])
        gm_hit <- x_sv[i, 'gM_hit']
        bm_hit <- x_sv[i, 'bM_hit']
        tmp <- strsplit(sc, '\\|')[[1]]
        if (length(tmp) == 1) {
            sc <- as.numeric(strsplit(tmp[1], ',')[[1]])
        } else {
            sc <- c(as.numeric(strsplit(tmp[1], ',')[[1]]), as.numeric(strsplit(tmp[2], ',')[[1]]))
        }
        if (!is.na(gm_hit) & gm_sc[gm_hit, 'frac1_A'] != 1) {
            x_sv[i, 'cnv_true'] <- 'subclonal'
        } 
        if (!is.na(bm_hit) & bm_sc[bm_hit, 'frac1_A'] != 1) {
            x_sv[i, 'cnv_true'] <- 'subclonal'
        }
        x_sv[i, 'cnv_true'] <- trimws(x_sv[i, 'cnv_true'], which = 'left')
        if (!is.na(gm_hit) & !is.na(bm_hit) & is.na(x_sv[i, 'cnv_true'])) {
            if (length(sc) == 3) {
                if (gm_sc[gm_hit, 'nMaj1_A'] == sc[1] & gm_sc[gm_hit, 'nMin1_A'] == sc[2]) {
                    if (bm_sc[bm_hit, 'nMaj1_A'] == sc[1] & bm_sc[bm_hit, 'nMin1_A'] == sc[2]) {
                        x_sv[i, 'cnv_true'] <- 'match'
                    }
                }
            } else if (gm_sc[gm_hit, 'nMaj1_A'] == sc[1] & gm_sc[gm_hit, 'nMin1_A'] == sc[2]) {
                if (bm_sc[bm_hit, 'nMaj1_A'] == sc[4] & bm_sc[bm_hit, 'nMin1_A'] == sc[5]) {
                    x_sv[i, 'cnv_true'] <- 'match'
                }
            } else if (gm_sc[gm_hit, 'nMaj1_A'] == sc[4] & gm_sc[gm_hit, 'nMin1_A'] == sc[5]) {
                if (bm_sc[bm_hit, 'nMaj1_A'] == sc[1] & bm_sc[bm_hit, 'nMin1_A'] == sc[2]) {
                    x_sv[i, 'cnv_true'] <- 'match'
                }
            }
        }
    }
    x_sv[is.na(x_sv$cnv_true), 'cnv_true'] <- 'no_match'
    return(x_sv[!is.na(x_sv$truth),])
}

get_best_cn <- function(combos, bx, bn) {
    lls <- NULL
    for (j in 1:nrow(combos)) {
        pv <- combos$pv[j]
        lls <- c(lls, bx * log(pv) + (bn - bx) * log(1 - pv))
    }
    cbest <- combos[which(max(lls,na.rm=T)==lls),]
    return(cbest)
}

estimate_best_ccfs <- function(x_sv, pur, p1, p2) {
    x_sv$best_ccf <- NA
    x_sv$num_cn_states <- NA
    x_sv$cn_v <- NA
    x_sv$cn_r <- NA
    x_sv$mu_v <- NA
    x_sv$pv <- NA
    x_sv$target <- NA
    x_sv$cn_frac <- NA
    for (i in 1:nrow(x_sv)) {
        bx <- as.numeric(x_sv$adjusted_support[i])
        bn <- as.numeric(x_sv$adjusted_depth[i])        
        target <- 1
        if (x_sv$truth[i]=='bM_only') {
            target <- p1/100
        } else if (x_sv$truth[i]=='gM_only') {
            target <- p2/100
        }
        x_sv$target[i] <- target
        combos <- get_combos(x_sv$gtype[i], target)
        x_sv$num_cn_states[i] <- nrow(combos)
        if (nrow(combos) == 1) {
            cbest <- combos[1,]
        } else {
            cbest <- get_best_cn(combos, bx, bn)
        }
        x_sv$most_likely_variant_copynumber[i] <- cbest$cn_v
        x_sv$most_likely_ref_copynumber[i] <- cbest$cn_r
        x_sv$prop_chrs_bearing_mutation[i] <- cbest$mu_v
        x_sv$cn_frac[i] <- cbest$frac
        x_sv$pv[i]   <- cbest$pv
    }
    x_sv$best_ccf <- adjust_vafs(x_sv, pur)
    x_sv <- x_sv[!duplicated(x_sv),]
    return(x_sv)
}

get_dat <- function(id, pur) {
    vaf_ccf_file <- paste(run_dir, 'SVs/',id,'_vaf_ccf.txt',sep='')
    ccert_file <- paste(run_dir, 'SVs/',id,'_cluster_certainty.txt',sep='')
    mlcn_file <- paste(run_dir, 'SVs/',id,'_most_likely_copynumbers.txt',sep='')
    svinfo_file <- paste(svinfos, id, '_svinfo.txt',sep='')
    
    if(!file.exists(vaf_ccf_file) | !file.exists(ccert_file) & 
       !file.exists(mlcn_file) & !file.exists(svinfo_file)){return(NULL)}
    
    dat <- read.delim(vaf_ccf_file, stringsAsFactors = F)[,2:9]
    ccert <- read.delim(ccert_file, stringsAsFactors = F)
    mlcn <- read.delim(mlcn_file, stringsAsFactors = F)
    svinfo <- read.delim(svinfo_file, stringsAsFactors = F)[,-1]
    
    dat <- merge(dat, ccert, by=c(1:6))
    dat <- merge(dat, mlcn, by=c(1:6))
    colnames(dat)[colnames(dat)=='support'] <- 'adjusted_support'
    dat <- merge(dat, svinfo, by=c(1:6))
    
    colnames(dat)[8] <- 'adjusted_vaf'
    dat$CCF <- (1 / (dat$pv / (dat$average_proportion/pur))) * dat$adjusted_vaf
    dat$CCF[dat$CCF>2] <- 2
    dat$cluster_CCF <- dat$average_proportion / pur
    
    return(dat)
}

get_dat_snvs <- function(id, pur) {
    vaf_ccf_file <- paste(run_dir, 'SNVs/',id,'_vaf_ccf.txt',sep='')
    ccert_file <- paste(run_dir, 'SNVs/',id,'_cluster_certainty.txt',sep='')
    mlcn_file <- paste(run_dir, 'SNVs/',id,'_most_likely_copynumbers.txt',sep='')
    
    if(!file.exists(vaf_ccf_file) | !file.exists(ccert_file) & !file.exists(mlcn_file)){return(NULL)}
    
    dat <- read.delim(vaf_ccf_file, stringsAsFactors = F)
    ccert <- read.delim(ccert_file, stringsAsFactors = F)
    mlcn <- read.delim(mlcn_file, stringsAsFactors = F)
    
    dat <- merge(dat, ccert, by=c(1:2))
    dat <- merge(dat, mlcn, by=c(1:2))
    
    dat$adjusted_vaf <- dat$raw_VAF
    dat$CCF <- (1 / (dat$pv / (dat$average_proportion/pur))) * dat$raw_VAF
    dat$CCF[dat$CCF>2] <- 2
    
    return(dat)
}

get_sc <- function(id) {
    sc <- read.delim(paste(run_dir, 'SVs/',id,'_subclonal_structure.txt',sep=''))
    return(sc)
}

get_cumvars <- function(x, pur) {
    vaf <- x$adjusted_vaf/pur
    vaf <- vaf[vaf > 0.12 & vaf < 0.24]
    if (length(vaf)==0) {
        return(NULL)            
    } else {
        vaf <- sort(1/vaf)
        histinfo <- hist(vaf, breaks=40)
        cumvars <- rep(0,length(vaf))
        for (i in 1:length(vaf)) {
            cumvars[i] <- sum(histinfo$breaks<vaf[i])
        }
        return(list(vaf,cumvars))
    }
}

draw_circos_clonal <- function(dat, sample, pur, bbf) {
    bb <- read.delim(bbf, sep='\t', stringsAsFactors = F)
    
    pdat <- data.frame(chr=bb$chromosome, startpos=bb$start, end=bb$end, value=bb$total_cn)
    pdat$chr <- paste('chr', pdat$chr,sep='')
    colnames(pdat) <- c('chr','start','end','value')
    pdat$value[pdat$value > 6] <- 6
    
    colours <- c('#0000FF80','#FF000080','darkgreen','#0000FF40','#FF000040','#00FF0040')
    sample <- strsplit(sample, '-')[[1]][1]
    # pdf(paste('scnr_plots/circos/', sample, '_circos.pdf', sep=''), height=12, width=12, pointsize = 30)
    par(mar = c(1, 1, 1, 1))
    circos.clear()
    circos.initializeWithIdeogram(plotType = c('axis','labels'), major.by = 1e100)
    circos.genomicTrackPlotRegion(pdat,ylim=c(0,6),track.height = 0.1,
                                  panel.fun=function(region,value,...){
                                      i=getI(...)
                                      circos.genomicLines(region,value,type='segment',lwd=3,col=colours[i],...)
                                  })
    for(j in 1:nrow(dat))
    {
        x <- dat[j,]
        ccf <- x$average_proportion / pur
        lcol <- colours[1]; if(ccf < subclonal_cutoff){lcol <- colours[2]}
        chr1 <- paste('chr',as.character(x[1]),sep='')
        chr2 <- paste('chr',as.character(x[4]),sep='')
        pos1 <- as.numeric(x[2]); pos2 <- as.numeric(x[5])
        circos.link(chr1, pos1, chr2, pos2, col=lcol, lwd=3)
    }
    # dev.off()
}

############################################################################################################################
# Plotting functions
############################################################################################################################

get_run_info <- function(wd, base_name, run, postfix = '', snvs = FALSE, renum = TRUE) {
    snv_pref <- ''
    if (snvs) {snv_pref <- 'snvs/'}
    scs_file <-  paste(wd, base_name, postfix, '/', run, '/', snv_pref, base_name, '_subclonal_structure.txt', sep = '')

    scs <- read.table(scs_file, sep = '\t', header = T)
#     scs <- scs[scs$n_ssms>3,]
    scs <- scs[order(scs$CCF, decreasing=T), ]
    
    sv_df <- ''
    if(length(grep('post_assign',run))>0){
        if (snvs) {
            sv_df <- read.table(paste(wd, base_name, postfix, '/', base_name, '_filtered_snvs_post_assign.tsv', sep=''), header=T, sep='\t', stringsAsFactors=F)
        } else {
            sv_df <- read.table(paste(wd, base_name, postfix, '/', base_name, '_filtered_svs_post_assign.tsv', sep=''), header=T, sep='\t', stringsAsFactors=F)
        }
    } else {
        if (snvs) {
            sv_df <- read.table(paste(wd, base_name, postfix, '/', base_name, '_filtered_snvs.tsv', sep=''), header=T, sep='\t', stringsAsFactors=F)
        } else {
            sv_df <- read.table(paste(wd, base_name, postfix, '/', base_name, '_filtered_svs.tsv', sep=''), header=T, sep='\t', stringsAsFactors=F)
        }
    }

    pur <- read.table(paste(wd, base_name, postfix, '/purity_ploidy.txt', sep=''), 
                      header=T, sep='\t', stringsAsFactors=F)$purity
    
    cc_file <-  paste(wd, base_name, postfix, '/', run, '/', snv_pref, base_name, '_cluster_certainty.txt', sep = '')
    cc <- read.table(cc_file, sep = '\t', stringsAsFactors = F, header = T)
    
    mlcn_file <- paste(wd, base_name, postfix, '/', run, '/', snv_pref, base_name, '_most_likely_copynumbers.txt', sep='')
    mlcn <- read.table(mlcn_file, header = T, sep = '\t', stringsAsFactors = F)
    
    merge_cols <- c('chr1', 'pos1', 'dir1', 'chr2', 'pos2', 'dir2')
    if (snvs) {
        colnames(sv_df)[1] <- 'chr'
        colnames(cc)[1] <- 'chr'
        merge_cols <- c('chr', 'pos')
        sv_df$adjusted_vaf <- sv_df$var / (sv_df$ref + sv_df$var)
    }
    dat <- merge(sv_df, mlcn, by.x=merge_cols, by.y=merge_cols)
    dat <- merge(dat, cc, by.x=merge_cols, by.y=merge_cols)
    dat$cn_frac <- apply(dat, 1, function(x){get_frac(x, snvs)})
    dat <- cbind(dat, CCF=adjust_vafs(dat, pur))
    
    if(renum) {
        scs$new_cluster <- 1:nrow(scs)
        dat$cluster <- NA
        for (j in 1:nrow(scs)) {
            dat[dat$most_likely_assignment==scs$cluster[j], 'cluster'] <- scs$new_cluster[j]
        }
    } 
    if (snvs) {
        dat$sv <- paste(dat$chr, dat$pos, sep='_')
    } else {
        dat$sv <- paste(dat$chr1, dat$pos1, dat$dir1, dat$chr2, dat$pos2, dat$dir2, sep=':')
    }
    dat <- dat[!duplicated(dat),]
    if(renum){scs$cluster <- scs$new_cluster}
    scs <- scs[,1:4]
    if('n_variants'%in%colnames(scs)){scs$variant_proportion <- scs$n_variants/sum(scs$n_variants)}
    else{scs$variant_proportion <- scs$n_ssms/sum(scs$n_ssms)}
    return(list(dat,scs))
}

plot_mix <- function(wd, wd_truth, p1, p2, postfix='', snvs=FALSE, sample_factor=TRUE, lim=2) {
    base_name <- paste('001bM_p', p1/100, '_001gM_p', p2/100, sep = '')
    mix <- paste(p1, p2, sep='-')
    
    runs <- list.files(paste(wd, base_name, postfix, sep=''))
    runs <- runs[grep('^run', runs)]
    
    all_runs_ccfs <- c()
    all_runs_scs <- c()
    
    var_id <- c('chr1','pos1','dir1','chr2','pos2','dir2')
    if (snvs) {var_id <- c('chr', 'pos')}
    
    for(run in runs) {
        runinf <- get_run_info(wd, base_name, run, postfix, snvs)
        dat <- runinf[[1]]
        if (snvs) {
            dat <- attach_001_ground_truth_snvs(wd_truth, dat)
        } else {
            dat <- attach_001_ground_truth_svs(wd_truth, dat)
        }
    
        scs <- runinf[[2]]
        scs$run <- run
        ccfs <- data.frame(sv = apply(dat[,var_id], 1, paste, collapse=':'), 
                           mix = paste(mix, 'svs', sep='_'), CCF = dat$CCF, cluster = dat$cluster,
                           run = run, sample = dat$sample)
        all_runs_ccfs <- rbind(all_runs_ccfs, ccfs)
        all_runs_scs <- rbind(all_runs_scs, scs)
    }
    
    all_runs_ccfs$run <- factor(all_runs_ccfs$run, levels=mixedsort(unique(as.character(all_runs_ccfs$run))))
    p <- ggplot(data = all_runs_scs) + theme_minimal() +
        scale_y_continuous(breaks = seq(0, lim, 0.2), limits = c(0,lim)) + xlab('')
    if (sample_factor) {
        p <- p + geom_line(data = all_runs_ccfs, aes(x = run, y = CCF, group = sv, 
                                                    colour = factor(sample)), alpha=0.2) +
            geom_jitter(data = all_runs_ccfs, aes(x = run, y = CCF, colour = factor(sample)), 
                        position=position_jitter(width=0.4), alpha = 0.5, size = 2) +
            scale_size(range = c(4,25), guide = FALSE) + 
            geom_point(aes(x = run, y = CCF, size=variant_proportion), colour = '#4d4d4d', alpha=0.8) +
            ggtitle(paste(mix, 'mix'))
    } else {
        p <- p + geom_line(data = all_runs_ccfs,aes(x = run, y = CCF, group = sv), alpha=0.05) + 
            geom_jitter(data = all_runs_ccfs, aes(x = run, y = CCF, colour = factor(cluster)), 
                        position=position_jitter(width=0.4), alpha = 0.5, size = 2) +
            scale_size(range = c(4,20), guide = FALSE) + 
            geom_point(aes(x = run, y = CCF, size=variant_proportion), colour = '#4d4d4d', alpha=0.8) +
            ggtitle(paste(mix, 'mix'))
    }
    return(p)
}

plot_ics <- function(wd, p1, p2, allowed_ics = c('BIC', 'AIC', 'AICc'), best_run_metric = 'BIC', ps='', snvs=FALSE) {
    snv_pref <- ''
    if (snvs) {snv_pref <- 'snvs/'}
    base_name <- paste('001bM_p', p1/100, '_001gM_p', p2/100, sep = '')
    mix <- paste(p1, p2, sep='-')
    
    runs <- get_runs(paste(wd, paste(base_name, ps, sep=''), sep='/'))
    ic_table <- get_ic_table(wd, base_name, runs, allowed_ics)
    
    best_run <- get_best_run(wd, base_name, ic_metric = best_run_metric)
    
    p <- ggplot(ic_table, aes(y=V2, x=run, group=V1, colour=factor(V1))) + geom_line() + theme_minimal() + ggtitle(best_run) + ylab('')
    return(p)
}

plot_nscs <- function(wd, p1, p2) {
    base_name <- paste('001bM_p', p1/100, '_001gM_p', p2/100, sep = '')
    mix <- paste(p1, p2, sep='-')
    samp_dir <- paste(wd, base_name, sep='/')
    
    nsc <- NULL
    runs <- get_runs(paste(wd, paste(base_name, ps, sep=''), sep='/'))
    for (run in runs) {
        sc <- read.table(paste(samp_dir, '/', run, '/', base_name, '_subclonal_structure.txt', sep=''), 
                         sep='\t', header=T, stringsAsFactors = F)
        nsc <- rbind(nsc, data.frame(run=run, n_clusts=nrow(sc)))
    }
    
    nsc$run <- factor(nsc$run, levels=mixedsort(unique(as.character(nsc$run))))
    p <- ggplot(nsc, aes(y=n_clusts, x=run, group=1)) + geom_line() + theme_minimal() + ylab('')
    return(p)
}

plot_lnls <- function(wd, p1, p2) {
    base_name <- paste('001bM_p', p1/100, '_001gM_p', p2/100, sep = '')
    mix <- paste(p1, p2, sep='-')
    
    runs <- get_runs(paste(wd, paste(base_name, ps, sep=''), sep='/'))
    ic_table <- get_ic_table(wd, base_name, runs, allowed_ics = c('lnL'))
    
    max_lnl <- as.character(ic_table[which(max(ic_table$V2)==ic_table$V2),'run'])
    max_lnl <- paste('Max lnl:', max_lnl, '; Mix:', mix)

    p <- ggplot(ic_table, aes(y=V2, x=run, group=V1, colour=factor(V1))) + geom_line() + theme_minimal() + ggtitle(max_lnl) + ylab('')
    return(p)
}

plot_pyc_hist <- function(sc_all, ccf_all, p1, p2) {
    mix <- paste(p1, '-', p2, '_pyc', sep='')
    
    clus_intercepts <- sc_all[sc_all$mix==mix,]$CCF
    pyc <- ccf_all[ccf_all$mix==mix,]
    
    p <- ggplot(pyc, aes(x=CCF, fill=factor(cluster), color=factor(cluster))) +
        theme_minimal() +
        xlim(0,2) + geom_histogram(alpha=0.3,position='identity',binwidth=0.05) + xlab('CCF') +
        geom_vline(xintercept=clus_intercepts, colour='blue', size=1) + ylab('') +
        scale_fill_brewer(palette = 'Set1', name = "Cluster") +
        scale_color_brewer(palette = 'Set1', name = "Cluster") + 
        theme(axis.title = element_text(size = 18),
              axis.text.x = element_text(size = 16),
              axis.text.y = element_text(size = 16)) +
        theme(legend.position="none") + ggtitle(paste(p1, '-', p2, 'mix'))
    
    return(p)
}

plot_pyc_mix <- function(sc_all, ccf_all, p1, p2, lim=2) {
    mix <- paste(p1, p2, sep='-')
    p_sc <- sc_all[grep(mix, sc_all$mix),]
    p_ccf <- ccf_all[grep(mix, ccf_all$mix),]
    
    p_sc$type <- 'sv'
    p_sc[grep('pyc',p_sc$mix),'type'] <- 'pyc'
    
    # collapse sv/pyc truth into one
    p_sc$mix <- as.character(p_sc$mix)
    cats <- unique(as.character(p_sc$mix))
    truth_label <- paste(strsplit(cats[grep('truth',cats)[1]], '_')[[1]][c(1,3)],collapse='_')
    p_sc[grep('truth', p_sc$mix), 'mix'] <- truth_label
    
    cats <- unique(as.character(p_sc$mix))
    p_sc$mix <- factor(p_sc$mix, levels = cats[c(2,3,1)])
    
    p_ccf$mix <- as.character(p_ccf$mix)
    p_ccf[grep('truth', p_ccf$mix), 'mix'] <- truth_label
    p_ccf$mix <- as.factor(p_ccf$mix)
    
    cols <- c('pyc' = '#4d4d4d', 
              'sv' = '#bf812d', 
              'bM_only' = '#1b9e77', 
              'gM_only' = '#d95f02', 
              'shared' = '#7570b3',
              'black' = '#000000')
    
    #rename p_sc and p_ccf to remove mix in name
    levels <- c('svs', 'truth', 'pyc')
    p_sc_split <- sapply(as.character(p_sc$mix), strsplit, split = '_')
    p_ccf_split <- sapply(as.character(p_ccf$mix), strsplit, split = '_')
    
    p_sc$mix <- factor(data.frame(p_sc_split, stringsAsFactors=F)[2,], levels = levels)
    p_ccf$mix <- factor(data.frame(p_ccf_split,stringsAsFactors = F)[2,], levels = levels)
    
    p_ccf$mix <- revalue(p_ccf$mix, c('svs'='SVclone', 'pyc'='PyClone', 'truth'='Truth'))
    p_sc$mix <- revalue(p_sc$mix, c('svs'='SVclone', 'pyc'='PyClone', 'truth'='Truth'))
    
    p_sc_var <- p_sc[grep('SV|Py',p_sc$mix),]
    p_sc_tr <- p_sc[grep('Truth',p_sc$mix),]
    
    # take maxes of truth sets
    p_sc_tr <- data.frame(data.table(p_sc_tr)[, j=list(max(n_ssms), 
                                                       max(CCF), max(variant_proportion)),by = list(mix, cluster, type)])
    colnames(p_sc_tr)[4:6] <- c('n_ssms', 'CCF', 'variant_proportion')
    p_sc_tr <- p_sc_tr[,c(1,2,4:6,3)]
    
    p <- ggplot(data = p_sc) + theme_minimal() +
        geom_point(aes(x = mix, y = CCF, size=variant_proportion), alpha=0) +
        #^ this geom point is just a hack to plot the order correctly
        scale_y_continuous(breaks = seq(0, lim, 0.2), limits = c(0,lim)) + xlab('') +
        geom_jitter(data = p_ccf[grep('Py',p_ccf$mix),], aes(x = mix, y = CCF, colour = factor(sample)), 
                    position=position_jitter(width=0.40), alpha = 0.2, size = 0.5) +
        geom_jitter(data = p_ccf[grep('SV',p_ccf$mix),], aes(x = mix, y = CCF, colour = factor(sample)), 
                    position=position_jitter(width=0.25), alpha = 0.7, size = 2) +
        geom_point(data = p_sc_tr, 
                   aes(x = mix, y = CCF, size=variant_proportion, fill = 'black'), colour = 'white', pch=21) +
        scale_color_manual(values = cols, guide = FALSE) + scale_size(range = c(4,25), guide = FALSE) + 
        scale_fill_manual(values = cols, guide = FALSE) +
        geom_line(data = p_ccf[grep('Py|Truth',p_ccf$mix),], 
                  aes(x = mix, y = CCF, group = sv, colour = factor(sample)), alpha=0.025) +
        geom_line(data = p_ccf[grep('SV|Truth',p_ccf$mix),], 
                  aes(x = mix, y = CCF, group = sv, colour = factor(sample)), alpha=0.2) +
        geom_point(data = p_sc_var, 
                   aes(x = mix, y = CCF, size=variant_proportion, fill = factor(type), colour='black'), pch = 21, alpha=0.6) +
        theme(axis.title = element_text(size = 16, face = "bold"),
              axis.text.x = element_text(size = 18, face = "bold"),
              axis.text.y = element_text(size = 18)) + ggtitle(paste(mix, 'Mix'))
    return(p)
}

plot_pyc_clusts_only <- function(sc_all, p1, p2, lim=1.2) {
    mix <- paste(p1, p2, sep='-')
    p_sc <- sc_all[grep(mix, sc_all$mix),]
    p_sc$type <- 'sv'
    p_sc[grep('pyc',p_sc$mix),'type'] <- 'pyc'
    
    # collapse sv/pyc truth into one
    p_sc$mix <- as.character(p_sc$mix)
    cats <- unique(as.character(p_sc$mix))
    truth_label <- paste(strsplit(cats[grep('truth',cats)[1]], '_')[[1]][c(1,3)],collapse='_')
    p_sc[grep('truth', p_sc$mix), 'mix'] <- truth_label
    
    cats <- unique(as.character(p_sc$mix))
    p_sc$mix <- factor(p_sc$mix, levels = cats[c(2,3,1)])
    
    cols <- c('pyc' = '#4d4d4d', 
              'sv' = '#bf812d', 
              'bM_only' = '#1b9e77', 
              'gM_only' = '#d95f02', 
              'shared' = '#7570b3',
              'black' = '#000000')
    
    # rename p_sc to remove mix in name
    levels <- c('svs', 'truth', 'pyc')
    p_sc_split <- sapply(as.character(p_sc$mix), strsplit, split = '_')
    p_sc$mix <- factor(data.frame(p_sc_split, stringsAsFactors=F)[2,], levels = levels)
    p_sc$mix <- revalue(p_sc$mix, c('svs'='SVclone', 'pyc'='PyClone', 'truth'='Truth'))
    p_sc_var <- p_sc[grep('SV|Py',p_sc$mix),]
    p_sc_tr <- p_sc[grep('Truth',p_sc$mix),]
    
    # take maxes of truth sets
    p_sc_tr <- data.frame(data.table(p_sc_tr)[, j=list(max(n_ssms), 
                                                       max(CCF), max(variant_proportion)),by = list(mix, cluster, type)])
    colnames(p_sc_tr)[4:6] <- c('n_ssms', 'CCF', 'variant_proportion')
    p_sc_tr <- p_sc_tr[,c(1,2,4:6,3)]
    
    p <- ggplot(data = p_sc) + theme_minimal() +
        geom_point(aes(x = mix, y = CCF, size=variant_proportion), alpha=0) +
        #^ this geom point is just a hack to plot the order correctly
        scale_y_continuous(breaks = seq(0, lim, 0.2), limits = c(0,lim)) + xlab('') +
        geom_point(data = p_sc_tr, 
                   aes(x = mix, y = CCF, size=variant_proportion, fill = 'black'), colour = 'white', pch=21) +
        scale_color_manual(values = cols, guide = FALSE) + scale_size(range = c(4,15), guide = FALSE) + 
        scale_fill_manual(values = cols, guide = FALSE) +
        geom_point(data = p_sc_var, 
                   aes(x = mix, y = CCF, size=variant_proportion, fill = factor(type), colour='black'), 
                   pch = 21, alpha=0.6) +
        theme(plot.title = element_text(size = 16, face = "bold"),
              axis.title = element_text(size = 16),
              axis.text.x = element_text(size = 18),
              axis.text.y = element_text(size = 18)) + 
        theme(legend.position="none")
    return(p)
}

plot_ccf_hist <- function(wd, base_name, p1, p2, pick_run='best', snvs = F, ylim = 25) {
    mix <- paste(p1, p2, sep='-')
    if (pick_run=='best') {
        pick_run <- get_best_run(wd, base_name, 'svc_IC', snvs = snvs)
    }

    x <- get_run_info(wd, base_name, pick_run, postfix = '', snvs = snvs, renum = F)
    pur <- read.delim(paste(wd, base_name, '/purity_ploidy.txt', sep=''))$purity
    dat <- x[[1]]
    sc <- x[[2]]
    
    clus_intercepts <- 1 / pur * as.numeric(sc$proportion)
    
    # make factor order consistent
    shared_clus <- sc$cluster[sc$CCF==max(sc$CCF)]
    sc_nonc <- sc[sc$cluster!=shared_clus,]
    if(p1 < p2) {
        bm_clus <- sc_nonc$cluster[sc_nonc$CCF==max(sc_nonc$CCF)]
        gm_clus <- sc_nonc$cluster[sc_nonc$CCF==min(sc_nonc$CCF)]
    } else {
        bm_clus <- sc_nonc$cluster[sc_nonc$CCF==min(sc_nonc$CCF)]
        gm_clus <- sc_nonc$cluster[sc_nonc$CCF==max(sc_nonc$CCF)]
    }
    dat$most_likely_assignment <- factor(dat$most_likely_assignment, levels=c(shared_clus, gm_clus, bm_clus))    
    
    variant <- 'SVs'
    if(snvs){variant <- 'SNVs'}
    ccf_hist <- ggplot(dat, aes(x=as.numeric(dat$CCF),
                                fill=factor(most_likely_assignment),color=factor(most_likely_assignment))) +
        theme_minimal() +
        xlim(0,2) + geom_histogram(alpha=0.3,position='identity',binwidth=0.05)+xlab('CCF') +
        geom_vline(xintercept=clus_intercepts, colour='blue', size=1) + ylab('') + ylim(0,ylim) +
        scale_fill_brewer(palette = 'Set1', name = "Cluster") +
        scale_color_brewer(palette = 'Set1', name = "Cluster") + 
        theme(axis.title = element_text(size = 22),
              axis.text.x = element_text(size = 20),
              axis.text.y = element_text(size = 20)) +
        theme(legend.position="none")
    
    return(ccf_hist)
}

plot_ccf_hist_snvs <- function(ccf_all, sc_all, p1, p2, ylim=450) {
    mix <- paste(p1, p2, sep='-')
    dat <- ccf_all[ccf_all$mix == paste(mix, 'pyc', sep='_'), ]
    sc <- sc_all[sc_all$mix == paste(mix, 'pyc', sep='_'), ]
    clus_intercepts <- as.numeric(sc$CCF)
    
    ccf_hist <- ggplot(dat, aes(x=as.numeric(dat$CCF),
                                fill=factor(cluster),color=factor(cluster))) +
        theme_minimal() +
        xlim(0,2) + geom_histogram(alpha=0.3,position='identity',binwidth=0.05) + xlab('CCF') +
        geom_vline(xintercept=clus_intercepts, colour='blue', size=1) + ylab('') + ylim(0,ylim) +
        scale_fill_brewer(palette = 'Set1', name = "Cluster") +
        scale_color_brewer(palette = 'Set1', name = "Cluster") + 
        theme(axis.title = element_text(size = 22),
              axis.text.x = element_text(size = 20),
              axis.text.y = element_text(size = 20)) +
        theme(legend.position="none")

    return(ccf_hist)
}

plot_truth_histogram <- function(p1, p2, truth_sv, run) {
    mix <- paste(p1, p2, sep='-')
    base_name <- paste('001bM_p', p1/100, '_001gM_p', p2/100, sep = '')
    x <- get_run_info(wd, base_name, paste('/', run, '/', sep=''), postfix = '', snvs = F)
    dat <- x[[1]]
    sc <- x[[2]]
    
    pp <- read.delim(paste(wd, base_name, '/purity_ploidy.txt', sep=''))
    pur <- pp$purity
    
    above_ssm_th <- sc$n_ssms / (sum(sc$n_ssms)) >= 0.04
    below_ssm_th <- sc$n_ssms / (sum(sc$n_ssms)) < 0.04
    clus_intercepts <- 1 / pur * as.numericget(sc$proportion[above_ssm_th & sc$n_ssms > 2])
    
    dat <- merge(dat, truth_sv, by='sv', all.x=T)
    ccf_hist <- ggplot(dat, aes(x=as.numeric(dat$CCF),
                              fill=factor(sample),color=factor(sample))) +
    theme_minimal() +
    xlim(0,2) + geom_histogram(alpha=0.3,position='identity',binwidth=0.05)+xlab('CCF') +
    geom_vline(xintercept=clus_intercepts, colour='blue', size=1) + ylab('') + ylim(0,20) +
    scale_fill_brewer(palette = 'Set1', name = "Cluster") +
    scale_color_brewer(palette = 'Set1', name = "Cluster") + 
    theme(plot.title = element_text(size = 16, face = "bold"),
          axis.title = element_text(size = 16),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18)) +
    ggtitle(paste(mix, 'Mix', run))
    
    return(ccf_hist)
}

plot_best_ccfs <- function(x_sv, p1, p2, ylim=20) {
    x_sv$truth <- factor(x_sv$truth, levels=c('shared', 'bM_only', 'gM_only'))
    p <- ggplot(x_sv, aes(x=as.numeric(x_sv$best_ccf),fill=factor(truth),color=factor(truth))) +
        theme_minimal() +
        geom_histogram(alpha=0.3,position='identity',binwidth=0.05) + xlab('CCF') + xlim(0,2) +
        scale_fill_brewer(palette = 'Set1', name = "Cluster") +
        scale_color_brewer(palette = 'Set1', name = "Cluster") + 
        theme(axis.title = element_text(size = 22),
              axis.text.x = element_text(size = 20),
              axis.text.y = element_text(size = 20)) +
        geom_vline(xintercept=p1/100, colour='blue', size=1) + ylim(0,ylim) +
        geom_vline(xintercept=p2/100, colour='blue', size=1) +
        geom_vline(xintercept=1, colour='blue', size=1) + ylab('') +
        theme(legend.position="none")
    return(p)
}

generate_circos <- function(wd, p1, p2, subclones_dir, colours, clonal=TRUE, subclonal=TRUE) {
    suppressMessages(library(circlize))    
    pur <- get_purity_mixes(wd, p1, p2)
    base_name <- paste('001bM_p', p1/100, '_001gM_p', p2/100, sep = '')
    bb <- read.table(paste(subclones_dir, '/', base_name,"_subclones.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
    clon <- bb[is.na(bb$nMaj2_A),]
    subclon <- bb[!is.na(bb$nMaj2_A),]    
    
    pdat<-clon[,c("chr","startpos","endpos")]
    pdat<-cbind(pdat,value=clon$nMaj1_A+clon$nMin1_A)
    pdat$chr<-paste("chr",pdat$chr,sep="")
    colnames(pdat)<-c("chr","start","end","value")
    pdat$value[pdat$value>6]<-6
    
    pdat2<-subclon[,c("chr","startpos","endpos")]
    pdat2<-cbind(pdat2,value=apply(cbind(subclon$nMaj2_A+subclon$nMin2_A,subclon$nMaj1_A+subclon$nMin1_A),1,mean))
    pdat2$chr<-paste("chr",pdat2$chr,sep="")
    pdat2$value[pdat2$value>6]<-6
    colnames(pdat2)<-c("chr","start","end","value")    
    pdat<-list(pdat,pdat2)
    
    dat <- get_run_info(wd, base_name, 'best_run_svs', postfix='', snvs=F)[[1]]
        
    par(mar = c(1, 1, 1, 1))
    circos.initializeWithIdeogram(plotType = c("axis"))
    circos.genomicTrackPlotRegion(pdat,ylim=c(0,6),
                                  panel.fun=function(region,value,...){
                                      i=getI(...)
                                      circos.genomicLines(region,value,type="segment",lwd=3,col=colours[i],...)
                                  })    
    subclones <- unique(dat[dat$average_proportion/pur < subclonal_cutoff, 'most_likely_assignment'])
    for(j in 1:nrow(dat))
    {
        x<-dat[j,]
        if(x$average_proportion/pur > subclonal_cutoff) {
            lcol=colours[1]
            if(!clonal) {next}
        } else if (x$most_likely_assignment == subclones[1]){
            lcol=colours[2]
            if (!subclonal) {next}
        } else {
            lcol=colours[3]
            if (!subclonal) {next}
        }
        circos.link(paste("chr",as.character(x[1]),sep=""),
                    as.numeric(x[2]),
                    paste("chr",as.character(x[4]),sep=""),
                    as.numeric(x[5]),col=lcol,lwd=3)
    }
}

plot_hist_func <- function(dat, sc, pur, varclass=FALSE, vaf=FALSE, clus=-1, title='', ylim=0) {
    var_types <- c('DEL', 'DUP', 'INTDUP', 'INTRX', 'INV', 'SNV')
    dat$classification <- factor(dat$classification, levels=var_types)
    set1 <- brewer.pal(7, 'Set1')[-6][6:1] #remove yellow
    type_cols <- set1; names(type_cols) <- var_types
    
    if (clus>=0) {
        dat <- dat[dat$most_likely_assignment==clus,]
        sc <- sc[sc$cluster==clus,]
    }    

    clus_intercepts <- sc$CCF
    dat$most_likely_assignment <- factor(dat$most_likely_assignment)
    
    plotvar <- 'Cluster'
    if (varclass) {
        plotvar <- 'Classification'
        if (!'classification'%in%colnames(dat)) {
            dat$classification <- rep('SNV', nrow(dat))
        }
        if (vaf) {
            var_hist <- ggplot(dat, aes(x=adjusted_vaf,fill=classification,color=classification)) +
                geom_histogram(alpha=0.3,binwidth=0.05) + xlab('VAF')
        } else {
            var_hist <- ggplot(dat, aes(x=CCF,fill=classification,color=classification)) +
                geom_histogram(alpha=0.3,binwidth=0.05) + xlab('CCF')
        }
    } else {
        if (vaf) {
            var_hist <- ggplot(dat, aes(x=adjusted_vaf,fill=most_likely_assignment,color=most_likely_assignment)) +
                geom_histogram(alpha=0.3,position='identity',binwidth=0.05) + xlab('VAF')
        } else {
            var_hist <- ggplot(dat, aes(x=CCF,fill=most_likely_assignment,color=most_likely_assignment)) +
                geom_histogram(alpha=0.3,position='identity',binwidth=0.05) + xlab('CCF')
        }
    }
    
    if (ylim > 0) {var_hist <- var_hist + ylim(0, ylim)}
    var_hist <- var_hist + theme_minimal() + xlim(0,2) + ggtitle(title) + ylab('') +        
        theme(plot.title = element_text(size = 20),
              axis.title = element_text(size = 18),
              axis.text.x = element_text(size = 18),
              axis.text.y = element_text(size = 18),
              legend.text = element_text(size = 18),
              legend.position = 'bottom')
    if (varclass) {
        var_hist <- var_hist + scale_fill_manual(values = type_cols) + 
            scale_colour_manual(values = type_cols)
    } else {
        var_hist <- var_hist + scale_fill_brewer(palette = 'Set1', name = plotvar) + 
            scale_color_brewer(palette = 'Set1', name = plotvar)
    }
    
    if (!vaf) {
        var_hist <- var_hist + geom_vline(xintercept=clus_intercepts, colour='blue', size=1)
    }
    return(var_hist)
}

plot_hist <- function(wd, id, are_snvs=FALSE, run='best', varclass=FALSE, vaf=FALSE, clus=-1, title='', ylim=0) {
    if (run=='best') {
        run <- get_best_run(wd, id, 'svc_IC', snvs = are_snvs)
    }

    x <- get_run_info(wd, id, run, renum = FALSE)
    dat <- x[[1]]
    sc <- x[[2]]
    
    if (are_snvs) {
        y <- get_run_info(wd, id, run, snvs = are_snvs)
        dat_y <- y[[1]][,c('sv', 'most_likely_assignment', 'CCF')]
        dat_y$adjusted_vaf <- y[[1]]$var / (y[[1]]$ref + y[[1]]$var)
        dat_y$classification <- 'SNV'

        dat <- dat[,c('sv', 'most_likely_assignment', 'CCF', 'adjusted_vaf', 'classification')]
        dat <- rbind(dat, dat_y)
    }
    pur <- get_purity(wd, id)
    plot_hist_func(dat, sc, pur, varclass, vaf, clus, title, ylim)
}