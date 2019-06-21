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
    scs <- scs[,c('cluster', 'n_variants', 'proportion', 'CCF')]
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
    
    # run_files <- list.files(paste(wd, paste(base_name, ps, sep=''), 'best_run_svs', sep='/'))
    # best_run <- strsplit(run_files[grep('^run',run_files)],'\\.')[[1]][1]
    # best_run <- paste('Best run:', best_run)
    # best_run <- ic_table[min(ic_table$V2)==ic_table$V2,'run']
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
        # sc <- sc[sc$n_ssm/sum(sc$n_ssm)>=0.05,]
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
    
    #p_sc <- p_sc[-c(grep('snv',p_sc$mix)),]
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
        #ggtitle(paste(mix, 'Mix')) + 
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
    # dat$most_likely_assignment <- factor(dat$most_likely_assignment, levels=c(shared_clus, gm_clus, bm_clus))    
    dat$most_likely_assignment <- factor(dat$most_likely_assignment)    
    
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
    
    # above_ssm_th <- sc$n_ssms / (sum(sc$n_ssms)) >= 0.04
    # below_ssm_th <- sc$n_ssms / (sum(sc$n_ssms)) < 0.04
    # clus_intercepts <- as.numeric(sc$CCF[above_ssm_th & sc$n_ssms > 2])
    clus_intercepts <- as.numeric(sc$CCF)
    # clus_intercepts_minor <- 1 / as.numeric(sc$CCF[below_ssm_th | sc$n_ssms<=2])
    
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
    #ggtitle(paste(mix, 'Mix'))
    
    #     if (length(clus_intercepts_minor) > 0) {
    #         ccf_hist <- ccf_hist + geom_vline(xintercept=clus_intercepts_minor,colour='red',lty=2)
    #     }
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
    clus_intercepts <- 1 / pur * as.numeric(sc$proportion[above_ssm_th & sc$n_ssms > 2])
    # clus_intercepts_minor <- 1 / pur * as.numeric(sc$proportion[below_ssm_th | sc$n_ssms<=2])
    
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
    
    #     if (length(clus_intercepts_minor) > 0) {
    #         ccf_hist <- ccf_hist + geom_vline(xintercept=clus_intercepts_minor,colour='red',lty=2)
    #     }
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
    #ggtitle(paste(p1, '-', p2, 'Mix'))
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
    
    # clus_intercepts <- 1 / pur * as.numeric(sc$proportion)
    # clus_intercepts <- unique(dat$average_proportion/pur)
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

plot_ccube_mix <- function(mydata, res, purity, mix) {
    par(mfrow=c(2,2))
    uniqLabels <- unique(res$label)
    plot(mydata$ccube_ccf, mydata$vaf, col = myColors[res$label], 
         xlab = "cancer cell fraction", ylab = "variant allele frequecy", 
         main = "ccf vs vaf (colored by cluster memebership)")
    mydata$total_cn =mydata$major_cn+mydata$minor_cn
    uniqueTotCn = unique(mydata$total_cn)
    xx = seq(0,2, length.out = 100)
    for (cn in uniqueTotCn) {
        for (i in 1:cn) {
            points(MapVaf2CcfPyClone(xx, purity, 2, cn, cn, i, constraint = F), xx, type = 'l')
        }
    }
    
    Emu <- res$full.model$ccfMean
    Esigma <- res$full.model$ccfCov
    Epi <- res$full.model$Epi
    params <- data.frame(Emu=as.numeric(Emu), Esigma=as.numeric(Esigma), Epi=as.numeric(Epi))
    xx <- seq(0, 2,  length.out = 1000)
    ll <- 0
    ll1 <- 0 
    
    for (j in seq_len(nrow(params))) {
        ll <- ll + params[j,]$Epi * dnorm(xx, mean = params[j,]$Emu, sd = sqrt(params[j,]$Esigma))
    }
    
    hist(mydata$ccube_ccf, density=20, breaks=20, prob=TRUE, 
         main = "ccf histogram +
         fitted marginal denstiy (red)",
         xlab = "cancer cell fraction")
    lines(xx,ll, lwd=2, col = "darkred")
    
    names(Epi) <- as.character(format(round(Emu, 2), nsmall = 2))
    barplot(Epi[sort(uniqLabels)], las = 2, col = myColors[sort(uniqLabels)], 
            xlab = "cluster mean", ylab="expected weight", 
            main = paste("cluster weights", paste('(', mix, ':', 1-mix, ')', sep='')))
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
        # lcol=colours[x$most_likely_assignment+1]
        chr1 <- paste('chr',as.character(x[1]),sep='')
        chr2 <- paste('chr',as.character(x[4]),sep='')
        pos1 <- as.numeric(x[2]); pos2 <- as.numeric(x[5])
        circos.link(chr1, pos1, chr2, pos2, col=lcol, lwd=3)
    }
    # dev.off()
}