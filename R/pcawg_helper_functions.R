get_sv_hits <- function(samples, pcawg_pp, svinfos, sv_dir, pcawg_drivers, drivers_only=T) {
    svs_list <- vector("list", length(samples))
    for(i in 1:length(samples)){
        sample <- samples[i]
        svinfo_file <- paste(svinfos,sample,'_svinfo.txt',sep='')
        ccf_file <- paste(sv_dir, sample, '/ccube_out/', sample, '_cluster_certainty.txt', sep='')
        svi <- paste(svinfos, sample, '_svinfo.txt', sep='')        
        
        if(!file.exists(svinfo_file) | !file.exists(ccf_file) | !file.exists(svi)) {
            next
        }
        
        pur <- pcawg_pp[pcawg_pp$samplename==sample,'purity']
        # x <- get_dat(sample, pur)
        x <- read.delim(ccf_file)
        x$average_proportion <- apply(x[,c('average_proportion1', 'average_proportion2')], 1, mean)
        x$CCF <- x$average_proportion / pur
        svinf <- read.delim(svi)
        x$chr1 <- as.character(x$chr1); x$chr2 <- as.character(x$chr2)
        svinf$chr1 <- as.character(svinf$chr1); svinf$chr2 <- as.character(svinf$chr2)
        x <- left_join(x, svinf, by=c('chr1', 'pos1', 'dir1', 'chr2', 'pos2', 'dir2'))
        
        svs_list[[i]] <- data.frame(x, sample=sample)
    }
    all_svs <- NULL
    all_svs <- rbind(all_svs, do.call(rbind, svs_list))
    
    grx_tmp <- grx
    grx_ex_tmp <- grx_ex
    if (drivers_only) {
        grx_tmp <- grx[grx$genes%in%pcawg_drivers$gene]
        grx_ex_tmp <- grx_ex[grx_ex$genes%in%pcawg_drivers$gene]
    }
    gene_df <- data.frame(grx_tmp@ranges)
    gene_df$genes <- grx_tmp$genes
    
    # inter-chromsomal translocations
    int_svs <- all_svs[all_svs$classification%in%c('INTRX'),]
    
    chroms <- c(int_svs$chr1, int_svs$chr2)
    starts <- c(int_svs$pos1, int_svs$pos2)
    ends <- c(int_svs$pos1, int_svs$pos2)
    
    int_x <- GRanges(seqnames = chroms, ranges = IRanges(start = starts, end = ends))
    olaps <- suppressWarnings(findOverlaps(int_x, grx_tmp))
    
    int_svs <- rbindlist(list(int_svs, int_svs), use.names = T)
    hits <- int_svs[queryHits(olaps)]
    hits$genes <- grx_tmp[subjectHits(olaps),]$genes
    if (nrow(hits)>0) {
        hits$exon_start <- NA
        hits$exon_end <- NA
    }
    hits <- merge(hits, gene_df, by='genes')
    
    # other
    other_svs <- all_svs[!all_svs$classification%in%c('INTRX'),]
    other_x <- GRanges(seqnames = other_svs$chr1, 
                       ranges = IRanges(start = other_svs$pos1, end = other_svs$pos2))
    
    olaps <- findOverlaps(other_x, grx_ex_tmp)
    hits2 <- data.table(other_svs)[queryHits(olaps)]
    hits2$genes <- grx_ex_tmp[subjectHits(olaps),]$genes
    hits2$exon_start <- start(grx_ex_tmp[subjectHits(olaps),])
    hits2$exon_end <- end(grx_ex_tmp[subjectHits(olaps),])
    
    hits2 <- inner_join(hits2, gene_df, by='genes')
    outside_gene <- hits2$pos1 < hits2$start & hits2$pos2 > hits2$end
    hits2 <- hits2[!outside_gene,]
    
    if(nrow(hits) + nrow(hits2) == 0) {return(NA)}
    
    # combine lists
    all_hits <- rbindlist(list(hits, hits2), use.names = T)
    all_hits$pcawg_driver <- all_hits$genes%in%pcawg_drivers$gene
    all_hits <- all_hits[,c('exon_start','exon_end'):=NULL]
    
    all_hits$id <- paste(all_hits$chr1, all_hits$pos1, all_hits$dir1, 
                         all_hits$chr2, all_hits$pos2, all_hits$dir2, all_hits$sample, sep=':')
    return(all_hits)
}

get_loh_hits <- function(samples, cn_dir, clin, pcawg_drivers) {
    all_loh <- NULL
    cn_files <- list.files(cn_dir)
    for (sample in samples) {
        sample <- as.character(sample)
        # print(sample)
        cnf <- cn_files[grep(sample, cn_files)]
        cn <- read.delim(paste(cn_dir, '/', cnf, sep=''), stringsAsFactors = F)
        sample <- strsplit(cnf, '\\.')[[1]][1]
        # id <- pcawg_list[pcawg_list$tumor_wgs_aliquot_id==sample,]$donor_unique_id
        # sex <- clin[clin$X..donor_unique_id==id,]$donor_sex
        sex <- clin[clin$sample%in%sample,]$inferred_sex
        cn <- cn[cn$minor_cn==0,]
        if (nrow(cn) > 0 & sex=='male') {
            xy_with_loss <- cn[(cn$chromosome=='X' | cn$chromosome=='Y') & cn$total_cn==0,]
            cn <- cn[!cn$chromosome%in%c('X','Y'),]
            cn <- rbind(xy_with_loss, cn)
        }
        if(nrow(cn) > 0){all_loh <- rbind(all_loh, data.frame(cn, sample=sample))}
    }
    all_loh <- all_loh[!is.na(all_loh$chromosome),]
    lohx <- GRanges(seqnames = all_loh$chromosome, ranges = IRanges(start = all_loh$start, end = all_loh$end))
    olaps <- suppressWarnings(findOverlaps(lohx, grx))
    loh_hits <- all_loh[queryHits(olaps),]
    loh_hits$genes <- grx[subjectHits(olaps),]$genes
    loh_hits$pcawg_driver <- loh_hits$genes%in%pcawg_drivers$gene
    
    return(unique(loh_hits))
}

get_cn_of_gene <- function(sample, cn_dir, gene, grx) {
    cn_files <- list.files(cn_dir)
    cnf <- cn_files[grep(sample, cn_files)]
    cn <- read.delim(paste(cn_dir, '/', cnf, sep=''), stringsAsFactors = F)
    if (nrow(cn) > 0) {
        cnx <- GRanges(seqnames = cn$chromosome, ranges = IRanges(start = cn$start, end = cn$end))
        olaps <- suppressWarnings(findOverlaps(cnx, grx))
        hits <- cn[queryHits(olaps),]
        hits$genes <- grx[subjectHits(olaps),]$genes
        return(unique(hits[hits$genes%in%gene,]))
    } else {return(NULL)}
}

get_snv_drivers <- function(driver_list, samples, snv_dir='', pcawg_pp='', get_CCF=FALSE) {
    sam_driv <- driver_list[driver_list$sample%in%as.character(samples),]
    sam_driv$id <- paste(sam_driv$chr, sam_driv$pos, sep=':')
    
    driver_snvs <- NULL
    if (get_CCF & snv_dir!='') {
        for(sample in unique(sam_driv$sample)) {
            pur <- pcawg_pp[pcawg_pp$samplename==sample,]$purity
            snv_cc <- paste(snv_dir, sample, '/ccube_out/snvs/', sample, '_cluster_certainty.txt', sep='')
            if(file.exists(snv_cc)){
                snvs <- read.delim(snv_cc)
                snvs$CCF <- snvs$average_proportion / pur
                if(nrow(snvs)==0){next}
                snvs$id <- paste(snvs$chr, snvs$pos, sep=':')
                dlist <- sam_driv[sam_driv$sample%in%sample,]
                snvs <- snvs[snvs$id%in%dlist$id,]
                if (nrow(snvs)>0) {
                    snvs$sample <- sample
                    driver_snvs <- rbind(driver_snvs, snvs)
                }
            }
            # dat_snv <- get_dat_snvs(sample, pur)
            # if (nrow(dat_snv)>0) {
            #     dat_snv$id <- paste(dat_snv$chr, dat_snv$pos, sep=':')
            #     dlist <- sam_driv[sam_driv$sample_id%in%sample,]
            #     snvs <- dat_snv[dat_snv$id%in%dlist$id,]
            #     if (nrow(snvs)>0) {
            #         snvs$sample <- sample
            #         driver_snvs <- rbind(driver_snvs, snvs)
            #     }
            # }
        }
        driver_snvs <- merge(driver_snvs, sam_driv, all.x=T)
    } else {driver_snvs <- sam_driv}
    
    return(driver_snvs)
}

get_coding_snv_hits <- function(samples, pcawg_pp, snvs_dir, pcawg_drivers) {
    snvs_list <- vector("list", length(samples))
    for(i in 1:length(samples)){
        sample <- samples[i]
        
        # dat <- get_dat_snvs(sample, pur)
        snv_cc <- paste(snv_dir, sample, '/ccube_out/snvs/', sample, '_cluster_certainty.txt', sep='')
        if(!file.exists(snv_cc)) {next}
        snvs <- read.delim(snv_cc, sep='\t')
        pur <- pcawg_pp[pcawg_pp$samplename==sample,]$purity
        snvs$CCF <- snvs$average_proportion / pur
        if(nrow(snvs)==0){next}
        
        snv_vcf_file <- paste(snvs_dir, sample, '.consensus.20160830.somatic.snv_mnv.vcf.gz', sep='')
        vcf <- read.delim(gzfile(snv_vcf_file), comment.char = '#', stringsAsFactors = F, header = F, sep='\t')
        vcf <- vcf[grep('Nonsense_Mutation|Missense_Mutation|Frameshift', vcf$V8),]
        if(nrow(vcf)==0){next}
        
        snvs_list[[i]] <- data.frame(merge(snvs, vcf, by=c(1,2)), sample=sample)
    }
    all_snvs <- NULL
    all_snvs <- rbind(all_snvs, do.call(rbind, snvs_list))
    
    snvx <- GRanges(seqnames = all_snvs$chr, ranges = IRanges(start = all_snvs$pos, end=all_snvs$pos))
    olaps <- suppressWarnings(findOverlaps(snvx, grx_ex))
    coding_hits <- all_snvs[queryHits(olaps),]
    coding_hits$gene <- grx_ex[subjectHits(olaps),]$genes
    
    var_clas <- as.character(sapply(coding_hits$V8, function(x){x<-strsplit(x,';')[[1]]; return(x[grep('Variant',x)])}))
    coding_hits$var_clas <- as.character(sapply(var_clas, function(x){strsplit(x,'=')[[1]][2]}))
    
    x <- coding_hits
    coding_hits <- data.frame(chrom=x$chr, pos=x$pos, CCF=x$CCF, 
                              ref=x$V4, alt=x$V5, variant=x$var_clas, gene=x$gene, sample=x$sample,
                              pcawg_driver=x$gene%in%pcawg_drivers$gene)
    
    return(coding_hits)
}

get_all_hits <- function(samples, drivers, drivers_only=T) {
    if(drivers_only) {
        snv_hits <- get_snv_drivers(driver_list, samples, pcawg_pp, snv_dir=snv_dir, get_CCF = T)
        snv_hits$pcawg_driver <- T
    } else {
        snv_hits <- get_coding_snv_hits(samples, pcawg_pp, snvs_dir, drivers)
    }
    sv_hits <- get_sv_hits(samples, pcawg_pp, svinfos, sv_dir, drivers, drivers_only)
    loh_hits <- get_loh_hits(samples, cn_dir, xclin, drivers)
    if(drivers_only){loh_hits <- loh_hits[loh_hits$pcawg_driver,]}
    
    clonality_sv <- sapply(sv_hits$CCF<cutoff, function(x){if(x){return('subclonal')}else{return('clonal')}})
    clonality_snv <- sapply(snv_hits$CCF<cutoff, function(x){if(x){return('subclonal')}else{return('clonal')}})
    all_hits <- data.frame(sample=sv_hits$sample, gene=sv_hits$genes, 
                           event=sv_hits$classification, clonality=clonality_sv,
                           driver=sv_hits$pcawg_driver)
    all_hits <- rbind(all_hits, data.frame(sample=snv_hits$sample, gene=snv_hits$gene, event='SNV',
                                           clonality=clonality_snv, driver=snv_hits$pcawg_driver))
    all_hits <- rbind(all_hits, data.frame(sample=loh_hits$sample, gene=loh_hits$genes, event='LoH', clonality='clonal',
                                           driver=loh_hits$pcawg_driver))
    all_hits <- distinct(all_hits)
    return(all_hits)
}

get_hit_proportions <- function(hz_events, cutoff=0.05) {
    hit_proportion <- table(as.character(unique(hz_events[,c('sample','gene')])$gene))/length(bal_enr_samples)
    rank <- order(hit_proportion, decreasing = T)
    hit_proportion <- hit_proportion[rank]
    return(hit_proportion[hit_proportion>cutoff])
}

get_homozygous_hits <- function(all_hits, subclonal_component=F) {
    events_pergene <- data.table(all_hits)[,length(event), by=c('gene','sample')]
    mult_events <- events_pergene[events_pergene$V1>1,]
    hz_events <- inner_join(mult_events, all_hits, by=c('sample','gene'))
    
    hz_events$id <- paste(hz_events$sample,hz_events$gene,sep=':')
    events_pergene <- data.table(hz_events)[,length(event), by=list(gene,sample,id)]
    mult_events <- events_pergene[events_pergene$V1>1,]$id
    hz_events <- hz_events[hz_events$id%in%mult_events | hz_events$event%in%'Loss',]
    
    hz_events <- unique(data.frame(hz_events))
    keep <- c()
    
    for (id in unique(hz_events$id)) {
        consider <- hz_events[hz_events$id==id,]
        events <- paste(consider$event, consider$clonality)
        if ('LoH clonal'%in%events & 'DEL clonal'%in%events) {
            consider <- consider[!events%in%'DEL clonal',]
        }
        if (nrow(consider)>1 | 'Loss clonal'%in%events) {
            if (subclonal_component) {
                if ('subclonal'%in%consider$clonality){keep <- c(keep, rownames(consider))}
            } else {keep <- c(keep, rownames(consider))}
        }
    }
    return(data.table(hz_events[keep,]))
}

get_homozygous_hits_quick <- function(hz_events) {
    # remove clonal deletions, assume they must be covered by LoH or loss
    hz_events <- unique(hz_events[!(hz_events$clonality%in%'clonal' & hz_events$event%in%'DEL'),])
    
    hz_events$id <- paste(hz_events$sample,hz_events$gene,sep=':')
    events_pergene <- data.table(hz_events)[,length(event), by=list(gene,sample,id)]
    mult_events <- events_pergene[events_pergene$V1>1,]$id
    hz_events <- hz_events[hz_events$id%in%mult_events | hz_events$event%in%'Loss',]
    
    hz_events <- unique(hz_events)
    return(hz_events)
}

post_assign <- function(topa, pur, clus_ccfs) {
    topa$purity <- pur
    topa$var_counts1 <- topa$support
    topa$var_counts2 <- topa$support
    topa$total_counts1 <- topa$support + topa$norm1_adjusted
    topa$total_counts2 <- topa$support + topa$norm2_adjusted
    topa$ccf_true <- NA
    
    cdev <- NULL
    cccf <- NULL
    valid_cn <- topa$total_cn1 > 0
    if(any(valid_cn)) {
        for(ccf in clus_ccfs) {
            topa[valid_cn,]$ccf_true <- ccf
            mults <- as.numeric(t(data.frame(apply(topa[valid_cn,], 1, get_true_cn_by_side, side=1))))
            pv1 <- as.numeric(get_pv(topa[valid_cn,]$purity, 1, topa[valid_cn,]$total_cn1, 
                                     topa[valid_cn,]$total_cn1, mults / topa[valid_cn,]$total_cn1))
            ccfs <- (1/pv1) * 
                topa[valid_cn,]$support / (topa[valid_cn,]$norm1_adjusted + topa[valid_cn,]$support)
            cdev <- cbind(cdev, (ccfs - ccf))
            cccf <- cbind(cccf, ccfs)
        }
        if(ncol(cccf) == 1) {
            best_ccf <- cccf[,1]
        } else {
            select_ccf <- apply(cdev, 1, function(x){which(min(abs(x))==abs(x))})
            best_ccf <- NULL
            for(i in 1:nrow(cccf)) {
                select <- select_ccf[i]
                best_ccf <- c(best_ccf, cccf[i, select])
            }
        }
        topa[valid_cn,]$average_proportion1 <- as.numeric(best_ccf) * pur
    }
    
    cdev <- NULL
    cccf <- NULL
    valid_cn <- topa$total_cn2 > 0
    if(any(valid_cn)) {
        for(ccf in clus_ccfs) {
            topa[valid_cn,]$ccf_true <- ccf
            mults <- as.numeric(t(data.frame(apply(topa[valid_cn,], 1, get_true_cn_by_side, side=2))))
            pv2 <- as.numeric(get_pv(topa[valid_cn,]$purity, 1, topa[valid_cn,]$total_cn2, 
                                     topa[valid_cn,]$total_cn2, mults / topa[valid_cn,]$total_cn2))
            ccfs <- (1/pv2) * 
                topa[valid_cn,]$support / (topa[valid_cn,]$norm2_adjusted + topa[valid_cn,]$support)
            cdev <- cbind(cdev, (ccfs - ccf))
            cccf <- cbind(cccf, ccfs)
        }
        if(ncol(cccf) == 1) {
            best_ccf <- cccf[,1]
        } else {
            select_ccf <- apply(cdev, 1, function(x){which(min(abs(x))==abs(x))})
            best_ccf <- NULL
            for(i in 1:nrow(cccf)) {
                select <- select_ccf[i]
                best_ccf <- c(best_ccf, cccf[i, select])
            }
        }
        topa[valid_cn,]$average_proportion2 <- as.numeric(best_ccf) * pur
    }
    return(list(topa$average_proportion1, topa$average_proportion2))
}

is_gtype_valid <- function(gtype) {
    if(is.na(gtype)) {
        return(FALSE)
    }
    cns <- strsplit(gtype, '\\|')[[1]]; 
    cn1 <- as.numeric(strsplit(cns[1], ',')[[1]][1:2]);
    cn2 <- as.numeric(strsplit(cns[2], ',')[[1]][1:2]);
    positive_cn <- sum(cn1) + sum(cn2, na.rm=T) > 0
    return(as.logical(positive_cn))
}

post_assign_subclonal <- function(topa, pur, sc) {
    clus_ccfs <- sc$proportion / pur
    clus_labs <- sc$cluster
    
    topa$purity <- pur
    topa$var_counts1 <- topa$support
    topa$var_counts2 <- topa$support
    topa$total_counts1 <- topa$support + topa$norm1_adjusted
    topa$total_counts2 <- topa$support + topa$norm2_adjusted
    topa$true_ccf <- NA
    topa$adjusted_support <- topa$support
    topa$mindev1 <- NA; topa$mindev2 <- NA
    topa$ass1 <- NA; topa$ass2 <- NA
    cdev <- NULL
    cccf <- NULL
    valid_cn1 <- as.logical(sapply(topa$gtype1, is_gtype_valid))
    if(any(valid_cn1)) {
        for(ccf in clus_ccfs) {
            topa[valid_cn1,]$true_ccf <- ccf
            ccfs <- as.numeric(t(data.frame(apply(topa[valid_cn1,], 1, 
                                                  get_true_sc_cn_by_side, pur=pur, side=1)))[,4])
            cdev <- cbind(cdev, (ccfs - ccf))
            cccf <- cbind(cccf, ccfs)
        }
        if(ncol(cccf) == 1) {
            best_ccf <- cccf[,1]
            select_ccf <- rep(1, nrow(cccf))
        } else {
            select_ccf <- apply(cdev, 1, function(x){which(min(abs(x))==abs(x))[1]})
            best_ccf <- NULL
            for(i in 1:nrow(cccf)) {
                select <- select_ccf[i]
                best_ccf <- c(best_ccf, cccf[i, select])
            }
        }
        topa[valid_cn1,]$average_proportion1 <- as.numeric(best_ccf) * pur
        topa[valid_cn1,]$mindev1 <- apply(cdev, 1, function(x){min(abs(x))}) 
        # ^ smallest deviation from clus CCF
        topa[valid_cn1,]$ass1 <- clus_labs[select_ccf] # assignments 1
    }
    
    cdev <- NULL
    cccf <- NULL
    valid_cn2 <- as.logical(sapply(topa$gtype2, is_gtype_valid))
    if(any(valid_cn2)) {
        for(ccf in clus_ccfs) {
            topa[valid_cn2,]$true_ccf <- ccf
            ccfs <- as.numeric(t(data.frame(apply(topa[valid_cn2,], 1, 
                                                  get_true_sc_cn_by_side, pur=pur, side=2)))[,4])
            cdev <- cbind(cdev, (ccfs - ccf))
            cccf <- cbind(cccf, ccfs)
        }
        if(ncol(cccf) == 1) {
            best_ccf <- cccf[,1]
            select_ccf <- rep(1, nrow(cccf))
        } else {
            select_ccf <- apply(cdev, 1, function(x){which(min(abs(x))==abs(x))[1]})
            best_ccf <- NULL
            for(i in 1:nrow(cccf)) {
                select <- select_ccf[i]
                best_ccf <- c(best_ccf, cccf[i, select])
            }
        }
        topa[valid_cn2,]$average_proportion2 <- as.numeric(best_ccf) * pur
        topa[valid_cn2,]$mindev2 <- apply(cdev, 1, function(x){min(abs(x))}) # smallest deviation from clus CCF
        topa[valid_cn2,]$ass2 <- clus_labs[select_ccf] # assignments 1
    }
    
    # bulk-assign (handles any with NA on one side)   
    topa[valid_cn1,]$cluster <- topa[valid_cn1,]$ass1
    topa[valid_cn2,]$cluster <- topa[valid_cn2,]$ass2
    
    # now for these we have to look at the mimimum deviation
    tmp <- topa[valid_cn1 & valid_cn2, c('mindev1', 'mindev2')]
    select <- as.numeric(apply(tmp, 1, function(x){which(min(abs(x),na.rm=T)==abs(x))[1]}))
    tmp <- topa[valid_cn1 & valid_cn2, c('ass1', 'ass2')]
    
    clus <- NULL
    for(i in 1:nrow(tmp)) {clus <- c(clus, tmp[i, select[i]])}
    topa[valid_cn1 & valid_cn2, 'cluster'] <- clus
    
    return(list(topa$average_proportion1, topa$average_proportion2, topa$cluster))
}

post_assign_all <- function(svs, sc, cna_file, pur, pl, all=F) {
    svs <- svs[svs$support>0, ]
    svs$norm1_adjusted <- svs$norm1
    svs$norm2_adjusted <- svs$norm2
    gains <- svs$classification%in%c('DUP','INTDUP')
    af <- 1 - (pur / pl)
    svs$norm1_adjusted[gains] <- svs$norm1_adjusted[gains] * af
    svs$norm2_adjusted[gains] <- svs$norm2_adjusted[gains] * af
    
    cna <- read.delim(cna_file, sep='\t')
    # svs <- assign_cnas(svs, cna)
    svs <- assign_cnas_subclonal(svs, cna)
    # svs <- svs[which((svs$total_cn1 + svs$total_cn2) > 0),]
    
    topa <- svs[is.na(svs$average_proportion1),]
    if(all){topa <- svs}
    if(nrow(topa)>0) {
        # props <- post_assign(topa, pur, clus_ccfs)
        props <- post_assign_subclonal(topa, pur, sc)
        if(all) {
            svs$average_proportion1 <- props[[1]]
            svs$average_proportion2 <- props[[2]]    
            svs$cluster <- props[[3]]   
        } else {
            svs[is.na(svs$average_proportion1),'average_proportion1'] <- props[[1]]
            svs[is.na(svs$average_proportion2),'average_proportion2'] <- props[[2]]    
            svs[is.na(svs$cluster),'cluster'] <- props[[3]]        
        }
    }
    return(svs)
}
get_svs_new <- function(sv_dir, svinfos, aid) {
    sv_ccert <- paste(sv_dir, aid, '/ccube_out/', aid, '_cluster_certainty.txt', sep='')
    sv_ss <- paste(sv_dir, aid, '/ccube_out/', aid, '_subclonal_structure.txt', sep='')
    sv_filt <- paste(sv_filt_dir, aid, '/', aid, '_filtered_svs.tsv', sep='')
    sv_info <- paste(svinfos, aid, '_svinfo.txt', sep='')
    cna_file <- paste('~/data/full_cna/', 
                      aid, '.full.consensus.20170119.somatic.cna.annotated.txt', sep='')
    pur <- pcawg_pp[pcawg_pp$samplename%in%aid,'purity']
    pl <- pcawg_pp[pcawg_pp$samplename%in%aid,'ploidy']
    svs <- NULL
    sc <- NULL
    if (file.exists(sv_ccert)) {
        svs <- read.delim(sv_ccert, sep='\t')
        svf <- read.delim(sv_filt, sep='\t')
        svi <- read.delim(sv_info, sep='\t')
        sc <- read.delim(sv_ss, sep='\t')
        
        svs$cluster <- svs$most_likely_assignment
        svs <- inner_join(svs, sc, by='cluster')
        
        ########################################################################
        # post assign
        #########################################################################
        svs$chr1 <- as.character(svs$chr1); svs$chr2 <- as.character(svs$chr2)
        svi$chr1 <- as.character(svi$chr1); svi$chr2 <- as.character(svi$chr2)
        
        svf$chr1 <- as.character(svf$chr1); svf$chr2 <- as.character(svf$chr2)
        svs <- left_join(svi, svs, by=c('chr1', 'pos1', 'dir1', 'chr2', 'pos2', 'dir2'))
        svf_cols <- c('chr1', 'pos1', 'dir1', 'chr2', 'pos2', 'dir2', 
                      colnames(svf[!colnames(svf)%in%colnames(svi)]))
        svs <- left_join(svs, svf[,svf_cols], by=c('chr1', 'pos1', 'dir1', 'chr2', 'pos2', 'dir2'))
        
        svs <- post_assign_all(svs, sc, cna_file, pur, pl)
        svs$post_assign <- is.na(svs$proportion)
        svs$proportion <- sc$proportion[svs$cluster]
        post_assign_frac <- sum(svs$post_assign) / nrow(svs)
        
        #########################################################################
        
        ccfs <- svs[,c('average_proportion1', 'average_proportion2')]
        ccfs <- ccfs / pur
        svs$ccf_max <- as.numeric(apply(ccfs, 1, max, na.rm=T))
        svs$ccf_mean <- as.numeric(apply(ccfs, 1, mean, na.rm=T))
        
        svs$true_ccf <- svs$proportion / pur
        valid_cn <- as.logical(sapply(svs$gtype1, is_gtype_valid))
        svs$adjusted_support <- svs$support
        
        tmp <- t(data.frame(apply(svs[valid_cn,], 1, get_true_sc_cn_by_side, side=1, pur=pur)))
        mult1 <- as.numeric(tmp[,1]); total_cn1 <- as.numeric(tmp[,2])
        svs$mult1 <- NA; svs[valid_cn, 'mult1'] <- mult1
        svs$total_cn1 <- NA; svs[valid_cn, 'total_cn1'] <- total_cn1
        
        valid_cn <- as.logical(sapply(svs$gtype2, is_gtype_valid))
        tmp <- t(data.frame(apply(svs[valid_cn,], 1, get_true_sc_cn_by_side, side=2, pur=pur)))
        mult2 <- as.numeric(tmp[,1]); total_cn2 <- as.numeric(tmp[,2])
        svs$mult2 <- NA; svs[valid_cn, 'mult2'] <- mult2
        svs$total_cn2 <- NA; svs[valid_cn, 'total_cn2'] <- total_cn2
    }
    return(list(svs, sc))
}

get_snvs_new <- function(snv_dir, snvs_dir, aid) {
    snv_mult <- paste(snv_dir, aid, '/', aid, '_multiplicity.txt.gz', sep='')
    snv_ss <- paste(snv_dir, aid, '/ccube_out/snvs/', aid, '_subclonal_structure.txt', sep='')
    snv_cc <- paste(snv_dir, aid, '/ccube_out/snvs/', aid, '_cluster_certainty.txt', sep='')
    snv_mu <- paste(snv_dir, aid, '/ccube_out/snvs/', aid, '_multiplicity.txt', sep='')
    pc_snvs <- paste0(snvs_dir, aid, '.consensus.20160830.somatic.snv_mnv.vcf.gz')
    
    if (file.exists(snv_cc) & file.exists(pc_snvs)) {
        pur <- pcawg_pp[pcawg_pp$samplename%in%aid,'purity']
        
        snvs <- read.delim(snv_cc, sep='\t')
        sc <- read.delim(snv_ss, sep='\t')
        sc$CCF <- sc$proportion / pur
        mult <- read.delim(snv_mu, sep='\t')
        pcs <- read.delim(gzfile(pc_snvs), comment='#', sep='\t', header=F)
        colnames(pcs)[1:2] <- c('chr', 'pos')
        pcs$vaf <- get_vafs(pcs)
        
        snvs$cluster <- snvs$most_likely_assignment
        snvs <- inner_join(snvs, sc, by='cluster')
        snvs <- inner_join(snvs, mult, by=c('chr', 'pos'))
        snvs$pv <- get_pv(pur, 1,
                          snvs$tumour_copynumber,
                          snvs$tumour_copynumber,
                          snvs$multiplicity/snvs$tumour_copynumber, 1)
        snvs$chr <- as.character(snvs$chr)
        snvs <- inner_join(snvs, pcs[,c('chr', 'pos', 'vaf')], by=c('chr', 'pos'))
        snvs$CCF <- (1/snvs$pv) * snvs$vaf
        return(list(snvs,sc))
    }
    return(NULL)
}

get_vafs <- function(vcf) {
    vafs <- sapply(as.character(vcf$V8), function(x) {
        ss <- strsplit(x, ';')[[1]]
        alt <- ss[grep('t_alt', ss)]
        ref <- ss[grep('t_ref', ss)]
        if(length(alt) > 0 & length(ref) > 0)  {
            alt <- as.numeric(strsplit(as.character(alt[1]), '=')[[1]][2])
            ref <- as.numeric(strsplit(as.character(ref[1]), '=')[[1]][2])
            return(alt / (ref + alt))
        } else {return(NA)}
    })
    return(as.numeric(unlist(vafs)))
}