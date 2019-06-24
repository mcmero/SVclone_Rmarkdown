get_sv_hits <- function(samples, pcawg_pp, svinfos, sv_dir, pcawg_driver_genes, drivers_only=T) {
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
        grx_tmp <- grx[grx$genes%in%pcawg_driver_genes$gene]
        grx_ex_tmp <- grx_ex[grx_ex$genes%in%pcawg_driver_genes$gene]
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
    all_hits$pcawg_driver <- all_hits$genes%in%pcawg_driver_genes$gene
    all_hits <- all_hits[,c('exon_start','exon_end'):=NULL]
    
    all_hits$id <- paste(all_hits$chr1, all_hits$pos1, all_hits$dir1, 
                         all_hits$chr2, all_hits$pos2, all_hits$dir2, all_hits$sample, sep=':')
    return(all_hits)
}

get_loh_hits <- function(samples, consensus_cna_dir, clin, pcawg_driver_genes) {
    all_loh <- NULL
    cn_files <- list.files(consensus_cna_dir)
    for (sample in samples) {
        sample <- as.character(sample)
        cnf <- cn_files[grep(sample, cn_files)]
        cn <- read.delim(paste(consensus_cna_dir, '/', cnf, sep=''), stringsAsFactors = F)
        sample <- strsplit(cnf, '\\.')[[1]][1]
        sex <- xclin[xclin$sample%in%sample,]$inferred_sex
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
    loh_hits$pcawg_driver <- loh_hits$genes%in%pcawg_driver_genes$gene
    
    return(unique(loh_hits))
}

get_snv_drivers <- function(pcawg_driver_snvs, samples, snv_dir='', pcawg_pp='', get_CCF=FALSE) {
    sam_driv <- pcawg_driver_snvs[pcawg_driver_snvs$sample%in%as.character(samples),]
    sam_driv$id <- paste(sam_driv$chr, sam_driv$pos, sep=':')
    
    driver_snvs <- NULL
    if (get_CCF & snv_dir!='' & pcawg_pp!='') {
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
        }
        driver_snvs <- merge(driver_snvs, sam_driv, all.x=T)
    } else {driver_snvs <- sam_driv}
    
    return(driver_snvs)
}

get_coding_snv_hits <- function(samples, pcawg_pp, consensus_snvs_dir, pcawg_driver_genes) {
    snvs_list <- vector("list", length(samples))
    for(i in 1:length(samples)){
        sample <- samples[i]
        
        snv_cc <- paste(snv_dir, sample, '/ccube_out/snvs/', sample, '_cluster_certainty.txt', sep='')
        if(!file.exists(snv_cc)) {next}
        snvs <- read.delim(snv_cc, sep='\t')
        pur <- pcawg_pp[pcawg_pp$samplename==sample,]$purity
        snvs$CCF <- snvs$average_proportion / pur
        if(nrow(snvs)==0){next}
        
        snv_vcf_file <- paste(consensus_snvs_dir, sample, '.consensus.20160830.somatic.snv_mnv.vcf.gz', sep='')
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
                              pcawg_driver=x$gene%in%pcawg_driver_genes$gene)
    
    return(coding_hits)
}

get_all_hits <- function(samples, pcawg_driver_genes,
                         pcawg_driver_snvs, pcawg_pp, 
                         sv_dir, snv_dir, xclin,
                         consensus_snvs_dir, consensus_cna_dir,
                         drivers_only=T, cutoff=0.7) {
    if(drivers_only) {
        snv_hits <- get_snv_drivers(pcawg_driver_snvs, samples, pcawg_pp, snv_dir=snv_dir, get_CCF = T)
        snv_hits$pcawg_driver <- T
    } else {
        snv_hits <- get_coding_snv_hits(samples, pcawg_pp, consensus_snvs_dir, pcawg_driver_genes)
    }
    sv_hits <- get_sv_hits(samples, pcawg_pp, svinfos, sv_dir, pcawg_driver_genes, drivers_only)
    loh_hits <- get_loh_hits(samples, consensus_cna_dir, xclin, pcawg_driver_genes)
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

calculate_driver_enrichment <- function(driver_hits, all_driver_hits, n, n_total, perc_cutoff=0.05, qval=0.05) {
    x <- NULL
    genes <- names(table(driver_hits$gene)[table(driver_hits$gene)/n>perc_cutoff])
    for(gene in genes) {
        enr_total <- length(unique(all_driver_hits[all_driver_hits$gene==gene,]$sample))
        n_enr_in_sample <- length(unique(driver_hits[driver_hits$gene==gene,]$sample))
        pval <- phyper(n_enr_in_sample, enr_total, 
                       (n_total - enr_total), n, lower.tail=FALSE)
        x <- rbind(x, data.frame(gene=gene, pval=pval))
    }
    x$qval <- p.adjust(x$pval, method='BH')
    x <- x[order(x$pval),]
    return(x[x$qval < qval,])
}