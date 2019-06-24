draw_circos_clonal <- function(dat, sample, pur, bbf) {
    bb <- read.delim(bbf, sep='\t', stringsAsFactors = F)
    
    pdat <- data.frame(chr=bb$chromosome, startpos=bb$start, end=bb$end, value=bb$total_cn)
    pdat$chr <- paste('chr', pdat$chr,sep='')
    colnames(pdat) <- c('chr','start','end','value')
    pdat$value[pdat$value > 6] <- 6
    
    colours <- c('#0000FF80','#FF000080','darkgreen','#0000FF40','#FF000040','#00FF0040')
    sample <- strsplit(sample, '-')[[1]][1]
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