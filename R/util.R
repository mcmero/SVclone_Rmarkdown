get_ccube_f <- function(t, phi, cn, mult, epi=1e-3) {
    t <- as.numeric(t); phi <- as.numeric(phi)
    cn <- as.numeric(cn); mult <- as.numeric(mult)
    
    w = t * (mult * (1 - epi) - cn * epi)
    w = w / ((1 - t) * cn + t * cn)
    f = w * phi + epi
    return(f)
}

get_true_sc_cn_by_side <- function(row, pur, side=1) {
    if(side==1){
        gtype <- row['gtype1']
        norm <- row['norm1_adjusted']
    } else {
        gtype <- row['gtype2']
        norm <- row['norm2_adjusted']
    }
    if (gtype == '') {
        return(c(NA,NA,NA,NA))
    }
    
    frac <- 1
    gtype <- strsplit(as.character(gtype), '\\|')[[1]]
    ccf <- as.numeric(row['true_ccf'])
    s <- as.numeric(row['adjusted_support'])
    n <- round(as.numeric(norm) + s)
    norm_cn <- as.numeric(row['norm_cn'])
    
    if (length(gtype)>1) {
        cn1 <- as.numeric(strsplit(gtype[1], ',')[[1]])
        cn2 <- as.numeric(strsplit(gtype[2], ',')[[1]])
        frac <- cn1[3]
        pvs <- c()
        
        total_cn1 = cn1[1] + cn1[2]
        total_cn2 = cn2[1] + cn2[2]
        total_cn = (total_cn1 * frac) + (total_cn2 * (1 - frac))
        
        pv <- as.numeric(get_ccube_f(pur, ccf, total_cn, total_cn))
        m <- min((s/n)/pv*total_cn, total_cn)
    } else {
        cn <- as.numeric(strsplit(gtype[1], ',')[[1]])
        pvs <- c()
        total_cn = cn[1] + cn[2]
        
        ms = 1:cn[1]
        for (i in 1:length(ms)) {
            pvs <- c(pvs, as.numeric(get_ccube_f(pur, ccf, total_cn, ms[i])))
        }
        probs <- pbinom(s, n, pvs) - pbinom(s-1, n, pvs)
        m = ms[which(probs==max(probs))[1]]
        pv <- pvs[which(probs==max(probs))[1]]
    }
    
    ccf <- MapVaf2CcfPyClone(s/n,
                             pur,
                             norm_cn,
                             total_cn,
                             total_cn,
                             m,
                             constraint=F)
    ccf <- min(2, as.numeric(ccf))
    return(c(m, total_cn, pv, as.numeric(ccf)))
}

get_true_sc_cn <- function(row) {
    frac <- 1
    pur <- as.numeric(row['purity'])
    gtype <- strsplit(as.character(row['gtype']), '\\|')[[1]]
    ccf <- as.numeric(row['true_ccf'])
    s <- as.numeric(row['adjusted_support'])
    n <- round(as.numeric(row['adjusted_depth']))
    norm_cn <- as.numeric(row['norm_cn'])
    
    if (length(gtype)>1) {
        cn1 <- as.numeric(strsplit(gtype[1], ',')[[1]])
        cn2 <- as.numeric(strsplit(gtype[2], ',')[[1]])
        frac <- cn1[3]
        pvs <- c()
        
        total_cn1 = cn1[1] + cn1[2]
        total_cn2 = cn2[1] + cn2[2]
        total_cn = (total_cn1 * frac) + (total_cn2 * (1 - frac))
        
        major_cn = total_cn
        pv <- as.numeric(get_ccube_f(pur, ccf, total_cn, total_cn))
        m <- min((s/n)/pv*total_cn, total_cn)
    } else {
        cn <- as.numeric(strsplit(gtype[1], ',')[[1]])
        pvs <- c()
        
        total_cn = cn[1] + cn[2]
        ms = 1:cn[1]
        for (i in 1:length(ms)) {
            pvs <- c(pvs, as.numeric(get_ccube_f(pur, ccf, total_cn, ms[i])))
        }
        probs <- pbinom(s, n, pvs) - pbinom(s-1, n, pvs)
        m = ms[which(probs==max(probs))[1]]
        pv <- pvs[which(probs==max(probs))[1]]
    }
    
    ccf <- MapVaf2CcfPyClone(s/n,
                             pur,
                             norm_cn,
                             total_cn,
                             total_cn,
                             m,
                             constraint=F)
    ccf <- min(2, as.numeric(ccf))
    return(c(m, total_cn, pv, ccf))
}

# matlab style helper
repmat = function(X,m,n){
    ##R equivalent of repmat (matlab)
    mx = dim(X)[1]
    nx = dim(X)[2]
    matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}

GetMultFromCcf <- function(bn, dn, ccf, major_cn, minor_cn, purity, epi = 1e-3) {
    total_cn = major_cn + minor_cn
    z = (1-purity)*2 + purity*total_cn
    k = max(major_cn)
    multPool <- seq(0, k, by = 1)
    aa = repmat(as.matrix(purity * ccf * (1 - epi) / z), 1, k+1)
    bb = repmat(as.matrix(epi * ( z - purity * (1-ccf) * total_cn)), 1, k+1)
    pp = bsxfun.se("*", aa, multPool) + bb
    ll = bsxfun.se("*", log(pp), bn) + bsxfun.se("*", log(1-pp), dn-bn)
    return(multPool[apply(ll, 1, which.max)])
}