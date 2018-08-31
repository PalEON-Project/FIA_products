## setup for composition, following code in 'composition' and 'composition_fia'
## repositories

library(dplyr)
library(tidyr)

counts <- fia %>%
    rename(statecd = STATECD, plt_cn = PLT_CN) %>% 
    filter(!is.na(level3s)) %>% 
    group_by(statecd, plt_cn, level3s, x, y, cell) %>%
    summarize(count = n()) %>%
    spread(level3s, count, fill = 0) %>% ungroup()

## note all NAs in level3s are Douglas fir; omitting these

list_of_states <- list(eastern = paleon_states_east,
                       western = paleon_states_west_ohio)

for(region in c('eastern', 'western')) {

    nonTaxaCols <- c('statecd', 'plt_cn', 'x', 'y', 'cell')
    tmp <- counts %>% filter(statecd %in% list_of_states[[region]])
    plotInfo <- tmp[ , nonTaxaCols]
    data <- tmp[ , !names(tmp) %in% nonTaxaCols]

    nTaxa <- ncol(data)
    taxaUse <- names(data)
    taxaUse <- gsub("/", ",", taxaUse)
    taxa <- data.frame(taxonID = 1:nTaxa, taxonName = taxaUse, stringsAsFactors = FALSE)

    domainX <- eval(as.name(paste0(region, 'DomainX')))
    domainY <- eval(as.name(paste0(region, 'DomainY')))
    m1 <- length(domainX)
    m2 <- length(domainY)
    nCells <- m1 * m2

    ## add ID that just indexes the Paleon grid cells in rectangular
    ## subdomain being fit in statistical model
    ## number cells from NW corner, travelling east, then by row southwards
    coord <- expand.grid(x = xGrid[domainX], y = rev(yGrid)[domainY])
    coord$ID <- seq_len(nrow(coord))

    tmp <- plotInfo$cell
    plotInfo <- plotInfo %>% left_join(coord, c('x' = 'x', 'y' = 'y')) %>%
        dplyr::select(plt_cn, x, y, ID, cell)
    if(!identical(tmp, plotInfo$cell))
        stop("plotInfo row order not the same as data row order")
    data <- data[order(plotInfo$ID), ]
    plotInfo <- plotInfo[order(plotInfo$ID), ]

    type <- nbhdStructure
    substring(type, 1 ,1) = toupper(substring(type, 1, 1))
    fns <- rep("", 2)
    fns[1] <- paste('graph', type, '-',  m1, 'x', m2, '.csv', sep='')
    fns[2] <- paste('graphCats', type, '-', m1, 'x', m2, '.csv', sep='')

    if(!file.exists(file.path(interim_results_dir, fns[1])) || (nbhdStructure != 'bin' && !file.exists(file.path(interim_results_dir, fns[2])))) {
        fns <- graphCreate(m1, m2, type = nbhdStructure, dir = interim_results_dir, fn = fns[1], fn_cats = fns[2])
    } 

    nbhd <- graphRead(fns[1], fns[2], m1, m2, type = nbhdStructure, dir = interim_results_dir)

    if(!nbhdStructure %in% c('bin', 'tps')) {  
                                        # remove boundary stuff for now while wait to hear from Finn about boundary correction
        nbhd@entries[nbhd@entries %in% c(-4, -6)] <- -8
        nbhd@entries[nbhd@entries %in% c(4, 10, 11, 18, 19)] <- 20
    }
    if(nbhdStructure == "lindgren_nu1" || nbhdStructure == "tps") {
        nbhdIndices <- list()
                                        # determine which elements correspond to what types of neighbors for fast filling in MCMC
        nbhdIndices$self <- which(nbhd@entries == 20)
        nbhdIndices$cardNbr <- which(nbhd@entries == -8)
        nbhdIndices$otherNbr <- which(nbhd@entries %in% c(1,2))
    } else nbhdIndices <- NULL


    ## create data objects for MCMC fitting 

    data <- as.matrix(data)

    nPlots <- nrow(plotInfo)

    plotInfo$n <- rowSums(data)
    plotInfo$plot <- 1:nrow(plotInfo)

    tmp <- plotInfo %>% group_by(ID) %>% summarize(n = sum(n))

    coord <- coord %>%
        left_join(tmp, c('ID' = 'ID')) %>%
        dplyr::select(ID, x, y, n)
    coord[is.na(coord)] <- 0

    taxon <- rep(rep(1:nTaxa, nPlots), times = c(t(data)))
    total <- rowSums(data)
    plot <- rep(which(total > 0), times = total[total > 0])
    cell <- plotInfo$ID[plot]

    dataTabular <- data
    data <- data.frame(taxon = taxon, cell = cell, plot = plot)

    outpath <- file.path(interim_results_dir, paste0('composition_raw_', region, '.Rda'))
    cat("Writing data file: ", outpath, "\n")
    save(nbhd, nbhdIndices, m1, m2, nTaxa, nCells, data, coord, taxa, plotInfo, dataTabular, file = outpath)
}

