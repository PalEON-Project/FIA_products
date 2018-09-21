## copied from composition_fia (which includes plot effect, unlike PLS composition code)

require(spam)

runMCMC <-function(y, cell = NULL, plot = NULL, plotInfo = NULL, gridInfo = NULL, C, Cindices = NULL, town = NULL, townCellOverlap = NULL, townCellIds = NULL,
                   S, thin, resumeRun, secondRun = FALSE, hyperpar = c(-0.5, 0),
                   adaptInterval = 100, adaptStartDecay = 3000, nbhdStructure = 'bin',                  
                   outputNcdfName, taxa, numCores = 1, runID = "",
                   dataDir, outputDir) {
    
    if(!identical(hyperpar, c(-0.5, 0)))
        stop("Error (runMCMC): joint sampling not set up to use any prior other than flat on sd scale.")
    
    exclude <- is.na(y)
    exclude <- exclude | is.na(cell) 
    
    y    <- y[!exclude]
    nTrees <-length(y)
    
    if(sum(exclude))
        warning("NAs found in 'y', 'cell', or 'town'; excluding these observations.")
    
    cell <- cell[!exclude]
    plot <- plot[!exclude]
    
    
    P <- length(taxa$taxonName)
    nCells <- nrow(C)
    nPlots <- nrow(plotInfo)
    
    if(!nbhdStructure %in% c('bin', 'tps')) {
        etaBounds <- c(log(.1), 5)
        muBounds <- c(-10, 10)
        mu_propSD <- rep(0.10, P)
        numAcceptMu <- rep(0, P)
        joint_propSD <- rep(1, P)
        adaptVals <- adaptedCov <- adaptedL <- list()
        adaptScale <- rep(.01, P)
        length(adaptScale) <- length(adaptVals) <- length(adaptedCov) <- length(adaptedL) <- P
        for(p in seq_len(P)) {
            adaptedCov[[p]] <- diag(rep(1, 2))
            adaptedL[[p]] <- t(chol(adaptedCov[[p]]))
            adaptVals[[p]] <- matrix(0, adaptInterval, 2)
        }
        numAcceptSigmaEta <- rep(0, P)
    } else {
        logSigma2_propSD <- rep(0.02, P)
        logTau2_propSD <- rep(0.02, P)
        numAcceptSigma2 <- rep(0, P)
        numAcceptTau2 <- rep(0, P)
    }
    
    infs <- rep(Inf, nTrees)
    negInfs <- rep(-Inf, nTrees)
    
    S_keep <- floor(S/thin) 
    S      <- S - S%%S_keep
    
    if(resumeRun || secondRun) {
        load(file.path(dataDir, paste0('lastState_', runID, '.Rda')))
        .Random.seed <<- .Random.seed
        if(secondRun) s <- 0
        sampleIterates <- (s+1):S

        alpha_current <- alpha_next
        eta_current <- eta_next
        h_current <- h_next
        mu_current <- mu_next
        sigma2_current <- sigma2_next
        tau2_current <- tau2_next
        if(secondRun) storeIndex <- 1
    } else {
        sampleIterates <- 1:S
        storeIndex <- 1
        
        W <- matrix(rnorm(nTrees*P), nTrees, P)
        
        
        sigma2store <- matrix(0, S_keep, P)
        tau2store <- matrix(0, S_keep, P)
        if(!nbhdStructure %in% c('bin', 'tps'))
            muStore <- etaStore <- sigma2store
        ##alpha_current <- alpha_next <- matrix(0, I, P)
        alpha_current <- alpha_next <- matrix(rnorm(nCells*P), nCells, P)
        
        sigma2_current <- sigma2_next <- runif(P, 0.1, 2)
        tau2_current <- tau2_next <- runif(P, 0.1, 2)
        cnt <- table(y)
        ## have less common taxa start with larger sigma2
        sigma2_current[cnt < 2000] <- sigma2_next[cnt < 2000] <- runif(sum(cnt < 2000), 0.5, 2)
        
        # within-cell heterogeneity (across FIA plots)
        h_current <- h_next <- matrix(rnorm(nPlots*P, 0, rep(sqrt(tau2_current), each = nPlots)), nPlots, P)
        
        if(!nbhdStructure %in% c('bin', 'tps')) {
            mu_current <- mu_next <- runif(P, -2, 2)
            eta_current <- eta_next <- runif(P, 0, 5) # rep(log(3), P)
        }
    }
    if(secondRun) runID <- paste0(runID, 'b') else runID <- paste0(runID, 'a')
    
    nTreesPerPlot <- plotInfo$n
    
    nTreesPerCell <- rep(0, nCells)
    tbl <- table(cell)
    nTreesPerCell[as.numeric(names(tbl))] <- tbl
    
    Wi.pbar_cell_centered <- compute_cell_sums_cpp(W,cell,nCells,P)/c(nTreesPerCell+1*(nTreesPerCell==0)) - compute_cell_sums_cpp(h_current*nTreesPerPlot,plotInfo$ID,nCells,P)/c(nTreesPerCell+1*(nTreesPerCell==0))
   
    U <- list()
    dfAdj <- 0
    if(nbhdStructure == "bin") dfAdj <- 1 # intercept
    if(nbhdStructure == "tps") dfAdj <- 3 # intercept and linear in cardinal directions
  
    if(!nbhdStructure %in% c('bin', 'tps')) {
        Uc <- UcProp <- list()
        Vinv <- list()
        for(p in seq_len(P)) {
            Vinv[[p]] <- C  # start with TPS
            a <- 4 + 1/exp(2*eta_current[p])  # a = 4 + kappa^2 = 4 + (1/rho)^2 = 4 + (1/exp(eta)^2)
            Vinv[[p]]@entries[Cindices$self] <- 4 + a*a
            Vinv[[p]]@entries[Cindices$cardNbr] <- -2 * a
            Vinv[[p]] <- Vinv[[p]] / (sigma2_current[p]*4*pi/exp(2*eta_current[p]))
                                        # do initial Cholesky based on sparseness pattern, with first taxon
            B <- Vinv[[p]]
            diag(B) <- diag(B) + nTreesPerCell
            if(!is.spam(B)) warning("B is not spam")
            U[[p]] <- chol(B)  
            Uc[[p]] <- chol(Vinv[[p]])
        }
        UcProp <- chol(B)
    } else {
        for(p in seq_len(P)) {
            B <- Vinv <- C / sigma2_current[p]
            diag(B) <- diag(B) + nTreesPerCell
            if(!is.spam(B)) warning("B is not spam")
            U[[p]] <- chol(B)  
        }
    }
    
    Uprop <- chol(B)


  # sample alphas given Ws and hs and sigmas (and mus/etas if relevant)
    if(!resumeRun && !secondRun) {
        for(p in 1:P){
            if(!nbhdStructure %in% c('bin', 'tps')) {
                means <- backsolve(U[[p]], forwardsolve(U[[p]], Wi.pbar_cell_centered[ , p] * nTreesPerCell + mu_current[p]*rowSums(Vinv[[p]])))
            } else {
                Vinv <- C / sigma2_current[p]
                means <- backsolve(U[[p]], forwardsolve(U[[p]], Wi.pbar_cell_centered[ , p] * nTreesPerCell))
            }
            alpha_next[,p] <- alpha_current[,p] <- means + backsolve(U[[p]], rnorm(nCells))
        }
    }

  # Gathering some indices outside the loop
  
    treeInd <- treeNonInd <- cellInd <- cellNonInd <- plotInd <- plotNonInd <- WtreeNonIndInfo <- list() 
    
    for(p in 1:P){
        treeInd[[p]]  <- which(y==p)
        treeNonInd[[p]] <- which(y!=p)
        cellInd[[p]]  <- cell[y==p]
        cellNonInd[[p]] <- cell[y!=p]
        plotInd[[p]]  <- plot[y==p]
        plotNonInd[[p]] <- plot[y!=p]
        WtreeNonIndInfo[[p]] <- cbind(treeNonInd[[p]], y[treeNonInd[[p]]])
    }
    
    nTreeInd  <- lapply(treeInd, length)
    nTreeNonInd <- lapply(treeNonInd, length)
  
    if(!identical(hyperpar, c(-0.5, 0)))
        stop("Error (runMCMC): joint sampling not set up to use any prior other than flat on sd scale.")

    count <- 0
    for(s in sampleIterates){
        count <- count + 1
                                        # sample the latent W's
        for(p in 1:P){
            if(nTreeInd[[p]] == 1){
                                        # need to use max, since don't have a matrix in this case; just use simple R version rtruncnorm
                W[treeInd[[p]], p] <- rtruncnorm(nTreeInd[[p]], alpha_current[cellInd[[p]], p] + h_current[plotInd[[p]], p],
                                                 max(W[treeInd[[p]], -p]), Inf) 
            } else {
                                        # since this calc involves smaller objects, do simple rowmax_cpp as rowmax2_cpp doesn't help here
                W[treeInd[[p]], p] <- rtruncnorm_cpp(nTreeInd[[p]], alpha_current[cellInd[[p]], p] + h_current[plotInd[[p]], p],
                                                     rowmax_cpp(W[treeInd[[p]],-p]), Inf)
            }
            if(numCores > 1) {
                                        # this saves another 0.5 s or so (for PLS, not directly relevant here for FIA), but means that RNG is not same as with single core operation
                
                W[treeNonInd[[p]], p] <- rtruncnorm_cpp_mp(nTreeNonInd[[p]], alpha_current[cellNonInd[[p]], p] + h_current[plotNonInd[[p]], p],
                                                           -Inf, rowmax2_cpp_mp(W, treeNonInd[[p]], p))
                                        # rowmax2_cpp_mp does not seem to help in terms of speed - 3.4sec/it vs 3.7sec per it for PLS code (not directly relevant here for FIA)
            } else {
                                        # non-MP version
                #W[treeNonInd[[p]], p] <- rtruncnorm_cpp(nTreeNonInd[[p]], alpha_current[cellNonInd[[p]], p] + h_current[plotNonInd[[p]], p],
                #                                       -Inf, rowmax2_cpp(W, treeNonInd[[p]], p))
                # recognizing I know what the max is simply from the y's
                W[treeNonInd[[p]], p] <- rtruncnorm_cpp(nTreeNonInd[[p]], alpha_current[cellNonInd[[p]], p] + h_current[plotNonInd[[p]], p],
                                                        -Inf, W[WtreeNonIndInfo[[p]]])
            }
            
        }
        
        # summarize the W's
        Wi.pbar_cell_centered <- compute_cell_sums_cpp(W,cell,nCells,P)/c(nTreesPerCell+1*(nTreesPerCell==0)) - compute_cell_sums_cpp(h_current*nTreesPerPlot,plotInfo$ID,nCells,P)/c(nTreesPerCell+1*(nTreesPerCell==0))

       # joint mu-alpha sample
        if(!nbhdStructure %in% c('bin', 'tps')) {
            for(p in 1:P){
                
                mu_next[p] <- rnorm(1, mu_current[p], mu_propSD[p])
                if(mu_next[p] < muBounds[1] || mu_next[p] > muBounds[2]) {
                    accept <- FALSE
                } else {
                    
                    VinvRowSums <- rowSums(Vinv[[p]])
                    VinvSum <- sum(Vinv[[p]])
                    
                    UtWi <- forwardsolve(U[[p]], Wi.pbar_cell_centered[ , p] * nTreesPerCell + mu_current[p] * VinvRowSums)
                    denominator <- 0.5*sum(UtWi^2) - 0.5*VinvSum*mu_current[p]^2
                    
                    UtWi <- forwardsolve(U[[p]], Wi.pbar_cell_centered[ , p] * nTreesPerCell + mu_next[p] * VinvRowSums)
                    numerator <- 0.5*sum(UtWi^2) - 0.5*VinvSum*mu_next[p]^2
                    
                    accept <- decide(numerator - denominator)
                }
                if(accept) {
                    numAcceptMu[p] <- numAcceptMu[p] + 1
                                        # sample alphas
                    means <- backsolve(U[[p]], UtWi)
                    alpha_next[,p] <- means + backsolve(U[[p]], rnorm(nCells))
                } else {
                    mu_next[p] <- mu_current[p]
                    alpha_next[,p] <- alpha_current[,p]
                }
            }
      
            mu_current <- mu_next
            alpha_current  <- alpha_next
    
            for(p in seq_len(P)) {
                                        # update alphas and etas/sigma2s    
                prop <- adaptScale[p] * adaptedL[[p]] %*% rnorm(2)
                sigma2_next[p] <- exp(log(sigma2_current[p]) + prop[1])
                eta_next[p] <- eta_current[p] + prop[2]
                if(sigma2_next[p] < 0.0001 || eta_next[p] < etaBounds[1] || eta_next[p] > etaBounds[2])  {
                    accept <- FALSE 
                } else { 
                                        #           B <- Vinv[[p]]
                                        #           diag(B) <- diag(B) + n
                                        #           update.spam.chol.NgPeyton(U[[p]], B)
                                        #           update.spam.chol.NgPeyton(Uc[[p]], Vinv[[p]])
                    denominator <-  - sum(log(diag(U[[p]]))) + sum(log(diag(Uc[[p]]))) + eta_current[p]
                    UtWi <- forwardsolve(U[[p]], Wi.pbar_cell_centered[ , p] * nTreesPerCell + mu_current[p] * rowSums(Vinv[[p]]))
                    denominator <- denominator + 0.5*sum(UtWi^2) - 0.5*sum(Vinv[[p]])*mu_current[p]^2 + 0.5*log(sigma2_current[p])
                    
                                        # terms for forward proposal
                    VinvProp <- C  # start with TPS
                    a <- 4 + 1/exp(2*eta_next[p])  # a = 4 + kappa^2 = 4 + (1/rho)^2 = 4 + (1/exp(eta)^2)
                    VinvProp@entries[Cindices$self] <- 4 + a*a
                    VinvProp@entries[Cindices$cardNbr] <- -2 * a
                    VinvProp <- VinvProp / (sigma2_next[p]*4*pi/exp(2*eta_next[p]))
                                        #          Vinv <- Vinv / sigma2_next[p]
                                        # simplify as Vinv from above *sigma2_current[p]/sigma2_next[p]
                    
                    B <- VinvProp
                    diag(B) <- diag(B) + nTreesPerCell
                    
                    Uprop <- update.spam.chol.NgPeyton(Uprop, B)
                    UcProp <- update.spam.chol.NgPeyton(UcProp, VinvProp)
                    
                    numerator <-  - sum(log(diag(Uprop))) + sum(log(diag(UcProp))) + eta_next[p]
                    UtWi <- forwardsolve(Uprop, Wi.pbar_cell_centered[ , p] * nTreesPerCell + mu_current[p] * rowSums(VinvProp))
                    numerator <- numerator + 0.5*sum(UtWi^2) - 0.5*sum(VinvProp)*mu_current[p]^2 + 0.5*log(sigma2_next[p])
                    
                    accept <- decide(numerator - denominator)
                }
                if(accept) {
                    numAcceptSigmaEta[p] <- numAcceptSigmaEta[p] + 1
                                        # sample alphas
                    means <- backsolve(Uprop, UtWi)
                    alpha_next[,p] <- means + backsolve(Uprop, rnorm(nCells))
                    
                    U[[p]] <- Uprop
                    Uc[[p]] <- UcProp
                                        # with update to spam 1.0.1 the machinations below no longer needed
                                        #U[[p]]@entries <- Uprop@entries
                                        #U[[p]]@entries[1] <- U[[p]]@entries[1] # force copy
                                        # since .Fortran in update.spam in 0.41-0 does not do a copy
                                        # and R 3.1 doesn't recognize that a copy needs to be made
                                        # only reason 3.0 does do copy is because it's in a list...!
                                        #Uc[[p]]@entries <- UcProp@entries
                                        #Uc[[p]]@entries[1] <- Uc[[p]]@entries[1]
                    Vinv[[p]] <- VinvProp
                } else {
                    sigma2_next[p] <- sigma2_current[p]
                    eta_next[p] <- eta_current[p]
                    alpha_next[,p] <- alpha_current[,p]
                }
                adaptVals[[p]][count, ] <- c(log(sigma2_next[p]), eta_next[p])
            }
            eta_current <- eta_next
            sigma2_current <- sigma2_next
            alpha_current  <- alpha_next
            warning("updating of h not incorporated into Lindgren")
            
        } else {  # nbhdStructure == 'bin' or == 'tps'
            
            for(p in 1:P){
                
                sigma2_next[p] <- exp(rnorm(1, log(sigma2_current[p]), logSigma2_propSD[p]))
                if(sigma2_next[p] < 0) {
                    accept <- FALSE 
                } else {
                    
                                        # terms for reverse proposal
                    denominator <- -(nCells-dfAdj)*log(sigma2_current[p])/2 - sum(log(diag(U[[p]]))) + 0.5*log(sigma2_current[p])
                    UtWi <- forwardsolve(U[[p]], Wi.pbar_cell_centered[ , p] * nTreesPerCell)
                    denominator <- denominator + 0.5*sum(UtWi^2)
                    
                                        # terms for forward proposal
                    Vinv <- C / sigma2_next[p]
                    B <- Vinv
                    diag(B) <- diag(B) + nTreesPerCell
                    
                    Uprop <- update.spam.chol.NgPeyton(Uprop, B)
                    
                    numerator <- -(nCells-dfAdj)*log(sigma2_next[p])/2 - sum(log(diag(Uprop))) + 0.5*log(sigma2_next[p])
                    UtWi <- forwardsolve(Uprop, Wi.pbar_cell_centered[ , p] * nTreesPerCell)
                    numerator <- numerator + 0.5*sum(UtWi^2)
                    
                    accept <- decide(numerator - denominator)
                }
                if(accept) {
                    numAcceptSigma2[p] <- numAcceptSigma2[p] + 1
                                        # sample alphas
                    means <- backsolve(Uprop, UtWi)
                    alpha_next[,p] <- means + backsolve(Uprop, rnorm(nCells))
                    U[[p]] <- Uprop
                                        # with update to spam 1.0.1 these steps not needed
                                        #U[[p]]@entries <- Uprop@entries
                                        #U[[p]]@entries[1] <- U[[p]]@entries[1] # force copy
                                        # since .Fortran in update.spam in 0.41-0 does not do a copy
                                        # and R 3.1 doesn't recognize that a copy needs to be made
                                        # only reason 3.0 does do copy is because it's in a list...!
                } else {
                    sigma2_next[p] <- sigma2_current[p]
                    alpha_next[,p] <- alpha_current[,p]
                }
            }
            sigma2_current <- sigma2_next
            alpha_current  <- alpha_next
            
                                        # alpha and sigma2 non-joint samples
            if(TRUE) {
                                        # for the moment include the non-joint sample as well so that alphas can move on their own; this adds about a second per iteration
                for(p in 1:P){
                    if(!nbhdStructure %in% c('bin', 'tps')) {
                        means <- backsolve(U[[p]], forwardsolve(U[[p]], Wi.pbar_cell_centered[ , p] * n + mu_current[p]*rowSums(Vinv[[p]])))
                    } else {
                        Vinv <- C / sigma2_current[p]
                        means <- backsolve(U[[p]], forwardsolve(U[[p]], Wi.pbar_cell_centered[ , p] * nTreesPerCell))
                    }
                    alpha_next[,p] <- means + backsolve(U[[p]], rnorm(nCells))
                    
                    if(!nbhdStructure %in% c('bin', 'tps')) {
                        ss <- sigma2_current[p] * t(alpha_next[,p] - mu_current[p]) %*% (Vinv[[p]] %*% (alpha_next[,p] - mu_current[p]))
                        sigma2_next[p] <- 1 / rgamma(1, shape = hyperpar[1] + (nCells)/2, scale = 1/(.5*ss + hyperpar[2])) 
                                        # need nCells-1 because of the zero eigenvalue in the precision matrix (but not for Lindgren)
                    } else {
                        ss <- sigma2_current[p] * t(alpha_next[,p]) %*% (Vinv %*% alpha_next[,p])
                        sigma2_next[p] <- 1 / rgamma(1, shape = hyperpar[1] + (nCells-dfAdj)/2, scale = 1/(.5*ss + hyperpar[2])) 
                                        # need nCells-1 because of the zero eigenvalue in the precision matrix
                    }
                    
                    if(is.na(sigma2_next[p])) { # sigma2 too small
                        sigma2_next[p] <- sigma2_current[p]
                        alpha_next[,p] <- alpha_current[,p]
                    }
                }
                sigma2_current <- sigma2_next
                alpha_current  <- alpha_next
                
                                        # update U[[p]] to reflect new sigma2 (if aggregated this will be done below)
                for(p in seq_len(P)) {
                    B <- Vinv <- C / sigma2_current[p]
                    diag(B) <- diag(B) + nTreesPerCell
                    if(!is.spam(B)) warning("B is not spam")
                    U[[p]] <- chol(B)  
                }
            }
        }
        
        
        Wi.pbar_plot_centered <-compute_cell_sums_cpp(W,plot,nPlots,P)/c(nTreesPerPlot +1*(nTreesPerPlot == 0)) - alpha_current[plotInfo$ID, ]
        
        # sample h and tau2
        Wcentered <- W - alpha_current[cell, ]
        for(p in 1:P){
            
            tau2_next[p] <- exp(rnorm(1, log(tau2_current[p]), logTau2_propSD[p]))
            if(tau2_next[p] < 0) {
                accept <- FALSE 
            } else {

                hV <- 1 / (nTreesPerPlot + 1/tau2_next[p])
                hMean <- nTreesPerPlot * Wi.pbar_plot_centered[ , p] * hV
                h_next[, p] <- rnorm(nPlots, hMean, sqrt(hV))
                denominator <- 0.5*log(tau2_current[p]) +
                    sum(dnorm(h_next[ , p], hMean, sqrt(hV), log = TRUE)) -
                    # next line faster than using dnorm
                    0.5*sum((Wcentered[ , p] - h_current[plot, p])^2) +
                    sum(dnorm(h_current[ , p], 0, sqrt(tau2_current[p]), log = TRUE))
                                    
                numerator <- 0.5*log(tau2_next[p])
                hV <- 1 / (nTreesPerPlot + 1/tau2_current[p])
                hMean <- nTreesPerPlot * Wi.pbar_plot_centered[ , p] * hV
                numerator <- 0.5*log(tau2_next[p]) +
                    sum(dnorm(h_current[ , p], hMean, sqrt(hV), log = TRUE)) -
                    0.5*sum((Wcentered[ , p] - h_next[plot, p])^2) +
                    sum(dnorm(h_next[ , p], 0, sqrt(tau2_next[p]), log = TRUE))
                
                accept <- decide(numerator - denominator)
            }
            if(accept) {
                numAcceptTau2[p] <- numAcceptTau2[p] + 1
            } else {
                tau2_next[p] <- tau2_current[p]
                h_next[,p] <- h_current[,p]
            }
        }
        tau2_current <- tau2_next
        h_current  <- h_next
        
      # h and tau2 non-joint samples
        if(TRUE) {
                                        # for the moment include the non-joint sample as well so that alphas can move on their own; this adds about a second per iteration
            for(p in 1:P){
                hV <- 1 / (nTreesPerPlot + 1/tau2_current[p])
                hMean <- nTreesPerPlot * Wi.pbar_plot_centered[ , p] * hV
                h_next[,p] <- rnorm(nPlots, hMean, sqrt(hV))
                
                ss <- sum(h_next[ , p]^2)
                tau2_next[p] <- 1 / rgamma(1, shape = hyperpar[1] + nPlots/2, scale = 1/(.5*ss + hyperpar[2])) 
                
                if(is.na(tau2_next[p])) { # sigma2 too small
                    tau2_next[p] <- tau2_current[p]
                    h_next[,p] <- h_current[,p]
                }
            }
            tau2_current <- tau2_next
            h_current  <- h_next
            
        }
    
        if (s%%thin==0){
            outputNcdfPtr <- nc_open(file.path(dataDir, outputNcdfName), write = TRUE)
            sigma2store[storeIndex, ] <- sigma2_next
            tau2store[storeIndex, ] <- tau2_next
            if(!nbhdStructure %in% c('bin', 'tps')) {
                muStore[storeIndex, ] <- mu_next
                etaStore[storeIndex, ] <- eta_next
            }
            for(p in 1:P)
                ncvar_put(outputNcdfPtr, taxa$taxonName[p], matrix(alpha_next[ , p], m1, m2), start = c(1, 1, storeIndex), count = c(-1, -1, 1))
            nc_close(outputNcdfPtr)
            storeIndex <- storeIndex + 1
        }
        
        
        if(s%%100 == 0)
            cat("Finished MCMC iteration ", s, " at ", date(), ".\n", sep = "")
        
        if(s%%adaptInterval == 0) {
            count <- 0
            if(!nbhdStructure %in% c('bin', 'tps')) {
                mu_propSD <- mu_propSD * adaptJump (n = rep(1, P), pjump = numAcceptMu / adaptInterval,
                                                    type = 'ben', i = s, K = adaptInterval)
                cat("Acceptance rate for mu/alpha joint proposals: ", round(numAcceptMu/adaptInterval, 2), ".\n", sep = " ")
                numAcceptMu <- rep(0, P)
                adaptScale <- adaptScale * adaptJump(n = rep(2, P),
                                                     pjump = numAcceptSigmaEta / adaptInterval,
                                                     type = 'ben',
                                                     i =s,
                                                     K = adaptInterval)
                cat("Acceptance rate for sigma/eta/alpha joint proposals: ", round(numAcceptSigmaEta/adaptInterval, 2), ".\n", sep = " ")
                numAcceptSigmaEta <- rep(0, P)
            } else {
                logSigma2_propSD <- logSigma2_propSD * adaptJump (n = rep(1, P), pjump = numAcceptSigma2 / adaptInterval,
                                                                  type = 'ben', i = s, K = adaptInterval)
                cat("Acceptance rate for sigma/alpha joint proposals: ", round(numAcceptSigma2/adaptInterval, 2), ".\n", sep = " ")
                numAcceptSigma2 <- rep(0, P)
                logTau2_propSD <- logTau2_propSD * adaptJump (n = rep(1, P), pjump = numAcceptTau2 / adaptInterval,
                                                              type = 'ben', i = s, K = adaptInterval)
                cat("Acceptance rate for tau/h joint proposals: ", round(numAcceptTau2/adaptInterval, 2), ".\n", sep = " ")
                numAcceptTau2 <- rep(0, P)
            }
            
            if(!nbhdStructure %in% c('bin', 'tps')) print(c(s, exp(eta_current))) else print(c(s, sigma2_current))
            
            if(!nbhdStructure %in% c('bin', 'tps')) {
                if(s > adaptStartDecay){
                    gamma2     <- 1 / (max(s - adaptStartDecay, 1)/adaptInterval + 3)^(.8) # only let degree of adaptation decay after burnin = 2*gStart
                } else{
                    gamma2     <- 1/ (      1        + 3)^(.8) #adapt every time as if it's the first time
                }
                
                                        # update prop covariance (loop thru taxa)
                for(p in seq_len(P)) {
                    newCov <- cov(adaptVals[[p]])
                    adaptVals[[p]] <- matrix(0, adaptInterval, 2)
                    adaptedCov[[p]] <- (1e-10)*diag(2) + adaptedCov[[p]] + gamma2*(newCov - adaptedCov[[p]])
                    adaptedL[[p]] <- t(chol(adaptedCov[[p]]))
                }
            }
        }
        
        if(s %% 250 == 0) {
            if(nbhdStructure %in% c('bin', 'tps')) 
                etaStore <- muStore <- eta_next <- mu_next <- adaptedCov <- adaptedL <- adaptScale <- NULL else logSigma2_propSD <- NULL
            save(logSigma2_propSD, logTau2_propSD, adaptedCov, adaptedL, adaptScale, alpha_next, h_next, sigma2_next, tau2_next, eta_next, mu_next, s, cell, .Random.seed, W, muStore, sigma2store, tau2store, etaStore, storeIndex, file = file.path(dataDir, paste0('lastState_', runID, '.Rda')))
        }
    } # end for s loop
    
    save(sigma2store, file = file.path(outputDir, paste0("sigma2_", runID, ".Rda")))
    save(tau2store, file = file.path(outputDir, paste0("tau2_", runID, ".Rda")))
    if(!nbhdStructure %in% c('bin', 'tps')) {
        save(muStore, file = file.path(outputDir, paste0("mu_", runID, ".Rda")))
        save(etaStore, file = file.path(outputDir, paste0("eta_", runID, ".Rda")))
    }
    invisible(NULL)
}

drawProportions <- function(latentNcdfPtr, outputNcdfPtr, numMCsamples = 1000, numInputSamples, secondThin = 1, nCells, taxa){
  require(doParallel)
  registerDoParallel(numCoresForProps)

  samples <- seq(1, numInputSamples, by = secondThin)
  P <- length(taxa)

  alphas <- phat <- array(0, c(nCells, P, length(samples)))
  if(secondThin == 1) {
    for(p in 1:P) 
      alphas[ , p, ] <- ncvar_get(latentNcdfPtr, varid = taxa[p], start = c(1, 1, 1), count = c(-1, -1, -1))
  } else {
    for(s in seq_along(samples)) {
      for(p in 1:P) 
        alphas[ , p, s] <- ncvar_get(latentNcdfPtr, varid = taxa[p], start = c(1, 1, samples[s]), count = c(-1, -1, 1))
    }
  }
  
  phat <- foreach(s = seq_along(samples)) %dopar% {
    set.seed(s)
    cat("Started sampling probabilities for (thinned) MCMC sample number ", samples[s], " at ", date(), ".\n", sep = "")
    compute_cell_probabilities_cpp(alphas[ , , s], numMCsamples, nCells, P)  
  }

  phat <- array(c(unlist(phat)), c(nCells, P, length(samples)))
  
  for(p in 1:P) 
    ncvar_put(outputNcdfPtr, taxa[p], phat[ , p, ], start = c(1, 1, 1), count = c(-1, -1, -1))
  
  invisible(NULL)
}



