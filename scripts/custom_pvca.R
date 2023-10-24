pvcaBatchAssess <- function (expr, annot, batch.factors, threshold, threads=1, skip.unique=F) {
    library(RhpcBLASctl)
    blas_p = blas_get_num_procs()
    blas_set_num_threads(1)
    omp_p = omp_get_num_procs()
    omp_set_num_threads(1)

    theDataMatrix <- as.matrix(expr)
    dataRowN <- nrow(theDataMatrix)
    dataColN <- ncol(theDataMatrix)

    ########## Center the data (center rows) ##########
    theDataMatrixCentered <- matrix(data = 0, nrow = dataRowN, ncol = dataColN)
    theDataMatrixCentered_transposed = apply(theDataMatrix, 1, scale, center = TRUE, scale = FALSE)
    theDataMatrixCentered = t(theDataMatrixCentered_transposed)

    ########## Compute correlation matrix &  Obtain eigenvalues ##########
    # this allows to get more PCs than using prcomp directly
    # since we are doing cor and not cov the scaling does not matter
    theDataCor <- cor(theDataMatrixCentered)
    eigenData <- eigen(theDataCor)
    eigenValues = eigenData$values
    ev_n <- length(eigenValues)
    eigenVectorsMatrix = eigenData$vectors
    eigenValuesSum = sum(eigenValues)
    percents_PCs = eigenValues / eigenValuesSum

    ##===========================================
    ##    Getting the experimental information
    ##===========================================
    #expInfo <- pData(abatch)[,batch.factors]
    exp_design <- as.data.frame(annot)[,batch.factors, drop=FALSE]
    expDesignRowN <- nrow(exp_design)
    expDesignColN <- ncol(exp_design)

    ########## Merge experimental file and eigenvectors for n components ##########
    my_counter_2 = 0
    my_sum_2 = 1
    for (i in ev_n:1){
        my_sum_2  = my_sum_2 - percents_PCs[i]
        if ((my_sum_2) <= threshold ){
            my_counter_2 = my_counter_2 + 1
        }

    }

    if (my_counter_2 < 3){
        pc_n  = 3
    }else {
        pc_n = my_counter_2
    }

    ## pc_n is the number of principal components to model
    pc_data_matrix <- matrix(data = 0, nrow = (expDesignRowN*pc_n), ncol = 1)
    mycounter = 0
    for (i in 1:pc_n){
        for (j in 1:expDesignRowN){
            mycounter <- mycounter + 1
            pc_data_matrix[mycounter,1] = eigenVectorsMatrix[j,i]
        }
    }

    AAA <- exp_design[rep(1:expDesignRowN,pc_n),]
    AAA[] <- lapply(AAA, factor)
    Data <- cbind(AAA,pc_data_matrix)

    ########## Mixed linear model ##########
    op <- options(warn = (-1))
    effects_n = expDesignColN  + choose(expDesignColN, 2) + 1
    randomEffectsMatrix <- matrix(data = 0, nrow = pc_n, ncol = effects_n)

    ##============================#
    ##    Get model functions
    ##============================#
    variables <-c (colnames(exp_design))
    model.func <- c()
    index <- 1

    ##    level-1
    for (i in 1:length(variables))
    {
        skip = F
        if(skip.unique) {
            if(length(unique(levels(exp_design[,variables[i]]))) == expDesignRowN) {
                warning(sprintf("Skipping the batch factor %s since levels are unique along the data.", variables[i]))
                skip = T
            }
        }
        if(!skip) {
            mod = paste("(1|", variables[i], ")",   sep="")
            model.func[index] = mod
            index = index + 1
        }
    }

    ##    two-way interaction
    if(length(variables) > 1){
      for (i in 1:(length(variables)-1))
      {
          for (j in (i+1):length(variables))
          {
              skip = F
              if(skip.unique) {
                  combs <- interaction(exp_design[,variables[i]], exp_design[,variables[j]], drop=T)
                  if(length(unique(levels(combs))) == expDesignRowN) {
                      warning(sprintf("Skipping the interaction term of %s and %s since their combinations are unique along the data.", variables[i], variables[j]))
                      skip = T
                  }
              }
              if(!skip) {
                  mod = paste("(1|", variables[i], ":", variables[j], ")",   sep="")
                  model.func[index] = mod
                  index = index + 1
              }
          }
      }
    }

    function.mods <- paste (model.func , collapse = " + ")

    ##============================#
    ##    Get random effects          #
    ##============================#

    res = pbapply::pblapply(1:pc_n, function(i){
      y = (((i-1)*expDesignRowN)+1)
      funct <- paste ("pc_data_matrix", function.mods, sep =" ~ ")
      Rm1ML <- lmer( funct ,
                     Data[y:(((i-1)*expDesignRowN)+expDesignRowN),,drop=FALSE],
                     REML = TRUE, verbose = FALSE, na.action = na.omit)
      return(c(unlist(VarCorr(Rm1ML)), resid=sigma(Rm1ML)^2))
    }, cl=threads)

    effectsNames = names(res[[1]])
    randomEffectsMatrix = matrix(unlist(res), nrow=length(res), ncol=length(res[[1]]), byrow = T)
    ########## Standardize Variance ##########
    randomEffectsMatrixStdze <- randomEffectsMatrix / rowSums(randomEffectsMatrix)

    ########## Compute Average Weighted Proportions ##########
    prop.var <- colSums(randomEffectsMatrixStdze*percents_PCs[1:pc_n])
    randomEffectsMatrixWtAveProp <- prop.var/sum(prop.var)

    # restore to original values
    blas_set_num_threads(blas_p)
    omp_set_num_threads(omp_p)

    return(list(dat=randomEffectsMatrixWtAveProp, label=effectsNames))
}
