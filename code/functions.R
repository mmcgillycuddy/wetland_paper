# Set up starting params to be the same as the estimates from model.fitted
# Obtain initial parameter estimates model.fitted is the model object from glmmTMB for which 
startParFitted <- function(model.fitted, model.noFit){
  ## Getting model estimates from model 1a in the same form as expected for start option
  fit.par <- model.fitted$obj$env$last.par.best
  fit.beta <- fit.par[names(fit.par) == "beta"] # fit model beta parameters
  fit.theta <- fit.par[names(fit.par) == "theta"] # fit model theta parameters
  fit.b <- fit.par[names(fit.par) == "b"] # fit model b parameters
  
  ## assign names to parameters so they can be matched to other model params
  # beta names (fixed params)
  names.beta <- colnames(getME(model.fitted, "X"))
  names(fit.beta) <- names.beta
  
  ## Get general structure of random effects in the model
  cnms <- model.fitted$modelInfo$reTrms[["cond"]]$cnms   ## list of (named) terms and X columns
  flist <- model.fitted$modelInfo$reTrms[["cond"]]$flist ## list of grouping variables
  nc <- vapply(model.fitted$modelInfo$reStruc$condReStruc, function(x) x$blockSize, numeric(1)) ## number of RE params per block
  nb <- vapply(model.fitted$modelInfo$reStruc$condReStruc, function(x) x$blockReps, numeric(1)) ## number of blocks per RE (may != nlevs in some cases)
  nt <- vapply(model.fitted$modelInfo$reStruc$condReStruc, function(x) x$blockNumTheta, numeric(1)) ## number of blocks per RE (may != nlevs in some cases)
  bc <- vapply(model.fitted$modelInfo$reStruc$condReStruc, function(x) x$blockCode, numeric(1)) ## number of blocks per RE (may != nlevs in some cases)
  levs <- lapply(flist, levels)
  diag.block <- which(bc == 0) # diag has 0 code
  rr.block <- which(bc == 9) # diag has 0 code
  # theta (random effects params) names 
  diag.theta.names <- model.fitted$modelInfo$reTrms$cond$cnms[[diag.block]]
  nfl <- nt[[rr.block]] # number of factor loadings (same for every model)
  names(fit.theta) <- c(diag.theta.names, paste0("fl", 1:nfl))
  ## Names of bs (latent variables)
  ## bs for diag (For now assuming diag is the first (1) re in the model)
  # need to calculate how many bs are in the diagonal structure (nc * nb)
  names(fit.b)[1:(nc[[diag.block]] * nb[[diag.block]])] <- paste0(cnms[[diag.block]], ".", rep(levs[[diag.block]], each = nc[[diag.block]]))
  
  ### Starting parameters for comp_time.water_fixed model 
  start.beta <- model.noFit$parameters$beta
  start.theta <- model.noFit$parameters$theta
  start.b <- model.noFit$parameters$b
  
  ### Assigning names for these starting parameters
  names(start.beta) <- colnames(model.noFit$data.tmb$X)  ## assign names to beta so we can match based on the old model with just water
  ### Assuming the model.noFit has the same order of random effects as before
  names(start.theta) <- c(model.noFit$condList$reTrms$cnms[[diag.block]], paste0("fl", 1:nfl))
  names(start.b) <- rep("b", length(start.b))
  
  cnms <- model.noFit$condList$reTrms$cnms[[diag.block]]  ## list of (named) terms and X columns
  flist <- model.noFit$condList$reTrms$flist[[diag.block]] ## list of grouping variables
  levs <- levels(flist)
  nc <- model.noFit$condReStruc[[diag.block]]$blockSize
  nb <- model.noFit$condReStruc[[diag.block]]$blockReps
  names(start.b)[1: (nc*nb)] <- paste0(cnms, ".", rep(levs, each = nc))
  
  ### Initialise the values to the previous model by matching names
  ### start.beta 
  start.beta[names(fit.beta)] <- fit.beta
  # start.theta for diag
  start.theta[intersect(names(start.theta), names(fit.theta))] <- fit.theta[ intersect(names(start.theta), names(fit.theta))]
  ## bs for diag 
  start.b[intersect(names(start.b), names(fit.b))] <- fit.b[ intersect(names(start.b), names(fit.b))]
  start.pars <- list(start.beta = start.beta, start.theta = start.theta, start.b = start.b)
  
  return(start.pars)
}

library(ochRe)
### For plots
clrs2 <- ochre_palettes$nolan_ned[c(3,5)]
clrs3 <- ochre_palettes$lorikeet[4:2]
trait_clrs <- ochre_palettes$jumping_frog[c(2, 3, 4)]