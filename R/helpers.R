
# panel cor plots taken from Zuur et al. somewhere

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y, use = "pairwise.complete.obs"),2)
  if(missing(cex.cor)) cex.cor <- 4
  text(0.5, 0.5, r, cex = cex.cor * abs(r))
}



# function to get median and 90%credible interval from Stan fit objects

get_ci <- function(object, pars = pars, probs = c(0.05, 0.5, 0.95), 
                   names = NULL){
  ci <- rstan::summary(object, pars = pars, 
                       probs = probs)$summary
  ci <- ci[ , -c(1:3, ncol(ci) ) ]
  ci <- round(ci, 2)
  ci[ , "n_eff"] <- round(ci[ , "n_eff"])
  ci <- as.data.frame(ci)
  names(ci)[ 1:length(probs) ] <- paste("p", probs, sep = "") 
  if ( !is.null(names) ){ 
    ci$parameter <- factor(names)
    ci <- ci[c("parameter", setdiff(names(ci), "parameter"))]
  }
  return(ci)
}



# annotation custom function to put inset plots into ggplot2 facet plots
# (found here: https://stackoverflow.com/questions/37867758/insetting-on-facet-grided-and-grid-arrangeed-plot)

annotation_custom2 <- 
  function (grob, xmin = -Inf, xmax = Inf, 
            ymin = -Inf, ymax = Inf, data){ 
    layer(data = data, stat = StatIdentity, position = PositionIdentity,
          geom = ggplot2:::GeomCustomAnn, inherit.aes = TRUE, 
          params = list(grob = grob, xmin = xmin, xmax = xmax, 
                        ymin = ymin, ymax = ymax))
  }




# bayesR2 method, follows: https://avehtari.github.io/bayes_R2/bayes_R2.html

bayesR2_res <- function(fit = fit, y = y, ypred = ypred) {
  ypred <- rstan::extract(fit, pars = ypred, permuted = TRUE)[[1]]
  ypred <- as.matrix(ypred)
  e <- -1 * sweep(ypred, 2, y)
  var_ypred <- matrixStats::rowVars(ypred)
  var_e <- matrixStats::rowVars(e)
  return(as.matrix(var_ypred / (var_ypred + var_e)))
}




bayesR2 <- function(fit = fit, ypred = ypred, error = error, family = c("normal", "Gamma", "poisson")) { 
  
  if(missing(family)) {
    stop("Warning: family not specified")
  } 
  
  family  <- match.arg(family)
  
   bayesR2_normal <- function(fit = fit, ypred = ypred, error = error) {
    ypred     <- rstan::extract(fit, pars = ypred, permuted = TRUE)[[1]]
    ypred     <- as.matrix(ypred)
    var_fit   <- apply(ypred, 1, var)
    var_res   <- extract(fit, pars = error, permuted = TRUE)[[1]]
    var_res   <- as.matrix(var_res)
    var_res   <- var_res^2
    R2        <- var_fit / (var_fit + var_res)
    names(R2) <- "bayes_R2"
    R2
  }
  
  bayesR2_Gamma <- function(fit = fit, ypred = ypred, error = error) {
    ypred     <- rstan::extract(fit, pars = ypred, permuted = TRUE)[[1]]
    ypred     <- as.matrix(ypred)
    var_fit   <- apply(ypred, 1, var)
    var_res   <- extract(fit, pars = error, permuted = TRUE)[[1]]
    var_res   <- as.matrix(var_res)
    R2        <- var_fit / (var_fit + var_res)
    names(R2) <- "bayes_R2"
    R2
  }
  
  
  bayesR2_poisson <- function(fit = fit, ypred = ypred, error = error) {
  ypred     <- rstan::extract(fit, pars = ypred, permuted = TRUE)[[1]]
  ypred     <- as.matrix(ypred)
  var_fit   <- apply(ypred, 1, var)
  var_res   <- extract(fit, pars = error, permuted = TRUE)[[1]]
  var_res   <- as.matrix(var_res)
  warning("no variance defined for poisson Bayes R2 here")
  warn <- paste("failed")
  warn
  #var_res   <- var_res^2
  #R2        <- var_fit / (var_fit + var_res)
  #names(R2) <- "bayes_R2"
  #R2
  }
  
  switch(family, 
         normal   =  bayesR2_normal( fit, ypred, error),
         Gamma    =  bayesR2_Gamma(  fit, ypred, error),
         poisson  =  bayesR2_poisson(fit, ypred, error)
  )
  
}


# posterior predictive checks using the bayesplo package, kernel density overlays
ppc_dens_overlay_stanfit <- function(object, pars = pars, y = y, n = 50){
  yrep <- as.matrix(object, pars = pars)
  yrep <- yrep[ sample(nrow(yrep), n) , ]
  bayesplot::ppc_dens_overlay(y, yrep) + theme(plot.title = element_text(hjust=0.3, size = 9),
                                    plot.subtitle = element_text(hjust=0.2, size = 9))
}

