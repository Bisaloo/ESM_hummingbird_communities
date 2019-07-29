rescale_hugo <- function(x, to = c(0, 1), from = range(x, na.rm = TRUE, finite = TRUE)) {
  if (zero_range(from) || zero_range(to)) {
    return(ifelse(is.na(x), NA, mean(to)))
  }

  (x - from[1]) / diff(from) * diff(to) + to[1]
}

# The following set of functions provides the full functionality to calculate
# the different distance matrices described and discussed in our study. Rather
# then writing a single function, the functionality is provided by a set of
# separate functions, some of which can be used to achieve subtasks of the
# overall task to calculate the desired distance matrices. The main function is
# decouple, which has as output the whole set of distance matrices. decouple
# depends on the other functions to be read in prior to using it. Of course one
# can also use the source function to load this entire script into an R
# workspace.

# Libraries required, install if needed
library(vegan)
library(ape)
library(scales)
library(multcomp)
library(FD)
library(hierfstat)
library(psych)

# ______________________________________________________________________________
# decouple function
# ______________________________________________________________________________
#
# decouple calculates a set of distance matrices from a given trait (or set of
# traits) and the accompanying phylogeny of the species. Trait(s) must be in a
# data.frame with species as rownames, the phylogeny in the phylo format of the
# ape package. Species names in the trait data and the phylogenetic tree need
# to be exactly the same and in the exact same order. The function returns a
# list, each element containing a distance matrix. The sequence of these
# matrices is as follows:
# Pdist, Fdist, dcPdist, dcFdist, covFPDist, where dc means decoupled, cov
# covariance, F functional, and P phylogenetic. keep_axes_method and ordi_adj_p
# are passed on to the keep_axes function (see above). When sqrt_tree = TRUE,
# the patristic distances of the phylogeny are squarerooted before further
# calculations as suggested by Letten and Cornwell (2014).

decouple <- function(traits,
                     tree,
                     keep_axes_method = c("lm", "ordi"),
                     ordi_adj_p = TRUE,
                     sqrt_tree = TRUE,
                     ...) {
  pfdist <- function(tree_dist, trait_dist, p, a) {
    ((a*tree_dist^p) + ((1-a)*trait_dist^p))^(1/p)
  }
  if(!is.ultrametric(tree, tol = 0.1)) {
    stop("The provided phylogeny needs to be ultrametric")
  }
  if(ncol(traits) == 1 & any(is.na(traits))) {
    warning("Only one trait was provided, which contained NAs. Species with NA were removed from trait data and the phylogenetic tree")
    traits <- na.omit(traits)
    int <- intersect(tree$tip.label, rownames(traits))
    tree <- drop.tip(tree, setdiff(tree$tip.label, int))
    traits <- traits[tree$tip.label, , drop = FALSE]
  }
  if(nrow(traits) != length(tree$tip.label)) {
    stop("The number of species in the traits data frame and the phylogenetic tree must be the same")
  }
  if(!all(rownames(traits) == tree$tip.label)) {
    stop("Species names in traits data frame and the phylogenetic tree must be identical and be in the same order")
  }
  if(sqrt_tree) {
    pdist <- rescale_hugo(sqrt(cophenetic(tree)))
  } else{
    pdist <- rescale_hugo(cophenetic(tree))
  }

  fdist <- rescale_hugo(trait2dist(traits, ...)[[1]])
  pcoa_tree <- hierfstat::pcoa(pdist, plotit = FALSE)
  cum_eig <- cumsum((pcoa_tree$valp)/sum(pcoa_tree$valp))
  max_axes <- which(cum_eig > 0.99)[1]-1
  pcoa_scores <- as.data.frame(pcoa_tree$vecp[, 1:max_axes])
  names(pcoa_scores) <- paste("a", as.character(1:ncol(pcoa_scores)), sep = "")
  rownames(pcoa_scores) <- rownames(pdist)
  kept <- keep_axes(traits, axes = pcoa_scores, method = keep_axes_method,
                    adj_p = ordi_adj_p)
  if(kept[1] == "none") {
    stop("No PCoA axes were selected when regressing them on the trait(s). Hence, no decoupled distances can be calculated")
  } else {
    if(any(unlist(lapply(traits, is.factor))) |
       any(unlist(lapply(traits, is.character)))) {
      sel <- c(which(unlist(lapply(traits, is.factor))),
               which(unlist(lapply(traits, is.character))))
      traits_dum <- as.data.frame(do.call("cbind", lapply(traits[sel],
                                                          dummy.code)))
      traits <- traits[-sel]
      traits <- cbind(traits, traits_dum)
    }
    rda_traits <- rda2dist(X = traits, Y = pcoa_scores[, kept], weight = TRUE)
    rda_tree <- rda2dist(Y = traits, X = pcoa_scores[, kept], weight = TRUE)
    dc_fdist <- rescale_hugo(rda_traits[[1]])
    dc_pdist <- rescale_hugo(rda_tree[[1]])
    covar_dist <- rescale_hugo(rda_traits[[2]])
    pfdist <- rescale_hugo(pfdist(as.matrix(pdist), as.matrix(fdist), 2, 0.5))
    dist_mat <- list(Pdist = as.matrix(pdist),
                     Fdist = as.matrix(fdist),
                     dcPdist = as.matrix(dc_pdist),
                     dcFdist = as.matrix(dc_fdist),
                     jointFPdist = as.matrix(covar_dist),
                     PFdist = as.matrix(pfdist))
    return(dist_mat)
  }
}
# end decouple

# ______________________________________________________________________________
# keep axes function
# ______________________________________________________________________________
#
# keep_axes determines wich axes to keep from a PCoA of the phylogenetic
# patristic distance matrix. It returns a character vector of axes names that
# are then used in the rda2dist function to subset the PCoA axes. This is done
# by the nested keep function for a single trait. keep_axes is thus a wrapper
# that runs keep over all traits. The "method" argument specifies if this
# selection is done from a multiple lm model with stepwise selection by AIC, or
# if from the ordistep function available in package vegan, which selects by
# p-values after permutation tests. Speed is highly dependent on the number of
# species (i.e. axes from the PCoA), even more so for "ordi" (because of
# permutations). If multiple traits are used, adj_p will decrease the default
# p-value of 0.05 for axis selection by vegan::ordistep to 0.05/number of
# traits. If the argument traits contains more than one trait, the selection
# procedure will be run for each trait separately, and the selected axes will be
# combined.
keep_axes <- function(traits, axes, method = c("lm", "ordi"), adj_p = TRUE) {
  method <- match.arg(method)
  keep <- function(trait, axes, method, p_cor) {
    na_trait <- names(trait)
    dat <- as.data.frame(cbind(trait, axes))
    if(is.factor(trait[, 1]) | is.character(trait[, 1])) {
      method = "ordi"
      trait <- dummy.code(trait[, 1])
      dat <- as.data.frame(cbind(trait, axes))
    }
    if(method == "lm"){
      form <- formula(paste(na_trait,
                            "~",
                            paste(colnames(axes),
                                  collapse = " + ")))
      step_lm <- step(lm(form, data = dat), trace = FALSE)
      select <- names(summary(step_lm)$coefficients[, 4])
      if(is.null(select) || length(select) == 0) {
        select <- NULL
      }else if (select[1] == "(Intercept)") {
        select <- select[-1]
      }
      return(select)
    } else if(method =="ordi") {
      assign("trait_tmp", trait, envir = .GlobalEnv)
      form <- formula(paste("trait_tmp", "~",
                            paste(colnames(axes), collapse = " + ")))
      step_ordi <- ordistep(rda(form, data = dat), Pin = p_cor, trace = FALSE)
      select <- dimnames(step_ordi$CCA$biplot)[[1]]
      rm(trait_tmp, envir = .GlobalEnv)
      return(select)
    }
  }
  n_traits <- ncol(traits)
  traits_l <- vector("list", n_traits)
  for(i in 1:n_traits) {
    traits_l[[i]] <- traits[, i, drop = F]
  }
  if(adj_p) {
    p_cor <- 0.05/n_traits
  } else{
    p_cor <- 0.05
  }
  kept <- lapply(traits_l, keep, axes, p_cor = p_cor, method = method)
  kept_uniq <- unique(unlist(kept))
  if(is.null(kept_uniq)) {
    kept_uniq = "none"
    }
  return(kept_uniq)
}
# end keep_axes
#
# ______________________________________________________________________________
# rda2dist function
# ______________________________________________________________________________
#
# rda2dist calculates RDA (redundancy analysis, from package vegan) on traits
# and phylogeny and provides the pairwise distances based on the scores of the
# constrained and unconstrained axes that are needed to caclulate the decoupled
# and covariance components in the decouple function. The function is wrapped
# around the dist_scores function that caclulates the distances from the scores,
# while the rest of rda2dist calculates the RDA using the rda function from
# vegan and extracts the relevant scores. Arguments X and Y are: (a) species
# (rows) x traits (columns) matrix or data frame, or (b) a species (rows) x
# selected phylogenetic axes (columns) matrix or data frame from a PCoA of the
# patristic distances of the phylogeny of the species. Which matrix is X and
# which Y determines the direction of the partial RDA. X is the dependent matrix
# on which the effect of Y will be removed, i.e. X~Y. Distances can be weighted
# by the eigenvalues of the RDA axes, or not, via the "weight" argument.
rda2dist<-function(X, Y, weight = TRUE){
  dist_scores <- function(rda.scores, w){
    # Calculates the (weighted) eucledian distances from the scores of the RDA
    # axes. Eigenvalues of the axes are used as weights.
    weight <- w/sum(w)
    scores <- rda.scores
    dissim <- matrix(0, dim(scores)[1], dim(scores)[1])
    colnames(dissim) <- rownames(scores)
    rownames(dissim) <- rownames(scores)
    for(i in 1:dim(scores)[2]){
      distai <- w[i]*(as.matrix(dist(scores[, i])))^2
      dissim <- distai + dissim
    }
    results <- sqrt(as.dist(dissim))
    return(results)
  }
  RDA <- rda(X = X, Y = Y)
  eig_unconstr <- RDA$CA$eig
  scores_unconstr<-RDA$CA$u
  eig_constr <- RDA$CCA$eig
  scores_constr<-RDA$CCA$u
  if (weight) {
    w_unconstr <- eig_unconstr/sum(eig_unconstr)
    w_constr <- eig_constr/sum(eig_constr)
  } else {
    w_unconstr <- rep(1, nrow(X))
    w_constr <- w_unconstr
  }
  result_unconstr <- dist_scores(scores_unconstr, w_unconstr)
  result_constr <- dist_scores(scores_constr, w_constr)
  results <- list(result_unconstr, result_constr)
  names(results) <- c("distances_from_unconstrained_axes",
                     "distances_from_constrained_axes")
  return(results)
}
# end rda2dist

# ______________________________________________________________________________
# trait2dist function
# ______________________________________________________________________________
#
# trait2dist is an excert from the dbFD function of the FD package by Etienne
# Laliberté. It applies the functionality of the dbFD function that is
# responsible for calculating a distance matrix from a single trait or a set of
# traits. Depending on the type of variables of the trait(s), different
# procedures are used to calculate the distance matrix. See the help for the
# dbFD and the gowdis functions for description of the arguments and further
# details. Note that the arguments of this function can be passed on from
# decouple by use of the "..." argument.

trait2dist <- function(
  x,
  w,
  stand.x = TRUE,
  ord = c("podani", "metric"),
  asym.bin = NULL,
  corr = c("sqrt", "cailliez", "lingoes", "none"),
  m = "max",
  dist.bin = 2,
  messages = TRUE) {
  require(ade4)
  tol <- .Machine$double.eps
  corr <- match.arg(corr)
  ord <- match.arg(ord)
  if (is.matrix(x) | is.data.frame(x)) {
    is.dist.x <- FALSE
    s.x <- dim(x)[1]
    t.x <- dim(x)[2]
    if (is.null(row.names(x)))
      stop("'x' must have row names.", "\n")
    else x.rn <- row.names(x)
  }
  if (is.vector(x) | is.factor(x)) {
    is.dist.x <- FALSE
    s.x <- length(x)
    t.x <- 1
    if (is.null(names(x)))
      stop("'x' must have names.", "\n")
    else x.rn <- names(x)
  }
  if (class(x)[1] == "dist" | class(x)[1] == "dissimilarity") {
    is.dist.x <- TRUE
    s.x <- attr(x, "Size")
    t.x <- 1
    if (is.null(attr(x, "Labels")))
      stop("'x' must have labels.", "\n")
    else x.rn <- attr(x, "Labels")
  }
  if (!missing(w) & is.dist.x)
    stop("When 'x' is a distance matrix, 'w' should be left missing.",
         "\n")
  if (!missing(w) & !is.dist.x) {
    if (!is.numeric(w) | length(w) != t.x)
      stop("'w' should be a numeric vector of length = number of traits.",
           "\n")
    else w <- w/sum(w)
  }
  if (missing(w))
    w <- rep(1, t.x)/sum(rep(1, t.x))
  if (is.matrix(x) | is.data.frame(x)) {
    x <- data.frame(x)
    if (t.x >= 2) {
      x.class <- sapply(x, data.class)
      if (any(x.class == "character"))
        x[, x.class == "character"] <- as.factor(x[,
                                                   x.class == "character"])
      else x <- x
      if (all(x.class == "numeric") & all(!is.na(x))) {
        if (length(unique(w)) == 1) {
          x.s <- apply(x, 2, scale, center = TRUE, scale = stand.x)
          x.dist <- dist(x.s)
        }
        else {
          x.dist <- gowdis(x, w = w, ord = ord, asym.bin = asym.bin)
        }
      }
      else {
        x.dist <- gowdis(x, w = w, ord = ord, asym.bin = asym.bin)
      }
    }
    if (t.x == 1) {
      if (is.numeric(x[, 1])) {
        if (all(!is.na(x))) {
          x.s <- apply(x, 2, scale, center = TRUE, scale = stand.x)
          x.dist <- dist(x.s)
        }
        if (any(is.na(x))) {
          pos.NA <- which(is.na(x), arr.ind = TRUE)
          x <- na.omit(x)
          x.s <- apply(x, 2, scale, center = TRUE, scale = stand.x)
          x.dist <- dist(x.s)
          row.excl.ab <- pos.NA[, 1]
          if (messages)
            cat("Warning: Species with missing trait values have been excluded.",
                "\n")
        }
      }
      if (is.factor(x[, 1]) | is.character(x[, 1])) {
        if (is.ordered(x[, 1]))
          x <- x
        else x[, 1] <- as.factor(x[, 1])
        if (any(is.na(x))) {
          pos.NA <- which(is.na(x), arr.ind = TRUE)
          x <- na.omit(x)
          x.rn <- x.rn[-pos.NA]
          if (messages)
            cat("Warning: Species with missing trait values have been excluded.",
                "\n")
        }
        if (is.ordered(x[, 1])) {
          x.s <- data.frame(rank(x[, 1]))
          names(x.s) <- x.rn
          x.dist <- dist(x.s)
        }
        else {
          x.f <- as.factor(x[, 1])
          x.dummy <- diag(nlevels(x.f))[x.f, ]
          x.dummy.df <- data.frame(x.dummy, row.names = x.rn)
          sequence <- 1:10
          if (all(dist.bin != sequence[any(sequence)]))
            stop("'dist.bin' must be an integer between 1 and 10.",
                 "\n")
          x.dist <- dist.binary(x.dummy.df, method = dist.bin)
        }
      }
    }
  }
  if (is.vector(x) & is.numeric(x)) {
    if (any(is.na(x))) {
      pos.NA <- which(is.na(x))
      x <- na.omit(x)
      x.rn <- x.rn[-pos.NA]
      if (messages)
        cat("Warning: Species with missing trait values have been excluded.",
            "\n")
    }
    else x <- x
    x.s <- scale(x, center = T, scale = stand.x)
    x.dist <- dist(x.s)
    x <- data.frame(x)
    dimnames(x) <- list(x.rn, "Trait")
  }
  if (is.vector(x) & is.character(x)) {
    x <- as.factor(x)
    if (any(is.na(x))) {
      pos.NA <- which(is.na(x))
      x <- na.omit(x)
      x.rn <- x.rn[-pos.NA]
      if (messages)
        cat("Warning: Species with missing trait values have been excluded.",
            "\n")
    }
    else x <- x
    dimnames(x) <- list(x.rn, "Trait")
    x.dummy <- diag(nlevels(x))[x, ]
    x.dummy.df <- data.frame(x.dummy, row.names = x.rn)
    sequence <- 1:10
    if (all(dist.bin != sequence[any(sequence)]))
      stop("'dist.bin' must be an integer between 1 and 10.",
           "\n")
    x <- data.frame(x)
    x.dist <- dist.binary(x.dummy.df, method = dist.bin)
  }
  if (is.ordered(x)) {
    if (any(is.na(x))) {
      pos.NA <- which(is.na(x))
      x <- na.omit(x)
      x.rn <- x.rn[-pos.NA]
      cat("Warning: Species with missing trait values have been excluded.",
          "\n")
    }
    else x <- x
    x <- data.frame(x)
    dimnames(x) <- list(x.rn, "Trait")
    x.dist <- gowdis(x, w = w, ord = ord, asym.bin = asym.bin)
  }
  if (is.factor(x) & !is.ordered(x)) {
    if (any(is.na(x))) {
      pos.NA <- which(is.na(x))
      x <- na.omit(x)
      x.rn <- x.rn[-pos.NA]
      if (messages)
        cat("Warning: Species with missing trait values have been excluded.",
            "\n")
    }
    else x <- x
    x.dummy <- diag(nlevels(x))[x, ]
    x.dummy.df <- data.frame(x.dummy, row.names = x.rn)
    sequence <- 1:10
    if (all(dist.bin != sequence[any(sequence)]))
      stop("'dist.bin' must be an integer between 1 and 10.",
           "\n")
    x.dist <- dist.binary(x.dummy.df, method = dist.bin)
    x <- data.frame(x)
    dimnames(x) <- list(x.rn, "Trait")
  }
  if (class(x)[1] == "dist" | class(x)[1] == "dissimilarity") {
    if (any(is.na(x)))
      stop("When 'x' is a distance matrix, it cannot have missing values (NA).",
           "\n")
    x.dist <- x
  }
  if (any(is.na(x.dist)))
    stop("NA's in the distance matrix.", "\n")
  if (!is.dist.x) {
    no.traits <- apply(x, 1, function(v) length(v[!is.na(v)]))
    if (any(no.traits == 0))
      stop("At least one species has no trait data.", "\n")
  }
  attr(x.dist, "Labels") <- x.rn
  if (is.euclid(x.dist))
    x.dist2 <- x.dist
  if (!is.euclid(x.dist)) {
    if (corr == "lingoes") {
      x.dist2 <- lingoes(x.dist)
      if (messages)
        cat("Species x species distance matrix was not Euclidean. Lingoes correction was applied.",
            "\n")
    }
    if (corr == "cailliez") {
      x.dist2 <- cailliez(x.dist)
      if (messages)
        cat("Species x species distance matrix was not Euclidean. Cailliez correction was applied.",
            "\n")
    }
    if (corr == "sqrt") {
      x.dist2 <- sqrt(x.dist)
      if (!is.euclid(x.dist2))
        stop("Species x species distance matrix was still is not Euclidean after 'sqrt' correction. Use another correction method.",
             "\n")
      if (is.euclid(x.dist2))
        if (messages)
          cat("Species x species distance matrix was not Euclidean. 'sqrt' correction was applied.",
              "\n")
    }
    if (corr == "none") {
      x.dist2 <- quasieuclid(x.dist)
      if (messages)
        cat("Species x species distance was not Euclidean, but no correction was applied. Only the PCoA axes with positive eigenvalues were kept.",
            "\n")
    }
  }
  x.pco <- dudi.pco(x.dist2, scannf = FALSE, full = TRUE)
  traits <- round(x.pco$li, .Machine$double.exponent)
  res <- list(x.dist2, traits, x.pco)
  return(res)
}
# end trait2dist
