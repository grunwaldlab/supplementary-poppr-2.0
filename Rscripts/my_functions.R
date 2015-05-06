#' Generate random sample names
#' 
#' Generates random sample names by sampling 10 letters with replacement.
#' 
#' @param n the number of random names to generate
#' @author Zhian N. Kamvar
getNames <- function(n){
  vapply(1:n, function(x) paste(sample(letters, 10, replace = TRUE), collapse = ""), character(1))
}

#' Mutate a single SNP.
#' 
#' One random SNP in one random chromosome. Used by sample_mutator.
#' 
#' @param chrom a vector of raw elements
#' @param mutation a single integer representing a bit mask for a random
#'   mutation.
#'   
#' @return a raw vector with one bit changed.
#' 
#' @details the XOR gate is utilized to flip a single bit. 
#'   
#' @author Zhian N. Kamvar
#' @keyword internal
snp_mutator <- function(chrom, mutation = sample(2^(0:7), 1)){
  posi        <- sample(length(chrom) - 1, 1) # sample chunk of 8 loci
  orig        <- as.integer(chrom[posi])      # convert to integer
  chrom[posi] <- as.raw(bitwXor(orig, mutation)) # Exclusive OR will flip the bit. 
  return(chrom)
}

#' Mutate a single SNPbin object. 
#' 
#' propogate mutations randomly across chromosomes of a SNPbin object.
#' 
#' @param snpbin a SNPbin object
#' @param mu the rate of mutation
#' @param nLoc the number of loci in the object.
#' @param rawchars an integer vector indicating the masking bits for mutation.
#' 
#' @note The number of mutations is chosen from a single draw of a poisson 
#'   distribution. Each mutation is propogated on a random chromosome in a
#'   random location.
#' 
#' @return a mutated SNPbin object
#' @author Zhian N. Kamvar
sample_mutator <- function(snpbin, mu, nLoc, rawchars = 2^(0:7)){
  nmutations <- rpois(1, lambda = round(nLoc*mu))
  for (i in seq(nmutations)){
    chrom_index               <- sample(length(snpbin@snp), 1)
    chrom                     <- snpbin@snp[[chrom_index]]
    snpbin@snp[[chrom_index]] <- snp_mutator(chrom, sample(rawchars, 1))
  }
  return(snpbin)
}

#' Insert mutation into a genlight object
#' 
#' This will randomly insert mutations at a given rate per sample. 
#' 
#' @param glt a genlight object
#' @param mu a mutation rate
#' @param samples a vector indicating which samples should be mutated
#' 
#' @return a mutated genlight object
#' @author Zhian N. Kamvar
pop_mutator <- function(glt, mu = 0.05, samples = TRUE){
  rawchrs <- 2^(0:7)
  glt@gen[samples] <- lapply(glt@gen[samples], sample_mutator, mu, nLoc(glt), rawchrs)
  return(glt)
}

#' Generate missing data in a single chromosome
#' 
#' internal function utilized by NA_generator
#' 
#' @param chrom a vector of type raw.
#' @param NA.posi the position of missing data across the chromosome. Note that 
#'   these positions refer to positions of individual bits, not the individual
#'   raw elements.
#' @param rawchars a vector representing the masking bits. Passed for
#'   convenience.
#'   
#' @details This is the trickiest part to generate missing data. Missing data in
#'   genlight objects are represented by an integer vector representing the
#'   missing position, but in the binary coding, all missing data are
#'   represented as "off" bits. When generating missing data, all bits in the
#'   missing positions must be flipped to the "off" position using an XOR gate.
#'   Given a single NA position, when challenged with an "on" bit in an OR gate,
#'   if the result is "on" the same as the input, then the bit must be flipped
#'   because the input was in the "on" position.
#'   
#' @author Zhian N. Kamvar
#' @keyword internal
NA_zeromancer <- function(chrom, NA.posi, rawchars = 2^(0:7)){
  nas <- ceiling(NA.posi/8) # Getting the location in the RAW vector
  zero_bits <- NA.posi %% 8 # Getting the location of the locus in a RAW element.
  zero_bits[zero_bits == 0] <- 8
  
  for (i in seq(length(nas))){
    eight_bit <- as.integer(chrom[nas[i]])
    the_locus <- rawchars[zero_bits[i]]
    # If the locus does not change in the OR, then there is a 1 and it needs to
    # be changed to a zero with XOR.
    if (eight_bit == bitwOr(eight_bit, the_locus)){
      chrom[nas[i]] <- as.raw(bitwXor(eight_bit, the_locus))
    }
  }
  return(chrom)
}

#' Generate missing data for a SNPbin object
#' 
#' internal function utilized by pop_NA
#' 
#' @param snpbin a SNPbin object
#' @param nloc the number of loci in the SNPbin object
#' @param na.perc the percentage of missing data to propogate
#' @param rawchars a vector representing masking bits. Passed for convenience.
#'   
#' @details the number of missing loci is generated here by drawing from a
#'   poisson distribution.
#'   
#' @author Zhian N. Kamvar
#' @keyword internal
NA_generator <- function(snpbin, nloc, na.perc = 0.01, rawchars = 2^(0:7)){
  nas <- rpois(1, lambda = round(nloc*na.perc))
  NA.posi <- sort(sample(nloc, nas))
  for (i in seq(length(snpbin@snp))){
    snpbin@snp[[i]] <- NA_zeromancer(snpbin@snp[[i]], NA.posi, rawchars)
  }
  snpbin@NA.posi <- NA.posi
  return(snpbin)
}

#' Propogate missing data in a genlight object
#' 
#' @param na.perc the percentage of missing data per sample
#' @param parallel a logical specifying whether or not parallel processing
#'   should be utilized
#' @param n.cores the number of cores to use with parallel processing.
#'   
#' @return a genlight object with missing data propogated.
#'   
#' @author Zhian N. Kamvar
pop_NA <- function(glt, na.perc = 0.01, parallel = require('parallel'), n.cores = 2L){
  rawchars <- 2^(0:7)
  if (parallel){
    glt@gen <- mclapply(glt@gen, NA_generator, nLoc(glt), na.perc, rawchars,
                        mc.cores = getOption("mc.cores", n.cores))
  } else {
    glt@gen <- lapply(glt@gen, NA_generator, nLoc(glt), na.perc, rawchars)    
  }

  return(glt)
}

#' Create crossover in a SNPbin object
#' 
#' @param snpbin a SNPbin object
#'   
#' @return a SNPbin object with its chromosomes crossed over
#'   
#' @details this will randomly choose a cut point in the chromosomes of a
#'   diploid snpbin object and create a crossover event.
#'   
#' @author Zhian N. Kamvar
#' @keyword internal
crossover <- function(snpbin){
  chr1 <- snpbin@snp[[1]]
  chr2 <- snpbin@snp[[2]]
  chrlen <- length(chr1)
  cutpoint <- sample(2:(chrlen - 1), 1)
  first <- 1:(cutpoint - 1)
  last  <- cutpoint:chrlen
  snpbin@snp[[1]] <- c(chr1[first], chr2[last])
  snpbin@snp[[2]] <- c(chr2[first], chr1[last])
  return(snpbin)
}

#' mate two individuals from a genlight object
#' 
#' @param snpbin a genlight object (yes, I know its mis-named).
#' @param ind1 the index of parent 1
#' @param ind2 the index of parent 2
#' @return a SNPbin object
#'   
#' @details this will take two parental samples (they do not have to be
#'   different), have them each undergo crossover, and then randomly sample one
#'   chromosome from each to create the offspring.
#' @author Zhian N. Kamvar
#' @keyword internal
mate <- function(snpbin, ind1, ind2){
  snpbin@gen[[ind1]] <- crossover(snpbin@gen[[ind1]])
  snpbin@gen[[ind2]] <- crossover(snpbin@gen[[ind2]])
  snpout <- snpbin@gen[[ind1]]
  snpout@snp[[1]] <- snpbin@gen[[ind1]]@snp[[sample(1:2, 1)]]
  snpout@snp[[2]] <- snpbin@gen[[ind2]]@snp[[sample(1:2, 1)]]
  return(snpout)
}

#' randomly mate a genlight object once
#' 
#' @param glt a genlight object
#' @param err the number of mutations per mating event.
#'   
#' @details Samples are chosen with replacement from the population twice to
#'   create each set of parents. They are then mated and mutated. See getSims
#'   for details
#'   
#' @return a genlight object with the same number of individuals as glt
#' @author Zhian N. Kamvar
#' @keyword internal
random_mate <- function(glt, err){
  mat_pair <- matrix(integer(1), nrow = nInd(glt), ncol = 2)
  mat_pair[, 1] <- sample(nInd(glt), replace = TRUE)
  mat_pair[, 2] <- sample(nInd(glt), replace = TRUE)
  res <- apply(mat_pair, 1, function(x) mate(glt, x[1], x[2]))
  glt@gen <- res
  if (err > 0) glt <- pop_mutator(glt, err)
  return(glt)
}

#' randomly mate a genlight object for n generations.
#' 
#' @param glt a genlight object
#' @param err the number of mutations per mating event
#' @param gen the number of generations to mate the same population
#' 
#' @return a genlight object
#' @author Zhian N. Kamvar
random_mate_gen <- function(glt, err = 5e-3, gen = 1){
  for (i in seq(gen)){
    glt <- random_mate(glt, err)
  }
  return(glt)
}

#' Generate simulated SNP data with error and missing
#' 
#' This function is a wrapper to glSim, but it gives the user the option to 
#' introduce error, missing data, and clonality into the simulation.
#' 
#' @param z a dummy parameter for use with lapply
#' @param n the number of samples to simulate
#' @param snps the total number of snps to simulate
#' @param strucrat the fraction of unstructured snps
#' @param clone if \code{TRUE}, the data will be randomly sampled with 
#'   replacement before error propogation
#' @parm err the percent error per sample. The number of mutations per sample is
#'   determined by a draw from a poisson distribution.
#' @param na.perc the percent missing per sample. The number of missing per
#'   sample is determined by a draw from a posson distribution. Missing is
#'   always propogated after error propogation.
#' @param n.cores a parameter passed on to glSim, specifying the number of
#'   parallel processes.
#' @param mate_gen the nubmer of generations to "randomly mate" the population
#'   before error propogation/resampling. See details.
#' @param mate_err the probability of mutation per mating event.
#'   
#' @return a genlight object
#'   
#' @details Mating here is a very simple method. The population is sampled with
#'   replacement twice to create each set of parents. The resulting population
#'   will always have the same number of individuals as the starting population.
#'   Before mating, each parent will undergo a crossover event. Mating consists
#'   of randomly selecting one chromosome from each parent and then mutating the
#'   resultant offspring at the rate specified by mate_err.
#'   
#' @author Zhian N. Kamvar
getSims <- function(z = 1, n = 10, snps = 1e6, strucrat = c(0.25, 0.75), 
                    clone = TRUE, err = 0.1, na.perc = 0.1, n.cores = 4, 
                    mate_gen = NULL, mate_err = 5e-3, ...){
#   lam <- sample(strucrat, 1)
#   if (lam == 1){
#     perc <- 1
#   } else {
#     perc <- ceiling(sample(rpois(1000, lambda = snps*lam), 1)/snps)    
#   }

  perc <- sample(strucrat, 1)

  # n <- sample(rpois(1000, lambda = n), 1)
  the_names <- getNames(n)
  res <- glSim(n, n.snp.nonstruc = snps*perc, n.snp.struc = snps*(1-perc), 
               n.cores = n.cores, ...)
  if (!is.null(mate_gen)){
    res <- random_mate_gen(res, mate_err, mate_gen)
  }
  if (clone){
    clones <- sample(n, replace = TRUE)
    the_names <- paste(the_names[clones], 1:n, sep = ".")
    res <- res[clones]
    clones <- duplicated(clones)
    # res <- pop_mutator(res, err, clones)
  }
  if (err > 0){
    res <- pop_mutator(res, err)    
  }
  if (na.perc > 0){
    res <- pop_NA(res, na.perc = na.perc, n.cores = n.cores)    
  }
  indNames(res) <- the_names
  return(res)
}

#' Helper function to display bitcode.
#' 
#' @param y a raw string
#' @return a character string of 1s and 0s
#' @example
#' # Bitshifters
#' shifts <- as.raw(2^(0:7))
#' binary_char_from_hex(shifts)
binary_char_from_hex <- function(y){
  vapply(y, function(x) paste(as.integer(rawToBits(x)), collapse=""), character(1))
}

#' Utilize all algorithms of mlg.filter
#' 
#' This function is a wrapper to mlg.filter. It will calculate all of the stats 
#' for mlg.filter utilizing all of the algorithms.
#' 
#' @param x a genind or genlight object
#' @param distance a distance function or matrix
#' @param threshold a threshold to be passed to mlg.filter
#' @param stats what statistics should be calculated.
#' @param missing how to treat missing data with mlg.filter
#' @param If the threshold is a maximum threshold, should the statistics be 
#'   plotted (Figure 2)
#' @param nclone the number of multilocus genotypes you expect for the data. 
#'   This will draw horizontal line on the graph at the value nclone and then
#'   vertical lines showing the cutoff thresholds for each algorithm.
#'   
#' @return a list of results from mlg.filter from the three algorithms.
#'   
#' @author Zhian N. Kamvar, Jonah C. Brooks
filter_stats <- function(x, distance = bitwise.dist, threshold = 1, 
                         stats = "All", missing = "ignore", plot = FALSE, 
                         nclone = NULL, ...){
  if (!"dist" %in% class(distance)){
    distmat <- distance(x, ...)
  }
  f <- mlg.filter(x, threshold, missing, algorithm = "f", distance = distmat, 
                  stats = stats, ...)
  a <- mlg.filter(x, threshold, missing, algorithm = "a", distance = distmat, 
                  stats = stats, ...)
  n <- mlg.filter(x, threshold, missing, algorithm = "n", distance = distmat, 
                  stats = stats, ...)
  fanlist <- list(farthest = f, average = a, nearest = n)
  if (stats == "All"){
    for (i in names(fanlist)){
      names(fanlist[[i]]) <- c("MLG", "thresholds", "mat", "size")
    }
    if (plot){
      plot_filter_stats(x, fanlist, distmat, nclone)
    }
  }
  return(fanlist)
}

#' Predict cutoff thresholds for use with mlg.filter
#' 
#' Given a series of thresholds for a data set that collapse it into one giant 
#' cluster, this will search the top fraction of threshold differences to find 
#' the largest difference. The average between the thresholds spanning that 
#' difference is the cutoff threshold defining the clonal lineage threshold.
#' 
#' @param thresholds a vector of numerics coming from mlg.filter where the
#'   threshold has been set to the maximum threshold theoretically possible.
#' @param fraction the fraction of the data to seek the threshold.
#'   
#' @return a numeric value representing the threshold at which multilocus
#'   lineages should be defined.
#'   
#' @note This is a bit of a blunt instrument.
#'   
#' @author Zhian N. Kamvar
threshold_predictor <- function(thresholds, fraction = 0.5){
  frac <- 1:round(length(thresholds)*fraction)
  diffs <- diff(thresholds[frac])
  diffmax <- which.max(diffs)
  mean(thresholds[diffmax:(diffmax + 1)])
}

#' Plot the results of filter_stats
#' 
#' @param x a genlight of genind object
#' @param fstats the list passed from filter_stats
#' @param distmat a distance matrix passed from filter_stats
#' @param nclone see filter_stats
#'   
#' @return a plot depicting how many MLLs are collapsed as the genetic distance 
#'   increases for each algorithm.
#'   
#' @author Zhian N. Kamvar
#' @keyword internal
plot_filter_stats <- function(x, fstats, distmat, nclone = NULL){
  upper <- round(max(distmat), digits = 1)
  ylims <- c(ifelse(is.genind(x), mlg(x, quiet = TRUE), nInd(x)), 1)
  plot(x = c(upper, 0), y = ylims, type = "n",
       ylab = "Number of Multilocus Lineages",
       xlab = "Genetic Distance Cutoff")
  a <- fstats$average$thresholds
  n <- fstats$nearest$thresholds
  f <- fstats$farthest$thresholds
  plotcols <- RColorBrewer::brewer.pal(3, "Set1")
  names(plotcols) <- c("f", "a", "n")
  points(x = rev(a), y = 1:length(a), col = plotcols["a"])
  points(x = rev(f), y = 1:length(f), col = plotcols["f"])
  points(x = rev(n), y = 1:length(n), col = plotcols["n"])
  if (!is.null(nclone)){
    abline(v = a[1 + length(a) - nclone] + .Machine$double.eps^0.5, lty = 2,
           col = plotcols["a"])
    abline(v = f[1 + length(f) - nclone] + .Machine$double.eps^0.5, lty = 2, 
           col = plotcols["f"])
    abline(v = n[1 + length(n) - nclone] + .Machine$double.eps^0.5, lty = 2, 
           col = plotcols["n"])
    abline(h = nclone)
    text(upper, nclone, labels = paste0("n = ", nclone), 
         adj = c(1, -0.5))
    legend("topright", 
           legend = c("Nearest Neighbor", "UPGMA", "Farthest Neighbor"), 
           col = plotcols[c("n", "a", "f")], 
           pch = 1, 
           lty = 2,
           title = "Clustering Method")
  } else {
    legend("topright", 
           legend = c("Nearest Neighbor", "UPGMA", "Farthest Neighbor"), 
           col = plotcols[c("n", "a", "f")], 
           pch = 1, 
           title = "Clustering Method")    
  }
}

#' Function to color the tips of a dendrogram belonging to the same lineage
#' 
#' This function will color the tips of a dendrogram from ape based on whether 
#' or not the tips are defined to be in the same multilocus lineage. Red tips 
#' indicate that two tips are defined as the same multilocus lineage, blue tips 
#' indicate that they are separate.
#' 
#' @param x a genind or genlight object.
#' @param tree a tree derived from that object.
#' @param newmlgs a vector assigning the multilocus genotypes from x.
#' @param ... methods passed on to plot.phylo
#'   
#' @return nothing. A graphic.
#'   
#' @author Zhian N. Kamvar
color_mlg_tree <- function(x, tree, newmlgs, ...){
  mlg_vector <- newmlgs
  all_edges <- length(ape::which.edge(tree, tree$tip.label))
  edge_cols <- rep("black", all_edges)
  for (i in tree$tip.label){
    edge_cols[ape::which.edge(tree, i)] <- "blue"
  }
  edge_widths <- rep(1, all_edges)
  for (i in unique(mlg_vector)){
    indices <- ape::which.edge(tree, tree$tip.label[mlg_vector == i])
    if (length(indices) > 1){
      edge_cols[indices] <- "red"      
      edge_widths[indices] <- 3
    }
  }
  tree$tip.label <- paste(tree$tip.label, pop(x), sep = "_")
  ape::plot.phylo(tree, edge.color = edge_cols, edge.width = edge_widths,
                  adj = 0, label.offset = 0.005, ...)
}

#' Create an HTML animation of trees collapsing multilocus genotypes
#' 
#' This function utilizes saveHTML from the animation package to color
#' dendrograms as they are serially collapsed given increasing thresholds to
#' mlg.filter.
#' 
#' @param measure one of the algorithm options for mlg.filter
#' @param x a genind or genlight object.
#' @param treefun a function used to generate a tree.
#' @param distfun a function to calculate distance
#' @param fstats a list from filter_stats where stats = "All"
#' @param destdir a destination directory
#' @param ... parameters passed on to mlg.filter
#'   
#' @return an animation
#'   
#' @author Zhian N. Kamvar
HTML_collapse <- function(measure, x, treefun, distfun, fstats, destdir = NULL, ...){
  desc <- "An exploration of collapsing MLGs using"
  descname <- switch(measure, nearest = "Nearest Neighbor",
                     farthest = "Farthest Neighbor",
                     average = "UPGMA")
  desc <- paste(desc, descname)
  TREEFUN <- match.fun(treefun)
  DISTFUN <- match.fun(distfun)
  cwd <- getwd()
  if (!is.null(destdir)){
    setwd(destdir)
  }
  the_tree <- TREEFUN(DISTFUN(x, ...))
  isnj <- length(grep("nj", substitute(treefun))) > 0
  if (isnj){
    the_tree <- ape::ladderize(poppr:::fix_negative_branch(the_tree))
  }
  animation::saveHTML({
    for (i in unique(fstats[[measure]]$thresholds)){
      z <- mlg.filter(x, threshold = i, missing = "ignore", algorithm = measure, 
                      distance = distfun, stats = "MLGs", ...)
      if (isnj){
        color_mlg_tree(x, the_tree, z, type = "u", lab4ut = "axial")
      } else {
        color_mlg_tree(x, the_tree, z)
      }
      
      if (isnj){
        add.scale.bar()
      } else {
        axisPhylo(1)
      }
      title(paste("Threshold:", round(i, 3), "MLG:", length(unique(z))))
    }
  }, img.name = paste0(treefun, "_", measure), imgdir = paste0(measure), 
  htmlfile = paste0(treefun, "_", measure, ".html"), 
  autobrowse = FALSE, title = paste("mlg.filter option", measure), 
  description = desc)
  setwd(cwd)
}