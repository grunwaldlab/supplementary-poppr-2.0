#' Demonstrating Graph Walking algorithms in R
#' =============
#' 
#' In version 2.0, poppr introduces the option to include reticulation in 
#' minimum spanning networks. This is essential for clonal organisms as mutation
#' has a much higher influence on genetic diversity as opposed to meiotic 
#' processes. 
#' 
#' In this examle, I will utilize *P. ramorum* data from Kamvar et al. (2015). The source
#' should be downloaded with the following command:
#+ eval = FALSE
devtools::install_github("zkamvar/Sudden_Oak_Death_in_Oregon_Forests/PramCurry")
#' 
#' After you've downloaded the package, you can load the library and the data.
library("PramCurry")
library("poppr")
library("igraph")
library("ggplot2")
library("reshape2")
data(for2nur)
data(comparePal)
#' 
#' Below are some functions we need to convert the data to the proper format as
#' well as some helper functions from the PramCurry package.
# Conversion from adegenet 1.4 to 2.0
data_replacer <- function(x){
  # manipulate your data here
  if (x@type == "codom"){
    xtab <- x@tab
    newtab <- as.integer(x@tab * x@ploidy)
    x@tab <- matrix(newtab, nrow = nrow(xtab), ncol = ncol(xtab),
                    dimnames = dimnames(xtab))
    x@ploidy <- rep(x@ploidy, nrow(x@tab))
    names(x@loc.names) <- NULL
    rownames(x@tab) <- x@ind.names
    names(x@all.names) <- x@loc.names
    colnames(x@tab) <- unlist(lapply(x@loc.names, 
                                     function(i) paste(i, x@all.names[[i]], 
                                                       sep = ".")), 
                              use.names = FALSE)
    levels(x@loc.fac) <- x@loc.names
    names(x@loc.nall) <- x@loc.names
  }
  
  if (!is.null(x@pop)){
    levels(x@pop) <- x@pop.names
  }
  if ("genclone" %in% class(x)){
    x@strata <- x@hierarchy
    x@hierarchy <- NULL
  } else {
    x@hierarchy <- NULL
    x@strata    <- NULL
  }
  return(x)
}

# These functions are from Kamvar et al. 2015
make_node_list <- function(clusts, pal, cutoff = 3){
  PAL <- match.fun(pal)
  clustPal  <- table(clusts$membership)
  theClusts <- clustPal > 3
  clustOrder <- order(clustPal, decreasing = TRUE)
  clustPal[clustOrder][theClusts[clustOrder]]   <- PAL(sum(theClusts))
  clustPal[clustOrder][!theClusts[clustOrder]]  <- gray.colors(sum(!theClusts))
  nodeList <- lapply(1:length(clustPal), function(x) which(clusts$membership == x))
  names(nodeList) <- clustPal
  return(nodeList)
}

clust_cutoff <- function(clusts, cutoff = 3){
  table(clusts$membership) > 3
}

match_communities <- function(dat, comm, the_graph){
  members <- membership(comm)
  mlgs <- mlg.vector(dat)
  graph_mlgs <- PramCurry::mlgFromString(V(the_graph)$label)
  members[match(mlgs, graph_mlgs)]
}

#' 
#' ### Data setup
#' Conversion
for2nur <- data_replacer(for2nur)
setPop(for2nur) <- ~SOURCE/STATE
for2nur
#' 
#' Repeat lengths are necessary:
newReps <- other(for2nur)$REPLEN
(newReps[3] <- 4) # Tetranucleotide repeat
(newReps <- fix_replen(for2nur, newReps))
#'
#' ### Calculating the minimum spanning networks
#' 
#' We will create two msns: one with no reticulations (`nf.msn.noties`) and
#' one with reticulations (`nf.msn.ties`). these will be analyzed with the
#' infomap community detection algorithm and then the memberships will be compared
#' to the observed populations. 
# No reticulations in the network.
nf.msn.noties <- bruvo.msn(for2nur, replen = newReps, showplot = FALSE)
# Adding reticulations
nf.msn.ties <- bruvo.msn(for2nur, replen = newReps, include.ties = TRUE,
                         showplot = FALSE)
#' 
#' Now, we are going to use the igraph function `infomap.community`. Basically,
#' it will send random walkers through the graph and see where they spend the
#' most time. It considers edge weights and vertex weights. For vertex weights,
#' we are setting them to be the number of samples in each vertex. 
#+ community_detecting, cache = TRUE
set.seed(20150427)
clusts.noties <- infomap.community(nf.msn.noties$graph, nb.trials = 1e4, 
                                   v.weights = V(nf.msn.noties$graph)$size)
set.seed(20150427)
clusts.ties <- infomap.community(nf.msn.ties$graph, nb.trials = 1e4, 
                                 v.weights = V(nf.msn.ties$graph)$size)
#'
#' The membership of each sample is matched here. 
members.ties <- match_communities(for2nur, clusts.ties, nf.msn.ties$graph)
members.noties <- match_communities(for2nur, clusts.noties, nf.msn.noties$graph)
#'
#' Now, we are comparing the memberships vs. the populations
ties_table   <- table(members.ties, pop(for2nur))
noties_table <- table(members.noties, pop(for2nur))
plot(sizes(clusts.ties), main = "with reticulation",
     xlab = "Community", ylab = "Size")
plot(sizes(clusts.noties), main = "without reticulation",
     xlab = "Community", ylab = "Size")
table.value(ties_table, col.labels = popNames(for2nur))
table.value(noties_table, col.labels = popNames(for2nur))
#' 
#' We can utilize the Shannon index to investigate entropy.
#' 
#+ entropy_graph_ramorum
ent_table <- data.frame(`No Reticulation` = vegan::diversity(noties_table, MARGIN = 2),
                        `With Reticulation` = vegan::diversity(ties_table, MARGIN = 2),
                        population = popNames(for2nur), check.names = FALSE)
ggplot(melt(ent_table), aes(x = population, y = value, color = variable)) +
  geom_point(position = "identity", size = 5) +
  ggtitle("P. ramorum") + 
  ylab("entropy") + 
  PramCurry:::myTheme +
  theme(legend.title = element_blank(), plot.title = element_text(face = "italic"))
#'
#' One other thing to look at is the sizes of the communities. With a clonal
#' pathogen, we expect to see larger communities. 
sizes(clusts.ties)
sizes(clusts.noties)
vegan::diversity(sizes(clusts.ties))
vegan::diversity(sizes(clusts.noties))
#'
#' What we are seeing is that the reticulations allow for more genotypes to be
#' clustered together because there are more meaningful connections. 
#' 
#' Below, we will see what the graphs look like with and without ties.
#' ### with ties
#+ ties, fig.width = 10, fig.height = 10
nf.msn <- nf.msn.ties
clusts <- clusts.ties
mypal <- ifelse(length(unique(clusts$membership)) < 14, 
                RColorBrewer::brewer.pal(x, "Paired"), 
                colorRampPalette("black"))
nodeList <- make_node_list(clusts, mypal, cutoff = 3)
theClusts <- clust_cutoff(clusts, 3)

goodSeed <- 6
thisPal <- function(x) comparePal[nf.msn$populations]

set.seed(goodSeed)
MASTER <- get_layout(nf.msn$graph, LAYOUT = layout.fruchterman.reingold)
plot_poppr_msn(for2nur, nf.msn, gad = 10, palette = thisPal, mlg = TRUE, 
               layfun = MASTER, nodebase = 1.75, vertex.label.font = 2,
               quantiles = FALSE, #inds = "none",
               vertex.label.color = "firebrick")#,
#                mark.groups = nodeList[theClusts], 
#                mark.border = names(nodeList)[theClusts],
#                mark.col = transp(names(nodeList)[theClusts], 0.05),
#                mark.expand = 2,
#                mark.shape = 0)
plot(clusts, nf.msn$graph, vertex.size = log(V(nf.msn$graph)$size, 1.75) + 3, 
     vertex.label = NA, layout = MASTER)

#' ### without ties
#+ noties, fig.width = 10, fig.height = 10
nf.msn <- nf.msn.noties
clusts <- clusts.noties
mypal <- ifelse(length(unique(clusts$membership)) < 14, 
                RColorBrewer::brewer.pal(x, "Paired"), 
                colorRampPalette("black"))
nodeList <- make_node_list(clusts, mypal, cutoff = 3)
theClusts <- clust_cutoff(clusts, 3)

goodSeed <- 6
thisPal <- function(x) comparePal[nf.msn$populations]
set.seed(goodSeed)
MASTER <- get_layout(nf.msn$graph, LAYOUT = layout.auto)
plot_poppr_msn(for2nur, nf.msn, gad = 10, palette = thisPal, mlg = TRUE, 
               layfun = MASTER, nodebase = 1.75, vertex.label.font = 2,
               quantiles = FALSE, #inds = "none",
               vertex.label.color = "firebrick")#,
#                mark.groups = nodeList[theClusts], 
#                mark.border = names(nodeList)[theClusts],
#                mark.col = transp(names(nodeList)[theClusts], 0.05),
#                mark.expand = 2,
#                mark.shape = 0)
plot(clusts, nf.msn$graph, vertex.size = log(V(nf.msn$graph)$size, 1.75) + 3, 
     vertex.label = NA, layout = MASTER)
#' 
#' H3N2 virus data 
#' ===============
#' 
#' For this example, we will use the H3N2 viral data set from the adegenet
#' package.
data("H3N2")
hc <- as.genclone(H3N2) # For faster calculation of MLGs
pop(hc) <- other(H3N2)$epid
virusPal <- colorRampPalette(RColorBrewer::brewer.pal(n = nPop(hc), "Set1"))
#' 
#' First the distance needs to be calculated. For this, we will simply use a 
#' dissimilarity distance.
#+ H3N2_distance, cache = TRUE
hdist <- diss.dist(hc, percent = TRUE)
#' 
#' Now we will calculate the two minimum spanning networks and walk the graphs.
#+ H3N2_graph, cache = TRUE
hmsn.noties <- poppr.msn(hc, distmat = hdist, showplot = FALSE, palette = virusPal)
hmsn.ties <- poppr.msn(hc, distmat = hdist, showplot = FALSE, palette = virusPal, 
                         include.ties = TRUE)
set.seed(20150428)
hmsn.communities.noties <- infomap.community(hmsn.noties$graph, 
                                             v.weights = V(hmsn.noties$graph), 
                                             nb.trials = 1e4)
set.seed(20150428)
hmsn.communities.ties <- infomap.community(hmsn.ties$graph, 
                                           v.weights = V(hmsn.ties$graph), 
                                           nb.trials = 1e4)
#' 
#' Let's compare these, shall we? First, let's look at the community sizes.
#' Since there will be a lot of them, we can visualize them as barplots:
plot(sizes(hmsn.communities.ties), main = "with reticulation",
     xlab = "Community", ylab = "Size")
vegan::diversity(sizes(hmsn.communities.ties))
max(sizes(hmsn.communities.ties))
plot(sizes(hmsn.communities.noties), main = "without reticulation",
     xlab = "Community", ylab = "Size")
vegan::diversity(sizes(hmsn.communities.noties))
max(sizes(hmsn.communities.noties))
#'
#' Now we can take a look at the contingency tables:
members.ties <- match_communities(hc, hmsn.communities.ties, hmsn.ties$graph)
ties.table <- table(members.ties, pop(hc))
table.value(ties.table, col.labels = popNames(hc))

members.noties <- match_communities(hc, hmsn.communities.noties, hmsn.noties$graph)
noties.table <- table(members.noties, pop(hc))
table.value(noties.table, col.labels = popNames(hc))
#'
#' A markded difference! Now onto the graphs!
#' 
#' We can examine the information contained within these tables:
#' 
#+ entropy_graph_virus
ent_table <- data.frame(`No Reticulation` = vegan::diversity(noties.table, MARGIN = 2),
                        `With Reticulation` = vegan::diversity(ties.table, MARGIN = 2),
                        year = popNames(hc), check.names = FALSE)
ent_plot <- ggplot(melt(ent_table), 
                   aes(x = year, y = value, shape = variable)) +
  geom_point(position = "identity") +
  ylab("entropy (H)") + 
  theme_linedraw() +
  ylim(0, NA) + # Set zero as lower limit
  scale_shape_manual(values = c(3, 5)) +  # shape reflects method
  theme(legend.title = element_blank(), legend.position = "top",
        panel.grid = element_blank())
ent_plot
ggsave(filename = "../main_article/poppr_frontiers_files/custom_figures/entropy.pdf",plot = ent_plot, width = 80, height = 50, units = "mm", scale = 1.5)
#' 
#+ H3N2_msns, cache = TRUE, fig.width = 10, fig.height = 10
set.seed(20150427)
HLAYOUT <- get_layout(hmsn.ties$graph)
plot_poppr_msn(hc, hmsn.ties, layout = HLAYOUT, inds = "none", nodebase = 1.9,
               gadj = 25, quantiles = FALSE, nodelab = 1e3, wscale = FALSE)

pdf("../main_article/poppr_frontiers_files/custom_figures/H3N2-ties.pdf")
plot_poppr_msn(hc, hmsn.ties, layout = HLAYOUT, inds = "none", nodebase = 1.9,
               gadj = 25, quantiles = FALSE, nodelab = 1e3, wscale = FALSE)
dev.off()

plot(hmsn.communities.ties, hmsn.ties$graph, 
     vertex.size = log(V(hmsn.ties$graph)$size, 1.9) + 3, 
     layout = HLAYOUT, vertex.label = NA)

set.seed(20150427)
HLAYOUT <- get_layout(hmsn.noties$graph)
plot_poppr_msn(hc, hmsn.noties, layout = HLAYOUT, inds = "none", nodebase = 1.9,
               gadj = 25, quantiles = FALSE, nodelab = 1e3, wscale = FALSE)

pdf("../main_article/poppr_frontiers_files/custom_figures/H3N2-noties.pdf")
plot_poppr_msn(hc, hmsn.noties, layout = HLAYOUT, inds = "none", nodebase = 1.9,
               gadj = 25, quantiles = FALSE, nodelab = 1e3, wscale = FALSE)
dev.off()

plot(hmsn.communities.noties, hmsn.noties$graph, 
     vertex.size = log(V(hmsn.noties$graph)$size, 1.9) + 3, 
     layout = HLAYOUT, vertex.label = NA)
#'
#' Now to show the dapc analysis for comparison.
#+ H3N2 dapc
dapc.H3N2 <- dapc(hc, var.contrib=TRUE, scale=FALSE, n.pca=30, n.da=10)
scatter(dapc.H3N2, pch=19, cstar=0, mstree=TRUE, lwd=2, lty=2, clabel = FALSE,
        col = virusPal(6), scree.da = FALSE)
points(dapc.H3N2$ind.coord[, 1:2], pch = 21, bg = NA)
legend("topleft", legend = popNames(hc), pt.bg = virusPal(6), pch = 21)
pdf("../main_article/poppr_frontiers_files/custom_figures/H3N2-scatter.pdf", 
    width = 4*3.14961,
    height = 3*3.14961)
scatter(dapc.H3N2, pch=19, cstar=0, mstree=TRUE, lwd=2, lty=2, clabel = FALSE,
        col = virusPal(6), scree.da = FALSE)
points(dapc.H3N2$ind.coord[, 1:2], pch = 21, bg = NA)
legend("topleft", legend = popNames(hc), pt.bg = virusPal(6), pch = 21)
dev.off()
#'
#' ## Session Info
options(width = 100)
devtools::session_info()