#' Figure 1
#' ==========
#' 
#' Short script to demonstrate the principle used for figure 1.
#' 
#' ## Constuction
library("igraph")
mat <- combn(1:5, 2)
g <- graph(mat, directed = FALSE)
my_layout <- matrix(c(0, 0, 4, 8, 8,
                      3, 0, 0, 0, 3), 
                    ncol = 2, dimnames = list(LETTERS[1:5], c("x", "y")))
V(g)$label <- LETTERS[1:5]
E(g)$weights <- as.vector(dist(my_layout))/10

plot.igraph(g, layout = my_layout, edge.label = E(g)$weights)
my_graph <- get.data.frame(g)

my_layout <- as.data.frame(list(my_layout, names = LETTERS[1:5]))
my_graph$fromx <- my_layout$x[my_graph$from]
my_graph$fromy <- my_layout$y[my_graph$from]
my_graph$tox <- my_layout$x[my_graph$to]
my_graph$toy <- my_layout$y[my_graph$to]

library("ggplot2")

graph_plot <- ggplot(my_layout) +
  geom_segment(aes(x = fromx, y = fromy, xend = tox, yend = toy, size = 1- weights), 
               data = my_graph) +
  geom_point(aes(x = x, y = y), size = 10, shape = 21, fill = "white") +
  geom_text(aes(x = x, y = y, label = names)) +
  theme(line = element_blank(),
        rect = element_blank(),
        text = element_blank(),
        title = element_blank(),
        legend.position = "none") +
  xlim(c(-0.5, 8.5)) +
  ylim(c(-0.5, 3.5))
#' 
#' Figure 1 base. The circles were done in inkscape.
#+ fig.width = 4.5, fig.height = 2.5
graph_plot
# ggsave(filename = "example_graph.pdf", plot = graph_plot, width = 4.5, height = 2.5)
#'
#' ## Explanation with data
#' 
#' Here, we'll construct an imaginary data set with the distances presented in 
#' the diagram above.
library("poppr")
dat <- as.genclone(df2genind(my_layout["names"], ploidy = 1))
pop(dat) <- indNames(dat)
(dat.dist <- dist(my_layout[1:2])/10)
#'
#' Showing the minimum spanning network.
p <- poppr.msn(dat, dat.dist, pal = rainbow, showplot = FALSE)
set.seed(99)
plot(p$graph, vertex.size = V(p$graph)$size * 5, vertex.label.dist = 0.5)
legend("top", horiz = TRUE, legend = p$population, fill = p$color, bty = "n", cex = 0.75)
mll(dat)
#'
#' ### Farthest neighbor
p <- poppr.msn(dat, dat.dist, pal = rainbow, threshold = 0.451, showplot = FALSE,
               clustering.algorithm = "f")
set.seed(99)
plot(p$graph, vertex.size = V(p$graph)$size * 5, vertex.label.dist = 0.5)
legend("top", horiz = TRUE, legend = p$population, fill = p$color, bty = "n", cex = 0.75)
mlg.filter(dat, threshold = 0.451, distance = dat.dist, algorithm = "f")
#'
#' ### Average neighbor
p <- poppr.msn(dat, dat.dist, pal = rainbow, threshold = 0.451, showplot = FALSE,
               clustering.algorithm = "a")
set.seed(99)
plot(p$graph, vertex.size = V(p$graph)$size * 5, vertex.label.dist = 0.5)
legend("top", horiz = TRUE, legend = p$population, fill = p$color, bty = "n", cex = 0.75)
mlg.filter(dat, threshold = 0.451, distance = dat.dist, algorithm = "a")
#'
#' ### Nearest neighbor
p <- poppr.msn(dat, dat.dist, pal = rainbow, threshold = 0.451, showplot = FALSE,
               clustering.algorithm = "n")
set.seed(99)
plot(p$graph, vertex.size = V(p$graph)$size * 5, vertex.label.dist = 0.5)
legend("top", horiz = TRUE, legend = p$population, fill = p$color, bty = "n", cex = 0.75)
mlg.filter(dat, threshold = 0.451, distance = dat.dist, algorithm = "n")
#'
#' ## Session Info
options(width = 100)
devtools::session_info()