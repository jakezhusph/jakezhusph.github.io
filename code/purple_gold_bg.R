rm(list=ls())
library(ggplot2) 
library(grid)
library(RColorBrewer)

make_gradient <- function(deg = 45, n = 100, cols = blues9) {
  cols <- colorRampPalette(cols)(n + 1)
  rad <- deg / (180 / pi)
  mat <- matrix(
    data = rep(seq(0, 1, length.out = n) * cos(rad), n),
    byrow = TRUE,
    ncol = n
  ) +
  matrix(
    data = rep(seq(0, 1, length.out = n) * sin(rad), n),
    byrow = FALSE,
    ncol = n
  )
  mat <- mat - min(mat)
  mat <- mat / max(mat)
  mat <- 1 + mat * n
  mat <- matrix(data = cols[round(mat)], ncol = n)
  grid::rasterGrob(
    image = mat,
    width = unit(1, "npc"),
    height = unit(1, "npc"), 
    interpolate = TRUE
  )
}


g <- make_gradient(
  deg = 0, n = 500, cols = c("#552583","#FDB927")
  # deg = 0, n = 500, cols = c("#FDB927","#552583","#FDB927","#552583","#FDB927")
)


png(file="/Users/jiaqiangzhu/Documents/Github/jakezhusph.github.io/images/purple_gold_gradient.png",width=2160,height=3840)
ggplot() +
  annotation_custom(
    grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) 
  dev.off()





