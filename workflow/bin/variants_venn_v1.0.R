suppressPackageStartupMessages(library(getopt))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(scales)) # for alpha()
#+--------------------
# get options
#+--------------------
spec <- matrix(c(
  'help', 'h', 0, "logical", "help",
  'verbose', 'v', 2, "integer", "verbose mode, default [1]",
  'input', 'i', 1, "character", "variants file, forced.",
  'outfile','o',1,"character","names of output file,forced"
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

#+--------------------
# check options
#+--------------------
if ( !is.null(opt$help) | is.null(opt$outfile) | is.null(opt$input) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

dat <- read.table(opt$input,sep="\t")
caller_var_num <- table(dat[,1])
category_names <- paste0(names(caller_var_num),":",caller_var_num)

dat_list <- split(dat[,2],dat[,1])
names(dat_list) <- category_names

caller_num <- length(names(caller_var_num))
venn_colors <- viridis::viridis(caller_num)

if (caller_num == 2) {
  cat_pos <- c(0, 180)
  cat_dist <- c(0.05, 0.05)
} else if (caller_num == 3) {
  cat_pos <- c(-45, 45, 180)
  cat_dist <- c(0.06, 0.06, 0.07)
} else if (caller_num == 4) {
  cat_pos <- c(-45, 45, 135, -135)
  cat_dist <- c(0.07, 0.07, 0.07, 0.07)
} else {
  # Fallback for 1 or >4 callers (rare)
  cat_pos <- rep(0, caller_num)
  cat_dist <- rep(0.05, caller_num)
}

a <- venn.diagram(
  x = dat_list,
  category.names = category_names,
  filename = opt$outfile,
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col = venn_colors,
  fill = alpha(venn_colors, 0.3),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = cat_pos,
  cat.dist = cat_dist,
  cat.fontfamily = "sans",
  cat.col = venn_colors,
  rotation = 1
)

