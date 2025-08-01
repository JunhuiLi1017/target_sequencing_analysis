library(getopt)
library(ggplot2)
#+--------------------
# get options
#+--------------------
spec <- matrix(c(
  'help', 'h', 0, "logical", "help",
  'verbose', 'v', 2, "integer", "verbose mode, default [1]",
  'input', 'i', 1, "character", "variants of mutect2, forced.",
  'outfile','o',1,"character","names of output file,forced"
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

#+--------------------
# check options
#+--------------------
if ( !is.null(opt$help) | is.null(opt$input) | is.null(opt$outfile) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}
a <- read.table(opt$input,sep="\t",header=F)
colnames(a) <- c("caller","id","depth","af")
a$caller <- as.factor(a$caller)
a$af <- as.numeric(a$af)
a$depth <- as.numeric(a$depth)

p <- ggplot(data=a,aes(x=depth,y=af)) + 
  geom_point(aes(color=caller))

p <- p + scale_x_log10() + theme_bw()

png(file = opt$outfile, width = 1000, height = 1000)
print(p)
dev.off()
