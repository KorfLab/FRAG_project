#################
## SET UP PLOT ##
#################

# change to coordinates instead of scientific notation
options(scipen=10)

# change style of point that is used
par(pch=20)

# set up to make two plots, one above the other
#par(mfrow=c(2,1))

# set up a matrix to describe graph layout and then send to layout function
# 0 values can be used for non-plotting, for when you want a blank region
#m <- matrix(c(1, 2), nrow = 2, ncol = 1)
#layout(m, heights=c(6, 2))

x_min <-  1
x_max <-  27700001

x_min <-  1
x_max <-  8000000; # all good through 1-8 Mbp

x_min <-  8000000; # possible small missing 3x block at 8.3 Mbp?
x_max <- 10000000; 

x_min <- 12000000; # possible erroneous edge at 11.5 Mbp
x_max <- 14000000; 

x_min <- 14000000; # possible problems with too much variance
x_max <- 16000000; 

x_min <- 16000000; # possible problems with too much variance
x_max <- 18000000; 

x_min <- 18000000; 
x_max <- 20000000; 

x_min <- 20000000; 
x_max <- 22000000; 

x_min <- 22000000; 
x_max <- 24000000; 

x_min <- 24000000; # possible problems with too much variance
x_max <- 26000000; 

x_min <- 26000000; 
x_max <- 28000000; 

x_min <- 18000000; # possible missing block
x_max <- 18600000; 


x_min <- 16200000; # possible problems with too much variance
x_max <- 17200000; 

x_min <- 7600000;
x_max <- 10000000;

x_min <-  11000000;
x_max <-  13000000; # all good through 1-8 Mbp

x_min <-  12500000;
x_max <-  14500000; # all good through 1-8 Mbp

################
## FIRST PLOT ##
################

# read data
frag_data = read.table("~/Work/Chan_lab/CGRs_1kb_plots_forkeith.tsv", header=TRUE, na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
frag_chr1 = subset(frag_data, frag_data$Chrom == "Chr1")

# set margins
#par(mar=c(2, 5, 2, 5))


# make main plot, type 'l' for line connections only (b for points + lines)
plot(x=frag_chr1$Strt, y=frag_chr1$FRAG00062.SuperF1diploid, type="l", ylim=c(1,5), xlab=NA, ylab="Copy", col="magenta", xlim=c(x_min,x_max), axes=FALSE)

# now overlay points in different color (default = black)
points(x=frag_chr1$Strt, y=frag_chr1$FRAG00062.SuperF1diploid)

box()
axis(side=1, at=seq(x_min, x_max, by=100000))
axis(side=2, at=c(0:5))


#################
## SECOND PLOT ##
#################

# read data
t_test_data = read.table("~/Work/Chan_lab/T-test_results/test20.tsv", header=TRUE, na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
t_test_data = read.table("~/Work/Chan_lab/T-test_results/test30.tsv", header=TRUE, na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
t_test_data = read.table("~/Work/Chan_lab/T-test_results/test40.tsv", header=TRUE, na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
t_test_data = read.table("~/Work/Chan_lab/T-test_results/test50.tsv", header=TRUE, na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)

chr1_1x_block = subset(t_test_data, t_test_data$Chr == "Chr1" & t_test_data$Copy_number == "1x" & t_test_data$Type == "block")
chr1_2x_block = subset(t_test_data, t_test_data$Chr == "Chr1" & t_test_data$Copy_number == "2x" & t_test_data$Type == "block")
chr1_3x_block = subset(t_test_data, t_test_data$Chr == "Chr1" & t_test_data$Copy_number == "3x" & t_test_data$Type == "block")

chr1_1x_L = subset(t_test_data, t_test_data$Chr == "Chr1" & t_test_data$Copy_number == "1x" & t_test_data$Type == "edge-L")
chr1_2x_L = subset(t_test_data, t_test_data$Chr == "Chr1" & t_test_data$Copy_number == "2x" & t_test_data$Type == "edge-L")
chr1_3x_L = subset(t_test_data, t_test_data$Chr == "Chr1" & t_test_data$Copy_number == "3x" & t_test_data$Type == "edge-L")

chr1_1x_R = subset(t_test_data, t_test_data$Chr == "Chr1" & t_test_data$Copy_number == "1x" & t_test_data$Type == "edge-R")
chr1_2x_R = subset(t_test_data, t_test_data$Chr == "Chr1" & t_test_data$Copy_number == "2x" & t_test_data$Type == "edge-R")
chr1_3x_R = subset(t_test_data, t_test_data$Chr == "Chr1" & t_test_data$Copy_number == "3x" & t_test_data$Type == "edge-R")

#block_chr1_3x = subset(block_data, block_data$V1 =="Chr1" & block_data$V6 == "3x")
#block_chr1_3x = subset(block_data, block_data$V1 =="Chr1" & block_data$V6 %in% "cat")

# set margins
#par(mar=c(5, 5, 1, 5))

# draw rectangle over plot
#plot(x=chr1$Strt, y=frag_chr1$FRAG00062.SuperF1diploid, type="n", ylim=c(0,3), xlab="Chromosome 1", ylab="Blocks",yaxt='n',xlim=c(x_min,x_max))
#rect(xleft=block_chr1_3x$V4, xright=block_chr1_3x$V5, ybottom=1.5, ytop=3, col="blue", border="black")


# first plot for edges
rect(xleft=chr1_1x_L$Start_coord, xright=chr1_1x_L$End_coord, ybottom=1.55, ytop=1.8, col="lightblue", border="black")
rect(xleft=chr1_2x_L$Start_coord, xright=chr1_2x_L$End_coord, ybottom=2.4, ytop=2.65, col="lightblue", border="black")
rect(xleft=chr1_3x_L$Start_coord, xright=chr1_3x_L$End_coord, ybottom=3.35, ytop=3.6, col="lightblue", border="black")

rect(xleft=chr1_1x_R$Start_coord, xright=chr1_1x_R$End_coord, ybottom=1.55, ytop=1.8, col="red", border="black")
rect(xleft=chr1_2x_R$Start_coord, xright=chr1_2x_R$End_coord, ybottom=2.4, ytop=2.65, col="red", border="black")
rect(xleft=chr1_3x_R$Start_coord, xright=chr1_3x_R$End_coord, ybottom=3.35, ytop=3.6, col="red", border="black")

# now plot for whole blocks
rect(xleft=chr1_1x_block$Start_coord, xright=chr1_1x_block$End_coord, ybottom=1.62, ytop=1.73, col="yellow", border="black")
rect(xleft=chr1_2x_block$Start_coord, xright=chr1_2x_block$End_coord, ybottom=2.47, ytop=2.58, col="yellow", border="black")
rect(xleft=chr1_3x_block$Start_coord, xright=chr1_3x_block$End_coord, ybottom=3.42, ytop=3.53, col="yellow", border="black")
