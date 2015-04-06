
# Plotting average distance of various genomic features to the nearest
# breakpoint (edge of duplicated or triplicated block), and then comparing
# this to a distribution of average distances that might occur if we shuffled 
# the locations of the breakpoints 1,000 times.

# make large bottom outer margin
par(mar=c(0,2.5,1.5,2), oma=c(6,3,1,1))


# set out a layout of select features
m = c(1, 2, 17, 3,4, 16, 5,6, 15, 7,8, 14, 9,10, 13, 11,12, 18)
pp <- layout(matrix(m, 6,3, byrow=F), height = c(4, 4, 4, 4, 4, 3))
layout.show(pp)


# 2x gene data
real_data = read.table("~/Work/Code/FRAG_project/Results/real_distances_2x_gene.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
real_mean = mean(real_data$V1)
shuffled_data = read.table("~/Work/Code/FRAG_project/Results/shuffled_distances_2x_gene.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(shuffled_data$V1, xlim=c(0,10000), ylim=c(0,80), breaks=50, col="red", main="Genes", xaxt='n', xlab=NA, ylab=NA)
abline(v=real_mean, lwd=2, col="black")

# 3x genes data
real_data = read.table("~/Work/Code/FRAG_project/Results/real_distances_3x_gene.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
real_mean = mean(real_data$V1)
shuffled_data = read.table("~/Work/Code/FRAG_project/Results/shuffled_distances_3x_gene.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(shuffled_data$V1, xlim=c(0,10000), ylim=c(0,80), breaks=50, col="blue", main=NA, xlab=NA, ylab=NA)
abline(v=real_mean, lwd=2, col="black")



# 2x state1 data
real_data = read.table("~/Work/Code/FRAG_project/Results/real_distances_2x_open_chromatin_state_1.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
real_mean = mean(real_data$V1)
shuffled_data = read.table("~/Work/Code/FRAG_project/Results/shuffled_distances_2x_open_chromatin_state_1.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(shuffled_data$V1, xlim=c(0,75000), ylim=c(0,80), breaks=50, col="red", main="Chromatin state 1", xaxt='n', xlab=NA, ylab=NA)
abline(v=real_mean, lwd=2, col="black")

# 3x state1 data
real_data = read.table("~/Work/Code/FRAG_project/Results/real_distances_3x_open_chromatin_state_1.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
real_mean = mean(real_data$V1)
shuffled_data = read.table("~/Work/Code/FRAG_project/Results/shuffled_distances_3x_open_chromatin_state_1.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(shuffled_data$V1, xlim=c(0,75000), ylim=c(0,80), breaks=50, col="blue", main=NA, xlab=NA, ylab=NA)
abline(v=real_mean, lwd=2, col="black")



# 2x state3 data
real_data = read.table("~/Work/Code/FRAG_project/Results/real_distances_2x_open_chromatin_state_3.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
real_mean = mean(real_data$V1)
shuffled_data = read.table("~/Work/Code/FRAG_project/Results/shuffled_distances_2x_open_chromatin_state_3.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(shuffled_data$V1, xlim=c(0,75000), ylim=c(0,80), breaks=50, col="red", main="Chromatin state 3", xaxt='n', xlab=NA, ylab=NA)
abline(v=real_mean, lwd=2, col="black")

# 3x state3 data
real_data = read.table("~/Work/Code/FRAG_project/Results/real_distances_3x_open_chromatin_state_3.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
real_mean = mean(real_data$V1)
shuffled_data = read.table("~/Work/Code/FRAG_project/Results/shuffled_distances_3x_open_chromatin_state_3.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(shuffled_data$V1, xlim=c(0,75000), ylim=c(0,80), breaks=50, col="blue", main=NA, xlab=NA, ylab=NA)
abline(v=real_mean, lwd=2, col="black")


# 2x state6 data
real_data = read.table("~/Work/Code/FRAG_project/Results/real_distances_2x_open_chromatin_state_6.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
real_mean = mean(real_data$V1)
shuffled_data = read.table("~/Work/Code/FRAG_project/Results/shuffled_distances_2x_open_chromatin_state_6.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(shuffled_data$V1, xlim=c(0,35000), ylim=c(0,80), breaks=50, col="red", main="Chromatin state 6", xaxt='n', xlab=NA, ylab=NA)
abline(v=real_mean, lwd=2, col="black")

# 3x state6 data
real_data = read.table("~/Work/Code/FRAG_project/Results/real_distances_3x_open_chromatin_state_6.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
real_mean = mean(real_data$V1)
shuffled_data = read.table("~/Work/Code/FRAG_project/Results/shuffled_distances_3x_open_chromatin_state_6.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(shuffled_data$V1, xlim=c(0,35000), ylim=c(0,80), breaks=50, col="blue", main=NA, xlab=NA, ylab=NA)
abline(v=real_mean, lwd=2, col="black")


# 2x transposable element data
real_data = read.table("~/Work/Code/FRAG_project/Results/real_distances_2x_transposable_element.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
real_mean = mean(real_data$V1)
shuffled_data = read.table("~/Work/Code/FRAG_project/Results/shuffled_distances_2x_transposable_element.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(shuffled_data$V1, xlim=c(0,12000), ylim=c(0,80), breaks=50, col="red", main="Transposable elements", xaxt='n', xlab=NA, ylab=NA)
abline(v=real_mean, lwd=2, col="black")

# 3x transposable elements data
real_data = read.table("~/Work/Code/FRAG_project/Results/real_distances_3x_transposable_element.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
real_mean = mean(real_data$V1)
shuffled_data = read.table("~/Work/Code/FRAG_project/Results/shuffled_distances_3x_tr", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(shuffled_data$V1, xlim=c(0,12000), ylim=c(0,80), breaks=50, col="blue", main=NA, xlab=NA, ylab=NA)
abline(v=real_mean, lwd=2, col="black")



# 2x origin data
real_data = read.table("~/Work/Code/FRAG_project/Results/real_distances_2x_replication_origins.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
real_mean = mean(real_data$V1)
shuffled_data = read.table("~/Work/Code/FRAG_project/Results/shuffled_distances_2x_replication_origins.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(shuffled_data$V1, xlim=c(0,70000), ylim=c(0,80), breaks=50, col="red", main="Replication origins", xaxt='n', xlab=NA, ylab=NA)
abline(v=real_mean, lwd=2, col="black")

# 3x origin data
real_data = read.table("~/Work/Code/FRAG_project/Results/real_distances_3x_replication_origins.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
real_mean = mean(real_data$V1)
shuffled_data = read.table("~/Work/Code/FRAG_project/Results/shuffled_distances_3x_replication_origins.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(shuffled_data$V1, xlim=c(0,70000), ylim=c(0,80), breaks=50, col="blue", main=NA, xlab=NA, ylab=NA)
abline(v=real_mean, lwd=2, col="black")

# Add series labels
mtext("Average distance of feature to nearest breakpoint", side=1, outer=T, at=0.5, cex=1.2)
mtext("Frequency", side=2, outer=T, at=0.5, cex=1.2)



par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("Duplicated shuffle data", "Triplicated shuffle data", "Real average distances"), xpd = TRUE, horiz = TRUE, inset = c(0, 0), bty = "n",  col = c("red", "blue", "black"), cex=1.2, lty = c(1, 1), lwd=3)



# xpd = TRUE tells R that it is OK to plot outside the region 
# horiz = TRUE tells R that I want a horizontal legend 
# inset = c(x,y) tells R how to move the legend relative to the 'bottom' location 
# bty = 'n' means that 'no' box will be drawn around it 
# pch and col are the types and colors of points 


