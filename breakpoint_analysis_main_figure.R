# Figure for Han's paper using data from FRAG00062 line
# plotting variation in frequency of various genomic features
# as we move through a breakpoint region (between duplicated and unduplicated block
# or between duplicated and triplicated blocks)
# these all represent flipped datasets (correcting to reflect inside & outside blocks)

# make large bottom outer margin
par(oma = c(7, 3, 1, 1))

# set layout of a 2x2 plot (and reduce margins of indvididual plots)
par(mfrow=c(2,1), mai=c(0.7, 0.8, 0.7, 0.2))


##############################
# Plot 1: Genes
##############################
genes_2x = read.table("Results/FRAG00062_2x_bias_data_genes_500_50_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
genes_3x = read.table("Results/FRAG00062_3x_bias_data_genes_500_50_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)


plot(x=genes_2x$V3,  y=genes_2x$V7, type="l",  xlim=c(-3000, 3000), ylim=c(45,70), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
plot(x=genes_3x$V3,  y=genes_3x$V7, type="l",  xlim=c(-3000, 3000), ylim=c(45,70), xlab=NA, lwd = 3, ylab=NA, col="blue",  axes=F)

abline(v=0, col=1, lty=3)
abline(h=50.47,col=1,lty=2)

title(main="Genes")



##############################
# Plot 2: Replication origins
##############################

origins_2x = read.table("Results/FRAG00062_2x_bias_data_replication_origins_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
origins_3x = read.table("Results/FRAG00062_3x_bias_data_replication_origins_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)


plot(x=origins_2x$V3, y=origins_2x$V7, type="l",  xlim=c(-20000,20000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
plot(x=origins_3x$V3, y=origins_3x$V7, type="l",  xlim=c(-20000,20000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="blue", axes=F)

abline(v=0,col=1,lty=3)
abline(h=4.17,col=1,lty=2)

title(main="Replication origins")



# Add series labels
mtext("Distance of windows from edge of block in bp", side=1, outer=T, at=0.5, cex=1.2)
mtext("% of window occupied by feature", side=2, outer=T, at=0.5, cex=1.2)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("Duplicated block", "Triplicated block"), xpd = TRUE, horiz = TRUE, inset = c(0, 0), bty = "n",  col = c("red", "blue"), cex=1.2, lty = c(1, 1), lwd=3)
# xpd = TRUE tells R that it is OK to plot outside the region 
# horiz = TRUE tells R that I want a horizontal legend 
# inset = c(x,y) tells R how to move the legend relative to the 'bottom' location 
# bty = 'n' means that 'no' box will be drawn around it 
# pch and col are the types and colors of points 
