# Supplemenal figure for Han's paper using data from FRAG00062 line
# plotting variation in frequency of various genomic features
# as we move through a breakpoint region (between duplicated and unduplicated block
# or between duplicated and triplicated blocks)
# these all represent flipped datasets (correcting to reflect inside & outside blocks)
# plotting using window size 2000 bp, and step size 200 bp

# make large bottom outer margin
par(oma = c(7, 3, 1, 1))

# set layout of a 2x2 plot (and reduce margins of indvididual plots)
par(mfrow=c(3,5), mai=c(0.7, 0.8, 0.7, 0.2))


##############################
# Plot 1: Genes
##############################
data_2x = read.table("Results/FRAG00062_2x_bias_data_genes_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data_3x = read.table("Results/FRAG00062_3x_bias_data_genes_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
plot(x=data_2x$V3,  y=data_2x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(45,70), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
plot(x=data_3x$V3,  y=data_3x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(45,70), xlab=NA, lwd = 3, ylab=NA, col="blue",  axes=F)
abline(v=0, col=1, lty=3)
abline(h=50.47,col=1,lty=2)
title(main="Genes")


##############################
# Plot 2: Pseudogenes
##############################
data_2x = read.table("Results/FRAG00062_2x_bias_data_pseudogene_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data_3x = read.table("Results/FRAG00062_3x_bias_data_pseudogene_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
plot(x=data_2x$V3,  y=data_2x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
plot(x=data_3x$V3,  y=data_3x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="blue",  axes=F)
abline(v=0, col=1, lty=3)
abline(h=0.9,col=1,lty=2)
title(main="Pseudogenes")


##############################
# Plot 3: Transposable elements
##############################
data_2x = read.table("Results/FRAG00062_2x_bias_data_transposable_element_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data_3x = read.table("Results/FRAG00062_3x_bias_data_transposable_element_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
plot(x=data_2x$V3,  y=data_2x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
plot(x=data_3x$V3,  y=data_3x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="blue",  axes=F)
abline(v=0, col=1, lty=3)
abline(h=20,col=1,lty=2)
title(main="Transposable elements")


##############################
# Plot 4: Satellite repeats
##############################
data_2x = read.table("Results/FRAG00062_2x_bias_data_satellite_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data_3x = read.table("Results/FRAG00062_3x_bias_data_satellite_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
plot(x=data_2x$V3,  y=data_2x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
plot(x=data_3x$V3,  y=data_3x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="blue",  axes=F)
abline(v=0, col=1, lty=3)
abline(h=6.2,col=1,lty=2)
title(main="Satellite repeats")


##############################
# Plot 5: Replication origins
##############################
data_2x = read.table("Results/FRAG00062_2x_bias_data_replication_origins_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data_3x = read.table("Results/FRAG00062_3x_bias_data_replication_origins_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
plot(x=data_2x$V3,  y=data_2x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
plot(x=data_3x$V3,  y=data_3x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="blue",  axes=F)
abline(v=0, col=1, lty=3)
abline(h=4.2,col=1,lty=2)
title(main="Replication origins")

##############################
# Plot 6: DHS
##############################
data_2x = read.table("Results/FRAG00062_2x_bias_data_DHS_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data_3x = read.table("Results/FRAG00062_3x_bias_data_DHS_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
plot(x=data_2x$V3,  y=data_2x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
plot(x=data_3x$V3,  y=data_3x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="blue",  axes=F)
abline(v=0, col=1, lty=3)
abline(h=12.5,col=1,lty=2)
title(main="DNAse I Hypersensitive sites")

##############################
# Plot 7: chromatin state 1
##############################
data_2x = read.table("Results/FRAG00062_2x_bias_data_state1_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data_3x = read.table("Results/FRAG00062_3x_bias_data_state1_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
plot(x=data_2x$V3,  y=data_2x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
plot(x=data_3x$V3,  y=data_3x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="blue",  axes=F)
abline(v=0, col=1, lty=3)
abline(h=11.1,col=1,lty=2)
title(main="Chromatin state 1")


##############################
# Plot 8: chromatin state 2
##############################
data_2x = read.table("Results/FRAG00062_2x_bias_data_state2_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data_3x = read.table("Results/FRAG00062_3x_bias_data_state2_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
plot(x=data_2x$V3,  y=data_2x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
plot(x=data_3x$V3,  y=data_3x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="blue",  axes=F)
abline(v=0, col=1, lty=3)
abline(h=9.6,col=1,lty=2)
title(main="Chromatin state 2")


##############################
# Plot 9: chromatin state 3
##############################
data_2x = read.table("Results/FRAG00062_2x_bias_data_state3_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data_3x = read.table("Results/FRAG00062_3x_bias_data_state3_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
plot(x=data_2x$V3,  y=data_2x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
plot(x=data_3x$V3,  y=data_3x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="blue",  axes=F)
abline(v=0, col=1, lty=3)
abline(h=10.1,col=1,lty=2)
title(main="Chromatin state 3")


##############################
# Plot 10: chromatin state 4
##############################
data_2x = read.table("Results/FRAG00062_2x_bias_data_state4_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data_3x = read.table("Results/FRAG00062_3x_bias_data_state4_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
plot(x=data_2x$V3,  y=data_2x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
plot(x=data_3x$V3,  y=data_3x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="blue",  axes=F)
abline(v=0, col=1, lty=3)
abline(h=14.8,col=1,lty=2)
title(main="Chromatin state 4")


##############################
# Plot 11: chromatin state 5
##############################
data_2x = read.table("Results/FRAG00062_2x_bias_data_state5_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data_3x = read.table("Results/FRAG00062_3x_bias_data_state5_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
plot(x=data_2x$V3,  y=data_2x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
plot(x=data_3x$V3,  y=data_3x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="blue",  axes=F)
abline(v=0, col=1, lty=3)
abline(h=13.4,col=1,lty=2)
title(main="Chromatin state 5")


##############################
# Plot 12: chromatin state 6
##############################
data_2x = read.table("Results/FRAG00062_2x_bias_data_state6_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data_3x = read.table("Results/FRAG00062_3x_bias_data_state6_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
plot(x=data_2x$V3,  y=data_2x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
plot(x=data_3x$V3,  y=data_3x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="blue",  axes=F)
abline(v=0, col=1, lty=3)
abline(h=9.9,col=1,lty=2)
title(main="Chromatin state 6")


##############################
# Plot 13: chromatin state 7
##############################
data_2x = read.table("Results/FRAG00062_2x_bias_data_state7_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data_3x = read.table("Results/FRAG00062_3x_bias_data_state7_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
plot(x=data_2x$V3,  y=data_2x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
plot(x=data_3x$V3,  y=data_3x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="blue",  axes=F)
abline(v=0, col=1, lty=3)
abline(h=9.9,col=1,lty=2)
title(main="Chromatin state 7")


##############################
# Plot 14: chromatin state 8
##############################
data_2x = read.table("Results/FRAG00062_2x_bias_data_state8_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data_3x = read.table("Results/FRAG00062_3x_bias_data_state8_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
plot(x=data_2x$V3,  y=data_2x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
plot(x=data_3x$V3,  y=data_3x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="blue",  axes=F)
abline(v=0, col=1, lty=3)
abline(h=9.0,col=1,lty=2)
title(main="Chromatin state 8")


##############################
# Plot 15: chromatin state 9
##############################
data_2x = read.table("Results/FRAG00062_2x_bias_data_state9_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data_3x = read.table("Results/FRAG00062_3x_bias_data_state9_2000_200_flipped.tsv", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
plot(x=data_2x$V3,  y=data_2x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="red")
par(new=TRUE)
plot(x=data_3x$V3,  y=data_3x$V7, type="l",  xlim=c(-10000, 10000), ylim=c(0,25), xlab=NA, lwd = 3, ylab=NA, col="blue",  axes=F)
abline(v=0, col=1, lty=3)
abline(h=12.3,col=1,lty=2)
title(main="Chromatin state 9")



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
