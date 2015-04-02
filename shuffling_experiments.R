

### FIXED DATA ###

# 2x state2 data
data1 = read.table("~/Work/Code/FRAG_project/shuffled_distances_2x_state2.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data2 = read.table("~/Work/Code/FRAG_project/real_distances_2x_state2.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(data1$V1, xlim=c(0,75000), breaks=50, col="red", main="Chromatin state 2 at duplicated breakpoints", xlab="Average distance to nearest breakpoint")
hist(data2$V1, add=T, breaks=50, col=rgb(0, 0, 1, 0.5))
abline(v=2778, lwd=3, col="blue")
legend(50000,45, c("Shuffled","Real"), lty=c(1,1), lwd=c(2.5,2.5),col=c("red","blue"))

# 3x state2 data
data1 = read.table("~/Work/Code/FRAG_project/shuffled_distances_3x_state2.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data2 = read.table("~/Work/Code/FRAG_project/real_distances_3x_state2.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(data1$V1, xlim=c(0,75000), breaks=50, col="red", main="Chromatin state 2 at triplicated breakpoints", xlab="Average distance to nearest breakpoint")
hist(data2$V1, add=T, breaks=50, col=rgb(0, 0, 1, 0.5))
abline(v=2064, lwd=3, col="blue")
legend(50000,45, c("Shuffled","Real"), lty=c(1,1), lwd=c(2.5,2.5),col=c("red","blue"))

# state2 data 2x vs 3x
data_2x = read.table("~/Work/Code/FRAG_project/shuffled_distances_2x_state2.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data_3x = read.table("~/Work/Code/FRAG_project/shuffled_distances_3x_state2.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(data_2x$V1, xlim=c(0,75000), breaks=50, col="red", main="Chromatin state 2", xlab="Average distance to nearest breakpoint")
hist(data_3x$V1, xlim=c(0,95000), breaks=50, col=rgb(0, 0, 1, 0.5), add=T)
abline(v=2778, lwd=3, col="red")
abline(v=2064, lwd=3, col="blue")
legend(50000,50, c("2x","3x"), lty=c(1,1), lwd=c(2.5,2.5),col=c("red","blue"))


# 2x state3 data
data1 = read.table("~/Work/Code/FRAG_project/shuffled_distances_2x_state3.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data2 = read.table("~/Work/Code/FRAG_project/real_distances_2x_state3.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(data1$V1, xlim=c(0,95000), breaks=50, col="red", main="Chromatin state 3 at duplicated breakpoints", xlab="Average distance to nearest breakpoint")
hist(data2$V1, add=T, breaks=50, col=rgb(0, 0, 1, 0.5))
abline(v=3176, lwd=3, col="blue")
legend(65000,60, c("Shuffled","Real"), lty=c(1,1), lwd=c(2.5,2.5),col=c("red","blue"))

# 3x state3 data
data1 = read.table("~/Work/Code/FRAG_project/shuffled_distances_3x_state3.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data2 = read.table("~/Work/Code/FRAG_project/real_distances_3x_state3.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(data1$V1, xlim=c(0,95000), breaks=50, col="red", main="Chromatin state 3 at triplicated breakpoints", xlab="Average distance to nearest breakpoint")
hist(data2$V1, add=T, breaks=50, col=rgb(0, 0, 1, 0.5))
abline(v=4242, lwd=3, col="blue")
legend(65000,60, c("Shuffled","Real"), lty=c(1,1), lwd=c(2.5,2.5),col=c("red","blue"))

# state3 data 2x vs 3x
data_2x = read.table("~/Work/Code/FRAG_project/shuffled_distances_2x_state3.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data_3x = read.table("~/Work/Code/FRAG_project/shuffled_distances_3x_state3.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(data_2x$V1, xlim=c(0,95000), breaks=50, col="red", main="Chromatin state 3", xlab="Average distance to nearest breakpoint")
hist(data_3x$V1, xlim=c(0,95000), breaks=50, col=rgb(0, 0, 1, 0.5), add=T)
abline(v=3176, lwd=3, col="red")
abline(v=4242, lwd=3, col="blue")
legend(75000,60, c("2x","3x"), lty=c(1,1), lwd=c(2.5,2.5),col=c("red","blue"))



# 2x DHS data
data1 = read.table("~/Work/Code/FRAG_project/shuffled_distances_2x_DHS.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data2 = read.table("~/Work/Code/FRAG_project/real_distances_2x_DHS.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(data1$V1, xlim=c(0,25000), breaks=50, col="red", main="DHS at duplicated breakpoints", xlab="Average distance to nearest breakpoint")
hist(data2$V1, add=T, breaks=50, col=rgb(0, 0, 1, 0.5))
abline(v=1313, lwd=3, col="blue")
legend(18000,70, c("Shuffled","Real"), lty=c(1,1), lwd=c(2.5,2.5),col=c("red","blue"))

# 3x DHS data
data1 = read.table("~/Work/Code/FRAG_project/shuffled_distances_3x_DHS.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data2 = read.table("~/Work/Code/FRAG_project/real_distances_3x_DHS.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(data1$V1, xlim=c(0,25000), breaks=50, col="red", main="DHS at triplicated breakpoints", xlab="Average distance to nearest breakpoint")
hist(data2$V1, add=T, breaks=50, col=rgb(0, 0, 1, 0.5))
abline(v=1556, lwd=3, col="blue")
legend(18000,70, c("Shuffled","Real"), lty=c(1,1), lwd=c(2.5,2.5),col=c("red","blue"))

# DHS data 2x vs 3x
data_2x = read.table("~/Work/Code/FRAG_project/shuffled_distances_2x_DHS.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data_3x = read.table("~/Work/Code/FRAG_project/shuffled_distances_3x_DHS.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(data_2x$V1, xlim=c(0,25000), breaks=50, col="red", main="DHS", xlab="Average distance to nearest breakpoint")
hist(data_3x$V1, xlim=c(0,25000), breaks=50, col=rgb(0, 0, 1, 0.5), add=T)
abline(v=1313, lwd=3, col="red")
abline(v=1556, lwd=3, col="blue")
legend(18000,70, c("Shuffled","Real"), lty=c(1,1), lwd=c(2.5,2.5),col=c("red","blue"))



# 2x origin data
data1 = read.table("~/Work/Code/FRAG_project/shuffled_distances_2x_origins.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data2 = read.table("~/Work/Code/FRAG_project/real_distances_2x_origins.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(data1$V1, xlim=c(0,80000), breaks=50, col="red", main="Replication origins at duplicated breakpoints", xlab="Average distance to nearest breakpoint")
hist(data2$V1, add=T, breaks=50, col=rgb(0, 0, 1, 0.5))
abline(v=51408, lwd=3, col="blue")
legend(58000,55, c("Shuffled","Real"), lty=c(1,1), lwd=c(2.5,2.5),col=c("red","blue"))

# 3x origins data
data1 = read.table("~/Work/Code/FRAG_project/shuffled_distances_3x_origins.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data2 = read.table("~/Work/Code/FRAG_project/real_distances_3x_origins.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(data1$V1, xlim=c(0,80000), breaks=50, col="red", main="Replication origins at triplicated breakpoints", xlab="Average distance to nearest breakpoint")
hist(data2$V1, add=T, breaks=50, col=rgb(0, 0, 1, 0.5))
abline(v=38386, lwd=3, col="blue")
legend(58000,55, c("Shuffled","Real"), lty=c(1,1), lwd=c(2.5,2.5),col=c("red","blue"))

# origins data 2x vs 3x
data_2x = read.table("~/Work/Code/FRAG_project/shuffled_distances_2x_origins.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data_3x = read.table("~/Work/Code/FRAG_project/shuffled_distances_3x_origins.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(data_2x$V1, xlim=c(0,80000), breaks=50, col="red", main="Replication origins", xlab="Average distance to nearest breakpoint")
hist(data_3x$V1, xlim=c(0,80000), breaks=50, col=rgb(0, 0, 1, 0.5), add=T)
abline(v=51408, lwd=3, col="red")
abline(v=38386, lwd=3, col="blue")
legend(58000,55, c("2x","3x"), lty=c(1,1), lwd=c(2.5,2.5),col=c("red","blue"))



# 2x pseudogene data
data1 = read.table("~/Work/Code/FRAG_project/shuffled_distances_2x_pseudogenes.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data2 = read.table("~/Work/Code/FRAG_project/real_distances_2x_pseudogenes.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(data1$V1, xlim=c(0,250000), breaks=50, col="red", main="Pseudogenes at duplicated breakpoints", xlab="Average distance to nearest breakpoint")
hist(data2$V1, add=T, breaks=50, col=rgb(0, 0, 1, 0.5))
abline(v=71393, lwd=3, col="blue")
legend(180000,65, c("Shuffled","Real"), lty=c(1,1), lwd=c(2.5,2.5),col=c("red","blue"))

# 3x pseudogenes data
data1 = read.table("~/Work/Code/FRAG_project/shuffled_distances_3x_pseudogenes.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data2 = read.table("~/Work/Code/FRAG_project/real_distances_3x_pseudogenes.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(data1$V1, xlim=c(0,250000), breaks=50, col="red", main="Pseudogenes at triplicated breakpoints", xlab="Average distance to nearest breakpoint")
hist(data2$V1, add=T, breaks=50, col=rgb(0, 0, 1, 0.5))
abline(v=146582, lwd=3, col="blue")
legend(180000,65, c("Shuffled","Real"), lty=c(1,1), lwd=c(2.5,2.5),col=c("red","blue"))

# pseudogenes data 2x vs 3x
data_2x = read.table("~/Work/Code/FRAG_project/shuffled_distances_2x_pseudogenes.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data_3x = read.table("~/Work/Code/FRAG_project/shuffled_distances_3x_pseudogenes.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(data_2x$V1, xlim=c(0,250000), breaks=50, col="red", main="Pseudogenes", xlab="Average distance to nearest breakpoint")
hist(data_3x$V1, xlim=c(0,250000), breaks=50, col=rgb(0, 0, 1, 0.5), add=T)
abline(v=71393, lwd=3, col="red")
abline(v=146582, lwd=3, col="blue")
legend(200000,65, c("2x","3x"), lty=c(1,1), lwd=c(2.5,2.5),col=c("red","blue"))





# 2x gene data
data1 = read.table("~/Work/Code/FRAG_project/shuffled_distances_2x_genes.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data2 = read.table("~/Work/Code/FRAG_project/real_distances_2x_genes.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(data1$V1, xlim=c(0,14000), breaks=50, col="red", main="Genes at duplicated breakpoints", xlab="Average distance to nearest breakpoint")
hist(data2$V1, add=T, breaks=50, col=rgb(0, 0, 1, 0.5))
abline(v=598, lwd=3, col="blue")
legend(10000,50, c("Shuffled","Real"), lty=c(1,1), lwd=c(2.5,2.5),col=c("red","blue"))

# 3x genes data
data1 = read.table("~/Work/Code/FRAG_project/shuffled_distances_3x_genes.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data2 = read.table("~/Work/Code/FRAG_project/real_distances_3x_genes.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(data1$V1, xlim=c(0,14000), breaks=50, col="red", main="Genes at triplicated breakpoints", xlab="Average distance to nearest breakpoint")
hist(data2$V1, add=T, breaks=50, col=rgb(0, 0, 1, 0.5))
abline(v=837, lwd=3, col="blue")
legend(10000,50, c("Shuffled","Real"), lty=c(1,1), lwd=c(2.5,2.5),col=c("red","blue"))



# pseudogenes data 2x vs 3x
data_2x = read.table("~/Work/Code/FRAG_project/shuffled_distances_2x_pseudogenes.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
data_3x = read.table("~/Work/Code/FRAG_project/shuffled_distances_3x_pseudogenes.txt", na.strings = "", blank.lines.skip=TRUE, fill=TRUE, stringsAsFactors=FALSE)
hist(data_2x$V1, xlim=c(0,14000), breaks=50, col="red", main="Pseudogenes", xlab="Average distance to nearest breakpoint")
hist(data_3x$V1, xlim=c(0,14000), breaks=50, col=rgb(0, 0, 1, 0.5), add=T)
abline(v=598, lwd=3, col="red")
abline(v=837, lwd=3, col="blue")
legend(10000,50, c("2x","3x"), lty=c(1,1), lwd=c(2.5,2.5),col=c("red","blue"))











