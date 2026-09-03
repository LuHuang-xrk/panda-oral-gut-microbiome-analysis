#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use FindBin qw/$Bin/;
my $VERSION = "1.1"; 
my $DATE = "2013-10-22";

my %opts;
GetOptions( \%opts,"otu=s","filter=s","top=i","color=s","tax=s","o=s","lw=s","w=f","h=f","lcex=f","scex=f","p=f","spa=f");
my $usage = <<"USAGE";
	Program : $0
	Discription : updated_20190619_add_sample_name_cex  
	Version : $0 $VERSION 
	Usage : perl $0 [options]
		-otu*	otu table
		-tax*	tax table
		-filter filter by p/top,p:num/all<=p push to Low_abundance,top:top abundance ,default:p
		-top     must > 0 integer ,invalid when -filter p default:20
		-color   color list file
		
		-lw 	the proportion of the three parts in the graph default :"0.82-1.8-1"
		-w  	the width of the graph default : 10
		-h  	the height of the graph default : 5
		-lcex	the cex of the legend default : 1.3
		-scex   the cex of the sample name default : 1.0
		-p  	push those numbers in all samples below this threshold into "Others" default : 0.01
		-spa	barplot option "space" default : 0.3
USAGE

die $usage if(!($opts{otu}&&$opts{tax}));
$opts{lw}=$opts{lw}?$opts{lw}:"0.82-1.8-1";
$opts{w}=$opts{w}?$opts{w}:10;
$opts{h}=$opts{h}?$opts{h}:5;
$opts{lcex}=$opts{lcex}?$opts{lcex}:1.3;
$opts{scex}=$opts{scex}?$opts{scex}:1;
$opts{p}=$opts{p}?$opts{p}:0.01;
$opts{spa}=$opts{spa}?$opts{spa}:0.3;

$opts{color}=$opts{color}?$opts{color}:"FALSE";
$opts{top}=defined $opts{top}?$opts{top}:20;###map
$opts{filter}=defined $opts{filter}?$opts{filter}:"p";###map

#Version: 3.5

open CMD,">cmd.r";
print CMD "

#mycol <- c(119, 132, 147, 454, 89, 404, 123, 529, 463, 104, 552, 28, 54, 84, 256, 100, 558, 43, 652, 31, 610, 477, 588, 99, 81, 503, 562, 76, 96, 495)
#mycol <- c(34, 51, 142, 23, 50, 27, 31, 75, 525, 62, 119, 46, 475, 554, 622, 483, 657, 545, 402, 477, 503, 40, 115, 5, 376,473,546,482)
mycol <- c(34, 51, 142, 26, 31, 371, 36, 7, 12, 30, 84, 88, 116, 121, 77, 56, 386, 373, 423, 435, 438, 471, 512, 130, 52, 47, 6, 11, 43, 54, 367, 382, 422, 4, 8, 375, 124, 448, 419, 614, 401, 403, 613, 583, 652, 628, 633, 496, 638, 655, 132, 503, 24)
mycol <- colors()[rep(mycol, 3)]

otu0 <- read.table(\"$opts{otu}\", sep = \"\\t\", head = T, check.names = F)
rownames(otu0) <- otu0[, 1]
otu1 <- otu0[, -1]
otu_r <- apply(otu1, 2, function(x) x/sum(x)*100)
library(vegan)
dist_bray <- vegdist(t(otu_r), method = \"bray\")
write.table(as.matrix(dist_bray), \"otu_bray.xls\", row.names = TRUE, sep = \"\\t\", quote = F)
library(ape,lib=\"/home/lianwen/R/x86_64-pc-linux-gnu-library/3.4/\")
hc <- hclust(dist_bray, method = \"average\")
hc_tre <- as.phylo(hc)
write.tree(hc_tre, \"otu_hc.tre\")


tax0 <- read.table(\"$opts{tax}\", sep = \"\\t\", head = T, check.names = F)
rownames(tax0) <- tax0[, 1]
tax1 <- tax0[, -1]
tax2 <- tax1[, hc\$labels[hc\$order]]
rowsum <- sapply(1:nrow(tax2), function(x) sum(tax2[x, ]))

perc <- $opts{p}
filter <- \"$opts{filter}\"
top <- $opts{top}




if((perc > 0) & (filter == \"p\")){
  tax3 <- tax2[order(rowsum, decreasing = TRUE), ]
  rowsum_perc <- sapply(1:nrow(tax3), function(x) sum(tax3[x, ])/sum(rowsum))
  tax_max <- tax3[which(rowsum_perc >= perc),]
  tax_min <- tax3[which(rowsum_perc < perc),]
  other <- sapply(1:ncol(tax_min),function(x) sum(tax_min[,x]))
  tax_new <- rbind(tax_max,other)
  rownames(tax_new)[nrow(tax_new)] <-  \"Others\"
  tax3 <- tax_new
}else if ((top >0) & (filter == \"top\")) {
     tax3 <- tax2[order(rowsum, decreasing = TRUE), ]
     otu <- tax3
     if ($opts{top} >=  dim(otu)[1]){
     top <- dim(otu)[1]}else {
     top <- top } 
     otu <- otu[names(sort(apply(otu,1,sum),decreasing=T)),]
     otu1<-otu[1:top,]
     if (top < dim(otu)[1]) {
     otu_xp<-otu[(top+1):nrow(otu),]
     other <- sapply(1:ncol(otu_xp),function(y) sum(otu_xp[,y]))
     otu <-rbind(otu[1:top,],other)
     rownames(otu)[nrow(otu)] <-\"Others\"
     tax3 <- otu }
     }else {
  tax3 <- tax2[order(rowsum, decreasing = TRUE), ]
}

tax4 <- as.matrix(sapply(1:ncol(tax3), function(x) tax3[, x]/sum(tax3[, x])))
colnames(tax4) <- colnames(tax3)
rownames(tax4) <- rownames(tax3)
tax5 <- as.data.frame(tax4)
        tax5\$ID <- rownames(tax5)
	        tax5 <- tax5[,c(\"ID\",colnames(tax5)[1:(ncol(tax5)-1)])]


write.table(tax5,\"all.bar\",sep=\"\\t\",eol=\"\\n\",quote=FALSE,row.names=F)
mycol<-c(\"#EE3B3B\",\"#458100\",\"#FFD700\",\"#0000FF\",\"#8A2BE2\",\"#8B3A62\",\"#912323\",\"#898380\",\"#458F74\",\"#00008B\",\"#8B008B\",\"#A2CD5A\",\"#FF147F\",\"#00BFFF\",\"#EEAD0E\",\"#8B4F13\",\"#89864E\",\"#FF6A6A\",\"#8B5F65\",\"#8470FF\",\"#B0C4DE\",\"#5D478B\",\"#8B4789\",\"#1C86EE\",\"#D2691E\",\"#7FFF00\",\"#CDB6B0\",\"#66CDAA\",\"#98F5FF\",\"#EE7611\",\"#FF69B4\",\"#F0E68C\",\"#CD8C95\",\"#FFF9DB\",\"#7FFFD4\",\"#CD5555\",\"#009ACD\",\"#32CD32\",\"#FFB6C1\",\"#008B45\",\"#B2DFF8\",\"#68858B\",\"#00CD66\",\"#8B8686\",\"#FFFF00\",\"#CDB5CD\",\"#CD4F39\",\"#9ACD32\",\"#00C5CD\",\"#CDCD00\")
        mycol <- rep(mycol,50)
if (\"$opts{color}\" == \"FALSE\") {mycol<-mycol[1:(dim(tax5)[1])]} else {mycol <- readLines(\"$opts{color}\")}







pdf(\"treebar.pdf\", width = $opts{w}, height = $opts{h})
lw <- unlist(strsplit(\"$opts{lw}\",\"-\"))
layout(matrix(c(1, 2, 3), 1, 3), widths = lw)

par(mar = c(3.5, 2, 4.5, 0))
source(\"$Bin/plot.phylo2.r\")
plot.phylo2(hc_tre, type = \"phylogram\", cex = $opts{scex}, adj = 1)
title(main = \"Similarity\", line = 0, cex.main = 2, font.main= 4)

par(mar = c(3, 0, 4, 0))
barplot(tax4, space =$opts{spa}, horiz = T, border = NA, offset = 0.5, col = mycol[1:nrow(tax4)], axes = FALSE, axisnames = FALSE)
title(main = \"Taxonomic composition\", line = 0, cex.main = 2, font.main= 4)

plot.new()
par(mar = c(4, 0, 4, 1))
leng <- rownames(tax4)
legend(0.1, 1, legend = leng, fill = mycol[1:nrow(tax4)], bty = \"n\", cex = $opts{lcex})
title(main = \"Taxon\", line = 0, cex.main = 2, font.main= 4)
dev.off()

svg(\"treebar.svg\", width = $opts{w}, height = $opts{h})
lw <- unlist(strsplit(\"$opts{lw}\",\"-\"))
layout(matrix(c(1, 2, 3), 1, 3), widths = lw)

par(mar = c(3.5, 2, 4.5, 0))
source(\"$Bin/plot.phylo2.r\")
plot.phylo2(hc_tre, type = \"phylogram\", cex = $opts{scex}, adj = 1)
title(main = \"Similarity\", line = 0, cex.main = 2, font.main= 4)

par(mar = c(3, 0, 4, 0))
barplot(tax4, space =$opts{spa}, horiz = T, border = NA, offset = 0.5, col = mycol[1:nrow(tax4)], axes = FALSE, axisnames = FALSE)
title(main = \"Taxonomic composition\", line = 0, cex.main = 2, font.main= 4)

plot.new()
par(mar = c(4, 0, 4, 1))
leng <- rownames(tax4)
legend(0.1, 1, legend = leng, fill = mycol[1:nrow(tax4)], bty = \"n\", cex = $opts{lcex})
title(main = \"Taxon\", line = 0, cex.main = 2, font.main= 4)
dev.off()


png(\"treebar.png\", width = $opts{w}*240,height=$opts{h}*200,res=300)
lw <- unlist(strsplit(\"$opts{lw}\",\"-\"))
layout(matrix(c(1, 2, 3), 1, 3), widths = lw)

par(mar = c(3.5, 2, 4.5, 0))
source(\"$Bin/plot.phylo2.r\")
plot.phylo2(hc_tre, type = \"phylogram\", cex = $opts{scex}, adj = 1)
title(main = \"Similarity\", line = 0, cex.main = 2, font.main= 4)

par(mar = c(3, 0, 4, 0))
barplot(tax4, space =$opts{spa}, horiz = T, border = NA, offset = 0.5, col = mycol[1:nrow(tax4)], axes = FALSE, axisnames = FALSE)
title(main = \"Taxonomic composition\", line = 0, cex.main = 2, font.main= 4)

plot.new()
par(mar = c(4, 0, 4, 1))
leng <- rownames(tax4)
legend(0.1, 1, legend = leng, fill = mycol[1:nrow(tax4)], bty = \"n\", cex = $opts{lcex})
title(main = \"Taxon\", line = 0, cex.main = 2, font.main= 4)
dev.off()


";
`R --restore --no-save < cmd.r`;


