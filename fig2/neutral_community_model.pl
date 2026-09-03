#!/usr/bin/perl -w
###################################################################
#####Author : GuoqiLiu                                             #
#####Date   : 2024-04-18                                           #
#####Copyright (C) 2019~ Mingke Biotechnology (Hangzhou) Co., Ltd. #
#####Contact: liuguoqi@mingkebio.com                               #
#####Suppose: Auto diversity report  program                       #
#####step :                                                        #
#####1. diver_pipe_V3.pl generate basic report                     #
#####2. modify html and data                                       #
#####3. html2pdf report                                            #
#####4. Auto-send report ? whether or not                          #
#####Platform :                                                    #
##############Ubuntu 18.04 & Windows10 ,Perl v5.26 #################
####################################################################
#Log:
##V1 version 2024-04-21 calculate neutral_community_model (NCM)

##
use strict;
use warnings;
use Getopt::Long;
use POSIX; 
use FindBin qw($Bin);
my %opts;
my $v = "v2025-04-21";

my $help;

my $method =0;
GetOptions (\%opts,"i=s","t=s","o=s","grid=s","w=f","h=f","m=s","g=s","color=s","rlc=f","point=s","method!"=>\$method,"clc=f","llc=f","lg=i","angle=i","test=s");#,"scale=s");
my $usage = <<"USAGE";
#######################################################################################################################################################
#                 Program : perl $0 
#                 Version :  $v
#######################################################################################################################################################
#                -i*  	   <str>            otu/ASV table  
#               -color     <str>            color  file
#
#                Example : perl $0  -i otu_table.xls 
#
#########################################################################################################################################################

USAGE
#die $usage if (!($opts{i}&&$opts{g}));
##                -lw    split plot in width,defalt (four numbers) :0.1:0.2:4:1
#                -lh    split plot in heigth,defalt (three numbers) :0.3:5.5:1.2


#                -method                     mean distance file prepared by customer themselves 
die $usage if (!($opts{i}));
#die $usage if (!($opts{t}));
#die $usage if (!($opts{m}));
#####                -llc     <float>            legend size
#die $usage if (!($opts{i2}));
#$opts{o}=defined$opts{o}?$opts{o}:"output";
$opts{w} ||= 6;
$opts{h} ||= 6;
$opts{g} ||= "F";
$opts{test} ||= "Wilcox-test";
$opts{grid} ||= "F" ;
$opts{color}=$opts{color}?$opts{color}:"FALSE";
$opts{rlc} ||= 12 ;
$opts{clc} ||= 12 ;  
$opts{llc} ||= 12 ;
$opts{angle} ||= 0 ;
$opts{point} ||= "F" ;
$opts{o} ||= "output" ;


#$opts{color}=defined$opts{color}?$opts{color}:"default";
#$opts{row} ||= 0 ;
#$opts{col} ||= 0 ;
#$opts{scale} ||= "none" ;
my $test = $opts{test} ;
my $group = $opts{m} ;

#$opts{method} ||= 0 ;
#if ($method == 0) {$method = "FALSE";} else {$method = "TRUE";}


my @mycol=("#EE3B3B","#458100","#FFD700","#0000FF","#8A2BE2","#8B3A62","#912323","#898380","#458F74","#00008B","#8B008B","#A2CD5A","#FF147F","#00BFFF","#EEAD0E","#8B4F13","#89864E","#FF6A6A","#8B5F65","#8470FF","#B0C4DE","#5D478B","#8B4789","#1C86EE","#D2691E","#7FFF00","#CDB6B0","#66CDAA","#98F5FF","#EE7611","#FF69B4","#F0E68C","#CD8C95","#FFF9DB","#7FFFD4","#CD5555","#009ACD","#32CD32","#FFB6C1","#008B45","#B2DFF8","#68858B","#00CD66","#8B8686","#FFFF00","#CDB5CD","#CD4F39","#9ACD32","#00C5CD","#CDCD00");
my $center="black";
my $up="#29A6A6";
my $down="#A52A2A";

if ($opts{color} eq "FALSE") {
	open T,">color.txt";
	print T "up\t$up\ncenter\t$center\ndown\t$down\n";
}
else {
	open C,$opts{color};
	while(<C>){
	chomp;
	next if(/^#/); 
	my @ab=split/\t/,$_;
        if ($ab[0]=~/center/i){$center=$ab[1];}
        elsif($ab[0]=~/up/i){$up=$ab[1];}
        elsif($ab[0]=~/down/i){$down=$ab[1];}
        else {
        	print "Error!! Cant not find $_\n";
        }
  }
 }













open RSCRIPT,">cmd.r";

print RSCRIPT "

##Fits the neutral model from Sloan et al. 2006 to an OTU table and returns several fitting statistics. Alternatively, will return predicted occurrence frequencies for each OTU based on their abundance in the metacommunity
#Install the following packages if they haven't been availabled in your computer yet 
library(Hmisc)
library(minpack.lm)
library(stats4)
spp <-  read.table(\"$opts{i}\",header = T,sep = \"\\t\",check.names=F,comment.char=\"\",quote=\"\",row.names=1)
spp<-t(spp)
N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit #get the m value
m.ci <- confint(m.fit, \"m\", level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method=\"wilson\", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr# get the R2 value
#Optional: write 3 files: p.csv, freq.csv and freq.pred.csv
#write.csv(p, file = \"j:/iue/aehg11/neutr_model/p.csv\")
#write.csv(freq, file = \"j:/iue/aehg11/neutr_model/freq.csv\")
#write.csv(freq.pred, file = \"j:/iue/aehg11/neutr_model/freq.pred.csv\")
write.table(p,file=\"p.xls\",sep=\"\\t\",quote=F,row.names=FALSE)
write.table(freq,file=\"freq.xls\",sep=\"\\t\",quote=F,row.names=FALSE)
write.table(freq.pred,file=\"freq.pred.xls\",sep=\"\\t\",quote=F,row.names=FALSE)
#Drawing the figure using grid package:
#p is the mean relative abundance
#freq is occurrence frequency
#freq.pred is predicted occurrence frequency
bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
#inter.col<-rep(\'black\',nrow(bacnlsALL))
#inter.col[bacnlsALL\$freq <= bacnlsALL\$Lower]<-\'#A52A2A\'#define the color of below points
#inter.col[bacnlsALL\$freq >= bacnlsALL\$Upper]<-\'#29A6A6\'#define the color of up points
inter.col<-rep(\"$center\",nrow(bacnlsALL))
inter.col[bacnlsALL\$freq <= bacnlsALL\$Lower]<-\"$down\" #define the color of below points
inter.col[bacnlsALL\$freq >= bacnlsALL\$Upper]<-\"$up\" #define the color of up points



pdf(\"neutral_community_model.pdf\")
library(grid)
#grid.layout(nrow = 1, ncol = 2)
grid.newpage()
vplay <- grid.layout(nrow=1,ncol = 1)
pushViewport(viewport(layout = vplay))
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
pushViewport(viewport(h=0.7,w=0.7))
pushViewport(dataViewport(xData=range(log10(bacnlsALL\$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL\$p), bacnlsALL\$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL\$p),bacnlsALL\$freq.pred,gp=gpar(col=\'blue\',lwd=2),default=\'native\')

grid.lines(log10(bacnlsALL\$p),bacnlsALL\$Lower ,gp=gpar(col=\'blue\',lwd=2,lty=2),default=\'native\') 
grid.lines(log10(bacnlsALL\$p),bacnlsALL\$Upper,gp=gpar(col=\'blue\',lwd=2,lty=2),default=\'native\')  
grid.text(y=unit(0,\'npc\')-unit(2.5,\'lines\'),label=\'Mean Relative Abundance (log10)\', gp=gpar(fontface=2)) 
grid.text(x=unit(0,\'npc\')-unit(3,\'lines\'),label=\'Frequency of Occurance\',gp=gpar(fontface=2),rot=90) 

draw.text <- function(just, i, j) {
  grid.text(paste(\"Rsqr=\",round(Rsqr,3),\"\\n\",\"Nm=\",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)}
x <- unit(1:4/5, \"npc\")
y <- unit(1:4/5, \"npc\")
draw.text(c(\"centre\", \"bottom\"), 4, 1)
grid.text(\"neutral community model \", y=unit(26, \"line\"),gp=gpar(fontface=2))
dev.off()
";
#/usr/bin/R 
#/home/guest/bin/R-4.2.0/bin/R ;
#/home/guest/bin/R-4.4.1/bin/R ;
#`R --no-save <cmd.r` ;
`/home/guest/bin/R-4.2.0/bin/R  --no-save <cmd.r` ;

