#!/usr/bin/perl -w
###################################################################
#####Author : GuoqiLiu                                             #
#####Date   : 2024-04-17                                           #
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
##V1 version 2024-06-12 
# https://mp.weixin.qq.com/s?__biz=MzUzMjYyMDE2OQ==&mid=2247485864&idx=2&sn=737c9bc9da9e9552d4c00b7620ac546f&chksm=fab13f4ecdc6b658bca4b3e30ba29a88e2ed1146dc61cbaea8f8370c5db3ed8710adb5f016b1&scene=27 
##
use strict;
use warnings;
use Getopt::Long;
use POSIX; 

my %opts;
my $v = "v2024-06-13";

my $help;

GetOptions (\%opts,"i=s","grid=s","o=s","rep=s","l=s","type=i","alpha=f","sz=f","method=s","cex=f",,"c=i","pc=s","std=s","w=f","h=f","m=s","g=s","color=s","rlc=f","point=s","clc=f","llc=f","lg=i","angle=i","test=s");#,"scale=s");
my $usage = <<"USAGE";
#######################################################################################################################################################
#                 Program : perl $0 
#                 Version :  $v
#######################################################################################################################################################
#                -i*  	   <str>             matrix file  
#                -method*   <str>             {braycurtis},{jaccard},{euclidean},{weighted_unifrac}(need -rep),{unweighted_unifrac}(need -rep)  and {distance}(distance file prepared by customer themselves) can choice,default:braycurtis 
#                -rep      <str>             rep fasta file,if -method weighted_unifrac and unweighted_unifrac must need
#                -m        <str>             map file
#                -g        <str>             group sort list
#                -c        <int>             whether add confidence ellipse or polygon; choose 0/1/2/3; 0 mean not add,1 mean add 90% confidence ellipse ,2 mean add 95% confidence ellipse,3 mean polygon default:0 
#                -color    <str>             group color and pch  
#                -l        <str>             T/F, T display point label ,F not display point label,default:F 
#                -grid     <str>             T/F, T display grid  line ,F not display display grid  line,default:F  
#                -sz      <float>            point size default:4
#                -cex     <float>            point label size default:4
#                -angle    <int>             xlab angle 0-360 default:0
#                -rlc     <float>            xlim size default:12
#                -clc     <float>            ylim size default:12
#                -llc     <float>            legend size default:12
#                -type     <int>             1/2,must -c 1 or 2,  1 mean border color ,2 mean fill color,default 1
#                -alpha   <float>            must -type 2 , fill color transparency 0-1 , default 0.5
#		 -w       <float>            the width of the figure,default:6
#		 -h	  <float>            the height of the figure,default:6
#
#                Example :                   Usage:perl $0  -i otu/ASV_table.xls or distance-matrix.txt 
#########################################################################################################################################################

USAGE
#die $usage if (!($opts{i}&&$opts{g}));
##                -lw    split plot in width,defalt (four numbers) :0.1:0.2:4:1
#                -lh    split plot in heigth,defalt (three numbers) :0.3:5.5:1.2



#die $usage if ( !(defined $opts{i} && $opts{m} && $opts{rep}) );

#die $usage if (!($opts{i} && $opts{method}));

die $usage if ( !$opts{i} ) ;
#die $usage if (!($opts{m}));
#die $usage if (!($opts{rep}));
#

#####                -llc     <float>            legend size
#die $usage if (!($opts{i2}));
#$opts{o}=defined$opts{o}?$opts{o}:"output";
$opts{w} ||= 6;
$opts{h} ||= 6;
$opts{m} ||= "F";
$opts{g} ||= "F";
$opts{rep} ||= "F";
$opts{method} ||= "braycurtis";
$opts{test} ||= "Wilcox-test";
#$opts{top} ||= 15 ;
$opts{color} ||= "F" ;
$opts{grid} ||= "F" ;
$opts{rlc} ||= 12 ;
$opts{clc} ||= 12 ;  
$opts{llc} ||= 12 ;
$opts{angle} ||= 0 ;
$opts{l} ||= "F" ;
$opts{std} ||= "F" ;
#                -std      <str>             T/F, T mean data scale,F not data scale ,default:F
$opts{c} ||= 0 ;
$opts{sz} ||= 4 ;
$opts{cex} ||= 4 ;

$opts{type} ||= 1 ;
$opts{alpha} ||= 0.5 ;

$opts{pc}=defined $opts{pc}?$opts{pc}:"pc1-pc2";
##                -pc       <str>             pc to be draw ,you can choice pc1-pc2,pc2-pc3,pc1-pc3,default: pc1-pc2
my $pc_raw;
if ($opts{pc} eq "pc1-pc2") {$pc_raw="1-2";}
elsif($opts{pc} eq "pc1-pc3") {$pc_raw="1-3";}
elsif($opts{pc} eq "pc2-pc3") {$pc_raw="2-3";}

#$opts{color}=defined$opts{color}?$opts{color}:"default";
#$opts{row} ||= 0 ;
#$opts{col} ||= 0 ;
#$opts{scale} ||= "none" ;

my $angle=$opts{angle};
my $rlc=$opts{rlc};
my $llc=$opts{llc};
my $clc=$opts{clc};


if (($opts{method} eq "weighted_unifrac") & ($opts{rep} eq "F")) {die $usage ;}
if (($opts{method} eq "unweighted_unifrac") & ($opts{rep} eq "F")) {die $usage;}


open II,$opts{i};
chomp (my $headtmp=<II> );
my @headtmp_array1=split/\t/,$headtmp;
@headtmp_array1 = @headtmp_array1[1..$#headtmp_array1] ;
my @headtmp_array2=(); 
while (<II>) {
	chomp;
	my @headtmparray=split/\t/,$_;
	push @headtmp_array2,$headtmparray[0] ;
}
close II;

my $tag = "Lebel"; 
foreach my $i (0..$#headtmp_array1) {
	if ($headtmp_array1[$i] ne $headtmp_array2[$i]) {
		$tag = "error";last; 
	}
}
















# print "$pathqiime2\n"; 
# 
#sub get_time {
#	 my $hhh=shift;
#	 my $gettime=strftime("%Y-%m-%d %H:%M:%S",localtime());
#	 print  "\@\@$gettime====>>>>>>>$hhh\n";
#}
#&get_time("$0 program start running...") ;


#my $qiime2_2="/root/soft/anaconda3/envs/qiime2-2022.2/bin/qiime";
my $qiime2_2="/mnt/sdb/lgq/bin/software/py38/envs/qiime2-2022.2/bin/qiime";
my @qiime2_arr=split/\//,$qiime2_2;
my $pathqiime2=join("/",@qiime2_arr[0..($#qiime2_arr-1)]);
# print "$pathqiime2\n"; 
sub get_time {
	 my $hhh=shift;
	 my $gettime=strftime("%Y-%m-%d %H:%M:%S",localtime());
	 print  "\@\@$gettime====>>>>>>>$hhh\n";
}
&get_time("$0 program start running...") ;


#print "$opts{angle} ------------------\n";

########################################################

use FindBin qw($Bin);
chomp (my $current_dir=`pwd`);


sub beta {
	my ($qiime2,$asvtable,$rep,$method) = @_;
	if (-e "feature_table.biom") {
		`rm feature_table.biom`;
	         `biom convert -i $asvtable   -o feature_table.biom --table-type "OTU table" --to-json`;
	 }
	 else {
		 `biom convert -i $asvtable   -o feature_table.biom --table-type "OTU table" --to-json`;
	 }

	 #`rm tmp_feature_table.xls` if (-e "tmp_feature_table.xls") ;
	 # `rm tree.nwk` if (-e "tree.nwk") ;
	 #open I,$asvtable;
	 #open O,">tmp_feature_table.xls";
	 #chomp (my $sam = <I>); 
	 #my @samgroup = split/\t/,$sam; 
	 #print O "ID\t",join("\t",@samgroup[1..$#samgroup]),"\n";
	 #while (<I>){
	 #	 print O $_;
	 #}
	 #close I;close O;
	 #my @alpha_arr = ("chao1","observed_features","goods_coverage","pielou_e","simpson","shannon","ace");
	 #my @alpha_arr = ('jensenshannon','russellrao', 'matching','dice', 'correlation', 'aitchison','sokalmichener', 'braycurtis', 'yule','sqeuclidean', 'cosine', 'jaccard', 'chebyshev', 'rogerstanimoto','sokalsneath', 'euclidean', 'cityblock', 'kulsinski', 'seuclidean','canberra', 'canberra_adkins', 'hamming', 'minkowski');
	 #if (-d "tmp_beta") {`rm -r tmp_beta`;}else {mkdir "tmp_beta" ;}
	 #if (-e "tree.nwk"){`rm tree.nwk`;}
	  mkdir "tmp_beta"  unless(-d "tmp_beta");
	 open Oo,">beta.sh";
	 #print Oo "export PATH=\"/mnt/sdb/lgq/bin/software/py38/envs/qiime2-2022.2/bin:\$PATH\"\n"; 
	 #$pathqiime2
	 print Oo "export PATH=\"$pathqiime2:\$PATH\"\n";
	 #close Oo;
	 #`sh config.sh`;
	 print Oo "$qiime2 tools  import --input-path feature_table.biom   --output-path tmp_beta/feature_table.qza --type 'FeatureTable[Frequency]'  --input-format BIOMV100Format\n";
	 if (($method ne "unweighted_unifrac") & ($method ne "weighted_unifrac") & ($method ne "weighted_normalized_unifrac") & ($method ne "generalized_unifrac")) {
		 #print Oo "$qiime2 tools  import --input-path feature_table.biom   --output-path tmp_beta/feature_table.qza --type 'FeatureTable[Frequency]'  --input-format BIOMV100Format\n";
	 print Oo "$qiime2    diversity beta --i-table tmp_beta/feature_table.qza --p-metric $method  --o-distance-matrix tmp_beta/$method\n"; 
	 print Oo "$qiime2 tools  extract --input-path tmp_beta/$method.qza --output-path tmp_beta/$method\n";
 }
	 #foreach my $i (@alpha_arr) {
		 #`$qiime2 tools  import --input-path feature_table.biom   --output-path tmp_alpha/feature_table.qza --type 'FeatureTable[Frequency]'  --input-format BIOMV100Format`;
		 #	 print Oo "$qiime2    diversity  alpha --i-table tmp_alpha/feature_table.qza --p-metric $i  --o-alpha-diversity tmp_alpha/$i.qza \n";
		 # print Oo "$qiime2 tools  extract --input-path tmp_alpha/$i.qza --output-path tmp_alpha/$i\n";
		 # }
         elsif (($rep ne "F") & ($method eq "unweighted_unifrac")) {
		 if (-e "tmp_beta/rooted-tree.qza"){
			 print Oo "$qiime2 diversity      beta-phylogenetic --i-table tmp_beta/feature_table.qza --i-phylogeny tmp_beta/rooted-tree.qza --p-metric $method --o-distance-matrix tmp_beta/$method\n";
		 }
		 else {
	 print Oo "$qiime2 tools import --input-path $rep  --output-path tmp_beta/rep-seqs.qza  --type 'FeatureData[Sequence]' \n";
	 print Oo "$qiime2 phylogeny align-to-tree-mafft-fasttree  --i-sequences tmp_beta/rep-seqs.qza --o-alignment tmp_beta/aligned-rep-seqs.qza --o-masked-alignment tmp_beta/masked-aligned-rep-seqs.qza --o-tree tmp_beta/unrooted-tree.qza  --o-rooted-tree tmp_beta/rooted-tree.qza\n";
	 print Oo "$qiime2 diversity      beta-phylogenetic --i-table tmp_beta/feature_table.qza --i-phylogeny tmp_beta/rooted-tree.qza --p-metric $method --o-distance-matrix tmp_beta/$method\n";
 }
	 print Oo "$qiime2 tools  extract --input-path tmp_beta/$method.qza --output-path tmp_beta/$method\n";
	 #print Oo "$qiime2 tools  extract --input-path tmp_alpha/rooted-tree.qza --output-path tmp_alpha/rooted-tree\n";
	   #print Oo "#ln -s tmp_alpha/rooted-tree/*/data/tree.nwk .\n";
	   # print Oo "cp  -a  tmp_alpha/rooted-tree/*/data/tree.nwk .\n"; 
   }
   elsif (($rep ne "F") & ($method eq "weighted_unifrac")) {
	   if (-e "tmp_beta/rooted-tree.qza"){
		   print Oo "$qiime2 diversity      beta-phylogenetic --i-table tmp_beta/feature_table.qza --i-phylogeny tmp_beta/rooted-tree.qza --p-metric $method --o-distance-matrix tmp_beta/$method\n";
	   }
	   else {
	   print Oo "$qiime2 tools import --input-path $rep  --output-path tmp_beta/rep-seqs.qza  --type 'FeatureData[Sequence]' \n";
	   print Oo "$qiime2 phylogeny align-to-tree-mafft-fasttree  --i-sequences tmp_beta/rep-seqs.qza --o-alignment tmp_beta/aligned-rep-seqs.qza --o-masked-alignment tmp_beta/masked-aligned-rep-seqs.qza --o-tree tmp_beta/unrooted-tree.qza  --o-rooted-tree tmp_beta/rooted-tree.qza\n";
	   print Oo "$qiime2 diversity      beta-phylogenetic --i-table tmp_beta/feature_table.qza --i-phylogeny tmp_beta/rooted-tree.qza --p-metric $method --o-distance-matrix tmp_beta/$method\n";
   }
	   print Oo "$qiime2 tools  extract --input-path tmp_beta/$method.qza --output-path tmp_beta/$method\n";
   }
   elsif (($rep ne "F") & ($method eq "weighted_normalized_unifrac")) {
	   if (-e "tmp_beta/rooted-tree.qza"){
		   print Oo "$qiime2 diversity      beta-phylogenetic --i-table tmp_beta/feature_table.qza --i-phylogeny tmp_beta/rooted-tree.qza --p-metric $method --o-distance-matrix tmp_beta/$method\n";
	   }
	   else {
		   print Oo "$qiime2 tools import --input-path $rep  --output-path tmp_beta/rep-seqs.qza  --type 'FeatureData[Sequence]' \n";
		   print Oo "$qiime2 phylogeny align-to-tree-mafft-fasttree  --i-sequences tmp_beta/rep-seqs.qza --o-alignment tmp_beta/aligned-rep-seqs.qza --o-masked-alignment tmp_beta/masked-aligned-rep-seqs.qza --o-tree tmp_beta/unrooted-tree.qza  --o-rooted-tree tmp_beta/rooted-tree.qza\n";
		   print Oo "$qiime2 diversity      beta-phylogenetic --i-table tmp_beta/feature_table.qza --i-phylogeny tmp_beta/rooted-tree.qza --p-metric $method --o-distance-matrix tmp_beta/$method\n";
	   }
	   print Oo "$qiime2 tools  extract --input-path tmp_beta/$method.qza --output-path tmp_beta/$method\n";
   }
   elsif (($rep ne "F") & ($method eq "generalized_unifrac")) {
	   if (-e "tmp_beta/rooted-tree.qza"){
		   print Oo "$qiime2 diversity      beta-phylogenetic --i-table tmp_beta/feature_table.qza --i-phylogeny tmp_beta/rooted-tree.qza --p-metric $method --o-distance-matrix tmp_beta/$method\n";
	   }
	   else {
		   print Oo "$qiime2 tools import --input-path $rep  --output-path tmp_beta/rep-seqs.qza  --type 'FeatureData[Sequence]' \n";
		   print Oo "$qiime2 phylogeny align-to-tree-mafft-fasttree  --i-sequences tmp_beta/rep-seqs.qza --o-alignment tmp_beta/aligned-rep-seqs.qza --o-masked-alignment tmp_beta/masked-aligned-rep-seqs.qza --o-tree tmp_beta/unrooted-tree.qza  --o-rooted-tree tmp_beta/rooted-tree.qza\n";
		   print Oo "$qiime2 diversity      beta-phylogenetic --i-table tmp_beta/feature_table.qza --i-phylogeny tmp_beta/rooted-tree.qza --p-metric $method --o-distance-matrix tmp_beta/$method\n";
	   }
	   print Oo "$qiime2 tools  extract --input-path tmp_beta/$method.qza --output-path tmp_beta/$method\n";
   }








	   #print Oo "cp $Bin/pd.py .;\n";
	   #print Oo "./pd.py \n";
	   #close Oo;
	   # `sh beta.sh`;
   }
   #&beta($qiime2_2,$opts{i},$opts{rep},$opts{method}); 

sub sub_piont {
########################################################
	my ($input) = @_; 
	print $input."\n";
my @pcc = split("-",$pc_raw);
my $PC11 = "MDS".$pcc[0];
my $PC22 = "MDS".$pcc[1];
#my @mycol = ("#df89ff","#0000cd","#00c4ff","#ff8805","#ff5584","#00bd94","#d3b3b0","#4b0082","#c0c0c0","#ffd700","#8b0000","#00ffff","#ff0000","#0000cd","#006400","#ffff00","#008080","#d8bfd8","#40e0d0","#00ff7f","#6a5acd","#adff2f","#00ffff","#ff00ff","#8b4513","#6495ed","#ff6347","#800080","#dc143c","#000000","#7fff00","#d2691e","#ff7f50","#6495ed","#fff8dc","#dc143c","#00ffff");
#my @mycol = ("#df89ff","#0000cd","#00c4ff","#ff8805","#ff5584","#00bd94","#d3b3b0","#4b0082","#c0c0c0","#ffd700","#8b0000","#00ffff","#ff0000","#006400","#ffff00","#008080","#d8bfd8","#40e0d0","#00ff7f","#6a5acd","#adff2f","#ff00ff","#8b4513","#6495ed","#ff6347","#800080","#dc143c","#000000","#7fff00","#d2691e");

my @mycol = ("#d89041","#5fd8a3","#b051c5","#4bb8d6","#c54e54","#4364d0","#df89ff","#0000cd","#00c4ff","#ff8805","#ff5584","#00bd94","#d3b3b0","#4b0082","#c0c0c0","#ffd700","#8b0000","#00ffff","#ff0000","#006400","#ffff00","#008080","#d8bfd8","#40e0d0","#00ff7f","#6a5acd","#adff2f","#ff00ff","#8b4513","#6495ed","#ff6347","#800080","#dc143c","#000000","#7fff00","#d2691e");


my $col="";
my $pch="";
my @ttmp;
my @ttmp2;
my %group2=();
if ($opts{m} ne "F") {
open GROUP,$opts{m} ;
<GROUP>;
while(<GROUP>){
	chomp;
	my @a=split/\t/,$_;
	$group2{$a[1]} =  0 ;  
}
close GROUP;
}

my $sort_group;
my @group_m ; 
my %grp_color_pch=(); 

if ($opts{m} ne "F") {  if ($opts{g} eq "F") {@group_m =  sort (keys %group2); }
else { open TMM,$opts{g};chomp(@group_m =  <TMM>); } }
$sort_group = join("\",\"",@group_m);
$sort_group = "\"".$sort_group."\"";

if ($opts{m} ne "F") {

if ($opts{color} eq "F") {
	open T,">color.txt";
	foreach my $i (0..$#group_m) {
		$grp_color_pch{$group_m[$i]} = $mycol[$i]."\t"."19";
		print T "$group_m[$i]\t$mycol[$i]\t19\n";
		push @ttmp,"\"$group_m[$i]\""."="."\"$mycol[$i]\"";
		push @ttmp2,"\"$group_m[$i]\""."="."19";
	}
	$col=join(",",@ttmp);
	$pch=join(",",@ttmp2);
}
else {
	open C,$opts{color};
	while(<C>){
		chomp;

		my @ab=split/\t/,$_;
		$grp_color_pch{$ab[0]} = $ab[1]."\t".$ab[2] ;
		push @ttmp,"\"$ab[0]\""."="."\"$ab[1]\"";
		push @ttmp2,"\"$ab[0]\""."="."$ab[2]";
	}
	$col=join(",",@ttmp);
	$pch=join(",",@ttmp2);
 }
 }

if ($opts{m} ne "F") {
open IMP,$opts{m};
<IMP>;
open TMPGROUP,">tmp_group_color_pch.tsv";
print TMPGROUP "sample\tgroup\tcolor\tpch\n";
while (<IMP>) {
	chomp;
	        my @a=split/\t/,$_;
		print TMPGROUP "$_\t$grp_color_pch{$a[1]}\n";
	}
}
close IMP;close TMPGROUP;


open CMD,">cmd.r";
print CMD "
warnings_file <- file(\"warnings.log\", \"w\")
sink(warnings_file, type = \"message\")

options(stringsAsFactors=FALSE)
library(\"ggplot2\")
library(\"ggrepel\")
library(\"vegan\")
library(\"plyr\")
basename=\"nmds\"
da <-read.table(\"$input\",sep=\"\\t\",header=T,check.names = FALSE,comment.char=\"\",quote=\"\")
rownames(da) <-as.character(da[,1])
da <-da[,-1]
da <-t(da)
pc_num =as.numeric(unlist(strsplit(\"$pc_raw\",\"-\")))
pc_x =pc_num[1]
pc_y =pc_num[2]




da <-as.dist(da)
set.seed(123)
nmds <-metaMDS(da,k=2)
pc12 <- nmds\$points



#if (\"$opts{std}\" == \"F\"){
#pca <- prcomp(da,scale=FALSE)}else {pca <- prcomp(da,scale=TRUE)}
#pc12 <- pca\$x[,pc_num]
#pc <-summary(pca)\$importance[2,]*100
sites=paste(basename,\"_sites.xls\",sep=\"\")
#impo=paste(basename,\"_importance.xls\",sep=\"\")
#rotat=paste(basename,\"_rotation.xls\",sep=\"\")
#write.table(pca\$x,sites,sep=\"\\t\")
#write.table(summary(pca)\$importance[2,],impo,sep=\"\\t\")
#write.table(summary(pca)\$rotation,rotat,sep=\"\\t\")

write.table(nmds\$points,sites,sep=\"\\t\")
tmpp <- nmds\$stress
pc12 <- data.frame(pc12)
if (\"$opts{m}\" != \"F\") {
group <- read.table(\"tmp_group_color_pch.tsv\",sep=\"\\t\",header = T,check.names=F,comment.char=\"\",quote=\"\",row.names=1)
group\$sample <- rownames(group)
 pc12\$sample <- rownames(pc12)

pc12 <- merge(pc12,group,by=\"sample\")
#head(pc12)
pc12\$group <- factor (pc12\$group,levels=c($sort_group))




df <- pc12[,2:4]
group_border <- ddply(df, \"group\", function(df) df[chull(df[[1]], df[[2]]), ])

if (($opts{c} == 3) && ($opts{type} == 1)){
p<-ggplot(data=pc12,aes(x=$PC11, y=$PC22,shape=group,color=group))+ geom_polygon(data = group_border,fill=\"white\", alpha = 0.2, show.legend = F)   + geom_point(size=$opts{sz}) + theme_bw()+labs(title = paste(\"NMDS\",\"(\",\"stress:\",signif(tmpp,3),\")\",sep=\"\"))+theme(plot.title = element_text(hjust = 0.5))+xlab(\"nmds1\")+ylab(\"nmds2\")
}else if (($opts{c} == 3) && ($opts{type} == 2)) {
p<-ggplot(data=pc12,aes(x=$PC11, y=$PC22,shape=group,color=group,fill=group))+ geom_polygon(data = group_border, alpha = 0.2, show.legend = F)   + geom_point(size=$opts{sz}) + theme_bw()+labs(title =  paste(\"NMDS\",\"(\",\"stress:\",signif(tmpp,3),\")\",sep=\"\"))+theme(plot.title = element_text(hjust = 0.5))+xlab(\"nmds1\")+ylab(\"nmds2\")
}else if ($opts{c} < 3 ) {
p <- ggplot(data=pc12,aes(x=$PC11, y=$PC22,shape=group)) + geom_point(aes(color=group),size=$opts{sz}) + theme_bw()+
labs(title = paste(\"NMDS\",\"(\",\"stress:\",signif(tmpp,3),\")\",sep=\"\"))+theme(
 plot.title = element_text(hjust = 0.5))+xlab(\"nmds1\")+ylab(\"nmds2\") }
 #xlab(paste(\"PC\",pc_x,\": \",round(pc[pc_x],2),\"%\",sep=\"\"))+ylab(paste(\"PC\",pc_y,\" :  \",round(pc[pc_y],2),\"%\",sep=\"\"))
#stat_ellipse(aes(fill = group), level = 0.9,geom = 'polygon', alpha=.1,show.legend = FALSE)

#if ($opts{c} == 0) {p <- p}else if ($opts{c} == 1){
#p <- p+stat_ellipse(aes(color = group), level = 0.9, show.legend = FALSE)
#}else if ($opts{c} == 2) {
#p <- p+stat_ellipse(aes(color = group), level = 0.95, show.legend = FALSE)
#}



if ($opts{c} == 0) {p <- p}else if (($opts{c} == 1) && ($opts{type} == 1)){
p <- p+stat_ellipse(aes(color = group), level = 0.9, show.legend = FALSE)
}else if (($opts{c} == 2) && ($opts{type} == 2)){
p <- p+stat_ellipse(aes(fill = group), level = 0.95,geom = 'polygon', alpha=$opts{alpha},show.legend = FALSE)
}else if (($opts{c} == 1) && ($opts{type} == 2)){
p<- p+stat_ellipse(aes(fill = group), level = 0.9,geom = 'polygon', alpha=$opts{alpha},show.legend = FALSE)
}else if (($opts{c} == 2) && ($opts{type} == 1)){
p <- p+stat_ellipse(aes(color = group), level = 0.95, show.legend = FALSE)
}






if (\"$opts{l}\" != \"F\"){p <- p + geom_text_repel(aes($PC11, $PC22, label=sample),size=$opts{cex},show.legend = FALSE) } #+xlim(c(-12,12))+ylim(c(0,5))+

#p2 <- p+scale_color_manual(name=\"\",values=c($col))+scale_shape_manual(name=\"\",values=c($pch)) # + guides(color = guide_legend(override.aes = list(alpha = 0)))



if ($opts{c} == 0){p2 <- p+scale_color_manual(name=\"\",values=c($col))+scale_shape_manual(name=\"\",values=c($pch)) # + guides(color = guide_legend(override.aes = list(alpha = 0)))
}else if (($opts{c} == 1) && ($opts{type} == 1)){
p2 <- p+scale_color_manual(name=\"\",values=c($col))+scale_shape_manual(name=\"\",values=c($pch))
}else if (($opts{c} == 2) && ($opts{type} == 2)){
p2 <- p+scale_fill_manual(name=\"\",values=c($col))+scale_color_manual(name=\"\",values=c($col))+scale_shape_manual(name=\"\",values=c($pch))
}else if (($opts{c} == 1) && ($opts{type} == 2)){
p2 <- p+scale_fill_manual(name=\"\",values=c($col))+scale_color_manual(name=\"\",values=c($col))+scale_shape_manual(name=\"\",values=c($pch))
}else if (($opts{c} == 2) && ($opts{type} == 1)){
p2<- p+scale_color_manual(name=\"\",values=c($col))+scale_shape_manual(name=\"\",values=c($pch))
}else if (($opts{c} == 3) && ($opts{type} == 1)){
p2<- p+scale_color_manual(name=\"\",values=c($col))+scale_shape_manual(name=\"\",values=c($pch))
}else if (($opts{c} == 3) && ($opts{type} == 2)){
p2<- p+scale_fill_manual(name=\"\",values=c($col))+scale_color_manual(name=\"\",values=c($col))+scale_shape_manual(name=\"\",values=c($pch))
}






if ($angle == 90) {      
p2 <- p2  +    theme(axis.text.x=element_text(angle=$angle,hjust=1,vjust=0.5,size=$rlc,color=\"black\"),axis.title.x=element_text(size=$rlc))}else if ($angle == 0) {p2 <- p2  +    theme(axis.text.x=element_text(angle=$angle,hjust=0.5,size=$rlc,color=\"black\"),axis.title.x=element_text(size=$rlc))
}else {
p2 <- p2  +    theme(axis.text.x=element_text(angle=$angle,hjust=1,vjust=1,size=$rlc,color=\"black\"),axis.title.x=element_text(size=$rlc))
}
p2 <- p2+theme( legend.title = element_text(size=$llc),legend.text = element_text(size=$llc),axis.text.y=element_text(size=$clc,color=\"black\"),axis.title.y=element_text(size=$clc,color=\"black\"))
if (\"$opts{grid}\" == \"F\"){p2 <- p2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())}

#ggsave(p2,file = paste(basename,\".\",\"pc\",pc_x,\"-\",pc_y,\".pdf\",sep=\"\"),width = $opts{w},height=$opts{h})
#png(file = paste(basename,\".\",\"pc\",pc_x,\"-\",pc_y,\".png\",sep=\"\"),width = $opts{w}*240,height=$opts{h}*200,res=300)
ggsave(p2,file = paste(basename,\".pdf\",sep=\"\"),width = $opts{w},height=$opts{h})
png(file = paste(basename,\".png\",sep=\"\"),width = $opts{w}*240,height=$opts{h}*200,res=300)

print(p2)
dev.off()
#svg(file = paste(basename,\".\",\"pc\",pc_x,\"-\",pc_y,\".svg\",sep=\"\"),width = $opts{w},height=$opts{h})
svg(file = paste(basename,\".svg\",sep=\"\"),width = $opts{w},height=$opts{h})
print(p2)
dev.off()

}

if (\"$opts{m}\" == \"F\") { 
 pc12\$sample <- rownames(pc12)
p2 <- ggplot(data=pc12,aes(x=$PC11, y=$PC22)) + geom_point(size=$opts{sz},color=\"#1E90FF\",shape=19) + theme_bw()+
labs(title  = paste(\"NMDS\",\"(\",\"stress:\",signif(tmpp,3),\")\",sep=\"\"))+theme(
 plot.title = element_text(hjust = 0.5))+ xlab(\"nmds1\")+ylab(\"nmds2\")
 #xlab(paste(\"PC\",pc_x,\": \",round(pc[pc_x],2),\"%\",sep=\"\"))+ylab(paste(\"PC\",pc_y,\" :  \",round(pc[pc_y],2),\"%\",sep=\"\"))
if (\"$opts{l}\" != \"F\"){p2 <- p2 + geom_text_repel(aes($PC11, $PC22, label=sample),size=$opts{cex},show.legend = FALSE) }
if ($angle == 90) {      
p2 <- p2  +    theme(axis.text.x=element_text(angle=$angle,hjust=1,vjust=0.5,size=$rlc,color=\"black\"),axis.title.x=element_text(size=$rlc))}else  if ($angle == 0) {p2 <- p2  +    theme(axis.text.x=element_text(angle=$angle,hjust=0.5,size=$rlc,color=\"black\"),axis.title.x=element_text(size=$rlc))
}else {
p2 <- p2  +    theme(axis.text.x=element_text(angle=$angle,hjust=1,vjust=1,size=$rlc,color=\"black\"),axis.title.x=element_text(size=$rlc))
}
p2 <- p2+theme( legend.title = element_text(size=$llc),legend.text = element_text(size=$llc),axis.text.y=element_text(size=$clc,color=\"black\"),
axis.title.y=element_text(size=$clc,color=\"black\"))
if (\"$opts{grid}\" == \"F\"){p2 <- p2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())}

#ggsave(p2,file = paste(basename,\".\",\"pc\",pc_x,\"-\",pc_y,\".pdf\",sep=\"\"),width = $opts{w},height=$opts{h})
ggsave(p2,file = paste(basename,\".pdf\",sep=\"\"),width = $opts{w},height=$opts{h})
#png(file = paste(basename,\".\",\"pc\",pc_x,\"-\",pc_y,\".png\",sep=\"\"),width = $opts{w}*240,height=$opts{h}*200,res=300)
png(file = paste(basename,\".png\",sep=\"\"),width = $opts{w}*240,height=$opts{h}*200,res=300)
print(p2)
dev.off()
#svg(file = paste(basename,\".\",\"pc\",pc_x,\"-\",pc_y,\".svg\",sep=\"\"),width = $opts{w},height=$opts{h})
svg(file = paste(basename,\".svg\",sep=\"\"),width = $opts{w},height=$opts{h})
print(p2)
dev.off()
}








## close warnings
sink()
";

`R --no-save < cmd.r`;
}


if ((my @file = glob  "tmp_beta/$opts{method}/*/data/distance-matrix.tsv")  && ($opts{method} ne  "distance") ) {
	#print "tmp_beta/$opts{method}/*/data/distance-matrix.tsv\n"; 
	`cp tmp_beta/$opts{method}/*/data/distance-matrix.tsv  $opts{method}.distance-matrix.txt`;
	&sub_piont("$opts{method}.distance-matrix.txt");

}
else {
	if (($opts{method} eq "distance") && ($tag eq "Lebel")){
		&sub_piont($opts{i});
	}

	elsif (($opts{method} eq "distance") && ($tag ne "Lebel")) {
	         die "you input not distance file !!!!!!!!!!!!!\n"; 
	}
	else {
            if (($opts{method} ne "distance") && ($tag ne "Lebel"))   {
	&beta($qiime2_2,$opts{i},$opts{rep},$opts{method});
	`sh beta.sh`;
	`cp tmp_beta/$opts{method}/*/data/distance-matrix.tsv  $opts{method}.distance-matrix.txt`;
	&sub_piont("$opts{method}.distance-matrix.txt");
	} 
	else {die "you should input -method distance\n";}
     }
}



