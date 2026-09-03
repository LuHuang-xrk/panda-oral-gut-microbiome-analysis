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
##V1 version 2024-04-17 

##
use strict;
use warnings;
use Getopt::Long;
use POSIX; 

my %opts;
my $v = "v2024-08-28";

my $help;

my $method =0;
GetOptions (\%opts,"i=s","o=s","grid=s","w=f","h=f","m=s","g=s","color=s","rlc=f","point=s","method!"=>\$method,"clc=f","llc=f","lg=i","angle=i","test=s");#,"scale=s");
my $usage = <<"USAGE";
#######################################################################################################################################################
#                 Program : perl $0 
#                 Version :  $v
#######################################################################################################################################################
#                -i*  	   <str>             distance matrix file  
#                -m*       <str>             map file
#                -test     <str>             two group statistic test only can choice Wilcox-test,T-test,None,default:Wilcox-test
#                -g        <str>             group sort list
#                -color    <str>             group color 
#                -point    <str>             T/F, T display point ,F not display point ,default F   
#                -angle    <int>             xlab angle 0-360 default  0
#                -rlc     <float>            xlim size default 12
#                -clc     <float>            ylim size default 12
#                -grid     <str>             T/F, T display grid  line ,F not display display grid  line,default:F                 
#
#		 -w        <float>           the width of the figure,default:6
#		 -h	   <float>           the height of the figure,default:6
#
#                Example :                Usage:perl $0  -i braycurtis.distance-matrix.txt   -m map.txt 
#
#########################################################################################################################################################

USAGE
#die $usage if (!($opts{i}&&$opts{g}));
##                -lw    split plot in width,defalt (four numbers) :0.1:0.2:4:1
#                -lh    split plot in heigth,defalt (three numbers) :0.3:5.5:1.2


#                -method                     mean distance file prepared by customer themselves 
die $usage if (!($opts{i}));
die $usage if (!($opts{m}));
#die $usage if (!($opts{m}));
#####                -llc     <float>            legend size
#die $usage if (!($opts{i2}));
#$opts{o}=defined$opts{o}?$opts{o}:"output";
$opts{w} ||= 6;
$opts{h} ||= 6;
$opts{g} ||= "F";
$opts{test} ||= "Wilcox-test";
$opts{grid} ||= "F" ;
$opts{color} ||= "F" ;
$opts{rlc} ||= 12 ;
$opts{clc} ||= 12 ;  
$opts{llc} ||= 12 ;
$opts{angle} ||= 0 ;
$opts{point} ||= "F" ;
#$opts{color}=defined$opts{color}?$opts{color}:"default";
#$opts{row} ||= 0 ;
#$opts{col} ||= 0 ;
#$opts{scale} ||= "none" ;
my $test = $opts{test} ;
my $group = $opts{m} ;

$opts{method} ||= 0 ;

if ($method == 0) {$method = "FALSE";} else {$method = "TRUE";}


if ($method eq  "FALSE"){
open II,$opts{i};
chomp(my @tmp1 = <II>);
open III,$opts{m};
chomp(my @tmp2 = <III>);
my $tg1 =scalar @tmp1;
my $tg2 = scalar @tmp2;
if ($tg1 != $tg2) {die "YOU input is not  distance matrix ,so should need distance matrix file!!!!!!!!!!! \n";}
#while(<II>) {

}

; 
#print "$method\n";
sub get_time {
	 my $hhh=shift;
	 my $gettime=strftime("%Y-%m-%d %H:%M:%S",localtime());
	 print  "\@\@$gettime====>>>>>>>$hhh\n";
}
&get_time("$0 program start running...") ;


#print "$opts{angle} ------------------\n";

########################################################

#my @mycol = ("#df89ff","#0000cd","#00c4ff","#ff8805","#ff5584","#00bd94","#d3b3b0","#4b0082","#c0c0c0","#ffd700","#8b0000","#00ffff","#ff0000","#0000cd","#006400","#ffff00","#008080","#d8bfd8","#40e0d0","#00ff7f","#6a5acd","#adff2f","#00ffff","#ff00ff","#8b4513","#6495ed","#ff6347","#800080","#dc143c","#000000","#7fff00","#d2691e","#ff7f50","#6495ed","#fff8dc","#dc143c","#00ffff");
my @mycol = ("#df89ff","#0000cd","#00c4ff","#ff8805","#ff5584","#00bd94","#d3b3b0","#4b0082","#c0c0c0","#ffd700","#8b0000","#00ffff","#ff0000","#006400","#ffff00","#008080","#d8bfd8","#40e0d0","#00ff7f","#6a5acd","#adff2f","#ff00ff","#8b4513","#6495ed","#ff6347","#800080","#dc143c","#000000","#7fff00","#d2691e");
my $col;
my @ttmp;
#=head;
#

my %group=();





my %group2=();     
open GROUP,$opts{m} ;
<GROUP>;
while(<GROUP>){
        chomp;
        my @a=split/\t/,$_;
       push @{$group{$a[1]}},$a[0];
       $group2{$a[1]} =  0 ;  
} 
close GROUP;



my $mygroupfile;
$mygroupfile = $opts{color} ;
my $sort_group;
my @group_m ; 
if ($opts{g} eq "F") {@group_m =  sort (keys %group2); }
else { open TMM,$opts{g};chomp(@group_m =  <TMM>); }

$sort_group = join("\",\"",@group_m);
$sort_group = "\"".$sort_group."\"";

if ($mygroupfile eq "F") {
      open T,">color.txt";
      #my @ttmp;
      # my @group_m ; #=  keys %group;

      #if ($opts{g} eq "F") {@group_m =  sort (keys %group2); }
      #    else { open TMM,$opts{g};chomp(@group_m =  <TMM>); } 

      # unshift @group_m,"Between";

      # $sort_group = join("\",\"",@group_m);
      #	$sort_group = "\"".$sort_group."\"";
      foreach my $i (0..$#group_m) {
              print T "$group_m[$i]\t$mycol[$i]\n";
              #print "\"$group_m[$i]\""."="."\"$mycol[$i]\""."\n";
              push @ttmp,"\"$group_m[$i]\""."="."\"$mycol[$i]\"";
         }
      $col=join(",",@ttmp);
}
#=head;
else {
       #my @ttmp;
       open C,$mygroupfile;
       while(<C>){
              chomp; 
              my @ab=split/\t/,$_;
              push @ttmp,"\"$ab[0]\""."="."\"$ab[1]\"";
	      #$col=join(",",@ttmp);
  }
  $col=join(",",@ttmp);
}
##scale_fill_manual(name=\"\",values=c($col))
########################################################
########delete ns compare#######
sub delete_ns_compare {
         my ($groupfile,$file,$statistictest) = @_;



 #my %group=();
#my %group=();
open I,$groupfile ;
#open IO,">","$groupfile"."2"  ; 

my $hh2 = <I> ;

#print IO $hh2; 
my %renamegroup=();  
my %renamegroup2=();
my %group=();
my $tagg=0; 
while (<I>) {
	chomp;
	#$tagg++;
	my @g = split/\t/,$_;
	#$renamegroup{"a".$tagg}   =  $g[1] ; 
	#print IO "a".$tagg."\t".$g[1]."\n";  
	#$renamegroup2{$g[1]}   =  "a".$tagg ;
	$group{$g[1]} = 0;

}
close I;

open IX,">RNAMEGROUP";

foreach my $i (keys %group) {$tagg++;  $renamegroup{"a".$tagg}   =  $i ; print IX "a".$tagg."\t".$i."\n";    $renamegroup2{$i}   =  "a".$tagg ;}

#open IX,">RNAMEGROUP";
close IX ; 



open FILE,$file;
open FILE2,">",$file.".2";
my $headd=<FILE>;
print FILE2 $headd;
while (<FILE>) { 
	chomp;
	my @filee = split/\t/,$_;
	if (exists $renamegroup2{$filee[1]})  {
		print  FILE2 "$filee[0]\t$renamegroup2{$filee[1]}\t$filee[2]\n"; 
	}
	else {
		print "rename error check..............$filee[1]!!\n"; 
	}
}
close FILE;close FILE2;


#my @group2 = keys (%group) ; 
my @group2 = keys (%renamegroup) ;
#
#
my @ok=();
my @ok2=();

foreach my $i (0..$#group2) {
	foreach my $j ($i+1..$#group2) {

                my $tmp=$group2[$i].",".$group2[$j] ; 
                push @ok2,$tmp;
		push @ok,"c(".'"'.$group2[$i].'"'.",".'"'.$group2[$j].'"'.")";
		#print $group2[$i],"==>>$group2[$j]\n";
		#open O,">$group2[$i]\_\_$group2[$j]\.sbgroup.tsv";
		#open II,$ARGV[0] ;
		#my $hh = <II>;
		#print O $hh;
		#while (<II>) {
		#	chomp;
		#	my @gg= split/\t/,$_;
		#	if ($gg[1] eq $group2[$i]) {
		#		print O join("\t",@gg),"\n" ;
		#	}
		#	if ($gg[1] eq $group2[$j]) {
		#		print O join("\t",@gg),"\n" ;
		#	}
		#}
		#close O;
		#close II;
		#print $group2[$i] ,"====>>>",$group2[$j],"\n";
	}
}

#print join(",",@ok),"\n"; 
########delete ns compare#######
#sub delete_ns_compare {
#           my $file = shift; 
           my $tmp = "";
           my @ok3 = (); 
           my $ok4 = ""; 
           open Rspt,">cmd0.r" ; 
           print Rspt "
           re  <- read.table(\"$file.2\",header=T,sep=\"\\t\",check.names = F)
	   colnames(re) <- c(\"sample\",\"id\",\"value\")
           ";
           my @group = ();
           foreach my $k (@ok2) {
                my @m = split/,/,$k;
                push @group,join("_",@m);
                 my $ttt = $m[0]."_".$m[1];
                  $ttt=~s/-/_/g;
                print Rspt "
           ##re  <- read.table(\"$file\",header=T,sep=\"\\t\",check.names = F)
           tmp  <- \"\"
	   if (\"$statistictest\" == \"t.test\"){
	       tmp <- t.test(re\$value[re\$id==\"$m[0]\"],re\$value[re\$id==\"$m[1]\"],alternative = \"two.sided\")
	       }else if (\"$statistictest\" == \"wilcox.test\") {
	        tmp <- wilcox.test(re\$value[re\$id==\"$m[0]\"],re\$value[re\$id==\"$m[1]\"],alternative = \"two.sided\")
	   }#else {stop(\"Please attention statistic test only can :  T-test or Wilcox-test !!!!!!!!!!!!!!!!!!!\")}
	   $ttt  <- c(\"$m[0]\",\"$m[1]\",tmp\$p.value)
           " ;
           }
           my $tmpp = join(",",@group) ;
                $tmpp=~s/-/_/g;
           print Rspt "
           tt <- as.data.frame(rbind($tmpp)) 
           write.table(tt,file = paste(\"$file\",\".tmp.txt\",sep=\"\"),sep=\"\\t\",quote=F)
           ";
           close Rspt;   
           `/usr/bin/R --no-save<cmd0.r` ;
	   # `rm cmd.r` ;
            open II,$file.".tmp.txt" ;
            <II> ;
            while (<II>) {
                   chomp;
                   my @a = split/\s+/,$_;
                   if($a[3] <= 0.05) {
                              #push @ok3,"c(".'"'.$a[1].'"'.",".'"'.$a[2].'"'.")";
                              push @ok3,"c(".'"'.$renamegroup{$a[1]}.'"'.",".'"'.$renamegroup{$a[2]}.'"'.")";
                                 }
           }
           close II;
           $ok4 = join(",",@ok3);
            print "diff group: $ok4\n";
           return $ok4;

}


my $mp;
#Wilcox-test/T-test

open RSCRIPT1,">cmd1.r";
print RSCRIPT1 "
options(stringsAsFactors=FALSE)
library(vegan)
library(ggplot2)
#library(ggprism)
library(ggpubr)
library(reshape2)
library(plyr)
otu_raw <- read.table(file=\"$opts{i}\",sep=\"\\t\",header=T,check.names=FALSE ,row.names=1)
df <- read.table(file=\"$opts{m}\",header = T,sep = \"\\t\",check.names=F,comment.char=\"\",quote=\"\",row.names=1)
if (\"$opts{g}\" != \"F\"){
      sort.group <- readLines(\"$opts{g}\")
}

site_dis <- otu_raw

site_dis[upper.tri(site_dis)] <- 0
site_dis <- reshape2::melt(site_dis)
site_dis <- subset(site_dis, value != 0)

 group <- df
group\$variable <- rownames(group)
site_dis <- merge(group,site_dis)
re <- site_dis 
write.table(re,file = \"tmp\",sep=\"\\t\",quote=FALSE,row.names=FALSE)
";
`R --no-save<cmd1.r`;

















if ($test eq "Wilcox-test") {$test="wilcox.test"; $mp = &delete_ns_compare($group,"tmp",$test) ;}
if ($test eq "T-test") {$test="t.test";$mp = &delete_ns_compare($group,"tmp",$test) ;}
if ($test eq "None") {$mp="" ;}





open RSCRIPT,">cmd.r";
print RSCRIPT "
options(stringsAsFactors=FALSE)
library(vegan)
library(ggplot2)
#library(ggprism)
library(ggpubr)
library(reshape2)
library(plyr)
otu_raw <- read.table(file=\"$opts{i}\",sep=\"\\t\",header=T,check.names=FALSE ,row.names=1)
df <- read.table(file=\"$opts{m}\",header = T,sep = \"\\t\",check.names=F,comment.char=\"\",quote=\"\",row.names=1)
if (\"$opts{g}\" != \"F\"){
      sort.group <- readLines(\"$opts{g}\")
}

site_dis <- otu_raw

site_dis[upper.tri(site_dis)] <- 0
site_dis <- reshape2::melt(site_dis)
site_dis <- subset(site_dis, value != 0)

 group <-df
group\$variable <- rownames(group)
site_dis <- merge(group,site_dis)
re <- site_dis 




if (\"$opts{g}\" != \"F\"){
re\$Group <- factor(re\$Group,levels=c(sort.group))
}else {
re\$Group <- factor(re\$Group ,levels=c($sort_group))
}



p2 <- ggplot(re,aes(x=Group,y=value))+
stat_boxplot(geom = \"errorbar\", width=0.1,size=0.8)+#添加误差线,注意位置，放到最后则这条先不会被箱体覆盖
geom_boxplot(aes(fill=Group), 
outlier.colour=\"white\",size=0.5)
if (\"$opts{point}\" != \"F\"){
p2 <- p2+geom_jitter()}

p2 <- p2+ xlab(\"\") +ylab(\"Distance\") + scale_color_discrete(name=\"\",guide=FALSE)+
            theme(axis.text.x=element_text(angle=$opts{angle},hjust=1,vjust=1,size=$opts{rlc},color=\"black\"))+
	          theme(legend.title = element_text(size=$opts{llc}),legend.text = element_text(size=$opts{llc}),axis.text.y=element_text(size=$opts{clc},color=\"black\"))+scale_fill_manual(name=\"Group\",values=c($col),guide=FALSE) + 
#labs(title = paste(\"Bray-Curtis Anosim \",\"R=\",round(df_anosim\$statistic,3),\", \",\"p=\", round(df_anosim\$signif,3))) + 
#theme(
#    plot.title = element_text(hjust = 0.5))+
theme_bw()+theme(
 plot.title = element_text(hjust = 0.5))

if (\"$test\" == \"None\") {
p2 <- p2}else {
p2 <- p2+stat_compare_means(method = \"$test\",method.args = list(\"two.sided\"),comparisons = list($mp),label= \"p.signif\")}



if ($opts{angle} == 90) {      
p2 <- p2  +    theme(axis.text.x=element_text(angle=$opts{angle},hjust=1,vjust=0.5,size=$opts{rlc},color=\"black\"))}else if ($opts{angle} == 0) {p2 <- p2  +    theme(axis.text.x=element_text(angle=$opts{angle},hjust=0.5,size=$opts{rlc},color=\"black\"))
}else {
p2 <- p2  +    theme(axis.text.x=element_text(angle=$opts{angle},hjust=1,vjust=1,size=$opts{rlc},color=\"black\"))
}




p2 <- p2+theme( legend.title = element_text(size=$opts{llc}),legend.text = element_text(size=$opts{llc}),axis.text.y=element_text(size=$opts{clc},color=\"black\"),axis.title=element_text(size=$opts{llc}))
 
 
if (\"$opts{grid}\" == \"F\"){p2 <- p2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())}


ggsave(p2,file = \"distance_box.pdf\",width = $opts{w},height=$opts{h})
png(file = \"distance_box.png\",width = $opts{w}*240,height=$opts{h}*200,res=300)
print(p2)
dev.off()
svg(file = \"distance_box.svg\",width = $opts{w},height=$opts{h})
print(p2)
dev.off()


";
`R --no-save < cmd.r`;
#`rm cmd.r` ; 
#
#












if ($opts{m} ne  "F") {

my %renamegroup=();
open RR,"RNAMEGROUP";
while (<RR>) {
	chomp;
	my @ac=split/\t/,$_;
	$renamegroup{$ac[0]} = $ac[1];
}close RR;

	my @file =glob ("*.tmp.txt") ;
	open OUT,">distance-test.txt";
	print OUT "#Compare\tGroup1\tGroup2\tpvalue\n";
	foreach my $i (@file) {
	       my @split=split/\./,$i;
	        open II,$i;
		<II>;
		while (<II>) {
			chomp;
				my @aa = split/\s+/,$_;   
		print OUT "$renamegroup{$aa[1]}\_$renamegroup{$aa[2]}\t$renamegroup{$aa[1]}\t$renamegroup{$aa[2]}\t$aa[3]\n";
		#	print OUT "$_\n";
		}
		close II;
	}
	close OUT;
	`rm *.tmp.txt`;
}
`rm *.r tmp *.2 RNAMEGROUP`;
