#!/usr/bin/perl -w 
####################################################################
#####Author   : GuoqiLiu                                           #
#####Date     : 2019-03-29                                         #
#####Copyright (C) 2019 Mingke Biotechnology (Hangzhou) Co., Ltd.  #
#####Contact  : liuguoqi@mingkebio.com                             #
#####Suppose  : format picrust results by picrust_pipe_v2.pl       # 
#####Platform :                                                    #
#################Ubuntu18.04 & Windows10 ,perl v5.26################

use strict;
use warnings;
use FindBin qw/$Bin/;
die "perl $0   PICRUSt2output  PICRUSt2_Results  \n" if @ARGV != 2 ; 

chomp (my $dir = `pwd`) ;
my $dir1 = $ARGV[0] ;
my $dir2 = $ARGV[1] ;
mkdir $dir2 unless (-d "$dir2");
mkdir "$dir2/PICRUSt2_COG" unless (-d "$dir2/PICRUSt2_COG"); 
mkdir "$dir2/PICRUSt2_KEGG" unless (-d "$dir2/PICRUSt2_KEGG");
mkdir "$dir2/PICRUSt2_COG/bar" unless (-d "$dir2/PICRUSt2_COG/bar");
if (-e "$dir1/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz") {
              `gzip -dc $dir1/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz > $dir2/PICRUSt2_KEGG/pred_KO_profile.xls`;
              `perl  /mnt/sdb/lgq/bin/tools/percent_abundance/getPercent_otuTable.pl $dir2/PICRUSt2_KEGG/pred_KO_profile.xls F $dir2/PICRUSt2_KEGG/pred_KO_relative_profile.xls` ;
}
else {
               print "error !!! $dir1/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz not find !\n";
}

if (-e "$dir1/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz" )  {
              open EN,"$Bin/KO2enzyme.lst";
              my %en=();
              while (<EN>) {
                      chomp;
                       my @b=split/\t/,$_;
                       my @c=split/\s+/,$b[2] ;
                       $en{"EC:".$b[2]}     = $b[1] ;
                        foreach my $i (@c) {
                       $en{"EC:".$i}     = $b[1] ;   }
              }
             close EN;
              open EC,"gzip -dc $dir1/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz|";
              chomp(my $h=<EC>);
              open ECO,">$dir2/PICRUSt2_KEGG/pred_enzyme_profile.xls";
              print ECO "$h\tDefinition\n";
              while (<EC>) {
                    chomp;
                    my @a=split/\t/,$_;
                    if (exists $en{$a[0]}) {
                               print ECO join("\t",@a),"\t$en{$a[0]}\n";
                          }
                    else {
                             print ECO join("\t",@a),"\tNone\n";
                         }
                    }
             close EC;close ECO;
            `perl  /mnt/sdb/lgq/bin/tools/percent_abundance/getPercent_otuTable.pl  $dir2/PICRUSt2_KEGG/pred_enzyme_profile.xls T $dir2/PICRUSt2_KEGG/pred_enzyme_relative_profile.xls`;
}
else {
           print "error !!!  $dir1/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz \n";
}
if (-e "$dir2/PICRUSt2_KEGG/pred_KO_profile.xls") {
            `perl $Bin/filter_matrix_KO2ko.pl $Bin/KO2allko.lst  1 $dir2/PICRUSt2_KEGG/pred_KO_profile.xls  0 > $dir2/PICRUSt2_KEGG/predictions_ko.xls.ko.tmp`; 
             `perl $Bin/splitko.pl  $dir2/PICRUSt2_KEGG/predictions_ko.xls.ko.tmp  $dir2/PICRUSt2_KEGG/predictions_ko.xls.ko.tmp.tmp`;
             `/usr/bin/python $Bin/calculate_same_rownames.py $dir2/PICRUSt2_KEGG/predictions_ko.xls.ko.tmp.tmp   $dir2/PICRUSt2_KEGG/kegg.pathway.profile.xls.2` ;
                `perl  /mnt/sdb/lgq/bin/tools/percent_abundance/getPercent_otuTable.pl $dir2/PICRUSt2_KEGG/kegg.pathway.profile.xls.2 F $dir2/PICRUSt2_KEGG/kegg.pathway.profile.xls.2.2`; 
              `$Bin/filter_matrix.backup.pl   $Bin/ko2descrition2html.lst 1 $dir2/PICRUSt2_KEGG/kegg.pathway.profile.xls.2 0  > $dir2/PICRUSt2_KEGG/pred_kegg_pathway_profile.xls` ;
               `$Bin/filter_matrix.backup.pl   $Bin/ko2descrition2html.lst 1  $dir2/PICRUSt2_KEGG/kegg.pathway.profile.xls.2.2 0 > $dir2/PICRUSt2_KEGG/pred_kegg_pathway_relative_profile.xls` ;
               #ko2Levels123.lst#
              `perl $Bin/ko2Levels123.pl $dir2/PICRUSt2_KEGG/kegg.pathway.profile.xls.2 $Bin/ko2Levels123.lst Levels3 > $dir2/PICRUSt2_KEGG/pred_kegg_pathwayL3_profile.xls.tmp`;
               `perl $Bin/ko2Levels123.pl $dir2/PICRUSt2_KEGG/kegg.pathway.profile.xls.2.2  $Bin/ko2Levels123.lst Levels3 > $dir2/PICRUSt2_KEGG/pred_kegg_pathwayL3_relative_profile.xls.tmp`;
                `perl $Bin/ko2Levels123.pl $dir2/PICRUSt2_KEGG/kegg.pathway.profile.xls.2 $Bin/ko2Levels123.lst Levels2 > $dir2/PICRUSt2_KEGG/pred_kegg_pathwayL2_profile.xls.tmp`;
               `perl $Bin/ko2Levels123.pl $dir2/PICRUSt2_KEGG/kegg.pathway.profile.xls.2.2  $Bin/ko2Levels123.lst Levels2 > $dir2/PICRUSt2_KEGG/pred_kegg_pathwayL2_relative_profile.xls.tmp`;
               
              `perl $Bin/ko2Levels123.pl $dir2/PICRUSt2_KEGG/kegg.pathway.profile.xls.2 $Bin/ko2Levels123.lst Levels1 > $dir2/PICRUSt2_KEGG/pred_kegg_pathwayL1_profile.xls.tmp`;
               `perl $Bin/ko2Levels123.pl $dir2/PICRUSt2_KEGG/kegg.pathway.profile.xls.2.2  $Bin/ko2Levels123.lst Levels1 > $dir2/PICRUSt2_KEGG/pred_kegg_pathwayL1_relative_profile.xls.tmp`;
        

                `/usr/bin/python $Bin/calculate_same_rownames.py $dir2/PICRUSt2_KEGG/pred_kegg_pathwayL3_profile.xls.tmp $dir2/PICRUSt2_KEGG/pred_kegg_pathwayL3_profile.xls `;
                 `perl $Bin/Levels123.pl $Bin/ALLLevels123.lst 1 $dir2/PICRUSt2_KEGG/pred_kegg_pathwayL3_profile.xls  0 > $dir2/PICRUSt2_KEGG/pred_kegg_pathwayL123_profile.xls`;
                  `/usr/bin/python $Bin/calculate_same_rownames.py $dir2/PICRUSt2_KEGG/pred_kegg_pathwayL2_profile.xls.tmp $dir2/PICRUSt2_KEGG/pred_kegg_pathwayL2_profile.xls `;
                   `/usr/bin/python $Bin/calculate_same_rownames.py $dir2/PICRUSt2_KEGG/pred_kegg_pathwayL1_profile.xls.tmp $dir2/PICRUSt2_KEGG/pred_kegg_pathwayL1_profile.xls `;

                   `/usr/bin/python $Bin/calculate_same_rownames_per.py $dir2/PICRUSt2_KEGG/pred_kegg_pathwayL3_relative_profile.xls.tmp $dir2/PICRUSt2_KEGG/pred_kegg_pathwayL3_relative_profile.xls`;
                     `perl $Bin/Levels123.pl $Bin/ALLLevels123.lst 1 $dir2/PICRUSt2_KEGG/pred_kegg_pathwayL3_relative_profile.xls 0 > $dir2/PICRUSt2_KEGG/pred_kegg_pathwayL123_relative_profile.xls`;
                      `/usr/bin/python $Bin/calculate_same_rownames_per.py $dir2/PICRUSt2_KEGG/pred_kegg_pathwayL2_relative_profile.xls.tmp $dir2/PICRUSt2_KEGG/pred_kegg_pathwayL2_relative_profile.xls`;
                       `/usr/bin/python $Bin/calculate_same_rownames_per.py $dir2/PICRUSt2_KEGG/pred_kegg_pathwayL1_relative_profile.xls.tmp $dir2/PICRUSt2_KEGG/pred_kegg_pathwayL1_relative_profile.xls`;           

                  
       #`/usr/bin/python $Bin/calculate_same_rownames.py  
               `rm $dir2/PICRUSt2_KEGG/predictions_ko.xls.ko.tm*;rm $dir2/PICRUSt2_KEGG/kegg.pathway.profile.xls.2*;rm $dir2/PICRUSt2_KEGG/*.tmp ;` ;
               `cd $dir2/PICRUSt2_KEGG/;bar_pie.pl  -i pred_kegg_pathwayL1_profile.xls -pie F;bar_pie.pl  -i pred_kegg_pathwayL2_profile.xls -pie F;rm ALL*.xls` 
}
else {
         print "error !!! $dir2/PICRUSt2_KEGG/pred_KO_profile.xls not find \n";
}
print "PICRUSt2 KEGG Results handle successfully !!\n";
#step2 cog ;
##
if (-e "$dir1/COG_metagenome_out/pred_metagenome_unstrat.tsv.gz") {
                open Iii,"$Bin/NOG.description.txt" ;
                my %cog=();my %hash = ();
               while (<Iii>) {
	                chomp;
		        my @a = split/\t/,$_;
                        if (scalar @a == 1) {
                          $cog{$a[0]} = "None";
                         }
                        else {
	        	$cog{$a[0]} = $a[1];#join("\t",@a[1..$#a]) ;
		}
               }
                open COG, "gzip -dc $dir1/COG_metagenome_out/pred_metagenome_unstrat.tsv.gz |";
                open COGO,">$dir2/PICRUSt2_COG/pred_cog_description_profile.xls";
                chomp(my $headd = <COG>) ;
                print COGO $headd,"\tcog_description\n";
                while (<COG>) { 
                             chomp;
                             my @b =  split/\t/,$_;
                             if (exists $cog{$b[0]}) {
                                  print COGO $_,"\t",$cog{$b[0]},"\n";
                             }
                             else {
                                  print COGO $_,"\tNone\n";
                            }
              }
            close COG;close COGO;
          `perl  /mnt/sdb/lgq/bin/tools/percent_abundance/getPercent_otuTable.pl  $dir2/PICRUSt2_COG/pred_cog_description_profile.xls T $dir2/PICRUSt2_COG/pred_cog_description_relative_profile.xls` ;
        #NOG.funccat.txt filter_matrix_cog.pl splitcog.pl 
           `perl $Bin/filter_matrix_cog.pl $Bin/NOG.funccat.txt 1 $dir2/PICRUSt2_COG/pred_cog_description_profile.xls 0 > $dir2/PICRUSt2_COG/pred_cog_description_profile.xls.tmp`;
           `perl $Bin/splitcog.pl  $dir2/PICRUSt2_COG/pred_cog_description_profile.xls.tmp $dir2/PICRUSt2_COG/pred_cog_description_profile.xls.tmp.tmp`;
           `/usr/bin/python $Bin/calculate_same_rownames.py $dir2/PICRUSt2_COG/pred_cog_description_profile.xls.tmp.tmp $dir2/PICRUSt2_COG/pred_cog_description_profile.xls.tmp.tmp2`; 
           `perl $Bin/cog.pl $dir2/PICRUSt2_COG/pred_cog_description_profile.xls.tmp.tmp2 > $dir2/PICRUSt2_COG/pred_cog_category_profile.xls`;            `perl $Bin/addanno.pl $Bin/cog.func2picture.lst  $dir2/PICRUSt2_COG/pred_cog_description_profile.xls.tmp.tmp2 > $dir2/PICRUSt2_COG/h`;  
             `perl  /mnt/sdb/lgq/bin/tools/percent_abundance/getPercent_otuTable.pl $dir2/PICRUSt2_COG/pred_cog_category_profile.xls T $dir2/PICRUSt2_COG/pred_cog_category_relative_profile.xls`;
           `rm $dir2/PICRUSt2_COG/pred_cog_description_profile.xls.tmp*`;
}

else {
             print "error!! $dir1/COG_metagenome_out/pred_metagenome_unstrat.tsv.gz \n";
}
#``
`perl $Bin/4-4.cog_bar_eachSam.pl $dir2/PICRUSt2_COG/pred_cog_category_profile.xls $dir2/PICRUSt2_COG/bar/bar`;
`gzip -dc $dir1/COG_metagenome_out/pred_metagenome_unstrat.tsv.gz > $dir2/PICRUSt2_COG/predictions_cog.xls`;
`R --no-save < $Bin/cog.box.r  $dir2/PICRUSt2_COG/h $dir2/PICRUSt2_COG/predictions_cog.xls $dir2/PICRUSt2_COG/bar/cog.box.pdf`;
`rm $dir2/PICRUSt2_COG/predictions_cog.xls ; rm $dir2/PICRUSt2_COG/h `;
print "PICRUSt2 COG Results handle successfully !!\n"
