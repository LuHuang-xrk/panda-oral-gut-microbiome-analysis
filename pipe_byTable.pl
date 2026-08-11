#!/usr/bin/perl -w
use strict;
use warnings;
use FindBin qw($Bin);

my $opt;
my ($table,$rep,$map,$raw,$clean,$G3flag,$flag,$Rnum,$outDir,$typeFlag);
while($opt = shift){
	if($opt eq "-i"){
		$table = shift;
	}elsif($opt eq "-r"){
		$rep = shift;
	}elsif($opt eq "-m"){
		$map = shift;
	}elsif($opt eq "-rs"){
		$raw = shift;
	}elsif($opt eq "-cs"){
		$clean = shift;
#	}elsif($opt eq "-gf"){
#		$Gflag = shift;
	}elsif($opt eq "-g3"){
		$G3flag = shift;
	}elsif($opt eq "-af"){
		$flag = shift;
	}elsif($opt eq "-t"){
		$typeFlag = shift;
	}elsif($opt eq "-n"){
		$Rnum = shift;
	}elsif($opt eq "-o"){
		$outDir = shift;
	}elsif($opt eq "-h"){
		&usage;exit;
	}
}

unless($table and $rep and $raw and $clean and $outDir){
	&usage;
	exit;
}

if($G3flag and $G3flag eq "T"){
	unless($map){
		print STDERR "must provide map file for group analysis!!!\n";
		&usage;
		exit;
	}
}

$flag = "ASV" unless($flag);
$G3flag = "F" unless($G3flag);
$typeFlag = "T" unless($typeFlag);
my $pre = lc($flag);

$outDir =~ s/\/$//;
`mkdir $outDir` unless(-e $outDir);
`mkdir $outDir/0.Data` unless(-e "$outDir/0.Data");
`mkdir $outDir/1.Taxa` unless(-e "$outDir/1.Taxa");
`mkdir $outDir/1.Taxa/original` unless(-e "$outDir/1.Taxa/original");
`mkdir $outDir/1.Taxa/normalize` unless(-e "$outDir/1.Taxa/normalize");
`mkdir $outDir/2.Alpha/` unless(-e "$outDir/2.Alpha");

my ($samNumber,$groupNumber);
#`mkdir $outDir/2.Alpha/`
if($map){
	my %samNum;
	open OUS,"> $outDir/sub.txt" or die "$!\n";
	open INM,$map or die "$!\n";
	print OUS "sub";
	while(<INM>){
		chomp;
		my @temp = split;
		next if($. == 1 or /^#/);
		$samNumber++;
		$samNum{$temp[1]}++;
		print OUS "\t$temp[0]";
	}
	print OUS "\n";
	close INM;
	close OUS;
	foreach my $g(sort keys %samNum){
		$groupNumber++;
		if($samNum{$g} < 3 and $G3flag eq "T"){
			print STDERR "error: group $g has less than 3 samples than can not be used for statistical analysis, please set -g3 F\n";
			exit();
		}
	}
	`perl $Bin/bin/subTable.pl $table $outDir/sub.txt $outDir/$pre\_taxa_table T T`;
	`perl $Bin/bin/OTU_table_sortBySam.pl $outDir/$pre\_taxa_table.sub.xls $map T $outDir/1.Taxa/original/$pre\_taxa_table.xls`;
	`perl $Bin/bin/checkID.pl $outDir/1.Taxa/original/$pre\_taxa_table.xls $flag\-ID`;
#	`mv $outDir/$pre\_taxa_table.sub.xls $outDir/1.Taxa/original/$pre\_taxa_table.xls`;
	`perl $Bin/bin/getSeq_byName.pl $outDir/1.Taxa/original/$pre\_taxa_table.xls $rep $outDir/temp.fasta`;
	`perl $Bin/bin/group_sort.pl $map $outDir/group.sort`;
}else{
#	`perl $Bin/bin/OTU_table_sortBySam.pl $table $map T $outDir/1.Taxa/original/$pre\_taxa_table.xls`;
	`cp $table $outDir/1.Taxa/original/$pre\_taxa_table.xls`;
	`perl $Bin/bin/checkID.pl $outDir/1.Taxa/original/$pre\_taxa_table.xls $flag\-ID`;
	`cp $rep $outDir/temp.fasta`;
}

open INF,"$outDir/temp.fasta" or die "$!\n";
open OUTF,"> $outDir/1.Taxa/original/$pre\_rep.fasta" or die "$!\n";
my (%name2seq,$tempNum);
while(<INF>){
	chomp;
	my @temp = split;
	if($temp[0] =~ /^>/){
		if($temp[0] =~ /(\d+)/){
			$name2seq{$1} = "$temp[0]\n";
			$tempNum = $1;
		}else{
			die "error: $rep format not be OTUNum or ASVNum >> $_\n";
		}
	}else{
		$name2seq{$tempNum} .= "$temp[0]\n";
	}
}
close INF;
foreach my $n(sort {$a <=> $b} keys %name2seq){
	print OUTF $name2seq{$n};
}
close OUTF;

#`cp $raw $clean $outDir`;
`cp $map $outDir` if($map);
#$raw = (split /\//,$raw)[-1];
#$clean = (split /\//,$clean)[-1];
$map = (split /\//,$map)[-1] if($map);

open OUP,"> $outDir/pipe.sh" or die "$!\n";

`perl $Bin/bin/getSeqNum_inTable.pl $outDir/1.Taxa/original/$pre\_taxa_table.xls T $outDir/seqNum.txt`;
`perl $Bin/bin/dataStat.pl $raw $clean $outDir/seqNum.txt $outDir/1.Taxa/original/$pre\_taxa_table.xls $outDir/0.Data/data.stat.xls`;

unless($samNumber){
	chomp(my $temp_mess = `wc -l $outDir/seqNum.txt`);
	$samNumber = (split /\s+/,$temp_mess)[0] - 1;
}

if($groupNumber){
	print "$samNumber sample,$groupNumber group...\n";
}else{
	print "$samNumber sample,no group\n";
}

my %num2level = (1 => 'domain',	2 => 'phylum', 3 => 'class', 4 => 'order', 5 => 'family', 6 => 'genus', 7 => 'specie');

print OUP "###taxonomy...\n";
print OUP "biom convert --to-json --table-type \"OTU table\" --process-obs-metadata taxonomy -i 1.Taxa/original/$pre\_taxa_table.xls -o 1.Taxa/original/$pre\_taxa_table.biom\n";
print OUP "summarize_taxa.py -i 1.Taxa/original/$pre\_taxa_table.biom -o 1.Taxa/original/tax_summary_a -L 1,2,3,4,5,6,7 -a\n";
print OUP "summarize_taxa.py -i 1.Taxa/original/$pre\_taxa_table.biom -o 1.Taxa/original/tax_summary_r -L 1,2,3,4,5,6,7\n";
for(my $i = 1; $i <= 7;$i++){
	print OUP "$Bin/bin/sum_tax.pl -i 1.Taxa/original/tax_summary_a/$pre\_taxa_table_L$i.txt -o 1.Taxa/original/tax_summary_a/$num2level{$i}.xls ; perl $Bin/bin/checkID.pl 1.Taxa/original/tax_summary_a/$num2level{$i}.xls $num2level{$i} ; $Bin/bin/sum_tax.pl -i 1.Taxa/original/tax_summary_r/$pre\_taxa_table_L$i.txt -o 1.Taxa/original/tax_summary_r/$num2level{$i}.percent.xls ; perl $Bin/bin/checkID.pl 1.Taxa/original/tax_summary_r/$num2level{$i}.percent.xls $num2level{$i}\n";
}
print OUP "perl $Bin/bin/remove_tax.pl 1.Taxa/original/$pre\_taxa_table.xls 1.Taxa/original/$pre\_table.xls\n";
print OUP "biom convert --to-json --table-type \"OTU table\" -i 1.Taxa/original/$pre\_table.xls -o 1.Taxa/original/$pre\_table.biom\n\n";

print OUP "\n##normalize...\n";
if($Rnum){
	print OUP "single_rarefaction.py -i 1.Taxa/original/$pre\_taxa_table.biom -o 1.Taxa/normalize/$pre\_taxa_table.biom -d $Rnum\n";
	print OUP "biom convert --to-tsv --table-type \"OTU table\" --header-key taxonomy -i 1.Taxa/normalize/$pre\_taxa_table.biom -o 1.Taxa/normalize/$pre\_taxa_table.txt\n";
	print OUP "perl $Bin/bin/txt2xls.pl 1.Taxa/normalize/$pre\_taxa_table.txt 1.Taxa/normalize/$pre\_taxa_table.xls ; perl $Bin/bin/checkID.pl 1.Taxa/normalize/$pre\_taxa_table.xls $flag\-ID ; rm 1.Taxa/normalize/$pre\_taxa_table.txt\n";
	print OUP "perl $Bin/bin/remove_tax.pl 1.Taxa/normalize/$pre\_taxa_table.xls 1.Taxa/normalize/$pre\_table.xls\n";
	print OUP "biom convert --to-tsv --table-type \"OTU table\" -i 1.Taxa/normalize/$pre\_table.xls -o 1.Taxa/normalize/$pre\_table.biom\n";
	if($map){
		print OUP "perl $Bin/bin/filterMap_byTable.pl $map 1.Taxa/normalize/$pre\_table.xls F $map.filter ; mv $map.filter $map\n";
		print OUP "perl $Bin/bin/group_sort.pl $map $outDir/group.sort\n";
		if($G3flag eq "T"){
			print OUP "perl $Bin/bin/checkMap.pl $map\n";
			print OUP "if [ -e $map.error ];\n";
			print OUP "\tthen\n";
			print OUP "\t\techo \"error: At least one group has less than 3 samples than can not be used for statistical analysis, please set -g3 F\"\n";
			print OUP "\texit\nfi\n";
		}
	}
}else{
	chomp(my $temp_mess = `head -2 $outDir/seqNum.txt | tail -1`);
	my $temp_num = (split /\s+/,$temp_mess)[-1];
	print OUP "single_rarefaction.py -i 1.Taxa/original/$pre\_taxa_table.biom -o 1.Taxa/normalize/$pre\_taxa_table.biom -d $temp_num\n";
	print OUP "biom convert --to-tsv --table-type \"OTU table\" --header-key taxonomy -i 1.Taxa/normalize/$pre\_taxa_table.biom -o 1.Taxa/normalize/$pre\_taxa_table.txt\n";
	print OUP "perl $Bin/bin/txt2xls.pl 1.Taxa/normalize/$pre\_taxa_table.txt 1.Taxa/normalize/$pre\_taxa_table.xls ; perl $Bin/bin/checkID.pl 1.Taxa/normalize/$pre\_taxa_table.xls $flag\-ID ; rm 1.Taxa/normalize/$pre\_taxa_table.txt\n";
	print OUP "perl $Bin/bin/remove_tax.pl 1.Taxa/normalize/$pre\_taxa_table.xls 1.Taxa/normalize/$pre\_table.xls\n";
	print OUP "biom convert --to-json --table-type \"OTU table\" -i 1.Taxa/normalize/$pre\_table.xls -o 1.Taxa/normalize/$pre\_table.biom\n";
}
print OUP "perl $Bin/bin/getSeq_byName.pl 1.Taxa/normalize/$pre\_table.xls 1.Taxa/original/$pre\_rep.fasta 1.Taxa/normalize/$pre\_rep.fasta\n";
print OUP "summarize_taxa.py -i 1.Taxa/normalize/$pre\_taxa_table.biom -o 1.Taxa/normalize/tax_summary_a -L 1,2,3,4,5,6,7 -a\n";
print OUP "summarize_taxa.py -i 1.Taxa/normalize/$pre\_taxa_table.biom -o 1.Taxa/normalize/tax_summary_r -L 1,2,3,4,5,6,7\n";
for(my $i = 1;$i <= 7;$i++){
	print OUP "$Bin/bin/sum_tax.pl -i 1.Taxa/normalize/tax_summary_a/$pre\_taxa_table_L$i.txt -o 1.Taxa/normalize/tax_summary_a/$num2level{$i}.xls ; perl $Bin/bin/checkID.pl 1.Taxa/normalize/tax_summary_a/$num2level{$i}.xls $num2level{$i} ; $Bin/bin/sum_tax.pl -i 1.Taxa/normalize/tax_summary_r/$pre\_taxa_table_L$i.txt -o 1.Taxa/normalize/tax_summary_r/$num2level{$i}.percent.xls ; perl $Bin/bin/checkID.pl 1.Taxa/normalize/tax_summary_r/$num2level{$i}.percent.xls $num2level{$i}\n";
}

my $samSortFile;
if($map){
	$samSortFile = $map;
}else{
	print OUP "perl $Bin/bin/createSortByTable.pl 1.Taxa/normalize/$pre\_table.xls sam.sort\n";
	$samSortFile = "sam.sort";
}

print OUP "\n##phylogeny tree...\n";
print OUP "cd 1.Taxa/original ; perl /mnt/sdb/lgq/bin/tools/clound/11.pcoa/repfasta2tre.pl $pre\_rep.fasta ; mv tree.nwk phy.tre ; cd ../../\n";
print OUP "cd 1.Taxa/normalize ; perl /mnt/sdb/lgq/bin/tools/clound/11.pcoa/repfasta2tre.pl $pre\_rep.fasta ; mv tree.nwk phy.tre ; cd ../../\n";

print OUP "\n\n###alpha diversity...\n";
`mkdir $outDir/2.Alpha/Alpha` unless(-e "$outDir/2.Alpha/Alpha");
`mkdir $outDir/2.Alpha/2-1.Alpha-index` unless(-e "$outDir/2.Alpha/2-1.Alpha-index");
if($G3flag eq "T"){`mkdir $outDir/2.Alpha/2-2.Boxplot` unless(-e "$outDir/2.Alpha/2-2.Boxplot")};
`mkdir $outDir/2.Alpha/2-3.Rarefaction` unless(-e "$outDir/2.Alpha/2-3.Rarefaction");
`mkdir $outDir/2.Alpha/2-4.Rank-abundance` unless(-e "$outDir/2.Alpha/2-4.Rank-abundance");
`mkdir $outDir/2.Alpha/2-5.Core-pan` unless(-e "$outDir/2.Alpha/2-5.Core-pan");

print OUP "\n##alpha index...\n";
my ($temp_w,$temp_h,$temp_size);
if($G3flag eq "T"){
	$temp_w = $groupNumber + 1;
	$temp_w += 0.6 if($groupNumber <= 3);
	$temp_h = $temp_w + 0.5;
	print OUP "cd 2.Alpha/Alpha ; perl /mnt/sdb/lgq/bin/tools/clound/07.alpha/alpha2024-09-23.pl -i ../../1.Taxa/normalize/$pre\_table.xls -m ../../$map -rep ../../1.Taxa/normalize/phy.tre -g ../../group.sort -w $temp_w -h $temp_h ; cd ../../\n";
	print OUP "rm 2.Alpha/Alpha/Reads.* ; cp 2.Alpha/Alpha/*pdf 2.Alpha/Alpha/*png 2.Alpha/Alpha/*svg 2.Alpha/Alpha/alpha-test.txt 2.Alpha/2-2.Boxplot\n";
}else{
	print OUP "cd 2.Alpha/Alpha ; perl /mnt/sdb/lgq/bin/tools/clound/07.alpha/alpha2024-09-23.pl -i ../../1.Taxa/normalize/$pre\_table.xls -rep ../../1.Taxa/normalize/phy.tre ; cd ../../\n";
}
print OUP "perl $Bin/bin/sortByTable.pl 2.Alpha/Alpha/alpha.txt 1.Taxa/normalize/$pre\_table.xls 2.Alpha/2-1.Alpha-index/alpha.txt\n";

print OUP "\n##rarefaction...\n";
if($samNumber >= 30){
	$temp_size = 10;
}else{
	$temp_size = 10;
}
print OUP "cd 2.Alpha/2-3.Rarefaction ; perl /mnt/sdb/lgq/bin/tools/clound/17.rarefaction/rarefaction2024-09-02.pl -i ../../1.Taxa/normalize/$pre\_table.xls -w $temp_size -h $temp_size ; rm cmd.r cmd2.r mothur.* otu_seqids.* tmp_alpha -r ; cd ../../\n";
if($samNumber <= 10){
	$temp_w = 5;$temp_h = 4;
}elsif($samNumber <= 30){
	$temp_w = 6;$temp_h = 5;
}elsif($samNumber <= 50){
	$temp_w = 7;$temp_h = 6;
}else{
	$temp_w = 7 * $samNumber/50;$temp_h = 6 * $samNumber/50;
}

print OUP "\n##rank abundance...\n";
print OUP "cd 2.Alpha/2-4.Rank-abundance ; perl /mnt/sdb/lgq/bin/tools/clound/21.rank_abundance/rank_abundance2024-09-05.pl -i ../../1.Taxa/normalize/$pre\_table.xls -w $temp_w -h $temp_h ; rm cmd.r ; cd ../../\n";
if($samNumber <= 40){
	$temp_w = 6;
}else{
	$temp_w = 6 * $samNumber/40;
}

print OUP "\n##core-pan...\n";
print OUP "cd 2.Alpha/2-5.Core-pan ; perl /mnt/sdb/lgq/bin/tools/clound/08.core-pan/core-pan2024-06-06.pl -i ../../1.Taxa/normalize/$pre\_table.xls -angle 90 -w $temp_w ; rm okk* color.txt cmd.r ; cd ../../\n";


print OUP "\n\n###Community\n";
`mkdir $outDir/3.Community` unless(-e "$outDir/3.Community");
`mkdir $outDir/3.Community/3-1.Taxonomy-stat` unless(-e "$outDir/3.Community/3-1.Taxonomy-stat");
`mkdir $outDir/3.Community/3-2.Venn` unless(-e "$outDir/3.Community/3-2.Venn");
`mkdir $outDir/3.Community/3-3.Bar` unless(-e "$outDir/3.Community/3-3.Bar");
`mkdir $outDir/3.Community/3-4.Bubble` unless(-e "$outDir/3.Community/3-4.Bubble");
`mkdir $outDir/3.Community/3-5.Heatmap` unless(-e "$outDir/3.Community/3-5.Heatmap");
`mkdir $outDir/3.Community/3-6.TreeBar` unless(-e "$outDir/3.Community/3-6.TreeBar");
`mkdir $outDir/3.Community/3-7.topGenus` unless(-e "$outDir/3.Community/3-7.topGenus");
if($map){
	`mkdir $outDir/3.Community/3-2.Venn.Group` unless(-e "$outDir/3.Community/3-2.Venn.Group");
	`mkdir $outDir/3.Community/3-3.Bar.Group` unless(-e "$outDir/3.Community/3-3.Bar.Group");
	`mkdir $outDir/3.Community/3-4.Bubble.Group` unless(-e "$outDir/3.Community/3-4.Bubble.Group");
	`mkdir $outDir/3.Community/3-5.Heatmap.Group` unless(-e "$outDir/3.Community/3-5.Heatmap.Group");
	`mkdir $outDir/3.Community/3-6.TreeBar.Group` unless(-e "$outDir/3.Community/3-6.TreeBar.Group");
}

print OUP "\n##taxonomy stat...\n";
print OUP "perl $Bin/bin/taxonomy_stat.pl 1.Taxa/normalize/$pre\_taxa_table.xls $flag 3.Community/3-1.Taxonomy-stat/tax.num.xls\n";

print OUP "\n##venn or flower...\n";
if($samNumber <= 40){
	$temp_w = 6;$temp_h = 6.5;
}else{
	$temp_w = 6 * $samNumber/40;
	$temp_w = $temp_w * 0.8 if($temp_w > 10);
	$temp_h = $temp_w + 0.5;
}
print OUP "cd 3.Community/3-2.Venn ; perl /mnt/sdb/lgq/bin/tools/clound/13.venn/venn2024-09-12.pl -i ../../1.Taxa/normalize/$pre\_table.xls -w $temp_w -h $temp_h ; rm color.txt ; cd ../../\n";
if($map){
	die "error: group number error for $map!!!\n" unless($groupNumber);
	if($groupNumber <= 40){
		$temp_w = 6;$temp_h = 6.5;
	}else{
		$temp_w = 6 * $groupNumber/40;
		$temp_w = $temp_w * 0.8 if($temp_w > 10);
		$temp_h = $temp_w + 0.5;
	}
	print OUP "cd 3.Community/3-2.Venn.Group ; perl /mnt/sdb/lgq/bin/tools/clound/13.venn/venn2024-09-12.pl -i ../../1.Taxa/normalize/$pre\_table.xls -m ../../$map -g ../../group.sort -w $temp_w -h $temp_h ; rm color.txt tmp.mean.txt VennDiagram*log ; cd ../../\n";
}
#print OUP "cd 3.Community/3-2.Venn.Group ; perl /mnt/sdb/lgq/bin/tools/clound/13.venn/venn2024-06-25.pl -i ../../1.Taxa/normalize/$pre\_table.xls -m ../../$map -g ../../group.sort ; cd ../../\n" if($map);

print OUP "\n##barplot...\n";
print OUP "cp 1.Taxa/normalize/tax_summary_a/*xls 3.Community/3-3.Bar ; cp 1.Taxa/normalize/tax_summary_a/*xls 3.Community/3-4.Bubble ; cp 1.Taxa/normalize/tax_summary_a/*xls 3.Community/3-5.Heatmap ; cp 1.Taxa/normalize/tax_summary_a/*xls 3.Community/3-6.TreeBar\n\n";
my $temp_legend;
if($samNumber <= 30){
	$temp_w = 7;$temp_legend = 3;
}elsif($samNumber <= 50){
	$temp_w = 8;$temp_legend = 4;
}else{
	$temp_w = 8 * $samNumber/50 + 1;$temp_legend = 4 * $samNumber/50 + 1;$temp_legend = int($temp_legend);
}
print OUP "cd 3.Community/3-3.Bar ; perl /mnt/sdb/lgq/bin/tools/clound/02.bar/bar_pie2024-08-02.pl -d ./ -i phylum.xls -filter top -top 30 -t T -b_blas 2 -pie F -tp 1 -b_ncol $temp_legend -tv 0.5 -b_lcex 1 -b_w $temp_w -b_h 7 ; mv bar_link.pdf bar.phylum.pdf ; mv bar_link.png bar.phylum.png ; mv bar_link.svg bar.phylum.svg\n";
print OUP "perl /mnt/sdb/lgq/bin/tools/clound/02.bar/bar_pie2024-08-02.pl -d ./ -i class.xls -filter top -top 30 -t T -b_blas 2 -pie F -tp 1 -b_ncol $temp_legend -tv 0.5 -b_lcex 1 -b_w $temp_w -b_h 7 ; mv bar_link.pdf bar.class.pdf ; mv bar_link.png bar.class.png ; mv bar_link.svg bar.class.svg\n";
$temp_w++ if($temp_w <= 8);
$temp_legend-- if($temp_legend > 4);
print OUP "perl /mnt/sdb/lgq/bin/tools/clound/02.bar/bar_pie2024-08-02.pl -d ./ -i order.xls -filter top -top 30 -t T -b_blas 2 -pie F -tp 1 -b_ncol $temp_legend -tv 0.5 -b_lcex 0.9 -b_w $temp_w -b_h 7 ; mv bar_link.pdf bar.order.pdf ; mv bar_link.png bar.order.png ; mv bar_link.svg bar.order.svg\n";
print OUP "perl /mnt/sdb/lgq/bin/tools/clound/02.bar/bar_pie2024-08-02.pl -d ./ -i family.xls -filter top -top 30 -t T -b_blas 2 -pie F -tp 1 -b_ncol $temp_legend -tv 0.5 -b_lcex 0.9 -b_w $temp_w -b_h 7 ; mv bar_link.pdf bar.family.pdf ; mv bar_link.png bar.family.png ; mv bar_link.svg bar.family.svg\n";
#$temp_legend-- if($temp_w >= 5);
print OUP "perl /mnt/sdb/lgq/bin/tools/clound/02.bar/bar_pie2024-08-02.pl -d ./ -i genus.xls -filter top -top 30 -t T -b_blas 2 -pie F -tp 1 -b_ncol $temp_legend -tv 0.5 -b_lcex 0.8 -b_w $temp_w -b_h 7 ; mv bar_link.pdf bar.genus.pdf ; mv bar_link.png bar.genus.png ; mv bar_link.svg bar.genus.svg ; rm cmd.r out.txt bar.pdf bar.svg bar.png all.bar ; cd ../../\n";
print OUP "#cd 3.Community/3-3.Bar ; perl /mnt/sdb/lgq/bin/tools/clound/02.bar/bar_pie2024-08-02.pl -d ./ -i specie.xls -filter top -top 30 -t T -b_blas 2 -pie F -tp 1 -b_ncol $temp_legend -tv 0.5 -b_lcex 0.8 -b_w $temp_w -b_h 7 ; mv bar_link.pdf bar.specie.pdf ; mv bar_link.png bar.specie.png ; mv bar_link.svg bar.specie.svg ; rm cmd.r out.txt bar.pdf bar.svg bar.png all.bar ; cd ../../\n\n";

print OUP "\n##bubble...\n";
if($samNumber <= 10){
	$temp_w = 6;
}elsif($samNumber <= 30){
	$temp_w = 8;
}elsif($samNumber <= 50){
	$temp_w = 12;
}else{
	$temp_w = 12 * $samNumber/50;$temp_w = int($temp_w);
}
my $familyGenus_w = $temp_w;
$familyGenus_w++;
print OUP "cd 3.Community/3-4.Bubble ; perl /mnt/sdb/lgq/bin/tools/clound/05.bubble/bubble2024-05-09.pl -i phylum.xls -angle 60 -w $temp_w -clc 11 ; mv bubble.pdf bubble.phylum.pdf ; mv bubble.png bubble.phylum.png ; mv bubble.svg bubble.phylum.svg\n";
print OUP "perl /mnt/sdb/lgq/bin/tools/clound/05.bubble/bubble2024-05-09.pl -i class.xls -angle 60 -w $temp_w -clc 11 ; mv bubble.pdf bubble.class.pdf ; mv bubble.png bubble.class.png ; mv bubble.svg bubble.class.svg\n";
print OUP "perl /mnt/sdb/lgq/bin/tools/clound/05.bubble/bubble2024-05-09.pl -i order.xls -angle 60 -w $temp_w -clc 11 ; mv bubble.pdf bubble.order.pdf ; mv bubble.png bubble.order.png ; mv bubble.svg bubble.order.svg\n";
print OUP "perl /mnt/sdb/lgq/bin/tools/clound/05.bubble/bubble2024-05-09.pl -i family.xls -angle 60 -w $familyGenus_w -clc 10 ; mv bubble.pdf bubble.family.pdf ; mv bubble.png bubble.family.png ; mv bubble.svg bubble.family.svg\n";
print OUP "perl /mnt/sdb/lgq/bin/tools/clound/05.bubble/bubble2024-05-09.pl -i genus.xls -angle 60 -w $familyGenus_w -clc 10 ; mv bubble.pdf bubble.genus.pdf ; mv bubble.png bubble.genus.png ; mv bubble.svg bubble.genus.svg ; rm cmd.r ; cd ../../\n";
print OUP "#cd 3.Community/3-4.Bubble ; perl /mnt/sdb/lgq/bin/tools/clound/05.bubble/bubble2024-05-09.pl -i specie.xls -angle 60 -w $familyGenus_w -clc 10 ; mv bubble.pdf bubble.specie.pdf ; mv bubble.png bubble.specie.png ; mv bubble.svg bubble.specie.svg ; rm cmd.r ; cd ../../\n\n";

print OUP "\n##heatmap...\n";
if($samNumber <= 10){
	$temp_w = 6;
}elsif($samNumber <= 30){
	$temp_w = 1.1 * $samNumber/10 + 5;#$temp_w = int($temp_w);
}else{
	$temp_w = 1.1 * $samNumber/10 + 6;#$temp_w = int($temp_w);
}
print OUP "cd 3.Community/3-5.Heatmap ; perl /mnt/sdb/lgq/bin/tools/clound/01.heatmap/plot-heatmap2024-04-17_3.pl -d ./ -i phylum.xls -o heatmap.phylum -ct 0 -rt 0 -lg log2 -cs 2 -rtop 50 -slas 2 -rlc 0.9 -w $temp_w -marble 4-1-0-6\n";
print OUP "perl /mnt/sdb/lgq/bin/tools/clound/01.heatmap/plot-heatmap2024-04-17_3.pl -d ./ -i class.xls -o heatmap.class -ct 0 -rt 0 -lg log2 -cs 2 -rtop 50 -slas 2 -rlc 0.9 -w $temp_w -marble 4-1-0-10\n";
$temp_w++;
print OUP "perl /mnt/sdb/lgq/bin/tools/clound/01.heatmap/plot-heatmap2024-04-17_3.pl -d ./ -i order.xls -o heatmap.order -ct 0 -rt 0 -lg log2 -cs 2 -rtop 50 -slas 2 -rlc 0.9 -w $temp_w -marble 4-1-0-10\n";
print OUP "perl /mnt/sdb/lgq/bin/tools/clound/01.heatmap/plot-heatmap2024-04-17_3.pl -d ./ -i family.xls -o heatmap.family -ct 0 -rt 0 -lg log2 -cs 2 -rtop 50 -slas 2 -rlc 0.8 -w $temp_w -marble 4-1-0-14\n";
print OUP "perl /mnt/sdb/lgq/bin/tools/clound/01.heatmap/plot-heatmap2024-04-17_3.pl -d ./ -i genus.xls -o heatmap.genus -ct 0 -rt 0 -lg log2 -cs 2 -rtop 50 -slas 2 -rlc 0.9 -w $temp_w -marble 4-1-0-14 ; rm cmd.r cmd.r.log final.txt warnings.log ; cd ../../\n";
print OUP "#cd 3.Community/3-5.Heatmap ; perl /mnt/sdb/lgq/bin/tools/clound/01.heatmap/plot-heatmap2024-04-17_3.pl -d ./ -i specie.xls -o heatmap.specie -ct 0 -rt 0 -lg log2 -cs 2 -rtop 50 -slas 2 -rlc 0.9 -w $temp_w -marble 4-1-0-15 ; rm cmd.r cmd.r.log final.txt warnings.log ; cd ../../\n\n";

print OUP "\n##treebar...\n";
if($samNumber <= 10){
	$temp_h = 5;
}else{
	$temp_h = $samNumber/10 + 5;
}
print OUP "cd 3.Community/3-6.TreeBar ; perl /mnt/sdb/lgq/bin/tools/clound/20.treebar/treebar2024-09-05.pl -otu ../../1.Taxa/normalize/$pre\_table.xls -tax phylum.xls -filter top -scex 1.2 -lcex 1.2 -h $temp_h -lw 1.3-1.8-1.3 ; mv treebar.pdf treebar.phylum.pdf ; mv treebar.png treebar.phylum.png ; mv treebar.svg treebar.phylum.svg\n";
print OUP "perl /mnt/sdb/lgq/bin/tools/clound/20.treebar/treebar2024-09-05.pl -otu ../../1.Taxa/normalize/$pre\_table.xls -tax class.xls -filter top -scex 1.2 -lcex 1.2 -h $temp_h -lw 1.3-1.8-1.3 ; mv treebar.pdf treebar.class.pdf ; mv treebar.png treebar.class.png ; mv treebar.svg treebar.class.svg\n";
print OUP "perl /mnt/sdb/lgq/bin/tools/clound/20.treebar/treebar2024-09-05.pl -otu ../../1.Taxa/normalize/$pre\_table.xls -tax order.xls -filter top -scex 1.2 -lcex 1.2 -w 12 -h $temp_h -lw 1.3-1.8-1.3 ; mv treebar.pdf treebar.order.pdf ; mv treebar.png treebar.order.png ; mv treebar.svg treebar.order.svg\n";
print OUP "perl /mnt/sdb/lgq/bin/tools/clound/20.treebar/treebar2024-09-05.pl -otu ../../1.Taxa/normalize/$pre\_table.xls -tax family.xls -filter top -scex 1.2 -lcex 1.1 -w 13 -h $temp_h -lw 1.3-1.8-1.3 ; mv treebar.pdf treebar.family.pdf ; mv treebar.png treebar.family.png ; mv treebar.svg treebar.family.svg\n";
print OUP "perl /mnt/sdb/lgq/bin/tools/clound/20.treebar/treebar2024-09-05.pl -otu ../../1.Taxa/normalize/$pre\_table.xls -tax genus.xls -filter top -scex 1.2 -lcex 1.1 -w 13 -h $temp_h -lw 1.3-1.8-1.3 ; mv treebar.pdf treebar.genus.pdf ; mv treebar.png treebar.genus.png ; mv treebar.svg treebar.genus.svg ; rm cmd.r all.bar otu_hc.tre otu_bray.xls ; cd ../../\n";
print OUP "#cd 3.Community/3-6.TreeBar ; perl /mnt/sdb/lgq/bin/tools/clound/20.treebar/treebar2024-09-05.pl -otu ../../1.Taxa/normalize/$pre\_table.xls -tax specie.xls -filter top -scex 1.2 -lcex 1.1 -w 13 -h $temp_h -lw 1.3-1.8-1.3 ; mv treebar.pdf treebar.specie.pdf ; mv treebar.png treebar.specie.png ; mv treebar.svg treebar.specie.svg ; rm cmd.r all.bar otu_hc.tre otu_bray.xls ; cd ../../\n\n";

print OUP "\n##top genus...\n";
print OUP "cd 3.Community/3-7.topGenus ; perl /mnt/sdb/lgq/bin/tools/clound/09.dominant_species_bar/dominant_species_bar2024-07-05.pl -i ../../1.Taxa/normalize/$pre\_taxa_table.xls -angle 90 ; mv dominant_genus_bar.pdf bar.top15genus.pdf ; mv dominant_genus_bar.svg bar.top15genus.svg ; mv dominant_genus_bar.png bar.top15genus.png ; rm cmd2.r color.txt genus2phylum_tmp.xls cmd1.r genus.xls g2p.xls genus.xls.tmp ; cd ../../\n\n";

if($map){
	print OUP "\n##group mean barplot...\n";
	print OUP "perl $Bin/bin/combine_otuTable_byGroup.pl 1.Taxa/normalize/tax_summary_a/domain.xls $map F M 3.Community/3-3.Bar.Group/domain.mean.xls ; perl $Bin/bin/sortByGroup.pl 3.Community/3-3.Bar.Group/domain.mean.xls group.sort F 3.Community/3-3.Bar.Group/domain.mean.xls.2 ; mv 3.Community/3-3.Bar.Group/domain.mean.xls.2 3.Community/3-3.Bar.Group/domain.mean.xls\n";
	if($groupNumber <= 30){
		$temp_w = 7;$temp_legend = 3;
	}elsif($groupNumber <= 50){
		$temp_w = 8;$temp_legend = 4;
	}else{
		$temp_w = 8 * $groupNumber/50 + 1;$temp_legend = 4 * $groupNumber/50 + 1;$temp_legend = int($temp_legend);
	}
	print OUP "perl $Bin/bin/combine_otuTable_byGroup.pl 1.Taxa/normalize/tax_summary_a/phylum.xls $map F M 3.Community/3-3.Bar.Group/phylum.mean.xls ; perl $Bin/bin/sortByGroup.pl 3.Community/3-3.Bar.Group/phylum.mean.xls group.sort F 3.Community/3-3.Bar.Group/phylum.mean.xls.2 ; mv 3.Community/3-3.Bar.Group/phylum.mean.xls.2 3.Community/3-3.Bar.Group/phylum.mean.xls ; cd 3.Community/3-3.Bar.Group ; perl /mnt/sdb/lgq/bin/tools/clound/02.bar/bar_pie2024-08-02.pl -d ./ -i phylum.mean.xls -filter top -top 30 -t T -b_blas 2 -pie F -b_ncol $temp_legend -b_lcex 1 -b_w $temp_w -b_h 7 ; mv bar_link.pdf bar.phylum.pdf ; mv bar_link.png bar.phylum.png ; mv bar_link.svg bar.phylum.svg ; cd ../../\n";
	print OUP "perl $Bin/bin/combine_otuTable_byGroup.pl 1.Taxa/normalize/tax_summary_a/class.xls $map F M 3.Community/3-3.Bar.Group/class.mean.xls ; perl $Bin/bin/sortByGroup.pl 3.Community/3-3.Bar.Group/class.mean.xls group.sort F 3.Community/3-3.Bar.Group/class.mean.xls.2 ; mv 3.Community/3-3.Bar.Group/class.mean.xls.2 3.Community/3-3.Bar.Group/class.mean.xls ; cd 3.Community/3-3.Bar.Group ; perl /mnt/sdb/lgq/bin/tools/clound/02.bar/bar_pie2024-08-02.pl -d ./ -i class.mean.xls -filter top -top 30 -t T -b_blas 2 -pie F -b_ncol $temp_legend -b_lcex 1 -b_w $temp_w -b_h 7 ; mv bar_link.pdf bar.class.pdf ; mv bar_link.png bar.class.png ; mv bar_link.svg bar.class.svg ; cd ../../\n";
$temp_w++ if($temp_w <= 8);
$temp_legend-- if($temp_legend > 4);
	print OUP "perl $Bin/bin/combine_otuTable_byGroup.pl 1.Taxa/normalize/tax_summary_a/order.xls $map F M 3.Community/3-3.Bar.Group/order.mean.xls ; perl $Bin/bin/sortByGroup.pl 3.Community/3-3.Bar.Group/order.mean.xls group.sort F 3.Community/3-3.Bar.Group/order.mean.xls.2 ; mv 3.Community/3-3.Bar.Group/order.mean.xls.2 3.Community/3-3.Bar.Group/order.mean.xls ; cd 3.Community/3-3.Bar.Group ; perl /mnt/sdb/lgq/bin/tools/clound/02.bar/bar_pie2024-08-02.pl -d ./ -i order.mean.xls -filter top -top 30 -t T -b_blas 2 -pie F -b_ncol $temp_legend -b_lcex 0.9 -b_w $temp_w -b_h 7 ; mv bar_link.pdf bar.order.pdf ; mv bar_link.png bar.order.png ; mv bar_link.svg bar.order.svg ; cd ../../\n";
	print OUP "perl $Bin/bin/combine_otuTable_byGroup.pl 1.Taxa/normalize/tax_summary_a/family.xls $map F M 3.Community/3-3.Bar.Group/family.mean.xls ; perl $Bin/bin/sortByGroup.pl 3.Community/3-3.Bar.Group/family.mean.xls group.sort F 3.Community/3-3.Bar.Group/family.mean.xls.2 ; mv 3.Community/3-3.Bar.Group/family.mean.xls.2 3.Community/3-3.Bar.Group/family.mean.xls ; cd 3.Community/3-3.Bar.Group ; perl /mnt/sdb/lgq/bin/tools/clound/02.bar/bar_pie2024-08-02.pl -d ./ -i family.mean.xls -filter top -top 30 -t T -b_blas 2 -pie F -b_ncol $temp_legend -b_lcex 0.9 -b_w $temp_w -b_h 7 ; mv bar_link.pdf bar.family.pdf ; mv bar_link.png bar.family.png ; mv bar_link.svg bar.family.svg ; cd ../../\n";
	print OUP "perl $Bin/bin/combine_otuTable_byGroup.pl 1.Taxa/normalize/tax_summary_a/genus.xls $map F M 3.Community/3-3.Bar.Group/genus.mean.xls ; perl $Bin/bin/sortByGroup.pl 3.Community/3-3.Bar.Group/genus.mean.xls group.sort F 3.Community/3-3.Bar.Group/genus.mean.xls.2 ; mv 3.Community/3-3.Bar.Group/genus.mean.xls.2 3.Community/3-3.Bar.Group/genus.mean.xls ; cd 3.Community/3-3.Bar.Group ; perl /mnt/sdb/lgq/bin/tools/clound/02.bar/bar_pie2024-08-02.pl -d ./ -i genus.mean.xls -filter top -top 30 -t T -b_blas 2 -pie F -b_ncol $temp_legend -b_lcex 0.8 -b_w $temp_w -b_h 7 ; mv bar_link.pdf bar.genus.pdf ; mv bar_link.png bar.genus.png ; mv bar_link.svg bar.genus.svg ; rm all.bar bar.pdf bar.png bar.svg out.txt cmd.r ; cd ../../\n";
	print OUP "perl $Bin/bin/combine_otuTable_byGroup.pl 1.Taxa/normalize/tax_summary_a/specie.xls $map F M 3.Community/3-3.Bar.Group/specie.mean.xls ; perl $Bin/bin/sortByGroup.pl 3.Community/3-3.Bar.Group/specie.mean.xls group.sort F 3.Community/3-3.Bar.Group/specie.mean.xls.2 ; mv 3.Community/3-3.Bar.Group/specie.mean.xls.2 3.Community/3-3.Bar.Group/specie.mean.xls ; #cd 3.Community/3-3.Bar.Group ; perl /mnt/sdb/lgq/bin/tools/clound/02.bar/bar_pie2024-08-02.pl -d ./ -i specie.mean.xls -filter top -top 30 -t T -b_blas 2 -pie F -b_ncol $temp_legend -b_lcex 0.8 -b_w $temp_w -b_h 7 ; mv bar_link.pdf bar.specie.pdf ; mv bar_link.png bar.specie.png ; mv bar_link.svg bar.specie.svg ; rm all.bar bar.pdf bar.png bar.svg out.txt cmd.r ; cd ../../\n\n";

	print OUP "\n##group mean bubble...\n";
	print OUP "cp 3.Community/3-3.Bar.Group/*.mean.xls 3.Community/3-4.Bubble.Group ; cp 3.Community/3-3.Bar.Group/*.mean.xls 3.Community/3-5.Heatmap.Group ; cp 3.Community/3-3.Bar.Group/*.mean.xls 3.Community/3-6.TreeBar.Group\n\n";

	if($groupNumber <= 10){
		$temp_w = 6;
	}elsif($groupNumber <= 30){
		$temp_w = 8;
	}elsif($groupNumber <= 50){
		$temp_w = 12;
	}else{
		$temp_w = 12 * $groupNumber/50;$temp_w = int($temp_w);
	}
	my $familyGenus_w = $temp_w;
	$familyGenus_w++;
	print OUP "cd 3.Community/3-4.Bubble.Group/ ; perl /mnt/sdb/lgq/bin/tools/clound/05.bubble/bubble2024-05-09.pl -i phylum.mean.xls -angle 60 -w $temp_w -clc 11 ; mv bubble.pdf bubble.phylum.pdf ; mv bubble.png bubble.phylum.png ; mv bubble.svg bubble.phylum.svg\n";
	print OUP "perl /mnt/sdb/lgq/bin/tools/clound/05.bubble/bubble2024-05-09.pl -i class.mean.xls -angle 60 -w $temp_w -clc 11 ; mv bubble.pdf bubble.class.pdf ; mv bubble.png bubble.class.png ; mv bubble.svg bubble.class.svg\n";
	print OUP "perl /mnt/sdb/lgq/bin/tools/clound/05.bubble/bubble2024-05-09.pl -i order.mean.xls -angle 60 -w $temp_w -clc 11 ; mv bubble.pdf bubble.order.pdf ; mv bubble.png bubble.order.png ; mv bubble.svg bubble.order.svg\n";
	print OUP "perl /mnt/sdb/lgq/bin/tools/clound/05.bubble/bubble2024-05-09.pl -i family.mean.xls -angle 60 -w $familyGenus_w -clc 10 ; mv bubble.pdf bubble.family.pdf ; mv bubble.png bubble.family.png ; mv bubble.svg bubble.family.svg\n";
	print OUP "perl /mnt/sdb/lgq/bin/tools/clound/05.bubble/bubble2024-05-09.pl -i genus.mean.xls -angle 60 -w $familyGenus_w -clc 10 ; mv bubble.pdf bubble.genus.pdf ; mv bubble.png bubble.genus.png ; mv bubble.svg bubble.genus.svg ; rm cmd.r ; cd ../../\n";
	print OUP "#cd 3.Community/3-4.Bubble.Group/ ; perl /mnt/sdb/lgq/bin/tools/clound/05.bubble/bubble2024-05-09.pl -i specie.mean.xls -angle 60 -w $familyGenus_w -clc 10 ; mv bubble.pdf bubble.specie.pdf ; mv bubble.png bubble.specie.png ; mv bubble.svg bubble.specie.svg ; rm cmd.r ; cd ../../\n\n";

	print OUP "\n##group mean heatmap...\n";
	if($groupNumber <= 10){
		$temp_w = 6;
	}elsif($groupNumber <= 30){
		$temp_w = 1.1 * $groupNumber/10 + 5;#$temp_w = int($temp_w);
	}else{
		$temp_w = 1.1 * $groupNumber/10 + 6;#$temp_w = int($temp_w);
	}
	print OUP "cd 3.Community/3-5.Heatmap.Group ; perl /mnt/sdb/lgq/bin/tools/clound/01.heatmap/plot-heatmap2024-04-17_3.pl -d ./ -i phylum.mean.xls -o heatmap.phylum -ct 0 -rt 0 -lg log2 -cs 2 -rtop 50 -slas 2 -rlc 0.9 -w $temp_w -marble 4-1-0-6\n";
	print OUP "perl /mnt/sdb/lgq/bin/tools/clound/01.heatmap/plot-heatmap2024-04-17_3.pl -d ./ -i class.mean.xls -o heatmap.class -ct 0 -rt 0 -lg log2 -cs 2 -rtop 50 -slas 2 -rlc 0.9 -w $temp_w -marble 4-1-0-10\n";
	$temp_w++;
	print OUP "perl /mnt/sdb/lgq/bin/tools/clound/01.heatmap/plot-heatmap2024-04-17_3.pl -d ./ -i order.mean.xls -o heatmap.order -ct 0 -rt 0 -lg log2 -cs 2 -rtop 50 -slas 2 -rlc 0.9 -w $temp_w -marble 4-1-0-10\n";
	print OUP "perl /mnt/sdb/lgq/bin/tools/clound/01.heatmap/plot-heatmap2024-04-17_3.pl -d ./ -i family.mean.xls -o heatmap.family -ct 0 -rt 0 -lg log2 -cs 2 -rtop 50 -slas 2 -rlc 0.8 -w $temp_w -marble 4-1-0-14\n";
	print OUP "perl /mnt/sdb/lgq/bin/tools/clound/01.heatmap/plot-heatmap2024-04-17_3.pl -d ./ -i genus.mean.xls -o heatmap.genus -ct 0 -rt 0 -lg log2 -cs 2 -rtop 50 -slas 2 -rlc 0.9 -w $temp_w -marble 4-1-0-14 ; rm cmd.r cmd.r.log final.txt warnings.log ; cd ../../ ;\n";
	print OUP "#cd 3.Community/3-5.Heatmap.Group ; perl /mnt/sdb/lgq/bin/tools/clound/01.heatmap/plot-heatmap2024-04-17_3.pl -d ./ -i specie.mean.xls -o heatmap.specie -ct 0 -rt 0 -lg log2 -cs 2 -rtop 50 -slas 2 -rlc 0.9 -w $temp_w -marble 4-1-0-14; rm cmd.r cmd.r.log final.txt warnings.log ; cd ../../\n\n";

	print OUP "\n##group mean treebar...\n";
	if($groupNumber <= 10){
		$temp_h = 5;
	}else{
		$temp_h = $groupNumber/10 + 4;
	}
	print OUP "perl $Bin/bin/combine_otuTable_byGroup.pl 1.Taxa/normalize/$pre\_table.xls $map F M 3.Community/3-6.TreeBar.Group/$pre\_table.xls ; perl $Bin/bin/sortByGroup.pl 3.Community/3-6.TreeBar.Group/$pre\_table.xls group.sort F 3.Community/3-6.TreeBar.Group/$pre\_table.xls.2 ; mv 3.Community/3-6.TreeBar.Group/$pre\_table.xls.2 3.Community/3-6.TreeBar.Group/$pre\_table.xls\n";
	print OUP "cd 3.Community/3-6.TreeBar.Group ; perl /mnt/sdb/lgq/bin/tools/clound/20.treebar/treebar2024-09-05.pl -otu $pre\_table.xls -tax phylum.mean.xls -filter top -scex 1.2 -lcex 1.2 -h $temp_h -lw 1.3-1.8-1.3 ; mv treebar.pdf treebar.phylum.pdf ; mv treebar.png treebar.phylum.png ; mv treebar.svg treebar.phylum.svg\n";
	print OUP "perl /mnt/sdb/lgq/bin/tools/clound/20.treebar/treebar2024-09-05.pl -otu $pre\_table.xls -tax class.mean.xls -filter top -scex 1.2 -lcex 1.2 -h $temp_h -lw 1.3-1.8-1.3 ; mv treebar.pdf treebar.class.pdf ; mv treebar.png treebar.class.png ; mv treebar.svg treebar.class.svg\n";
	print OUP "perl /mnt/sdb/lgq/bin/tools/clound/20.treebar/treebar2024-09-05.pl -otu $pre\_table.xls -tax order.mean.xls -filter top -scex 1.2 -lcex 1.2 -w 12 -h $temp_h -lw 1.3-1.8-1.3 ; mv treebar.pdf treebar.order.pdf ; mv treebar.png treebar.order.png ; mv treebar.svg treebar.order.svg\n";
	print OUP "perl /mnt/sdb/lgq/bin/tools/clound/20.treebar/treebar2024-09-05.pl -otu $pre\_table.xls -tax family.mean.xls -filter top -scex 1.2 -lcex 1.1 -w 13  -h $temp_h -lw 1.3-1.8-1.3 ; mv treebar.pdf treebar.family.pdf ; mv treebar.png treebar.family.png ; mv treebar.svg treebar.family.svg\n";
	print OUP "perl /mnt/sdb/lgq/bin/tools/clound/20.treebar/treebar2024-09-05.pl -otu $pre\_table.xls -tax genus.mean.xls -filter top -scex 1.2 -lcex 1.1 -w 13 -h $temp_h -lw 1.3-1.8-1.3 ; mv treebar.pdf treebar.genus.pdf ; mv treebar.png treebar.genus.png ; mv treebar.svg treebar.genus.svg ; rm cmd.r all.bar otu_hc.tre otu_bray.xls $pre\_table.xls ; cd ../../\n";
	print OUP "#perl $Bin/bin/combine_otuTable_byGroup.pl 1.Taxa/normalize/$pre\_table.xls $map F M 3.Community/3-6.TreeBar.Group/$pre\_table.xls ; cd 3.Community/3-6.TreeBar.Group ; perl /mnt/sdb/lgq/bin/tools/clound/20.treebar/treebar2024-09-05.pl -otu $pre\_table.xls -tax specie.mean.xls -filter top -scex 1.2 -lcex 1.1 -w 13 -h $temp_h -lw 1.3-1.8-1.3 ; mv treebar.pdf treebar.specie.pdf ; mv treebar.png treebar.specie.png ; mv treebar.svg treebar.specie.svg ; rm cmd.r all.bar otu_hc.tre otu_bray.xls $pre\_table.xls ; cd ../../\n\n";
}


print OUP "\n\n###beta diversity...\n";
`mkdir $outDir/4.Beta` unless(-e "$outDir/4.Beta");
`mkdir $outDir/4.Beta/4-1.Distance` unless(-e "$outDir/4.Beta/4-1.Distance");
`mkdir $outDir/4.Beta/4-1.Distance/bray_curtis` unless(-e "$outDir/4.Beta/4-1.Distance/bray_curtis");
`mkdir $outDir/4.Beta/4-1.Distance/unweighted_unifrac` unless(-e "$outDir/4.Beta/4-1.Distance/unweighted_unifrac");
`mkdir $outDir/4.Beta/4-1.Distance/weighted_unifrac` unless(-e "$outDir/4.Beta/4-1.Distance/weighted_unifrac");
`mkdir $outDir/4.Beta/4-2.PCA` unless(-e "$outDir/4.Beta/4-2.PCA");
`mkdir $outDir/4.Beta/4-3.PCOA` unless(-e "$outDir/4.Beta/4-3.PCOA");
`mkdir $outDir/4.Beta/4-3.PCOA/bray_curtis` unless(-e "$outDir/4.Beta/4-3.PCOA/bray_curtis");
`mkdir $outDir/4.Beta/4-3.PCOA/unweighted_unifrac` unless(-e "$outDir/4.Beta/4-3.PCOA/unweighted_unifrac");
`mkdir $outDir/4.Beta/4-3.PCOA/weighted_unifrac` unless(-e "$outDir/4.Beta/4-3.PCOA/weighted_unifrac");
`mkdir $outDir/4.Beta/4-4.NMDS` unless(-e "$outDir/4.Beta/4-4.NMDS");
`mkdir $outDir/4.Beta/4-4.NMDS/bray_curtis` unless(-e "$outDir/4.Beta/4-4.NMDS/bray_curtis");
`mkdir $outDir/4.Beta/4-4.NMDS/unweighted_unifrac` unless(-e "$outDir/4.Beta/4-4.NMDS/unweighted_unifrac");
`mkdir $outDir/4.Beta/4-4.NMDS/weighted_unifrac` unless(-e "$outDir/4.Beta/4-4.NMDS/weighted_unifrac");
`mkdir $outDir/4.Beta/4-5.Hcluster_tree` unless(-e "$outDir/4.Beta/4-5.Hcluster_tree");
`mkdir $outDir/4.Beta/4-5.Hcluster_tree/bray_curtis` unless(-e "$outDir/4.Beta/4-5.Hcluster_tree/bray_curtis");
`mkdir $outDir/4.Beta/4-5.Hcluster_tree/unweighted_unifrac` unless(-e "$outDir/4.Beta/4-5.Hcluster_tree/unweighted_unifrac");
`mkdir $outDir/4.Beta/4-5.Hcluster_tree/weighted_unifrac` unless(-e "$outDir/4.Beta/4-5.Hcluster_tree/weighted_unifrac");
if($G3flag eq "T"){
	`mkdir $outDir/4.Beta/4-6.Anosim` unless(-e "$outDir/4.Beta/4-6.Anosim");
	`mkdir $outDir/4.Beta/4-6.Anosim/bray_curtis` unless(-e "$outDir/4.Beta/4-6.Anosim/bray_curtis");
	`mkdir $outDir/4.Beta/4-6.Anosim/unweighted_unifrac` unless(-e "$outDir/4.Beta/4-6.Anosim/unweighted_unifrac");
	`mkdir $outDir/4.Beta/4-6.Anosim/weighted_unifrac` unless(-e "$outDir/4.Beta/4-6.Anosim/weighted_unifrac");
	`mkdir $outDir/4.Beta/4-7.Adonis` unless(-e "$outDir/4.Beta/4-7.Adonis");
	`mkdir $outDir/4.Beta/4-7.Adonis/bray_curtis` unless(-e "$outDir/4.Beta/4-7.Adonis/bray_curtis");
	`mkdir $outDir/4.Beta/4-7.Adonis/unweighted_unifrac` unless(-e "$outDir/4.Beta/4-7.Adonis/unweighted_unifrac");
	`mkdir $outDir/4.Beta/4-7.Adonis/weighted_unifrac` unless(-e "$outDir/4.Beta/4-7.Adonis/weighted_unifrac");
}

print OUP "beta_diversity.py -i 1.Taxa/normalize/$pre\_table.biom -t 1.Taxa/normalize/phy.tre -m bray_curtis,unweighted_unifrac,weighted_unifrac -o 4.Beta/4-1.Distance\n";
print OUP "mv 4.Beta/4-1.Distance/bray_curtis_$pre\_table.txt 4.Beta/4-1.Distance/bray_curtis/bray_curtis.txt ; mv 4.Beta/4-1.Distance/unweighted_unifrac_$pre\_table.txt 4.Beta/4-1.Distance/unweighted_unifrac/unweighted_unifrac.txt ; mv 4.Beta/4-1.Distance/weighted_unifrac_$pre\_table.txt 4.Beta/4-1.Distance/weighted_unifrac/weighted_unifrac.txt\n";

print OUP "\n##distance heatmap...\n";
if($samNumber <= 40){
	$temp_w = 7;$temp_h = $temp_w + 1;
}else{
	$temp_w = 7 * $samNumber/40;$temp_h = $temp_w + 1;
}
print OUP "cd 4.Beta/4-1.Distance/bray_curtis ; perl /mnt/sdb/lgq/bin/tools/clound/16.sample_dist_heatmap/sample_dist_heatmap2024-08-28.pl -i bray_curtis.txt -w $temp_w -h $temp_h ; mv sample_dist_heatmap.pdf heatmap.distance.pdf ; mv sample_dist_heatmap.svg heatmap.distance.svg ; mv sample_dist_heatmap.png heatmap.distance.png ; cd ../../../\n";
print OUP "cd 4.Beta/4-1.Distance/unweighted_unifrac ; perl /mnt/sdb/lgq/bin/tools/clound/16.sample_dist_heatmap/sample_dist_heatmap2024-08-28.pl -i unweighted_unifrac.txt -w $temp_w -h $temp_h ; mv sample_dist_heatmap.pdf heatmap.distance.pdf ; mv sample_dist_heatmap.svg heatmap.distance.svg ; mv sample_dist_heatmap.png heatmap.distance.png ; cd ../../../\n";
print OUP "cd 4.Beta/4-1.Distance/weighted_unifrac ; perl /mnt/sdb/lgq/bin/tools/clound/16.sample_dist_heatmap/sample_dist_heatmap2024-08-28.pl -i weighted_unifrac.txt -w $temp_w -h $temp_h ; mv sample_dist_heatmap.pdf heatmap.distance.pdf ; mv sample_dist_heatmap.svg heatmap.distance.svg ; mv sample_dist_heatmap.png heatmap.distance.png ; cd ../../../\n";

if($G3flag eq "T"){
	print OUP "\n##distance boxplot...\n";
	if($groupNumber <= 3){
		$temp_w = 3;$temp_h = 4;
	}elsif($groupNumber <= 5){
		$temp_w = 4;$temp_h = 5;
	}elsif($groupNumber <= 10){
		$temp_w = 6;$temp_h = 6;
	}else{
		$temp_w = 6 * $groupNumber/10;$temp_h = $temp_w;
	}
	print OUP "cd 4.Beta/4-1.Distance/bray_curtis ; perl /mnt/sdb/lgq/bin/tools/clound/19.distance/distance_boxplot2024-08-28.pl -i bray_curtis.txt -g ../../../group.sort -m ../../../$map -w $temp_w -h $temp_h ; mv distance_box.pdf box.distance.pdf ; mv distance_box.svg box.distance.svg ; mv distance_box.png box.distance.png ; rm color.txt ; cd ../../../\n";
	print OUP "cd 4.Beta/4-1.Distance/unweighted_unifrac ; perl /mnt/sdb/lgq/bin/tools/clound/19.distance/distance_boxplot2024-08-28.pl -i unweighted_unifrac.txt -m ../../../$map -g ../../../group.sort -w $temp_w -h $temp_h ; mv distance_box.pdf box.distance.pdf ; mv distance_box.svg box.distance.svg ; mv distance_box.png box.distance.png ; rm color.txt ; cd ../../../\n";
	print OUP "cd 4.Beta/4-1.Distance/weighted_unifrac ; perl /mnt/sdb/lgq/bin/tools/clound/19.distance/distance_boxplot2024-08-28.pl -i weighted_unifrac.txt -m ../../../$map -g ../../../group.sort -w $temp_w -h $temp_h ; mv distance_box.pdf box.distance.pdf ; mv distance_box.svg box.distance.svg ; mv distance_box.png box.distance.png ; rm color.txt ; cd ../../../\n\n";
}

print OUP "\n##PCA...\n";
if($samNumber <= 50){
	$temp_w = 6;$temp_h = 6;
}else{
	$temp_w = 7;$temp_h = 7;
}
$temp_w += 1 if($map);
my $para = "-w $temp_w -h $temp_h";
$para .= " -m ../../$map -g ../../group.sort -c 2 -type 2" if($map);
print OUP "cd 4.Beta/4-2.PCA ; perl /mnt/sdb/lgq/bin/tools/clound/10.pca/pca2024-08-02.pl -i ../../1.Taxa/normalize/$pre\_table.xls $para -l T ; mv pca.pc1-2.pdf pca.pc1-2.label.pdf ; mv pca.pc1-2.svg pca.pc1-2.label.svg ; mv pca.pc1-2.png pca.pc1-2.label.png\n";
print OUP "perl /mnt/sdb/lgq/bin/tools/clound/10.pca/pca2024-08-02.pl -i ../../1.Taxa/normalize/$pre\_table.xls $para -l T -pc pc1-pc3 ; mv pca.pc1-3.pdf pca.pc1-3.label.pdf ; mv pca.pc1-3.svg pca.pc1-3.label.svg ; mv pca.pc1-3.png pca.pc1-3.label.png\n";
print OUP "perl /mnt/sdb/lgq/bin/tools/clound/10.pca/pca2024-08-02.pl -i ../../1.Taxa/normalize/$pre\_table.xls $para -l T -pc pc2-pc3 ; mv pca.pc2-3.pdf pca.pc2-3.label.pdf ; mv pca.pc2-3.svg pca.pc2-3.label.svg ; mv pca.pc2-3.png pca.pc2-3.label.png ; rm cmd.r color.txt tmp_group_color_pch.tsv pca_rotation.xls ; cd ../../\n";
if($map){
	print OUP "cd 4.Beta/4-2.PCA ; perl /mnt/sdb/lgq/bin/tools/clound/10.pca/pca2024-08-02.pl -i ../../1.Taxa/normalize/$pre\_table.xls $para\n";
	print OUP "perl /mnt/sdb/lgq/bin/tools/clound/10.pca/pca2024-08-02.pl -i ../../1.Taxa/normalize/$pre\_table.xls $para -pc pc1-pc3\n";
	print OUP "perl /mnt/sdb/lgq/bin/tools/clound/10.pca/pca2024-08-02.pl -i ../../1.Taxa/normalize/$pre\_table.xls $para -pc pc2-pc3 ; rm cmd.r color.txt tmp_group_color_pch.tsv pca_rotation.xls ; cd ../../\n";

}

print OUP "\n##PCOA...\n";
$para = "-w $temp_w -h $temp_h";
$para .= " -m ../../../$map -g ../../../group.sort -c 2 -type 2" if($map);
print OUP "\ncd 4.Beta/4-3.PCOA/bray_curtis ; perl /mnt/sdb/lgq/bin/tools/clound/11.pcoa/pcoa2024-09-25.pl -i ../../4-1.Distance/bray_curtis/bray_curtis.txt -method distance $para -l T ; mv pcoa.pc1-2.pdf pcoa.pc1-2.label.pdf ; mv pcoa.pc1-2.svg pcoa.pc1-2.label.svg ; mv pcoa.pc1-2.png pcoa.pc1-2.label.png\n";
print OUP "perl /mnt/sdb/lgq/bin/tools/clound/11.pcoa/pcoa2024-09-25.pl -i ../../4-1.Distance/bray_curtis/bray_curtis.txt -method distance $para -l T -pc pc1-pc3 ; mv pcoa.pc1-3.pdf pcoa.pc1-3.label.pdf ; mv pcoa.pc1-3.svg pcoa.pc1-3.label.svg ; mv pcoa.pc1-3.png pcoa.pc1-3.label.png\n";
print OUP "perl /mnt/sdb/lgq/bin/tools/clound/11.pcoa/pcoa2024-09-25.pl -i ../../4-1.Distance/bray_curtis/bray_curtis.txt -method distance $para -l T -pc pc2-pc3 ; mv pcoa.pc2-3.pdf pcoa.pc2-3.label.pdf ; mv pcoa.pc2-3.svg pcoa.pc2-3.label.svg ; mv pcoa.pc2-3.png pcoa.pc2-3.label.png ; rm warnings.log color.txt tmp_group_color_pch.tsv pcoa_rotation.xls ; cd ../../../\n";
if($map){
	print OUP "cd 4.Beta/4-3.PCOA/bray_curtis ; perl /mnt/sdb/lgq/bin/tools/clound/11.pcoa/pcoa2024-09-25.pl -i ../../4-1.Distance/bray_curtis/bray_curtis.txt -method distance $para\n";
	print OUP "perl /mnt/sdb/lgq/bin/tools/clound/11.pcoa/pcoa2024-09-25.pl -i ../../4-1.Distance/bray_curtis/bray_curtis.txt -method distance $para -pc pc1-pc3\n";
	print OUP "perl /mnt/sdb/lgq/bin/tools/clound/11.pcoa/pcoa2024-09-25.pl -i ../../4-1.Distance/bray_curtis/bray_curtis.txt -method distance $para -pc pc2-pc3 ; rm warnings.log color.txt tmp_group_color_pch.tsv pcoa_rotation.xls ; cd ../../../\n";
}

print OUP "cd 4.Beta/4-3.PCOA/unweighted_unifrac ; perl /mnt/sdb/lgq/bin/tools/clound/11.pcoa/pcoa2024-09-25.pl -i ../../4-1.Distance/unweighted_unifrac/unweighted_unifrac.txt -method distance $para -l T ; mv pcoa.pc1-2.pdf pcoa.pc1-2.label.pdf ; mv pcoa.pc1-2.svg pcoa.pc1-2.label.svg ; mv pcoa.pc1-2.png pcoa.pc1-2.label.png\n";
print OUP "perl /mnt/sdb/lgq/bin/tools/clound/11.pcoa/pcoa2024-09-25.pl -i ../../4-1.Distance/unweighted_unifrac/unweighted_unifrac.txt -method distance $para -l T -pc pc1-pc3 ; mv pcoa.pc1-3.pdf pcoa.pc1-3.label.pdf ; mv pcoa.pc1-3.svg pcoa.pc1-3.label.svg ; mv pcoa.pc1-3.png pcoa.pc1-3.label.png\n";
print OUP "perl /mnt/sdb/lgq/bin/tools/clound/11.pcoa/pcoa2024-09-25.pl -i ../../4-1.Distance/unweighted_unifrac/unweighted_unifrac.txt -method distance $para -l T -pc pc2-pc3 ; mv pcoa.pc2-3.pdf pcoa.pc2-3.label.pdf ; mv pcoa.pc2-3.svg pcoa.pc2-3.label.svg ; mv pcoa.pc2-3.png pcoa.pc2-3.label.png ; rm warnings.log color.txt tmp_group_color_pch.tsv pcoa_rotation.xls ; cd ../../../\n";
if($map){
	print OUP "cd 4.Beta/4-3.PCOA/unweighted_unifrac ; perl /mnt/sdb/lgq/bin/tools/clound/11.pcoa/pcoa2024-09-25.pl -i ../../4-1.Distance/unweighted_unifrac/unweighted_unifrac.txt -method distance $para\n";
	print OUP "perl /mnt/sdb/lgq/bin/tools/clound/11.pcoa/pcoa2024-09-25.pl -i ../../4-1.Distance/unweighted_unifrac/unweighted_unifrac.txt -method distance $para -pc pc1-pc3\n";
	print OUP "perl /mnt/sdb/lgq/bin/tools/clound/11.pcoa/pcoa2024-09-25.pl -i ../../4-1.Distance/unweighted_unifrac/unweighted_unifrac.txt -method distance $para -pc pc2-pc3 ; rm warnings.log color.txt tmp_group_color_pch.tsv pcoa_rotation.xls ; cd ../../../\n";
}

print OUP "cd 4.Beta/4-3.PCOA/weighted_unifrac ; perl /mnt/sdb/lgq/bin/tools/clound/11.pcoa/pcoa2024-09-25.pl -i ../../4-1.Distance/weighted_unifrac/weighted_unifrac.txt -method distance $para -l T ; mv pcoa.pc1-2.pdf pcoa.pc1-2.label.pdf ; mv pcoa.pc1-2.svg pcoa.pc1-2.label.svg ; mv pcoa.pc1-2.png pcoa.pc1-2.label.png\n";
print OUP "perl /mnt/sdb/lgq/bin/tools/clound/11.pcoa/pcoa2024-09-25.pl -i ../../4-1.Distance/weighted_unifrac/weighted_unifrac.txt -method distance $para -l T -pc pc1-pc3 ; mv pcoa.pc1-3.pdf pcoa.pc1-3.label.pdf ; mv pcoa.pc1-3.svg pcoa.pc1-3.label.svg ; mv pcoa.pc1-3.png pcoa.pc1-3.label.png\n";
print OUP "perl /mnt/sdb/lgq/bin/tools/clound/11.pcoa/pcoa2024-09-25.pl -i ../../4-1.Distance/weighted_unifrac/weighted_unifrac.txt -method distance $para -l T -pc pc2-pc3 ; mv pcoa.pc2-3.pdf pcoa.pc2-3.label.pdf ; mv pcoa.pc2-3.svg pcoa.pc2-3.label.svg ; mv pcoa.pc2-3.png pcoa.pc2-3.label.png ; rm warnings.log color.txt tmp_group_color_pch.tsv pcoa_rotation.xls ; cd ../../../\n";
if($map){
	print OUP "cd 4.Beta/4-3.PCOA/weighted_unifrac ; perl /mnt/sdb/lgq/bin/tools/clound/11.pcoa/pcoa2024-09-25.pl -i ../../4-1.Distance/weighted_unifrac/weighted_unifrac.txt -method distance $para\n";
	print OUP "perl /mnt/sdb/lgq/bin/tools/clound/11.pcoa/pcoa2024-09-25.pl -i ../../4-1.Distance/weighted_unifrac/weighted_unifrac.txt -method distance $para -pc pc1-pc3\n";
	print OUP "perl /mnt/sdb/lgq/bin/tools/clound/11.pcoa/pcoa2024-09-25.pl -i ../../4-1.Distance/weighted_unifrac/weighted_unifrac.txt -method distance $para -pc pc2-pc3 ; rm warnings.log color.txt tmp_group_color_pch.tsv pcoa_rotation.xls ; cd ../../../\n";
}

print OUP "\n##NMDS...\n";
$para = "-w $temp_w -h $temp_h";
$para .= " -m ../../../$map -g ../../../group.sort -c 2 -type 2" if($map);
print OUP "cd 4.Beta/4-4.NMDS/bray_curtis ; perl /mnt/sdb/lgq/bin/tools/clound/12.nmds/nmds2024-08-05.pl -i ../../4-1.Distance/bray_curtis/bray_curtis.txt -method distance $para -l T ; mv nmds.pdf nmds.label.pdf ; mv nmds.svg nmds.label.svg ; mv nmds.png nmds.label.png ; rm cmd.r color.txt warnings.log tmp_group_color_pch.tsv ; cd ../../../\n";
print OUP "cd 4.Beta/4-4.NMDS/bray_curtis ; perl /mnt/sdb/lgq/bin/tools/clound/12.nmds/nmds2024-08-05.pl -i ../../4-1.Distance/bray_curtis/bray_curtis.txt -method distance $para ; rm cmd.r color.txt warnings.log tmp_group_color_pch.tsv ; cd ../../../\n" if($map);
print OUP "cd 4.Beta/4-4.NMDS/unweighted_unifrac ; perl /mnt/sdb/lgq/bin/tools/clound/12.nmds/nmds2024-08-05.pl -i ../../4-1.Distance/unweighted_unifrac/unweighted_unifrac.txt -method distance $para -l T ; mv nmds.pdf nmds.label.pdf ; mv nmds.svg nmds.label.svg ; mv nmds.png nmds.label.png ; rm cmd.r color.txt warnings.log tmp_group_color_pch.tsv ; cd ../../../\n";
print OUP "cd 4.Beta/4-4.NMDS/unweighted_unifrac ; perl /mnt/sdb/lgq/bin/tools/clound/12.nmds/nmds2024-08-05.pl -i ../../4-1.Distance/unweighted_unifrac/unweighted_unifrac.txt -method distance $para ; rm cmd.r color.txt warnings.log tmp_group_color_pch.tsv ; cd ../../../\n" if($map);
print OUP "cd 4.Beta/4-4.NMDS/weighted_unifrac ; perl /mnt/sdb/lgq/bin/tools/clound/12.nmds/nmds2024-08-05.pl -i ../../4-1.Distance/weighted_unifrac/weighted_unifrac.txt -method distance $para -l T ; mv nmds.pdf nmds.label.pdf ; mv nmds.svg nmds.label.svg ; mv nmds.png nmds.label.png ; rm cmd.r color.txt warnings.log tmp_group_color_pch.tsv ; cd ../../../\n";
print OUP "cd 4.Beta/4-4.NMDS/weighted_unifrac ; perl /mnt/sdb/lgq/bin/tools/clound/12.nmds/nmds2024-08-05.pl -i ../../4-1.Distance/weighted_unifrac/weighted_unifrac.txt -method distance $para ; rm cmd.r color.txt warnings.log tmp_group_color_pch.tsv ; cd ../../../\n" if($map);

print OUP "\n##Hcluster tree...\n";
if($samNumber <= 25){
	$temp_w = 5;$temp_h = 6;
}elsif($samNumber <= 50){
	$temp_h = $samNumber/5 + 1;$temp_w = $temp_h * 0.7;
}else{
	$temp_h = $samNumber/5;$temp_w = $temp_h * 0.7;
}
$para = "-w $temp_w -h $temp_h";
$para .= " -s ../../../group.sort -d ../../../$map" if($map);
print OUP "cd 4.Beta/4-5.Hcluster_tree/bray_curtis ; perl /mnt/sdb/lgq/bin/tools/clound/15.hcluster_tree/sample_hcluster_tree2024-08-14.pl -i ../../4-1.Distance/bray_curtis/bray_curtis.txt -method distance $para ; mv sample.tre hcluster_tree.tre ; rm trecmd.r color.txt ; cd ../../../\n";
print OUP "cd 4.Beta/4-5.Hcluster_tree/unweighted_unifrac ; perl /mnt/sdb/lgq/bin/tools/clound/15.hcluster_tree/sample_hcluster_tree2024-08-14.pl -i ../../4-1.Distance/unweighted_unifrac/unweighted_unifrac.txt -method distance $para ; mv sample.tre hcluster_tree.tre ; rm trecmd.r color.txt ; cd ../../../\n";
print OUP "cd 4.Beta/4-5.Hcluster_tree/weighted_unifrac ; perl /mnt/sdb/lgq/bin/tools/clound/15.hcluster_tree/sample_hcluster_tree2024-08-14.pl -i ../../4-1.Distance/weighted_unifrac/weighted_unifrac.txt -method distance $para ; mv sample.tre hcluster_tree.tre ; rm trecmd.r color.txt ; cd ../../../\n";

if($G3flag eq "T"){
	print OUP "\n##Anosim...\n";
	if($groupNumber < 5){
		$temp_w = 4;$temp_h = 4;
	}elsif($groupNumber < 10){
		$temp_w = 5;$temp_h = 5;
	}else{
		$temp_w = 6;$temp_h = 6;
	}
	print OUP "cd 4.Beta/4-6.Anosim/bray_curtis ; perl /mnt/sdb/lgq/bin/tools/clound/06.ANOSIM/ANOSIM2024-09-26.pl -i ../../4-1.Distance/bray_curtis/bray_curtis.txt -m ../../../$map -g ../../../group.sort -method -angle 60 -w $temp_w -h $temp_h ; rm cmd.r color.txt data.txt ; cd ../../../\n";
	print OUP "cd 4.Beta/4-6.Anosim/unweighted_unifrac ; perl /mnt/sdb/lgq/bin/tools/clound/06.ANOSIM/ANOSIM2024-09-26.pl -i ../../4-1.Distance/unweighted_unifrac/unweighted_unifrac.txt -m ../../../$map -g ../../../group.sort -method -angle 60 -w $temp_w -h $temp_h ; rm cmd.r color.txt data.txt ; cd ../../../\n";
	print OUP "cd 4.Beta/4-6.Anosim/weighted_unifrac ; perl /mnt/sdb/lgq/bin/tools/clound/06.ANOSIM/ANOSIM2024-09-26.pl -i ../../4-1.Distance/weighted_unifrac/weighted_unifrac.txt -m ../../../$map -g ../../../group.sort -method -angle 60 -w $temp_w -h $temp_h ; rm cmd.r color.txt data.txt ; cd ../../../\n";

	print OUP "\n##Adonis...\n";
	print OUP "compare_categories.py --method adonis -i 4.Beta/4-1.Distance/bray_curtis/bray_curtis.txt -m $map -c Group -o 4.Beta/4-7.Adonis/bray_curtis\n";
	print OUP "compare_categories.py --method adonis -i 4.Beta/4-1.Distance/unweighted_unifrac/unweighted_unifrac.txt -m $map -c Group -o 4.Beta/4-7.Adonis/unweighted_unifrac\n";
	print OUP "compare_categories.py --method adonis -i 4.Beta/4-1.Distance/weighted_unifrac/weighted_unifrac.txt -m $map -c Group -o 4.Beta/4-7.Adonis/weighted_unifrac\n";
}


if($G3flag eq "T"){
	print OUP "\n\n###diferent analysis...\n";

	`mkdir $outDir/5.Diff` unless(-e "$outDir/5.Diff");
	if($groupNumber == 2){
		`mkdir $outDir/5.Diff/5-1.T-test` unless(-e "$outDir/5.Diff/5-1.T-test");
		`mkdir $outDir/5.Diff/5-2.Wilcoxon` unless(-e "$outDir/5.Diff/5-2.Wilcoxon");
	}else{
		`mkdir $outDir/5.Diff/5-1.ANOVA` unless(-e "$outDir/5.Diff/5-1.ANOVA");
		`mkdir $outDir/5.Diff/5-2.KW` unless(-e "$outDir/5.Diff/5-2.KW");
	}
	`mkdir $outDir/5.Diff/5-3.Lefse` unless(-e "$outDir/5.Diff/5-3.Lefse");
	if($groupNumber < 5){
		$temp_w = 7;$temp_h = 7;
	}elsif($groupNumber < 8){
		$temp_w = 8;$temp_h = 8;
	}else{
		$temp_w = $groupNumber;$temp_h = $groupNumber;
	}
	my $familyGenus_w = $temp_w;
	$familyGenus_w++ if($temp_w < 8);
	if($groupNumber == 2){
		print OUP "\n##T test...\n";
		print OUP "cd 5.Diff/5-1.T-test ; perl /mnt/sdb/lgq/bin/tools/clound/04.ANOVA/ANOVA-T.test2024-09-24.pl -i ../../1.Taxa/normalize/tax_summary_a/phylum.xls -m ../../$map -g ../../group.sort -w $temp_w -h $temp_h ; mv T.test.pdf phylum.T-test.pdf ; mv T.test.svg phylum.T-test.svg ; mv T.test.png phylum.T-test.png ; mv T.test.txt phylum.T-test.txt\n";
		print OUP "perl /mnt/sdb/lgq/bin/tools/clound/04.ANOVA/ANOVA-T.test2024-09-24.pl -i ../../1.Taxa/normalize/tax_summary_a/class.xls -m ../../$map -g ../../group.sort -w $temp_w -h $temp_h ; mv T.test.pdf class.T-test.pdf ; mv T.test.svg class.T-test.svg ; mv T.test.png class.T-test.png ; mv T.test.txt class.T-test.txt\n";
		print OUP "perl /mnt/sdb/lgq/bin/tools/clound/04.ANOVA/ANOVA-T.test2024-09-24.pl -i ../../1.Taxa/normalize/tax_summary_a/order.xls -m ../../$map -g ../../group.sort -w $temp_w -h $temp_h ; mv T.test.pdf order.T-test.pdf ; mv T.test.svg order.T-test.svg ; mv T.test.png order.T-test.png ; mv T.test.txt order.T-test.txt\n";
		print OUP "perl /mnt/sdb/lgq/bin/tools/clound/04.ANOVA/ANOVA-T.test2024-09-24.pl -i ../../1.Taxa/normalize/tax_summary_a/family.xls -m ../../$map -g ../../group.sort -w $familyGenus_w -h $temp_h ; mv T.test.pdf family.T-test.pdf ; mv T.test.svg family.T-test.svg ; mv T.test.png family.T-test.png ; mv T.test.txt family.T-test.txt\n";
		print OUP "perl /mnt/sdb/lgq/bin/tools/clound/04.ANOVA/ANOVA-T.test2024-09-24.pl -i ../../1.Taxa/normalize/tax_summary_a/genus.xls -m ../../$map -g ../../group.sort -w $familyGenus_w -h $temp_h ; mv T.test.pdf genus.T-test.pdf ; mv T.test.svg genus.T-test.svg ; mv T.test.png genus.T-test.png ; mv T.test.txt genus.T-test.txt ; rm z2.r zz.txt zz1.txt z.r color.txt rowsort ; cd ../../\n";

		print OUP "\n##Wilcoxon...\n";
		print OUP "cd 5.Diff/5-2.Wilcoxon ; perl /mnt/sdb/lgq/bin/tools/clound/03.KW/Kruskal-Wallis_Wilcoxon_test2024-09-24.pl -i ../../1.Taxa/normalize/tax_summary_a/phylum.xls -m ../../$map -g ../../group.sort -w $temp_w -h $temp_h ; mv wilcox.pdf phylum.wilcoxon.pdf ; mv wilcox.svg phylum.wilcoxon.svg ; mv wilcox.png phylum.wilcoxon.png ; mv wilcox.txt phylum.wilcoxon.txt\n";
		print OUP "perl /mnt/sdb/lgq/bin/tools/clound/03.KW/Kruskal-Wallis_Wilcoxon_test2024-09-24.pl -i ../../1.Taxa/normalize/tax_summary_a/class.xls -m ../../$map -g ../../group.sort -w $temp_w -h $temp_h ; mv wilcox.pdf class.wilcoxon.pdf ; mv wilcox.svg class.wilcoxon.svg ; mv wilcox.png class.wilcoxon.png ; mv wilcox.txt class.wilcoxon.txt\n";
		print OUP "perl /mnt/sdb/lgq/bin/tools/clound/03.KW/Kruskal-Wallis_Wilcoxon_test2024-09-24.pl -i ../../1.Taxa/normalize/tax_summary_a/order.xls -m ../../$map -g ../../group.sort -w $temp_w -h $temp_h ; mv wilcox.pdf order.wilcoxon.pdf ; mv wilcox.svg order.wilcoxon.svg ; mv wilcox.png order.wilcoxon.png ; mv wilcox.txt order.wilcoxon.txt\n";
		print OUP "perl /mnt/sdb/lgq/bin/tools/clound/03.KW/Kruskal-Wallis_Wilcoxon_test2024-09-24.pl -i ../../1.Taxa/normalize/tax_summary_a/family.xls -m ../../$map -g ../../group.sort -w $familyGenus_w -h $temp_h ; mv wilcox.pdf family.wilcoxon.pdf ; mv wilcox.svg family.wilcoxon.svg ; mv wilcox.png family.wilcoxon.png ; mv wilcox.txt family.wilcoxon.txt\n";
		print OUP "perl /mnt/sdb/lgq/bin/tools/clound/03.KW/Kruskal-Wallis_Wilcoxon_test2024-09-24.pl -i ../../1.Taxa/normalize/tax_summary_a/genus.xls -m ../../$map -g ../../group.sort -w $familyGenus_w -h $temp_h ; mv wilcox.pdf genus.wilcoxon.pdf ; mv wilcox.svg genus.wilcoxon.svg ; mv wilcox.png genus.wilcoxon.png ; mv wilcox.txt genus.wilcoxon.txt ; rm z2.r zz.txt zz1.txt z.r color.txt rowsort ; cd ../../\n";
	}else{
		print OUP "\n##ANOVA...\n";
		print OUP "cd 5.Diff/5-1.ANOVA ; perl /mnt/sdb/lgq/bin/tools/clound/04.ANOVA/ANOVA-T.test2024-09-24.pl -i ../../1.Taxa/normalize/tax_summary_a/phylum.xls -m ../../$map -g ../../group.sort -w $temp_w -h $temp_h ; mv ANOVA.pdf phylum.anova.pdf ; mv ANOVA.svg phylum.anova.svg ; mv ANOVA.png phylum.anova.png ; mv ANOVA.txt phylum.anova.txt\n";
		print OUP "perl /mnt/sdb/lgq/bin/tools/clound/04.ANOVA/ANOVA-T.test2024-09-24.pl -i ../../1.Taxa/normalize/tax_summary_a/class.xls -m ../../$map -g ../../group.sort -w $temp_w -h $temp_h ; mv ANOVA.pdf class.anova.pdf ; mv ANOVA.svg class.anova.svg ; mv ANOVA.png class.anova.png ; mv ANOVA.txt class.anova.txt\n";
		print OUP "perl /mnt/sdb/lgq/bin/tools/clound/04.ANOVA/ANOVA-T.test2024-09-24.pl -i ../../1.Taxa/normalize/tax_summary_a/order.xls -m ../../$map -g ../../group.sort -w $temp_w -h $temp_h ; mv ANOVA.pdf order.anova.pdf ; mv ANOVA.svg order.anova.svg ; mv ANOVA.png order.anova.png ; mv ANOVA.txt order.anova.txt\n";
		print OUP "perl /mnt/sdb/lgq/bin/tools/clound/04.ANOVA/ANOVA-T.test2024-09-24.pl -i ../../1.Taxa/normalize/tax_summary_a/family.xls -m ../../$map -g ../../group.sort -w $familyGenus_w -h $temp_h ; mv ANOVA.pdf family.anova.pdf ; mv ANOVA.svg family.anova.svg ; mv ANOVA.png family.anova.png ; mv ANOVA.txt family.anova.txt\n";
		print OUP "perl /mnt/sdb/lgq/bin/tools/clound/04.ANOVA/ANOVA-T.test2024-09-24.pl -i ../../1.Taxa/normalize/tax_summary_a/genus.xls -m ../../$map -g ../../group.sort -w $familyGenus_w -h $temp_h ; mv ANOVA.pdf genus.anova.pdf ; mv ANOVA.svg genus.anova.svg ; mv ANOVA.png genus.anova.png ; mv ANOVA.txt genus.anova.txt ; rm z2.r zz.txt zz1.txt z.r color.txt rowsort ; cd ../../\n";

		print OUP "\n##kruskal wallis...\n";
		print OUP "cd 5.Diff/5-2.KW ; perl /mnt/sdb/lgq/bin/tools/clound/03.KW/Kruskal-Wallis_Wilcoxon_test2024-09-24.pl -i ../../1.Taxa/normalize/tax_summary_a/phylum.xls -m ../../$map -g ../../group.sort -w $temp_w -h $temp_h ; mv kw.pdf phylum.kw.pdf ; mv kw.svg phylum.kw.svg ; mv kw.png phylum.kw.png ; mv kw.txt phylum.kw.txt\n";
		print OUP "perl /mnt/sdb/lgq/bin/tools/clound/03.KW/Kruskal-Wallis_Wilcoxon_test2024-09-24.pl -i ../../1.Taxa/normalize/tax_summary_a/class.xls -m ../../$map -g ../../group.sort -w $temp_w -h $temp_h ; mv kw.pdf class.kw.pdf ; mv kw.svg class.kw.svg ; mv kw.png class.kw.png ; mv kw.txt class.kw.txt\n";
		print OUP "perl /mnt/sdb/lgq/bin/tools/clound/03.KW/Kruskal-Wallis_Wilcoxon_test2024-09-24.pl -i ../../1.Taxa/normalize/tax_summary_a/order.xls -m ../../$map -g ../../group.sort -w $temp_w -h $temp_h ; mv kw.pdf order.kw.pdf ; mv kw.svg order.kw.svg ; mv kw.png order.kw.png ; mv kw.txt order.kw.txt\n";
		print OUP "perl /mnt/sdb/lgq/bin/tools/clound/03.KW/Kruskal-Wallis_Wilcoxon_test2024-09-24.pl -i ../../1.Taxa/normalize/tax_summary_a/family.xls -m ../../$map -g ../../group.sort -w $familyGenus_w -h $temp_h ; mv kw.pdf family.kw.pdf ; mv kw.svg family.kw.svg ; mv kw.png family.kw.png ; mv kw.txt family.kw.txt\n";
		print OUP "perl /mnt/sdb/lgq/bin/tools/clound/03.KW/Kruskal-Wallis_Wilcoxon_test2024-09-24.pl -i ../../1.Taxa/normalize/tax_summary_a/genus.xls -m ../../$map -g ../../group.sort -w $familyGenus_w -h $temp_h ; mv kw.pdf genus.kw.pdf ; mv kw.svg genus.kw.svg ; mv kw.png genus.kw.png ; mv kw.txt genus.kw.txt ; rm z2.r zz.txt zz1.txt z.r color.txt rowsort ; cd ../../\n";
	}

	print OUP "\n##Lefse...\n";
	print OUP "perl $Bin/bin/format_tax_forLefse.pl -i 1.Taxa/normalize/tax_summary_a -m $map -g Group -o 5.Diff/5-3.Lefse -t $pre\n";
	print OUP "perl $Bin/bin/filter_lefse.pl 5.Diff/5-3.Lefse/lefse_LDA.xls 5.Diff/5-3.Lefse/lefse_LDA.known.xls ; mv 5.Diff/5-3.Lefse/lefse_LDA.known.xls 5.Diff/5-3.Lefse/lefse_LDA.xls\n";
	print OUP "#perl $Bin/bin/filter_diffLefse.pl 5.Diff/5-3.Lefse/lefse_LDA.xls 5.Diff/5-3.Lefse/lefse_LDA.xls2 ; mv 5.Diff/5-3.Lefse/lefse_LDA.xls2 5.Diff/5-3.Lefse/lefse_LDA.xls\n";
	print OUP "cd 5.Diff/5-3.Lefse ; plot_res.py lefse_LDA.xls lefse_LDA.pdf --dpi 300 --format pdf --width 8\n";
	print OUP "plot_cladogram.py lefse_LDA.xls lefse_LDA.cladogram.pdf --format pdf ; cd ../../\n";


	print OUP "\n\n#module prediction...\n";
	`mkdir $outDir/6.Module` unless(-e "$outDir/6.Module");
	`mkdir $outDir/6.Module/6-1.Random_forest` unless(-e "$outDir/6.Module/6-1.Random_forest");
	print OUP "\n##random forest...\n";
	if($groupNumber == 2){
		$temp_h = 7;
	}elsif($groupNumber == 3){
		$temp_h = 8;
	}elsif($groupNumber == 4){
		$temp_h = 9;
	}elsif($groupNumber < 10){
		$temp_h = 10;
	}else{
		$temp_h = $groupNumber * 1.1;
	}
	print OUP "cd 6.Module/6-1.Random_forest ; perl /mnt/sdb/lgq/bin/tools/clound/14.randomForest/randomForest2024-07-11.pl -i ../../1.Taxa/normalize/tax_summary_a/phylum.xls -m ../../$map -g ../../group.sort -test kw -rlc 12 -clc 12 -llc 12 -h $temp_h ; mv randomForest.MeanDecreaseAccuracy.all.xls phylum.MDA.xls ; mv MeanDecreaseAccuracy_top_20.pdf phylum.top20MDA.pdf ; mv MeanDecreaseAccuracy_top_20.svg phylum.top20MDA.svg ; mv MeanDecreaseAccuracy_top_20.png phylum.top20MDA.png ; mv diff_top_20_MeanDecreaseAccuracy.pdf phylum.diff.top20MDA.pdf ; mv diff_top_20_MeanDecreaseAccuracy.svg phylum.diff.top20MDA.svg ; mv diff_top_20_MeanDecreaseAccuracy.png phylum.diff.top20MDA.png\n";
	print OUP "perl /mnt/sdb/lgq/bin/tools/clound/14.randomForest/randomForest2024-07-11.pl -i ../../1.Taxa/normalize/tax_summary_a/class.xls -m ../../$map -g ../../group.sort -test kw -rlc 12 -clc 12 -llc 12 -h $temp_h ; mv randomForest.MeanDecreaseAccuracy.all.xls class.MDA.xls ; mv MeanDecreaseAccuracy_top_20.pdf class.top20MDA.pdf ; mv MeanDecreaseAccuracy_top_20.svg class.top20MDA.svg ; mv MeanDecreaseAccuracy_top_20.png class.top20MDA.png ; mv diff_top_20_MeanDecreaseAccuracy.pdf class.diff.top20MDA.pdf ; mv diff_top_20_MeanDecreaseAccuracy.svg class.diff.top20MDA.svg ; mv diff_top_20_MeanDecreaseAccuracy.png class.diff.top20MDA.png\n";
	print OUP "perl /mnt/sdb/lgq/bin/tools/clound/14.randomForest/randomForest2024-07-11.pl -i ../../1.Taxa/normalize/tax_summary_a/order.xls -m ../../$map -g ../../group.sort -test kw -rlc 12 -clc 12 -llc 12 -w 13 -h $temp_h ; mv randomForest.MeanDecreaseAccuracy.all.xls order.MDA.xls ; mv MeanDecreaseAccuracy_top_20.pdf order.top20MDA.pdf ; mv MeanDecreaseAccuracy_top_20.svg order.top20MDA.svg ; mv MeanDecreaseAccuracy_top_20.png order.top20MDA.png ; mv diff_top_20_MeanDecreaseAccuracy.pdf order.diff.top20MDA.pdf ; mv diff_top_20_MeanDecreaseAccuracy.svg order.diff.top20MDA.svg ; mv diff_top_20_MeanDecreaseAccuracy.png order.diff.top20MDA.png\n";
	print OUP "perl /mnt/sdb/lgq/bin/tools/clound/14.randomForest/randomForest2024-07-11.pl -i ../../1.Taxa/normalize/tax_summary_a/family.xls -m ../../$map -g ../../group.sort -test kw -rlc 12 -clc 12 -llc 12 -w 13 -h $temp_h ; mv randomForest.MeanDecreaseAccuracy.all.xls family.MDA.xls ; mv MeanDecreaseAccuracy_top_20.pdf family.top20MDA.pdf ; mv MeanDecreaseAccuracy_top_20.svg family.top20MDA.svg ; mv MeanDecreaseAccuracy_top_20.png family.top20MDA.png ; mv diff_top_20_MeanDecreaseAccuracy.pdf family.diff.top20MDA.pdf ; mv diff_top_20_MeanDecreaseAccuracy.svg family.diff.top20MDA.svg ; mv diff_top_20_MeanDecreaseAccuracy.png family.diff.top20MDA.png\n";
	print OUP "perl /mnt/sdb/lgq/bin/tools/clound/14.randomForest/randomForest2024-07-11.pl -i ../../1.Taxa/normalize/tax_summary_a/genus.xls -m ../../$map -g ../../group.sort -test kw -rlc 12 -clc 12 -llc 12 -w 13 -h $temp_h ; mv randomForest.MeanDecreaseAccuracy.all.xls genus.MDA.xls ; mv MeanDecreaseAccuracy_top_20.pdf genus.top20MDA.pdf ; mv MeanDecreaseAccuracy_top_20.svg genus.top20MDA.svg ; mv MeanDecreaseAccuracy_top_20.png genus.top20MDA.png ; mv diff_top_20_MeanDecreaseAccuracy.pdf genus.diff.top20MDA.pdf ; mv diff_top_20_MeanDecreaseAccuracy.svg genus.diff.top20MDA.svg ; mv diff_top_20_MeanDecreaseAccuracy.png genus.diff.top20MDA.png ; rm diff_top_20_MeanDecreaseAccuracy.xls MeanDecreaseAccuracy_top_20.xls color.txt ; cd ../../\n";
}

if($typeFlag eq "T"){
	print OUP "\n\n###picrust2 function predict...\n";
	`mkdir $outDir/7.Picrust2` unless(-e "$outDir/7.Picrust2");
	`mkdir $outDir/7.Picrust2/7-1.KEGG` unless(-e "$outDir/7.Picrust2/7-1.KEGG");
	`mkdir $outDir/7.Picrust2/7-2.COG` unless(-e "$outDir/7.Picrust2/7-2.COG");
	`mkdir $outDir/7.Picrust2/7-3.Bar` unless(-e "$outDir/7.Picrust2/7-3.Bar");
	`mkdir $outDir/7.Picrust2/7-3.Bar/COG` unless(-e "$outDir/7.Picrust2/7-3.Bar/COG");
	`mkdir $outDir/7.Picrust2/7-3.Bar/KEGG` unless(-e "$outDir/7.Picrust2/7-3.Bar/KEGG");
	`mkdir $outDir/7.Picrust2/7-4.Heatmap` unless(-e "$outDir/7.Picrust2/7-4.Heatmap");
	`mkdir $outDir/7.Picrust2/7-4.Heatmap/COG` unless(-e "$outDir/7.Picrust2/7-4.Heatmap/COG");
	`mkdir $outDir/7.Picrust2/7-4.Heatmap/KEGG` unless(-e "$outDir/7.Picrust2/7-4.Heatmap/KEGG");
	print OUP "cd 7.Picrust2 ; sh /mnt/sdb/lgq/bin/tools/picrust/PICRUSt2/PICRUSt2.pipe.sh -i ../1.Taxa/normalize/$pre\_table.xls -s ../1.Taxa/normalize/$pre\_rep.fasta\n";
	print OUP "perl /mnt/sdb/lgq/bin/tools/picrust/PICRUSt2/PICRUSt2_Result_format.pl picrust2_out_pipeline ./ ; cd ../\n";
	print OUP "perl $Bin/bin/sortBySam.pl 7.Picrust2/PICRUSt2_KEGG/pred_KO_profile.xls $samSortFile 0 7.Picrust2/7-1.KEGG/orthology.xls\n";
	print OUP "perl $Bin/bin/sortBySam.pl 7.Picrust2/PICRUSt2_KEGG/pred_KO_relative_profile.xls $samSortFile 0 7.Picrust2/7-1.KEGG/orthology.percent.xls\n";
	print OUP "perl $Bin/bin/sortBySam.pl 7.Picrust2/PICRUSt2_KEGG/pred_enzyme_profile.xls $samSortFile 1 7.Picrust2/7-1.KEGG/enzyme.xls\n";
	print OUP "perl $Bin/bin/sortBySam.pl 7.Picrust2/PICRUSt2_KEGG/pred_enzyme_relative_profile.xls $samSortFile 1 7.Picrust2/7-1.KEGG/enzyme.percent.xls\n";
	print OUP "perl $Bin/bin/sortBySam.pl 7.Picrust2/PICRUSt2_KEGG/pred_kegg_pathway_profile.xls $samSortFile 2 7.Picrust2/7-1.KEGG/pathway.xls\n";
	print OUP "perl $Bin/bin/sortBySam.pl 7.Picrust2/PICRUSt2_KEGG/pred_kegg_pathway_relative_profile.xls $samSortFile 2 7.Picrust2/7-1.KEGG/pathway.percent.xls\n";
	print OUP "perl $Bin/bin/sortBySam.pl 7.Picrust2/PICRUSt2_KEGG/pred_kegg_pathwayL3_profile.xls $samSortFile 0 7.Picrust2/7-1.KEGG/L3.xls\n";
	print OUP "perl $Bin/bin/sortBySam.pl 7.Picrust2/PICRUSt2_KEGG/pred_kegg_pathwayL3_relative_profile.xls $samSortFile 0 7.Picrust2/7-1.KEGG/L3.percent.xls\n";
	print OUP "perl $Bin/bin/sortBySam.pl 7.Picrust2/PICRUSt2_KEGG/pred_kegg_pathwayL2_profile.xls $samSortFile 0 7.Picrust2/7-1.KEGG/L2.xls\n";
	print OUP "perl $Bin/bin/sortBySam.pl 7.Picrust2/PICRUSt2_KEGG/pred_kegg_pathwayL2_relative_profile.xls $samSortFile 0 7.Picrust2/7-1.KEGG/L2.percent.xls\n";
	print OUP "perl $Bin/bin/sortBySam.pl 7.Picrust2/PICRUSt2_KEGG/pred_kegg_pathwayL1_profile.xls $samSortFile 0 7.Picrust2/7-1.KEGG/L1.xls\n";
	print OUP "perl $Bin/bin/sortBySam.pl 7.Picrust2/PICRUSt2_KEGG/pred_kegg_pathwayL1_relative_profile.xls $samSortFile 0 7.Picrust2/7-1.KEGG/L1.percent.xls\n";
	print OUP "perl $Bin/bin/sortBySam.pl 7.Picrust2/PICRUSt2_KEGG/pred_kegg_pathwayL123_profile.xls $samSortFile 2 7.Picrust2/7-1.KEGG/L123.xls\n";
	print OUP "perl $Bin/bin/sortBySam.pl 7.Picrust2/PICRUSt2_KEGG/pred_kegg_pathwayL123_relative_profile.xls $samSortFile 2 7.Picrust2/7-1.KEGG/L123.percent.xls\n";

	print OUP "perl $Bin/bin/sortBySam.pl 7.Picrust2/PICRUSt2_COG/pred_cog_description_profile.xls $samSortFile 1 7.Picrust2/7-2.COG/COG.xls\n";
	print OUP "perl $Bin/bin/sortBySam.pl 7.Picrust2/PICRUSt2_COG/pred_cog_description_relative_profile.xls $samSortFile 1 7.Picrust2/7-2.COG/COG.percent.xls\n";
	print OUP "perl $Bin/bin/sortBySam.pl 7.Picrust2/PICRUSt2_COG/pred_cog_category_profile.xls $samSortFile 1 7.Picrust2/7-2.COG/category.xls\n";
	print OUP "perl $Bin/bin/sortBySam.pl 7.Picrust2/PICRUSt2_COG/pred_cog_category_relative_profile.xls $samSortFile 1 7.Picrust2/7-2.COG/category.percent.xls\n";

	print OUP "\n##barplot...\n";

	print OUP "cp 7.Picrust2/PICRUSt2_COG/bar/*pdf 7.Picrust2/7-3.Bar/COG ; rename cog.pdf category.pdf 7.Picrust2/7-3.Bar/COG/*pdf ; mv 7.Picrust2/7-3.Bar/COG/cog.box.pdf 7.Picrust2/7-3.Bar/COG/box.pdf ; rm 7.Picrust2/cmd.r 7.Picrust2/Rplots.pdf 7.Picrust2/pst2.sh\n";
	if($samNumber < 10){
		$temp_w = 5;$temp_h = 6;$temp_legend = 5;
	}elsif($samNumber < 30){
		$temp_w = 6;$temp_h = 6;$temp_legend = 6;
	}elsif($samNumber < 50){
		$temp_w = 8;$temp_h = 8;$temp_legend = 8;
	}else{
		$temp_w = $samNumber/50 * 8;$temp_h = 8;$temp_legend = int($temp_w);
	}
	print OUP "cd 7.Picrust2/7-3.Bar/KEGG ; perl /mnt/sdb/lgq/bin/tools/clound/02.bar/bar_pie2024-08-02.pl -d ./ -i ../../7-1.KEGG/orthology.xls -filter top -top 30 -t T -b_blas 2 -pie F -tp 1 -tv 0.5 -b_lcex 1 -b_ncol $temp_legend -b_w $temp_w -b_h $temp_h ; mv bar_link.pdf bar.orthology.pdf ; mv bar_link.svg bar.orthology.svg ; mv bar_link.png bar.orthology.png\n";
	if($samNumber < 10){
		$temp_w = 6;$temp_h = 6;$temp_legend = 4;
	}elsif($samNumber < 30){
		$temp_w = 7;$temp_h = 7;$temp_legend = 5;
	}elsif($samNumber < 50){
		$temp_w = 8;$temp_h = 8;$temp_legend = 6;
	}else{
		$temp_w = $samNumber/50*8;$temp_h = 8;$temp_legend = int($temp_w) - 2;
	}
	print OUP "perl $Bin/bin/remove_tax.pl ../../7-1.KEGG/enzyme.xls enzyme.xls ; perl /mnt/sdb/lgq/bin/tools/clound/02.bar/bar_pie2024-08-02.pl -d ./ -i enzyme.xls -filter top -top 30 -t T -b_blas 2 -pie F -tp 1 -tv 0.5 -b_lcex 1 -b_ncol $temp_legend -b_w $temp_w -b_h $temp_h ; mv bar_link.pdf bar.enzyme.pdf ; mv bar_link.svg bar.enzyme.svg ; mv bar_link.png bar.enzyme.png ; rm cmd.r out.txt all.bar bar.pdf bar.svg bar.png ; cd ../../../\n";

	print OUP "\n##heatmap...\n";

	if($samNumber <= 10){
		$temp_w = 10;
	}elsif($samNumber <= 30){
		$temp_w = $samNumber/10 * 1.1 + 5;$temp_w = int($temp_w);
	}else{
		$temp_w = $samNumber/10 * 1.1 + 6;$temp_w = int($temp_w);
	}

	print OUP "cd 7.Picrust2/7-4.Heatmap/KEGG ; perl /mnt/sdb/lgq/bin/tools/clound/01.heatmap/plot-heatmap2024-04-17_3.pl -d ./ -i ../../7-1.KEGG/orthology.xls -o heatmap.orthology -ct 0 -rt 0 -lg log2 -cs 2 -rtop 50 -slas 2 -rlc 0.9 -w $temp_w -marble 4-1-0-3\n";
	print OUP "perl /mnt/sdb/lgq/bin/tools/clound/01.heatmap/plot-heatmap2024-04-17_3.pl -d ./ -i ../../7-3.Bar/KEGG/enzyme.xls -o heatmap.enzyme -ct 0 -rt 0 -lg log2 -cs 2 -rtop 50 -slas 2 -rlc 0.9 -w $temp_w -marble 4-1-0-4 ; rm ../../7-3.Bar/KEGG/enzyme.xls cmd.r cmd.r.log warnings.log final.txt ; cd ../../../\n";
	print OUP "cd 7.Picrust2/7-4.Heatmap/COG ; perl $Bin/bin/remove_tax.pl ../../7-2.COG/COG.xls COG.xls ; perl $Bin/bin/remove_tax.pl ../../7-2.COG/category.xls category.xls ; perl /mnt/sdb/lgq/bin/tools/clound/01.heatmap/plot-heatmap2024-04-17_3.pl -d ./ -i COG.xls -o heatmap.COG -ct 0 -rt 0 -lg log2 -cs 2 -rtop 50 -slas 2 -rlc 0.9 -w $temp_w -marble 4-1-0-3\n";
	print OUP "perl /mnt/sdb/lgq/bin/tools/clound/01.heatmap/plot-heatmap2024-04-17_3.pl -d ./ -i category.xls -o heatmap.category -ct 0 -rt 0 -lg log2 -cs 2 -rtop 50 -slas 2 -rlc 0.9 -w $temp_w -marble 4-1-0-3 ; rm COG.xls category.xls cmd.r cmd.r.log warnings.log final.txt ; cd ../../../\n";
}
	close OUP;


	open OUR,"> $outDir/rm.sh" or die "$!\n";
	print OUR "rm sub.txt temp.fasta\n";
#print OUR "rm rawData.stat.xls cleanData.stat.xls\n";
	print OUR "rm $pre\_taxa_table.sub.xls\n";
	print OUR "rm -r 2.Alpha/Alpha\n";
	print OUR "rm -r 7.Picrust2/picrust2_out_pipeline/ 7.Picrust2/PICRUSt2_COG/ 7.Picrust2/PICRUSt2_KEGG/\n" if($typeFlag eq "T");
	print OUR "rm 2.Alpha/2-2.Boxplot/*svg 2.Alpha/2-2.Boxplot/*png\n" if($G3flag eq "T");
	print OUR "rm 2.Alpha/2-3.Rarefaction/*svg 2.Alpha/2-3.Rarefaction/*png\n";
	print OUR "rm 2.Alpha/2-4.Rank-abundance/*svg 2.Alpha/2-4.Rank-abundance/*png\n";
	print OUR "rm 2.Alpha/2-5.Core-pan/*svg 2.Alpha/2-5.Core-pan/*png\n";
	print OUR "rm 3.Community/3-2.Venn/*svg 3.Community/3-2.Venn/*png\n";
	print OUR "rm 3.Community/3-2.Venn.Group/*svg 3.Community/3-2.Venn.Group/*png\n" if($map);
	print OUR "rm 3.Community/3-3.Bar/*svg 3.Community/3-3.Bar/*png\n";
	print OUR "rm 3.Community/3-3.Bar.Group/*svg 3.Community/3-3.Bar.Group/*png\n" if($map);
	print OUR "rm 3.Community/3-4.Bubble/*svg 3.Community/3-4.Bubble/*png\n";
	print OUR "rm 3.Community/3-4.Bubble.Group/*svg 3.Community/3-4.Bubble.Group/*png\n" if($map);
	print OUR "rm 3.Community/3-5.Heatmap/*svg 3.Community/3-5.Heatmap/*png\n";
	print OUR "rm 3.Community/3-5.Heatmap.Group/*svg 3.Community/3-5.Heatmap.Group/*png\n" if($map);
	print OUR "rm 3.Community/3-6.TreeBar/*svg 3.Community/3-6.TreeBar/*png\n";
	print OUR "rm 3.Community/3-6.TreeBar.Group/*svg 3.Community/3-6.TreeBar.Group/*png\n" if($map);
	print OUR "rm 3.Community/3-7.topGenus/*svg 3.Community/3-7.topGenus/*png\n";
	print OUR "rm 4.Beta/4-1.Distance/*/*svg 4.Beta/4-1.Distance/*/*png\n";
	print OUR "rm 4.Beta/4-2.PCA/*svg 4.Beta/4-2.PCA/*png\n";
	print OUR "rm 4.Beta/4-3.PCOA/*/*svg 4.Beta/4-3.PCOA/*/*png\n";
	print OUR "rm 4.Beta/4-4.NMDS/*/*svg 4.Beta/4-4.NMDS/*/*png\n";
	print OUR "rm 4.Beta/4-5.Hcluster_tree/*/*svg 4.Beta/4-5.Hcluster_tree/*/*png\n";
	print OUR "rm 4.Beta/4-6.Anosim/*/*svg 4.Beta/4-6.Anosim/*/*png\n" if($G3flag eq "T");
	print OUR "rm 5.Diff/5-1.*/*svg 5.Diff/5-1.*/*png\n" if($G3flag eq "T");
	print OUR "rm 5.Diff/5-2.*/*svg 5.Diff/5-2.*/*png\n" if($G3flag eq "T");
	print OUR "rm 6.Module/6-1.Random_forest/*svg 6.Module/6-1.Random_forest/*png\n" if($G3flag eq "T");
	print OUR "rm 7.Picrust2/7-3.Bar/KEGG/*svg 7.Picrust2/7-3.Bar/KEGG/*png\n" if($typeFlag eq "T");
	print OUR "rm 7.Picrust2/7-4.Heatmap/*/*svg 7.Picrust2/7-4.Heatmap/*/*png\n" if($typeFlag eq "T");
	close OUR;

sub usage{
	print <<EOD
		usage: perl $0 -i otu/asv_taxa_table.xls -rep otu/asv_rep.fasta -m map.txt -rs rawData.stat.xls -cs cleanData.stat.xls -af ASV -gf T -g3f T -n number -o outDir
			-i*	asv(otu)_taxa_table.xls,required
			-r*	asv(otu)_rep.fasta,required
			-rs*	raw data stat,required
			-cs*	clean data stat,required
			-o*	out direction,required
			-m	map file,The header should be "#SampleID\tGroup"
			-af	analysis flag[OTU|ASV],default ASV
			-g3	flag of At least 3 samples per group,default F
			-t	16S data type,default T
			-n	number of sequences to subsample per sample,default Minimum sequences number of all samples
			-h	show this help
EOD
}
