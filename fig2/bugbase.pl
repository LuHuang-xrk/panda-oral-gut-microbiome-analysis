#!/usr/bin/perl -w
use FindBin qw/$Bin/;
die "perl $0 <meta.fasta>  <outdir> <map.txt>\n" if @ARGV != 3;
print "####get gg_13_5_fasta annotaion ...\n" ;
`pick_closed_reference_otus.py -i $ARGV[0] -r $Bin/gg_13_5.fasta -o $ARGV[1]/otus_w_tax/ -t $Bin/gg_13_5_taxonomy.txt `;
print "###biom to table ...\n";
`biom convert -i $ARGV[1]/otus_w_tax/otu_table.biom -o $ARGV[1]/otus_w_tax/otu_table.txt --table-type "OTU table" --to-tsv`;
print "###run bugbase workflow ...\n" ;
print "run.bugbase.r -i $ARGV[1]/otus_w_tax/otu_table.txt   -m $ARGV[2]  -c Group -o $ARGV[1]/bugbase_results\n"; 
#`run.bugbase.r -i $ARGV[1]/otus_w_tax/otu_table.txt   -m $ARGV[2]  -c Group -o $ARGV[1]/bugbase_results` ;
