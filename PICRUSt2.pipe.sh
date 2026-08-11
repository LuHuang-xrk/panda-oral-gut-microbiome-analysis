#!/usr/bin/env sh
################################################################
#Author : GuoqiLiu                                             #
#Date   : 2019-01-28                                           #
#Copyright (C) 2018~2019 precisiongenes.com.cn                 #
#Contact: liuguoqi@hmzkjy.cn                                   #
#Github  : https://github.com/liuguoqi/                        #
#Suppose:                                                      #
# 1. This shell script is module of cancer   annotation        #
# Platform : Centos 7.0 Perl v5.16.3 R 3.5.1 mysql 5.7.25      #
# sql and shell            requested                           #
# perl module :                                                #
#--------------------------------------------------------------# 
# DBI
#--------------------------------------------------------------#
################################################################
usage() {
      echo ""
      echo "Decription:"
      echo "usage : /mnt/sdb/lgq/bin/tools/picrust/PICRUSt2/PICRUSt2.pipe.sh"
      echo "This shell Script is PICRUSt2  pred   pipeline !!"
      echo ""
      echo ""
      echo -e "Option:"
      #echo "================================================================================================="
      echo "        USAGE:`basename $0` -i <otu> -s  <rep.seq>" 
      echo ""
      echo ""
      echo "       -i*    <str>     input otu table file"
      echo "       -s*    <str>     rep seq fasta"
      echo -e "\n\n" 
      echo "Author : liuguoqi"
      echo "Date   : 2018-12-11"
      echo "Contact: liuguoqi@hmzkjy.cn"
      echo ""
      exit 1
}
if [[ $# -lt 1 ]]
  then
     usage
  elif [[ $# != 4 ]] 
  then 
     usage

fi 
while getopts :i:s: OPTION;do
    case $OPTION in
    i)input=$OPTARG
    ;; 
    s)output=$OPTARG
    ;; 
    ?)usage
    ;; 
    esac
done

source  /mnt/sdb/lgq/bin/tools/picrust/PICRUSt2/config.sh

/mnt/sdb/lgq/bin/software/Mamba/bin/bin/picrust2_pipeline.py    -s ${output} -i ${input}  -o picrust2_out_pipeline --processes 64 --per_sequence_contrib  --coverage --skip_norm --in_traits COG,EC,KO,PFAM,TIGRFAM   --verbose 

# perl /mnt/sdb/lgq/bin/tools/picrust/PICRUSt2/PICRUSt2_Result_format.pl   picrust2_out_pipeline    PICRUSt2_Results 
cp /mnt/sdb/lgq/bin/tools/picrust/PICRUSt2/pst2.sh  .  
#sh pst2.sh 
##
#rm pst2.sh 
