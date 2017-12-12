#########################################################################
# File Name: run.sh
# Author: Wang Xinliang
# mail: clarencewxl@gmail.com
# Created Time: Thu 29 Sep 2016 02:41:40 PM CST
#########################################################################
#!/bin/bash

count=4
if [ "$#" -eq "$count" ];then
bsub $1 -p -b -q q_sw_expr -n 1 -np 1 -cgsp 64 -J SpTRSV -o ./results/result_$2_$3_$4 -host_stack 256 -share_size 6500 -cross_size 0 ./SwSpTRSV_R$2B$3  ../triangular_files/$4.cscu
else
bsub -p -b -q q_sw_expr -n 1 -np 1 -cgsp 64 -J SpTRSV -o ./results/result_$1_$2_$3 -host_stack 256 -share_size 6500 -cross_size 0 ./SwSpTRSV_R$1B$2  ../triangular_files/$3.cscu
fi
