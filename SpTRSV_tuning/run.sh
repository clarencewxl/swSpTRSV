#########################################################################
# File Name: run.sh
# Author: Wang Xinliang
# mail: clarencewxl@gmail.com
# Created Time: Thu 29 Sep 2016 02:41:40 PM CST
#########################################################################
#!/bin/bash

count=2
if [ "$#" -eq "$count" ];then
bsub $1 -p -b -q q_sw_expr -n 1 -np 1 -cgsp 64 -o ./results/result_$2 -J SpTRSV -host_stack 256 -share_size 6500 -cross_size 0 ./SwSpTRSV  ../triangular_files/$2.cscu
else
bsub -p -b -q q_sw_expr -n 1 -np 1 -cgsp 64 -o ./results/result_$1 -J SpTRSV -host_stack 256 -share_size 6500 -cross_size 0 ./SwSpTRSV  ../triangular_files/$1.cscu
fi
