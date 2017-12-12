#########################################################################
# File Name: test.sh
# Author: Wang Xinliang
# mail: clarencewxl@gmail.com
# Created Time: Thu 29 Sep 2016 02:41:40 PM CST
#########################################################################
#!/bin/bash

count=2
if [ "$#" -eq "$count" ];then
./run.sh $1 8 1024 $2
./run.sh $1 8 512  $2
./run.sh $1 8 256  $2
./run.sh $1 8 128  $2
./run.sh $1 8 64   $2
./run.sh $1 8 32   $2
./run.sh $1 4 1024 $2
./run.sh $1 4 512  $2
./run.sh $1 4 256  $2
./run.sh $1 4 128  $2
./run.sh $1 4 64   $2
./run.sh $1 4 32   $2
./run.sh $1 2 1024 $2
./run.sh $1 2 512  $2
./run.sh $1 2 256  $2
./run.sh $1 2 128  $2
./run.sh $1 2 64   $2
./run.sh $1 2 32   $2
./run.sh $1 1 1024 $2
./run.sh $1 1 512  $2
./run.sh $1 1 256  $2
./run.sh $1 1 128  $2
./run.sh $1 1 64   $2
./run.sh $1 1 32   $2
else
./run.sh 8 1024 $1
./run.sh 8 512  $1
./run.sh 8 256  $1
./run.sh 8 128  $1
./run.sh 8 64   $1
./run.sh 8 32   $1
./run.sh 4 1024 $1
./run.sh 4 512  $1
./run.sh 4 256  $1
./run.sh 4 128  $1
./run.sh 4 64   $1
./run.sh 4 32   $1
./run.sh 2 1024 $1
./run.sh 2 512  $1
./run.sh 2 256  $1
./run.sh 2 128  $1
./run.sh 2 64   $1
./run.sh 2 32   $1
./run.sh 1 1024 $1
./run.sh 1 512  $1
./run.sh 1 256  $1
./run.sh 1 128  $1
./run.sh 1 64   $1
./run.sh 1 32   $1
fi
