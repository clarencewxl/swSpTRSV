#########################################################################
# File Name: total_run.sh
# Author: Wang Xinliang
# mail: clarencewxl@gmail.com
# Created Time: 2017年12月12日 星期二 20时29分24秒
#########################################################################
#!/bin/bash

count=0
while read name
do
    ./run.sh $name
    ((count=$count+1))
    if [ "$count" -eq "80" ]; then
        ./wxl_sleep.sh 10 ##### 
        count=0
    fi
done < ./benchmarks_list
