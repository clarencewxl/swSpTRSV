#########################################################################
# File Name: wxl_sleep.sh
# Author: Wang Xinliang
# mail: clarencetc@163.com
# Created Time: 2017年06月09日 星期五 12时29分28秒
#########################################################################
#!/bin/bash
while true
do
    size=`bjobs | grep SpTRSV | wc -l`
    if [ "$size" -lt "$1" ];then
        break
    fi
    sleep 15s
    echo "Sleep 15s"
done


