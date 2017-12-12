#########################################################################
# File Name: makeall.sh
# Author: Wang Xinliang
# mail: clarencewxl@gmail.com
# Created Time: 2017年12月12日 星期二 10时49分31秒
#########################################################################
#!/bin/bash

rm -rf *.o; make PCR=8 PCRL=3 CXV4=1024 CXV4L=10
rm -rf *.o; make PCR=8 PCRL=3 CXV4=512  CXV4L=9
rm -rf *.o; make PCR=8 PCRL=3 CXV4=256  CXV4L=8
rm -rf *.o; make PCR=8 PCRL=3 CXV4=128  CXV4L=7
rm -rf *.o; make PCR=8 PCRL=3 CXV4=64   CXV4L=6
rm -rf *.o; make PCR=8 PCRL=3 CXV4=32   CXV4L=5
rm -rf *.o; make PCR=4 PCRL=2 CXV4=1024 CXV4L=10
rm -rf *.o; make PCR=4 PCRL=2 CXV4=512  CXV4L=9
rm -rf *.o; make PCR=4 PCRL=2 CXV4=256  CXV4L=8
rm -rf *.o; make PCR=4 PCRL=2 CXV4=128  CXV4L=7
rm -rf *.o; make PCR=4 PCRL=2 CXV4=64   CXV4L=6
rm -rf *.o; make PCR=4 PCRL=2 CXV4=32   CXV4L=5
rm -rf *.o; make PCR=2 PCRL=1 CXV4=1024 CXV4L=10
rm -rf *.o; make PCR=2 PCRL=1 CXV4=512  CXV4L=9
rm -rf *.o; make PCR=2 PCRL=1 CXV4=256  CXV4L=8
rm -rf *.o; make PCR=2 PCRL=1 CXV4=128  CXV4L=7
rm -rf *.o; make PCR=2 PCRL=1 CXV4=64   CXV4L=6
rm -rf *.o; make PCR=2 PCRL=1 CXV4=32   CXV4L=5
rm -rf *.o; make PCR=1 PCRL=0 CXV4=1024 CXV4L=10
rm -rf *.o; make PCR=1 PCRL=0 CXV4=512  CXV4L=9
rm -rf *.o; make PCR=1 PCRL=0 CXV4=256  CXV4L=8
rm -rf *.o; make PCR=1 PCRL=0 CXV4=128  CXV4L=7
rm -rf *.o; make PCR=1 PCRL=0 CXV4=64   CXV4L=6
rm -rf *.o; make PCR=1 PCRL=0 CXV4=32   CXV4L=5

