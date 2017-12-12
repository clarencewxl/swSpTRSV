MC = sw5cc
#CC = sw5cc.old
CC = sw5cc

CFLAGS = -DPRODUCER_CONSUMER_ROWS=${PCR}
CFLAGS+= -DPRODUCER_CONSUMER_ROWS_LOG2=${PCRL}
CFLAGS+= -DCACHE_X_V4=${CXV4}
CFLAGS+= -DCACHE_X_V4_LOG2=${CXV4L}
MCFLAGS = ${CFLAGS}
CCFLAGS = ${CFLAGS}
MCFLAGS+= -O2 -msimd -OPT:ieee_arith=1 -g -host
CCFLAGS+= -O2 -msimd -g -slave
LIBS    = -hybrid -lslave

MOBJS = main.o functions.o
COBJS = functions_slave.o

SwSpTRSV:$(MOBJS) $(COBJS)
	$(MC) $(LIBS) $(MOBJS) $(COBJS) -o SwSpTRSV_R${PCR}B${CXV4} 
	#ctags -R ./*.c
functions_slave.o:functions_slave.c my_solver.h
	$(CC) $(CCFLAGS) -c $< -o $@
%.o:%.c
	$(MC) $(MCFLAGS) -c $< -o $@

clean:
	rm -rf *.o SwSpTRSV_*
clear:
	rm -rf *.o SwSpTRSV_* ./results/* nohup.out

run: 
	./test.sh -I atmosmodd

