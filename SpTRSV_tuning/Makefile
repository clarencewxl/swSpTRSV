MC = mpicc
#CC = sw5cc.old
CC = sw5cc

MCFLAGS = ${CFLAGS}
CCFLAGS = ${CFLAGS}
MCFLAGS+= -O2 -msimd -OPT:ieee_arith=1 -g -host
CCFLAGS+= -O2 -msimd -g -slave
LIBS    = -hybrid -lslave

MOBJS = main.o functions_R8B1024.o functions_R4B512.o functions_R2B64.o
COBJS = functions_levelsets_slave.o functions_R8B1024_slave.o functions_R4B512_slave.o functions_R2B64_slave.o

SwSpTRSV:$(MOBJS) $(COBJS)
	$(MC) $(LIBS) $(MOBJS) $(COBJS) -o SwSpTRSV 
	ctags -R ./*.c
functions_levelsets_slave.o:functions_levelsets_slave.c my_solver.h
	$(CC) $(CCFLAGS) -c $< -o $@
functions_slave.o:functions_slave.c my_solver.h
	$(CC) $(CCFLAGS) -c $< -o $@
functions_R8B1024_slave.o:functions_R8B1024_slave.c my_solver.h
	$(CC) $(CCFLAGS) -c $< -o $@
functions_R4B512_slave.o:functions_R4B512_slave.c my_solver.h
	$(CC) $(CCFLAGS) -c $< -o $@
functions_R2B64_slave.o:functions_R2B64_slave.c my_solver.h
	$(CC) $(CCFLAGS) -c $< -o $@
%.o:%.c
	$(MC) $(MCFLAGS) -c $< -o $@

clean:
	rm -rf *.o SwSpTRSV
clear:
	rm -rf *.o SwSpTRSV ./results/*

run: SwSpTRSV
	./run.sh -I boyd1
	#./run.sh -I Ga3As3H12


