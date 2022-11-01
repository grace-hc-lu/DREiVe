# Description file for the "patser-v3" program.

###################################################
# Explicitly define the compiling rules assumed
COMPILE.c=$(CC) $(CFLAGS) $(CPPFLAGS) -c
.c.o:
	$(COMPILE.c) $(OUTPUT_OPTION) $<
###################################################

CFLAGS = -O2         # The basic version (optimized version)
#CFLAGS = -O2 -p     # Version for profiling function use (optimized version)
#CFLAGS = -g         # The version for debugging (limited optimization)

CPPFLAGS = -DOPTIONS   # Causes certain files to include options.h
#CPPFLAGS = -DTEST -DOPTIONS   # Defined only when patser-v3d.test is compiled
SGI = cc -O3 -cckr

#CC = cc
CC = gcc

OBJ = options.o main.o command-line.o parse-line.o parse-line-help.o\
	text-directions.o alpha.o process-matrix.o read-matrix.o\
	det-marginal-prob.o read-sequence.o process-sequence.o\
	heap.o alloc-error.o
OBJ_TEST = cutoff-exact.o cutoff-analytical.o gamma-prob.o root-find.o

SRC = options.c main.c command-line.c parse-line.c parse-line-help.c\
	text-directions.c alpha.c process-matrix.c read-matrix.c\
	det-marginal-prob.c read-sequence.c process-sequence.c\
	heap.c alloc-error.c
SRC_TEST = cutoff-exact.c cutoff-analytical.c gamma-prob.c root-find.c

SRC_H = definitions.h alloc-error.h options.h parse-line.h alpha.h heap.h


patser-v3d : $(OBJ)
	$(CC) $(CFLAGS) -o patser-v3d $(OBJ) -lm
lint :
	lint $(CPPFLAGS) $(SRC)
tar :
	tar -chvf patser-v3d.tar directions makefile $(SRC_H) $(SRC)
patser-v3d.sgi : $(SRC)
	$(SGI) -o patser-v3d $(SRC) -lm


patser-v3d.test : $(OBJ) $(OBJ_TEST)
	$(CC) $(CFLAGS) -o patser-v3d.test $(OBJ) $(OBJ_TEST) -lm
lint.test :
	lint $(CPPFLAGS) $(SRC) $(SRC_TEST)
tar.test :
	tar -chvf patser-v3d.test.tar\;
		directions makefile $(SRC_H) $(SRC) $(SRC_TEST)


options.o : definitions.h
main.o : definitions.h alpha.h alloc-error.h options.h
command-line.o : definitions.h alpha.h alloc-error.h options.h parse-line.h
parse-line.o : parse-line.h
parse-line-help.o : parse-line.h
text-directions.o :
alpha.o : definitions.h alpha.h alloc-error.h options.h
process-matrix.o : definitions.h alpha.h alloc-error.h options.h
read-matrix.o : definitions.h alpha.h alloc-error.h options.h
det-marginal-prob.o : definitions.h alpha.h alloc-error.h options.h
cutoff-exact.o : definitions.h alpha.h alloc-error.h options.h
cutoff-analytical.o : definitions.h alpha.h alloc-error.h options.h
read-sequence.o : definitions.h alpha.h alloc-error.h options.h
process-sequence.o : definitions.h alpha.h alloc-error.h options.h heap.h
gamma-prob.o : definitions.h alpha.h alloc-error.h options.h
heap.o : heap.h alloc-error.h
root-find.o : 
alloc-error.o :
