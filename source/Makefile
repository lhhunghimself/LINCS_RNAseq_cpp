CC=gcc
CPP=g++
LD=g++ -L /usr/local/lib
ifdef STATIC 
 STATIC_FLAG = -static
endif
ifdef DEBUG
CFLAGS= -I./ -I./optparse -ggdb -I -Wall -fopenmp -std=c++11 ${STATIC_FLAG}
else 
 ifdef AVX
  CFLAGS= -I./ -I./optparse -Ofast -std=c++11 -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-functions -fopenmp -m64 -mavx -msse3 -msse4 -frerun-loop-opt ${STATIC_FLAG} $(HDRS) $(DEFINES) 
 else
  ifdef SSE2
   ifdef 32Bit
    CFLAGS= -I./ -I./optparse -Ofast -std=c++11 -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-functions -fopenmp -frerun-loop-opt ${STATIC_FLAG} $(HDRS) $(DEFINES)
   else
    CFLAGS= -I./ -I./optparse -Ofast -std=c++11 -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-functions -fopenmp -m64 -frerun-loop-opt ${STATIC_FLAG} $(HDRS) $(DEFINES)
   endif
  else
   ifdef 32Bit
    CFLAGS= -I./ -I./optparse -Ofast -std=c++11 -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-functions -fopenmp -msse3 -frerun-loop-opt ${STATIC_FLAG} $(HDRS) $(DEFINES)
   else
    CFLAGS= -I./ -I./optparse -Ofast -std=c++11 -fexpensive-optimizations -funroll-all-loops -ffast-math -finline-functions -fopenmp -m64 -msse3 -frerun-loop-opt ${STATIC_FLAG} $(HDRS) $(DEFINES)
   endif 
  endif
 endif
endif
LIBS= -lboost_filesystem -lboost_system -lgomp -lm -lc -lgcc -lrt -lz
MKDIR_P = mkdir -p
#ifdef ALPINE
 ALPINE_FLAG= -fPIC -static
#endif
OCFLAGS = -std=c99 -Wall -Wextra -g3 ${ALPINE_FLAG}


.PHONY: clean all all96 all384 

default: all
all: umimerge_parallel umisplit cleanDeps
all96: clean directories96 umimerge_parallel umisplit
all384: clean directories384 umimerge_parallel umisplit 

OBJ_umisplit = umisplit.o optparse.o 
OBJ_umimerge_parallel = umimerge_parallel.o optparse.o 

directories96: 
	$(eval DIR+=w96/) 
	${MKDIR_P} ${DIR} 
directories384: 
	$(eval DIR+=w384/) 
	$(eval NWELLS+=w384/) 
	${MKDIR_P} ${DIR} 
umifilter : $(OBJ_umifilter)
	$(LD) -o ${DIR}umifilter $(LFLAGS) $(OBJ_umifilter) $(LIBS)
umifilter.o : umifilter.cpp 
	$(CPP) $(CFLAGS) -c umifilter.cpp	

umimerge_parallel : $(OBJ_umimerge_parallel)
	$(LD) -o ${DIR}umimerge_parallel $(LFLAGS) $(OBJ_umimerge_parallel) $(LIBS)
umimerge_parallel.o : umimerge_parallel.cpp 
	$(CPP) $(CFLAGS) -c umimerge_parallel.cpp			

umisplit : $(OBJ_umisplit)
	$(LD) -o ${DIR}umisplit $(LFLAGS) $(OBJ_umisplit) $(LIBS)
umisplit.o : umisplit.cpp 
	$(CPP) $(CFLAGS) -c umisplit.cpp			
	
optparse.o : optparse/optparse.c
	$(CC) $(OCFLAGS) -c optparse/optparse.c
#endif
clean: cleanBins cleanDeps

cleanDeps:
	-@rm -rf *.o 2>/dev/null || true
	-@rm -rf optparse/*.o 2>/dev/null || true
	-@rm -rf core.* 2>/dev/null || true

cleanBins: 
	-@rm  -f umimerge_parallel umisplit
