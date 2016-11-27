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
all: umimerge_parallel umimerge umimerge_single_pass umisplit_orig umisample umisplit umisplit_sam umifilter 
all96: clean directories96 umimerge_parallel umimerge umimerge_single_pass umisplit_orig umisample umisplit umisplit_sam umifilter
all384: clean directories384 umimerge_parallel umimerge umimerge_single_pass umisplit_orig umisample umisplit umisplit_sam umifilter

OBJ_umifilter = umifilter.o optparse.o 
OBJ_umisplit = umisplit.o optparse.o 
OBJ_umisplit_sam = umisplit_sam.o optparse.o 
OBJ_umimerge = umimerge.o optparse.o 
OBJ_umimerge_single_pass = umimerge_single_pass.o optparse.o  
OBJ_umisplit_orig = umisplit_orig.o optparse.o 
OBJ_umisample = umisample.o optparse.o 
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
umimerge : $(OBJ_umimerge)
	$(LD) -o ${DIR}umimerge $(LFLAGS) $(OBJ_umimerge) $(LIBS)
umimerge.o : umimerge.cpp 
	$(CPP) $(CFLAGS) -c umimerge.cpp		
umimerge_single_pass : $(OBJ_umimerge_single_pass)
	$(LD) -o ${DIR}umimerge_single_pass $(LFLAGS) $(OBJ_umimerge_single_pass) $(LIBS)
umimerge_single_pass.o : umimerge_single_pass.cpp 
	$(CPP) $(CFLAGS) -c umimerge_single_pass.cpp			

umisplit : $(OBJ_umisplit)
	$(LD) -o ${DIR}umisplit $(LFLAGS) $(OBJ_umisplit) $(LIBS)
umisplit.o : umisplit.cpp 
	$(CPP) $(CFLAGS) -c umisplit.cpp			

umisplit_sam : $(OBJ_umisplit_sam)
	$(LD) -o ${DIR}umisplit_sam $(LFLAGS) $(OBJ_umisplit_sam) $(LIBS)
umisplit_sam.o : umisplit_sam.cpp 
	$(CPP) $(CFLAGS) -c umisplit_sam.cpp			
	
umisplit_orig : $(OBJ_umisplit_orig)
	$(LD) -o ${DIR}umisplit_orig $(LFLAGS) $(OBJ_umisplit_orig) $(LIBS)
umisplit_orig.o : umisplit_orig.cpp 
	$(CPP) $(CFLAGS) -c umisplit_orig.cpp		
	
umisample : $(OBJ_umisample)
	$(LD) -o ${DIR}umisample $(LFLAGS) $(OBJ_umisample) $(LIBS)
umisample.o : umisample.cpp 
	$(CPP) $(CFLAGS) -c umisample.cpp
optparse.o : optparse/optparse.c
	$(CC) $(OCFLAGS) -c optparse/optparse.c
#endif
clean:
	-@rm  umimerge_parallel umimerge umimerge_single_pass umisplit_orig umisample umisplit umisplit_sam
	-@rm -rf *.o 2>/dev/null || true
	-@rm -rf optparse/*.o 2>/dev/null || true
	-@rm -rf core.* 2>/dev/null || true
