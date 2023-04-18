CC  = gcc
CXX = g++

SRW_DIR=      /home/erkn/Repositories/SRW
SRW_BASE_DIR= $(SRW_DIR)
SRW_CPP_DIR=	$(SRW_BASE_DIR)/cpp
SRW_SRC_DIR=	$(SRW_CPP_DIR)/src
SRW_SRC_GEN_DIR=	$(SRW_SRC_DIR)/core
SRW_SRC_LIB_DIR=	$(SRW_SRC_DIR)/lib
SRW_SRC_GENESIS_DIR=	$(SRW_SRC_DIR)/ext/genesis/genesis_july08
SH_SRC_PARSE_DIR=	$(SRW_SRC_DIR)/ext/auxparse
SH_SRC_GEN_MATH_DIR=	$(SRW_SRC_DIR)/ext/genmath
LIB_DIR=	$(SRW_BASE_DIR)/ext_lib

MCPL_LIB_DIR= /usr/share/mcstas/2.5/libs/mcpl

SRW_LIB_DIR=  $(SRW_DIR)/src/cpp/lib
SRW_EXT_LIB_DIR= $(SRW_DIR)/ext_lib
SRW_INC_DIR=  $(SRW_DIR)/src/cpp/lib

SRW_SRC_DEF=	-D_GNU_SOURCE -D__USE_XOPEN2K8 -DFFTW_ENABLE_FLOAT -D_GM_WITHOUT_BASE -DSRWLIB_STATIC -DNO_TIMER -DANSI_DECLARATORS -DTRILIBRARY -DLINUX
CFLAGS=	-O2 -fPIC -I$(SRW_SRC_GEN_DIR) -I$(SRW_SRC_LIB_DIR) -I$(SH_SRC_PARSE_DIR) -I$(SH_SRC_GEN_MATH_DIR) $(SRW_SRC_DEF) 

LDFLAGS=-L$(SRW_LIB_DIR) -L$(MCPL_LIB_DIR) -L$(LIB_DIR) -lm -lfftw 

OBJ=	auxparse.o gmfft.o gmfit.o gminterp.o gmmeth.o gmtrans.o srclcuti.o srcradint.o srctrjdt.o sremitpr.o srgsnbm.o srgtrjdt.o srisosrc.o srmagcnt.o srmagfld.o srmatsta.o sroptapt.o sroptcnt.o sroptdrf.o sroptel2.o sroptel3.o sroptelm.o sroptfoc.o sroptgrat.o sroptgtr.o sropthck.o sroptmat.o sroptpsh.o sroptshp.o sroptsmr.o sroptwgr.o sroptzp.o sroptzps.o srpersto.o srpowden.o srprdint.o srprgind.o srpropme.o srptrjdt.o srradinc.o srradint.o srradmnp.o srradstr.o srremflp.o srsase.o srsend.o srstowig.o srsysuti.o srthckbm.o srthckbm2.o srtrjaux.o srtrjdat.o srtrjdat3d.o all_com.o check.o diagno.o esource.o field.o incoherent.o initrun.o input.o loadbeam.o loadrad.o magfield.o main.o math.o mpi.o output.o partsim.o pushp.o rpos.o scan.o source.o stepz.o string.o tdepend.o timerec.o track.o	srerror.o srwlib.o 

PRG=	libsrw.a

libsrw.a: $(OBJ)
	ar -cvq $(PRG) *.o
	#cp $(PRG) $(LIB_DIR)/
	#rm -f *.o
	

lib:    libsrw.a

%.o: $(SRW_SRC_LIB_DIR)/%.cpp
	$(CC) $(CFLAGS) -c $<  

%.o: $(SH_SRC_PARSE_DIR)/%.cpp
	$(CC) $(CFLAGS) -c $<  

%.o: $(SH_SRC_GEN_MATH_DIR)/%.cpp
	$(CC) $(CFLAGS) -c $<  

%.o: $(SRW_SRC_GEN_DIR)/%.cpp
	$(CC) $(CFLAGS) -c $<  

%.o: $(SRW_SRC_GENESIS_DIR)/%.c
	$(CC) $(CFLAGS) -c $<  

srwlclient: srwlclient.cpp
	$(CXX) $(CFLAGS) -O3 -o srwlclient srwlclient.cpp libsrw.a $(LDFLAGS)

srwl2mcpl: srwl2mcpl.cpp
	$(CXX) $(CFLAGS) -I$(MCPL_LIB_DIR) -g -o srw2mcl srwl2mcpl.cpp libsrw.a $(LDFLAGS) -L$(MCPL_LIB_DIR) -lmcpl

clean:
	rm -f *.o *.so *.a srwlclient

all: lib srwl2mcpl srwclient 
