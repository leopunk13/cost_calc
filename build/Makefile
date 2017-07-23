# Makefile for imageLib

# you can compile versions for different architectures, and with and without debug (-g) info
# using "make clean; make" on different machines and with DBG commented in/out

#SRC = Convert.cpp Convolve.cpp Image.cpp ImageIO.cpp ImageIOpng.cpp RefCntMem.cpp
SRC_CPP = main.cpp SSCA.cpp GetMethod.cpp CC/CenCC.cpp CC/CGCC.cpp CC/GrdCC.cpp CC/AECensusCC.cpp CC/ADDCensusCC.cpp PPWM/WMPP.cpp\
CAFilter/BoxCA.cpp CAFilter/BFCA.cpp CAFilter/BilateralFilter.cpp CAFilter/GuidedFilter.cpp CAFilter/GFCA.cpp\
CAST/STCA.cpp CAST/SegmentTree.cpp CAST/StereoDisparity.cpp CAST/StereoHelper.cpp CAST/Toolkit.cpp\
PPSG/SGPP.cpp
DBG = -g
CC = g++
WARN = -Wall -Wconversion
OPT ?= -O3
CPPFLAGS = $(OPT) $(WARN) $(DBG) `pkg-config --cflags opencv`
#CPPFLAGS = $(OPT) $(WARN) `pkg-config --cflags opencv`
LIBS = `pkg-config --libs opencv`
IMLIB = ssca.a
OBJ = $(SRC_CPP:.cpp=.o)
all: $(IMLIB)

$(IMLIB):$(OBJ)
complie:
	$(CC) $(CPPFLAGS) -c $(SRC_CPP) $(LIBS)

link:
	$(CC) $(CPPFLAGS) $(OBJ) -o $(IMLIB) $(LIBS)

clean: 
	rm -f $(OBJ) $(IMLIB)

depend:
	@makedepend -Y -- $(CPPFLAGS) -- $(SRC_CPP) 2>> /dev/null

# DO NOT DELETE THIS LINE -- make depend depends on it.

main.o: CommFunc.h SSCA.h GetMethod.h CC/GrdCC.h CAFilter/GFCA.h 
SSCA.o: CommFunc.h SSCA.h
GetMethod.o: GetMethod.h
CC/CenCC.o: CC/CenCC.h
CC/CGCC.o: CC/CGCC.h CAFilter/GuidedFilter.h
CC/GrdCC.o: CC/GrdCC.h
CC/AECensusCC.o: CC/AECensusCC.h
CC/ADDCensusCC.o: CC/ADDCensusCC.h
PPWM/WMPP.o: PPWM/WMPP.h
CAFilter/BoxCA.o: CAFilter/BoxCA.h CAFilter/GuidedFilter.h
CAFilter/BFCA.o: CAFilter/BFCA.h CAFilter/BilateralFilter.h
CAFilter/BilateralFilter.o: CAFilter/BilateralFilter.h
CAFilter/GuidedFilter.o: CommFunc.h CAFilter/GuidedFilter.h
CAFilter/GFCA.o: CAFilter/GFCA.h CAFilter/GuidedFilter.h
#CAST/ctmf.o: CAST/ctmf.h
CAST/STCA.o: CAST/STCA.h CAST/StereoDisparity.h CAST/StereoHelper.h CAST/SegmentTree.h
CAST/SegmentTree.o: CAST/disjoint-set.h CAST/segment-graph.h CAST/Toolkit.h
CAST/StereoDisparity.o: CAST/StereoDisparity.h CAST/StereoHelper.h CAST/Toolkit.h CAST/SegmentTree.h
CAST/StereoHelper.o: CAST/StereoHelper.h CAST/Toolkit.h
CAST/Toolkit.o: CAST/Toolkit.h CAST/ctmf.h
PPSG/SGPP.o: PPSG/SGPP.h PPSG/image.h PPSG/misc.h PPSG/pnmfile.h PPSG/segment-image.h



