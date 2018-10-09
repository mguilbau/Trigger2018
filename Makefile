#std-c++11
CXX = g++
CXXFLAGS = -Wall -O2 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

ROOT = `root-config --cflags --glibs`


MKDIR_BIN = mkdir -p $(PWD)/bin
MKDIR_PDFDIR = mkdir -p $(PWD)/pdfDir
MKDIR_OUTPUT = mkdir -p $(PWD)/output

all: mkdirBin mkdirOutput mkdirPdfdir bin/l1Comp.exe bin/l1OfflineSubtract.exe bin/l1FiringFraction.exe bin/l1CaloTowerDist.exe bin/l1OfflineEvtDisp.exe bin/l1ToTTree.exe bin/hltFiringFraction.exe

mkdirBin:
	$(MKDIR_BIN)

mkdirOutput:
	$(MKDIR_OUTPUT)

mkdirPdfdir:
	$(MKDIR_PDFDIR)

bin/l1Comp.exe: src/l1Comp.C
	$(CXX) $(CXXFLAGS) src/l1Comp.C $(ROOT) -I $(PWD) -o bin/l1Comp.exe

bin/l1OfflineSubtract.exe: src/l1OfflineSubtract.C
	$(CXX) $(CXXFLAGS) src/l1OfflineSubtract.C $(ROOT) -I $(PWD) -o bin/l1OfflineSubtract.exe

bin/l1FiringFraction.exe: src/l1FiringFraction.C
	$(CXX) $(CXXFLAGS) src/l1FiringFraction.C $(ROOT) -I $(PWD) -o bin/l1FiringFraction.exe

bin/l1CaloTowerDist.exe: src/l1CaloTowerDist.C
	$(CXX) $(CXXFLAGS) src/l1CaloTowerDist.C $(ROOT) -I $(PWD) -o bin/l1CaloTowerDist.exe

bin/l1OfflineEvtDisp.exe: src/l1OfflineEvtDisp.C
	$(CXX) $(CXXFLAGS) src/l1OfflineEvtDisp.C $(ROOT) -I $(PWD) -o bin/l1OfflineEvtDisp.exe

bin/l1ToTTree.exe: src/l1ToTTree.C
	$(CXX) $(CXXFLAGS) src/l1ToTTree.C $(ROOT) -I $(PWD) -o bin/l1ToTTree.exe

bin/hltFiringFraction.exe: src/hltFiringFraction.C
	$(CXX) $(CXXFLAGS) src/hltFiringFraction.C $(ROOT) -I $(PWD) -o bin/hltFiringFraction.exe 


clean:
	rm *~ || true
	rm *# || true
	rm include/*~ || true
	rm include/#*# || true
	rm src/*~ || true
	rm src/#*# || true
	rm bin/*.exe || true
	rmdir bin || true
