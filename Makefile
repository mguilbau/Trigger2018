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

all: mkdirBin mkdirOutput mkdirPdfdir bin/l1Comp.exe bin/l1OfflineSubtract.exe bin/l1FiringFraction.exe bin/l1CaloTowerDist.exe bin/l1OfflineEvtDisp.exe

mkdirBin:
	$(MKDIR_BIN)

mkdirOutput:
	$(MKDIR_OUTPUT)

mkdirPdfdir:
	$(MKDIR_PDFDIR)

bin/l1Comp.exe: src/l1Comp.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -o bin/l1Comp.exe src/l1Comp.C

bin/l1OfflineSubtract.exe: src/l1OfflineSubtract.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -o bin/l1OfflineSubtract.exe src/l1OfflineSubtract.C

bin/l1FiringFraction.exe: src/l1FiringFraction.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -o bin/l1FiringFraction.exe src/l1FiringFraction.C

bin/l1CaloTowerDist.exe: src/l1CaloTowerDist.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -o bin/l1CaloTowerDist.exe src/l1CaloTowerDist.C

bin/l1OfflineEvtDisp.exe: src/l1OfflineEvtDisp.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -o bin/l1OfflineEvtDisp.exe src/l1OfflineEvtDisp.C

clean:
	rm *~ || true
	rm *# || true
	rm include/*~ || true
	rm include/#*# || true
	rm src/*~ || true
	rm src/#*# || true
	rm bin/*.exe || true
	rmdir bin || true