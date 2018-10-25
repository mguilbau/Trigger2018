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

all: mkdirBin mkdirOutput mkdirPdfdir bin/l1Comp.exe bin/l1OfflineSubtract.exe bin/l1FiringFraction.exe bin/l1CaloTowerDist.exe bin/l1OfflineEvtDisp.exe bin/l1ToTTree.exe bin/hltFiringFraction.exe bin/checkHLTAgainstL1Xml.exe bin/checkHLTAgainstPrescaleInput.exe bin/listHLT.exe bin/doPrescaling.exe  bin/doPrescalingL1.exe bin/quickHiBin.exe bin/emptyBunches.exe bin/createJsonFromL1.exe

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

bin/checkHLTAgainstL1Xml.exe: src/checkHLTAgainstL1Xml.C
	$(CXX) $(CXXFLAGS) src/checkHLTAgainstL1Xml.C $(ROOT) -I $(PWD) -o bin/checkHLTAgainstL1Xml.exe 

bin/checkHLTAgainstPrescaleInput.exe: src/checkHLTAgainstPrescaleInput.C
	$(CXX) $(CXXFLAGS) src/checkHLTAgainstPrescaleInput.C $(ROOT) -I $(PWD) -o bin/checkHLTAgainstPrescaleInput.exe 

bin/listHLT.exe: src/listHLT.C
	$(CXX) $(CXXFLAGS) src/listHLT.C $(ROOT) -I $(PWD) -o bin/listHLT.exe 

bin/doPrescaling.exe: src/doPrescaling.C
	$(CXX) $(CXXFLAGS) src/doPrescaling.C $(ROOT) -I $(PWD) -o bin/doPrescaling.exe 

bin/doPrescalingL1.exe: src/doPrescalingL1.C
	$(CXX) $(CXXFLAGS) src/doPrescalingL1.C $(ROOT) -I $(PWD) -o bin/doPrescalingL1.exe 

bin/quickHiBin.exe: src/quickHiBin.C
	$(CXX) $(CXXFLAGS) src/quickHiBin.C $(ROOT) -I $(PWD) -o bin/quickHiBin.exe 

bin/emptyBunches.exe: src/emptyBunches.C
	$(CXX) $(CXXFLAGS) src/emptyBunches.C $(ROOT) -I $(PWD) -o bin/emptyBunches.exe 

bin/createJsonFromL1.exe: src/createJsonFromL1.C
	$(CXX) $(CXXFLAGS) src/createJsonFromL1.C $(ROOT) -I $(PWD) -o bin/createJsonFromL1.exe 


clean:
	rm *~ || true
	rm *# || true
	rm include/*~ || true
	rm include/#*# || true
	rm src/*~ || true
	rm src/#*# || true
	rm bin/*.exe || true
	rmdir bin || true
