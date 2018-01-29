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

all: mkdirBin mkdirOutput mkdirPdfdir l1Comp l1FiringFraction l1CaloTowerDist

mkdirBin:
	$(MKDIR_BIN)

mkdirOutput:
	$(MKDIR_OUTPUT)

mkdirPdfdir:
	$(MKDIR_PDFDIR)

l1Comp: src/l1Comp.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -o bin/l1Comp.exe src/l1Comp.C

l1FiringFraction: src/l1FiringFraction.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -o bin/l1FiringFraction.exe src/l1FiringFraction.C

l1CaloTowerDist: src/l1CaloTowerDist.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -o bin/l1CaloTowerDist.exe src/l1CaloTowerDist.C

clean:
	rm *~ || true
	rm *# || true
	rm include/*~ || true
	rm include/#*# || true
	rm src/*~ || true
	rm src/#*# || true
	rm bin/*.exe || true
	rmdir bin || true