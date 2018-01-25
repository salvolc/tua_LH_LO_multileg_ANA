#CopyPaste Sachen damit hoffentlich root einigermaßen zufrieden läuft
#--------------------------------------------------------------------
# Root variables (SYSTEM)
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := -lMinuit $(shell root-config --libs)
ROOTTEST	 := -I$(shell root-config --prefix)/test
EXROOT = -L/home/salv/Dokumente/Masterarbeit/MG1/ExRootAnalysis/ -lExRootAnalysis

DELPHES = -L/home/salv/Dokumente/Masterarbeit/Delphes/ -lDelphes

# Root variables (.local) 
# ROOTCFLAGS :=
# ROOTLIBS := 

CXXFLAGS  =  -g -Wall -fPIC -Wno-deprecated -ggdb
CXXFLAGS    += $(ROOTCFLAGS) $(ROOTTEST) 
CXXFLAGS    += -I. -I./include -I/home/salv/Dokumente/Masterarbeit/MG1/ExRootAnalysis/ -I/home/salv/Dokumente/Masterarbeit/MG1/ExRootAnalysis/ExRootAnalysis -I/home/salv/Dokumente/Masterarbeit/Delphes/
LIBS 		= $(ROOTLIBS) -L./ -lTreePlayer -lMinuit $(EXROOT) $(DELPHES)
#--------------------------------------------------------------------
plots = plots/test.pdf
pres = Presentation.pdf
SRC = $(wildcard *.cxx)
#$(info $$plots is [${plots}])
OOBJ = $(SRC:%.cxx=%.o)
#$(info $$OOBJ is [${OOBJ}])

EX = $(wildcard *.cpp)
PRO = $(EX:%.cpp=%)
#$(info $$PRO is [${PRO}])

all: $(PRO) $(OOBJ) $(plots) 
# $(pres)

%.o: %.cxx %.h
	g++ -c $< -o $@ $(CXXFLAGS)

$(PRO): *.cpp $(OOBJ)
	g++ $^ -o $@ -fopenmp $(CXXFLAGS) $(LIBS)

$(plots): $(PRO) aufgabe.py
	./$(PRO)
	python3.5 aufgabe.py >> pyout

py: aufgabe.py
	python3.5 aufgabe.py >> pyout

light: $(PRO)
	./$(PRO)

vertex: 
	g++ vertex.cpp -o vertex -fopenmp $(CXXFLAGS) $(LIBS)
	./vertex

pres:
	cd Vortrag/ && make
	cp Vortrag/build/Presentation.pdf ./

pp:
	python3.5 aufgabe.py
	cd Vortrag/ && make
	cp Vortrag/build/Presentation.pdf ./

#aufgabe: aufgabe.cpp

#$(data): aufgabe
#	./aufgabe

#$(plots): $(data) aufgabe.py
#	python3.5 aufgabe.py


clean:
	@rm ./data/*
	@rm ./plots/*/*

