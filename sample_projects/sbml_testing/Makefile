VERSION := $(shell grep . VERSION.txt | cut -f1 -d:)
PROGRAM_NAME:= sbmltest1

CC := g++

# Check for environment definitions of compiler
# e.g., on CC = g++-7 on OSX
ifdef PHYSICELL_CPP
	CC := $(PHYSICELL_CPP)
endif

# override with clang++ due to libRoadrunner macOS version relying on it
#CC := clang++

ARCH := native # best auto-tuning

# For debugging with gdb just add: -g3 -ggdb -O3
# -DCLP_BUILD is an additional flags to link againt clp

CFLAGS := -march=$(ARCH) -g  -fomit-frame-pointer -Xpreprocessor -fopenmp -m64 -std=c++11 -D LIBROADRUNNER 
#CFLAGS := -march=$(ARCH) -fomit-frame-pointer -Xpreprocessor -fopenmp -m64 -std=c++11 -U LIBROADRUNNER 

# OSX
PATH_TO_OMP := -L/usr/local/opt/libomp/lib
OMP_LIB := $(PATH_TO_OMP) -lomp

#LIBRR_DIR := /Users/heiland/libroadrunner/roadrunner-osx-10.9-cp36m
#LIBRR_LIBS := /home/heiland/libroadrunner/cpplibroadrunner-1.3.0-linux_x86_64/lib

OSFLAG 	:=
ifeq ($(OS),Windows_NT)
	OSFLAG += -D WIN32
	ifeq ($(PROCESSOR_ARCHITECTURE),AMD64)
		OSFLAG += -D AMD64
	endif
	ifeq ($(PROCESSOR_ARCHITECTURE),x86)
		OSFLAG += -D IA32
	endif
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Linux)
		OSFLAG += -D LINUX
		OMP_LIB :=
		LIBRR_DIR := $(HOME)/libroadrunner/libroadrunner
		LIBRR_CFLAGS := -I${LIBRR_DIR}/include/rr/C
		LIBRR_LIBS := ${LIBRR_DIR}/lib
		CFLAGS := -march=$(ARCH) -g  -fomit-frame-pointer -fopenmp -m64 -std=c++11 -D LIBROADRUNNER 
	endif
	ifeq ($(UNAME_S),Darwin)
		CC := clang++
		OSFLAG += -D OSX
#		LIBRR_DIR := /Users/heiland/libroadrunner/roadrunner-osx-10.9-cp36m
		LIBRR_CFLAGS := -I${LIBRR_DIR}/include/rr/C
		LIBRR_LIBS := ${LIBRR_DIR}/lib
	endif
		UNAME_P := $(shell uname -p)
	ifeq ($(UNAME_P),x86_64)
		OSFLAG += -D AMD64
		LIBRR_DIR := $(HOME)/libroadrunner/libroadrunner
		LIBRR_CFLAGS := -I${LIBRR_DIR}/include/rr/C
		CFLAGS := -march=$(ARCH) -g  -fomit-frame-pointer -fopenmp -m64 -std=c++11 -D LIBROADRUNNER 
	endif
		ifneq ($(filter %86,$(UNAME_P)),)
	OSFLAG += -D IA32
		endif
	ifneq ($(filter arm%,$(UNAME_P)),)
		OSFLAG += -D ARM
	endif
endif

#all:
#	@echo $(OSFLAG)

COMPILE_COMMAND := $(CC) $(CFLAGS) $(LIBRR_CFLAGS) 

BioFVM_OBJECTS := BioFVM_vector.o BioFVM_mesh.o BioFVM_microenvironment.o BioFVM_solvers.o BioFVM_matlab.o \
BioFVM_utilities.o BioFVM_basic_agent.o BioFVM_MultiCellDS.o BioFVM_agent_container.o

PhysiCell_core_OBJECTS := PhysiCell_phenotype.o PhysiCell_cell_container.o PhysiCell_standard_models.o PhysiCell_cell.o PhysiCell_custom.o PhysiCell_utilities.o PhysiCell_constants.o

PhysiCell_module_OBJECTS := PhysiCell_SVG.o PhysiCell_pathology.o PhysiCell_MultiCellDS.o PhysiCell_various_outputs.o \
PhysiCell_pugixml.o PhysiCell_settings.o

# put your custom objects here (they should be in the custom_modules directory)
PhysiCell_custom_module_OBJECTS := sbml_test1.o

pugixml_OBJECTS := pugixml.o

PhysiCell_OBJECTS := $(BioFVM_OBJECTS)  $(pugixml_OBJECTS) $(PhysiCell_core_OBJECTS) $(PhysiCell_module_OBJECTS)
ALL_OBJECTS := $(PhysiCell_OBJECTS)  $(PhysiCell_custom_module_OBJECTS)


# compile the project

all: main.cpp $(ALL_OBJECTS)
	@echo Your OS= $(OSFLAG)
	@echo LIBRR_CFLAGS= $(LIBRR_CFLAGS)
	@echo LIBRR_LIBS= $(LIBRR_LIBS)
	@echo 
#	$(COMPILE_COMMAND) -L$(LIBRR_LIBS) -lroadrunner_c_api $(PATH_TO_OMP) -lomp -o $(PROGRAM_NAME) $(ALL_OBJECTS) main.cpp
	$(COMPILE_COMMAND) -L$(LIBRR_LIBS) -lroadrunner_c_api $(OMP_LIB) -o $(PROGRAM_NAME) $(ALL_OBJECTS) main.cpp
	@echo
	@echo created $(PROGRAM_NAME)
	@echo


# PhysiCell core components

PhysiCell_phenotype.o: ./core/PhysiCell_phenotype.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_phenotype.cpp

PhysiCell_digital_cell_line.o: ./core/PhysiCell_digital_cell_line.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_digital_cell_line.cpp

PhysiCell_cell.o: ./core/PhysiCell_cell.cpp ./core/PhysiCell_cell.h
	$(COMPILE_COMMAND) -c ./core/PhysiCell_cell.cpp

PhysiCell_cell_container.o: ./core/PhysiCell_cell_container.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_cell_container.cpp

PhysiCell_standard_models.o: ./core/PhysiCell_standard_models.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_standard_models.cpp

PhysiCell_utilities.o: ./core/PhysiCell_utilities.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_utilities.cpp

PhysiCell_custom.o: ./core/PhysiCell_custom.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_custom.cpp

PhysiCell_constants.o: ./core/PhysiCell_constants.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_constants.cpp 

# BioFVM core components (needed by PhysiCell)

BioFVM_vector.o: ./BioFVM/BioFVM_vector.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_vector.cpp

BioFVM_agent_container.o: ./BioFVM/BioFVM_agent_container.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_agent_container.cpp

BioFVM_mesh.o: ./BioFVM/BioFVM_mesh.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_mesh.cpp

BioFVM_microenvironment.o: ./BioFVM/BioFVM_microenvironment.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_microenvironment.cpp

BioFVM_solvers.o: ./BioFVM/BioFVM_solvers.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_solvers.cpp

BioFVM_utilities.o: ./BioFVM/BioFVM_utilities.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_utilities.cpp

BioFVM_basic_agent.o: ./BioFVM/BioFVM_basic_agent.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_basic_agent.cpp

BioFVM_matlab.o: ./BioFVM/BioFVM_matlab.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_matlab.cpp

BioFVM_MultiCellDS.o: ./BioFVM/BioFVM_MultiCellDS.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_MultiCellDS.cpp

pugixml.o: ./BioFVM/pugixml.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/pugixml.cpp

# standard PhysiCell modules

PhysiCell_SVG.o: ./modules/PhysiCell_SVG.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_SVG.cpp

PhysiCell_pathology.o: ./modules/PhysiCell_pathology.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_pathology.cpp

PhysiCell_MultiCellDS.o: ./modules/PhysiCell_MultiCellDS.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_MultiCellDS.cpp

PhysiCell_various_outputs.o: ./modules/PhysiCell_various_outputs.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_various_outputs.cpp

PhysiCell_pugixml.o: ./modules/PhysiCell_pugixml.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_pugixml.cpp

PhysiCell_settings.o: ./modules/PhysiCell_settings.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_settings.cpp

# PhysiCell core components
sbml_test1.o: ./custom_modules/sbml_test1.cpp ./intracellular/PhysiCell_intracellular.h
	$(COMPILE_COMMAND) -c ./custom_modules/sbml_test1.cpp

# cleanup

reset:
	rm -f *.cpp
	cp ./sample_projects/Makefile-default Makefile
	rm -f ./custom_modules/*
	touch ./custom_modules/empty.txt
	touch ALL_CITATIONS.txt
	rm ALL_CITATIONS.txt
	cp ./config/PhysiCell_settings_default.xml ./config/PhysiCell_settings.xml
	touch ./config/empty.txt

# next 2 targets for nanoHUB
install: all
	cp $(PROGRAM_NAME) ../bin

distclean: clean
	rm -f ../bin/$(PROGRAM_NAME)

clean:
	rm -f *.o
	rm -f $(PROGRAM_NAME)

data-cleanup:
	rm -f *.mat
	rm -f *.xml
	rm -f *.svg
	rm -f ./output/*
	touch ./config/empty.txt

# archival

zip:
	zip -r latest.zip Makefile* *.cpp *.h BioFVM/* config/* core/* custom_modules/* matlab/* modules/* sample_projects/*
	cp latest.zip $$(date +%b_%d_%Y_%H%M).zip
	cp latest.zip VERSION_$(VERSION).zip
	mv *.zip archives/

tar:
	tar --ignore-failed-read -czf latest.tar Makefile* *.cpp *.h BioFVM/* config/* core/* custom_modules/* matlab/* modules/* sample_projects/*
	cp latest.tar $$(date +%b_%d_%Y_%H%M).tar
	cp latest.tar VERSION_$(VERSION).tar
	mv *.tar archives/

unzip:
	cp ./archives/latest.zip .
	unzip latest.zip

untar:
	cp ./archives/latest.tar .
	tar -xzf latest.tar
