# Variables to be updated 
VERSION = 1_0_master
#    bin, debug
EXECUTION = debug
#    -O1 (bin), -g -O0 -fbacktrace -C (gfortran debug),
#    -g -O0 -traceback -C -check bounds -check noarg_temp_created -debug all (ifort debug)
COMPILATION_FLAGS = -g -O0 -fbacktrace -C
#    gfortran, ifort
COMPILER = gfortran
#    -fopenmp, -openmp, (blank for debug)  
OMP_FLAG = 
EXE_DIR = ../$(EXECUTION)/
SOURCES = \
$(wildcard ./*/*.f90)
# Other variables
MAIN_FILE_ROOT = ./Main/vtu2FEHM
CODE = $(EXE_DIR)vtu2FEHM_v_$(VERSION)_$(COMPILER)_$(EXECUTION)
OBJECTS = $(patsubst %.f90,%.o,$(SOURCES))
# Targets
#    Correct sequence for compile+link:
#       make touch
#       make 
all: compile_main 
compile_main:
	touch $(CODE).x
	$(COMPILER) $(MAIN_FILE_ROOT).f90 -o $(CODE).x $(OMP_FLAG) $(COMPILATION_FLAGS)
	rm -f $(MAIN_FILE_ROOT).o
compile: $(SOURCES)
%.f90: %.o
	$(COMPILER) $@ -o $< $(OMP_FLAG) $(COMPILATION_FLAGS) -c 
link: $(OBJECTS)
	touch $(CODE).x
	$(COMPILER) $(MAIN_FILE_ROOT).f90 $^ -o $(CODE).x $(OMP_FLAG) $(COMPILATION_FLAGS)
	rm -f $(MAIN_FILE_ROOT).o
remove: $(OBJECTS)
	rm -f $^
	rm -f *.mod
touch:
	touch $(OBJECTS)
clean: 
	rm -f $(CODE).x 

