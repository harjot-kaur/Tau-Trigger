#makefile 


CC   =   g++

#UCFLAGS = -O0 -g3 -Wall -gstabs+  
UCFLAGS = -O3 -Wall -gstabs+ -std=c++0x 


RUCFLAGS := $(shell root-config --cflags) -I${ROOFITSYS}/include/ -I./include/
LIBS :=  $(shell root-config --libs) -lTreePlayer -L${ROOFITSYS}/lib/ -lRooFit  -lRooFitCore ./obj/FuncCB_cdf_cpp.so 
GLIBS := $(shell root-config --glibs)

VPATH = ./src/

SRCPP = main.cpp\
	Utilities.cpp\
	TurnonFit.cpp\
	TurnonManager.cpp



         
#OBJCPP = $(SRCPP:.cpp=.o)
OBJCPP = $(patsubst %.cpp,obj/%.o,$(SRCPP))


all : obj/FuncCB_cdf_cpp.so fit.exe 

obj/%.o : %.cpp
	@mkdir -p obj/
	@echo compiling $*
	@$(CC) -c $< $(UCFLAGS) $(RUCFLAGS) -o $@

fit.exe : $(OBJCPP) 
	@echo linking
	@$(CC) $^ $(ACLIBS) $(LIBS) $(GLIBS) -o $@

obj/FuncCB_cdf_cpp.so: include/FuncCB_cdf.cpp createCBLib.sh
	sh createCBLib.sh

clean:
	@rm  -f obj/*.o
	@rm -f fit.exe

cleanall: clean
	@rm obj/FuncCB_cdf_cpp.so

