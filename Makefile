CNPY_DIR := /home/lucas/PADVANTG/build/cnpy
CNPY_FLAGS := -I$(CNPY_DIR)/include -L$(CNPY_DIR)/lib -Wl,-rpath,$(CNPY_DIR)/lib -lcnpy

DEBUG_FLAGS := -g -O0
OPT_FLAGS := -O3

EXEC_1 := g++ -o calc_dR_angular calc_dR_angular.cpp $(CNPY_FLAGS) -Wno-narrowing
EXEC_2 := g++ -o calc_dR_scalar  calc_dR_scalar.cpp  $(CNPY_FLAGS) -Wno-narrowing

all:
	$(EXEC_1) $(OPT_FLAGS)
	$(EXEC_2) $(OPT_FLAGS)

debug:
	$(EXEC_1) $(DEBUG_FLAGS)
	$(EXEC_2) $(DEBUG_FLAGS)

format:
	clang-format -i *.cpp
