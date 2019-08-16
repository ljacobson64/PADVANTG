CNPY_DIR := /home/lucas/PADVANTG/build/cnpy
CNPY_FLAGS := -I$(CNPY_DIR)/include -L$(CNPY_DIR)/lib -Wl,-rpath,$(CNPY_DIR)/lib -lcnpy

DEBUG_FLAGS := -g -O0
OPT_FLAGS := -O3

EXEC := g++ -o calc_dR calc_dR.cpp $(CNPY_FLAGS) -Wno-narrowing

all:
	$(EXEC) $(OPT_FLAGS)

debug:
	$(EXEC) $(DEBUG_FLAGS)

format:
	clang-format -i *.cpp
