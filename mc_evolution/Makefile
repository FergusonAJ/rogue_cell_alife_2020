EMP_DIR := /home/austin/research/tools/empirical_clean/Empirical/source

CXX := g++-8 

CFLAGS := -Wall -Wno-unused-function -iquote $(EMP_DIR)/ -std=c++17

OFLAGS_optim := -O3 -DNDEBUG
OFLAGS_debug := -g -pedantic -DEMP_TRACK_MEM  -Wnon-virtual-dtor -Wcast-align -Woverloaded-virtual

all: ./main.cc rogue_cell_mc_config.h
	$(CXX) ./main.cc $(CFLAGS) $(OFLAGS_optim) -o multicell_evolution

debug: ./main.cc rogue_cell_mc_config.h
	$(CXX) ./main.cc $(CFLAGS) $(OFLAGS_debug) -o multicell_evolution

clean:
	rm ./multicell_evolution
