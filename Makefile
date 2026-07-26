CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++17 -pthread -I.
BIN       = fungi_graphsv_tol_bin

.PHONY: all clean test

all: $(BIN)

$(BIN): main.cpp $(wildcard *.hpp)
	$(CXX) $(CXXFLAGS) main.cpp -o $@

test: $(BIN)
	python3 test_golden_smoke.py

clean:
	rm -f $(BIN)
