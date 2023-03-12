# Makefile for "memory.cpp" and additional function files
#
# 2016-10-19 by Jonas Nußdorfer (jonas.nussdorfer@gmail.com)
#
# =======================================================================================

CXX = g++
CXXFLAGS = -O2 -Wall -Wextra -Weffc++ -std=c++0x
OBJ = MCTdynamics.o Utils.o Mesh.o Correlator.o

all:  MCTdynamics

MCTdynamics: $(OBJ)
	$(CXX) $(CXXFLAGS) -o MCTdynamics $(OBJ)

Utils.o: Utils.cpp Utils.h
	$(CXX) $(CXXFLAGS) -c Utils.cpp
Mesh.o: Mesh.cpp Mesh.h
	$(CXX) $(CXXFLAGS) -c Mesh.cpp 
Correlator.o: Correlator.cpp Correlator.h 
	$(CXX) $(CXXFLAGS) -c Correlator.cpp
MCTdynamics.o: MCTdynamics.cpp
	$(CXX) $(CXXFLAGS) -c MCTdynamics.cpp 

.PHONY: clean
clean:
	rm *.o
