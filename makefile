CC=g++
CXXFLAGS=-O2 -Wall -Wextra -fopenmp
LIBS=-larmadillo -llapack -lblas -lfftw3 -lconfig++ -lgsl -lgslcblas
DEPS = Eliashberg2D.hpp 
OBJ = Eliashberg2D.o

/%.o: %cpp $(DEPS)
	$(CC) -c $@ $< $(CXXFLAGS) 

eliashberg2D: $(OBJ)
	$(CC) -o $@ $^ $(CXXFLAGS) $(LIBS)