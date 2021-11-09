CC=g++
CXXFLAGS=-O4 -Wall -Wextra -fopenmp
LIBS=-larmadillo -llapack -lblas -lfftw3 -lconfig++ -lgsl -lgslcblas -lalglib
DEPS = Eliashberg2D.hpp Eliashberg2DSettings.hpp
OBJ = Eliashberg2D.o Eliashberg2DSettings.o Eliashberg2DMain.o 

/%.o: %cpp $(DEPS)
	$(CC) -c $@ $< $(CXXFLAGS) 

eliashberg2D: $(OBJ)
	$(CC) -o $@ $^ $(CXXFLAGS) $(LIBS)