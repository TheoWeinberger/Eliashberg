CXX=icpx
CXXFLAGS=-Ofast -Wall -Wextra -fopenmp
LIBS=-larmadillo -lmkl_rt -llapack -lblas -lfftw3 -lconfig++ -lgsl -lgslcblas -lalglib
DEPS = Eliashberg2D.hpp Eliashberg2DSettings.hpp
OBJ = Eliashberg2D.o Eliashberg2DSettings.o Eliashberg2DMain.o 

/%.o: %cpp $(DEPS)
	$(CXX) -c $@ $< $(CXXFLAGS) 

eliashberg2D: $(OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)