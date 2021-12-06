CXX=icpx
CXXFLAGS=-Ofast -Wall -Wextra -std=c++14 -fopenmp
LIBS= -larmadillo -mkl=parallel -lconfig++ -lgsl -lgslcblas -lalglib
DEPS = Eliashberg2D.hpp Eliashberg2DSettings.hpp
OBJ = Eliashberg2D.o Eliashberg2DSettings.o Eliashberg2DMain.o 

/%.o: %cpp $(DEPS)
	$(CXX) -c $@ $< $(CXXFLAGS) 

eliashberg2D: $(OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)