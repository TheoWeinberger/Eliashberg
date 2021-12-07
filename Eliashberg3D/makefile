CXX=icpx
CXXFLAGS=-Ofast -Wall -Wextra -std=c++14 -fopenmp -g
LIBS= -larmadillo -mkl=parallel -lconfig++ -lgsl -lgslcblas -lalglib
DEPS = Eliashberg3D.hpp Eliashberg3DSettings.hpp
OBJ = Eliashberg3D.o Eliashberg3DSettings.o Eliashberg3DMain.o 

/%.o: %cpp $(DEPS)
	$(CXX) -c $@ $< $(CXXFLAGS) 

eliashberg3D: $(OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)