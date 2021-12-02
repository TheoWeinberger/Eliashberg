CXX=icpx
CXXFLAGS=-Ofast -Wall -Wextra -std=c++14 -fopenmp -fPIC 
LIBS= -larmadillo -mkl=parallel -lconfig++ -lgsl -lgslcblas -lalglib 
DEPS = Eliashberg2D/Eliashberg2D.hpp Eliashberg2D/Eliashberg2DSettings.hpp
OBJ = Eliashberg2D/Eliashberg2D.o Eliashberg2D/Eliashberg2DSettings.o

/%.o: %cpp $(DEPS)
	$(CXX) -c $@ $< $(CXXFLAGS) 

libeliashberg2D.so: $(OBJ)
	$(CXX) -shared -Wl,-soname,libeliashberg2D.so -o $@ $^ $(CXXFLAGS) $(LIBS)