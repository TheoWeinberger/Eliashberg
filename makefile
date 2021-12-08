CXX=icpx
CXXFLAGS=-Ofast -Wall -Wextra -std=c++14 -fopenmp -fPIC 
LIBS= -larmadillo -mkl=parallel -lconfig++ -lgsl -lgslcblas -lalglib 

#2D compilation
DEPS2D = Eliashberg2D/Eliashberg2D.hpp Eliashberg2D/Eliashberg2DSettings.hpp
OBJ2D = Eliashberg2D/Eliashberg2D.o Eliashberg2D/Eliashberg2DSettings.o

#specify both libraries
all: libeliashberg2D.so libeliashberg3D.so

/%.o: %cpp $(DEPS2D)
	$(CXX) -c $@ $< $(CXXFLAGS) 

libeliashberg2D.so: $(OBJ2D)
	$(CXX) -shared -Wl,-soname,libeliashberg2D.so -o $@ $^ $(CXXFLAGS) $(LIBS)

#3D compilation
DEPS3D = Eliashberg3D/Eliashberg3D.hpp Eliashberg3D/Eliashberg3DSettings.hpp
OBJ3D = Eliashberg3D/Eliashberg3D.o Eliashberg3D/Eliashberg3DSettings.o

#3D compilation
/%.o: %cpp $(DEPS3D)
	$(CXX) -c $@ $< $(CXXFLAGS) 

libeliashberg3D.so: $(OBJ3D)
	$(CXX) -shared -Wl,-soname,libeliashberg3D.so -o $@ $^ $(CXXFLAGS) $(LIBS)

#clean
.PHONY: clean

clean:
	-rm $(OBJ2D)
	-rm $(OBJ3D)   