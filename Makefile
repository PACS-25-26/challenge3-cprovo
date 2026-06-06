CXX = mpicxx
CXXFLAGS = -std=c++17 -Wall -O3 -fopenmp -Iinclude -I/usr/include/eigen3 $(mkEigenInc)
LDFLAGS = -fopenmp

SRCS = src/main.cpp src/Jacobi_solver.cpp src/BlockJacobi_solver.cpp
OBJS = $(patsubst src/%.cpp,obj/%.o,$(SRCS))
TARGET = jacobi_solver

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

obj/%.o: src/%.cpp | obj
	$(CXX) $(CXXFLAGS) -c -o $@ $<

obj:
	mkdir -p obj

clean:
	rm -rf obj $(TARGET)
