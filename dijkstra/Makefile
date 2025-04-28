CXX = g++
CXXFLAGS = -std=c++17 -fopenmp
TARGET = planner.out
SRC = planner.cpp

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

clean:
	rm -f $(TARGET)