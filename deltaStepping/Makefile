CXX = g++
CXXFLAGS = -std=c++17 -fopenmp -DTBB_INTERFACE_NEW=1
LDFLAGS = -ltbb

TARGET = planner.out
SRC = planner.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(TARGET)