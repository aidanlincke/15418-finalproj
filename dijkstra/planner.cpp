#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <limits>
#include <unordered_set>
#include <optional>

using Map = std::vector<std::vector<int>>;
struct Position {
    int row;
    int col;

    Position operator+(const Position& other) const {
        return {row + other.row, col + other.col};
    }
    bool operator==(const Position& other) const {
        return row == other.row && col == other.col;
    }
    bool valid(const Map& map) {
        int rows = map.size();
        int cols = map[0].size();
        return 0 <= row && row < rows && 0 <= col && col < cols;
    }
};
struct PositionHash {
    std::size_t operator()(const Position& position) const {
        return std::hash<int>()(position.row) ^ (std::hash<int>()(position.col) << 1);
    }
};
struct DijkstraNode {
    Position position;
    int cost;
    std::vector<Position> history;
    
    bool operator>(const DijkstraNode& other) const {
        return cost > other.cost;
    }
};
struct DijkstraResult {
    int totalCost;
    std::vector<Position> positions;
};

constexpr std::array<Position, 9> movements = {
    Position{0, 0},
    Position{1, 0},
    Position{-1, 0},
    Position{0, 1},
    Position{0, -1},
    Position{1, 1},
    Position{-1, -1},
    Position{1, -1},
    Position{-1, 1}
};


Map readGridFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open the map file.");
    }

    int rows, cols;
    char comma;
    file >> rows >> comma >> cols;
    if (file.fail() || comma != ',') {
        throw std::runtime_error("Invalid format for grid dimensions.");
    }

    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    Map grid(rows, std::vector<int>(cols));
    std::string line;
    for (int row = 0; row < rows; row++) {
        if (!std::getline(file, line)) {
            throw std::runtime_error("Unexpected end of file while reading grid data.");
        }

        std::stringstream ss(line);
        for (int col = 0; col < cols; col++) {
            int value;
            ss >> value;
            if (ss.fail()) {
                throw std::runtime_error("Invalid data at row " + std::to_string(row) + ", col " + std::to_string(col) + ".");
            }
            grid[row][col] = value;

            if (col < cols - 1) ss.ignore();
        }
    }

    return grid;
}

std::optional<DijkstraResult> dijkstras(const Map& map, Position start, Position end) {
    std::unordered_set<Position, PositionHash> closed;
    std::priority_queue<DijkstraNode, std::vector<DijkstraNode>, std::greater<DijkstraNode>> open;

    open.push({start, 0, {start}});
    while (open.size() > 0 && closed.find(end) == closed.end()) {
        DijkstraNode node = open.top();
        open.pop();
        if (closed.find(node.position) != closed.end()) { 
            continue; 
        }
        closed.insert(node.position);

        if (node.position == end) {
            return DijkstraResult{node.cost, node.history};
        }

        for (const Position& move : movements) {
            Position newPos = node.position + move;
            if (!newPos.valid(map) || closed.find(newPos) != closed.end()) { continue; }

            DijkstraNode newNode = DijkstraNode{newPos, node.cost + map[newPos.row][newPos.col], node.history};
            newNode.history.push_back(newPos);
            open.push(newNode);
        }
    }
    return std::nullopt;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <map_file>" << std::endl;
        return -1;
    }

    Map map = readGridFromFile(argv[1]);
    dijkstras(map, {0, 0}, {100, 100});

    return 0;
}