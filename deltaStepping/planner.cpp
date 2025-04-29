#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <limits>
#include <unordered_set>
#include <optional>
#include <queue>
#include <unordered_map>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <omp.h>
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_unordered_set.h>
#include <tbb/concurrent_unordered_map.h>

using Map = std::vector<std::vector<float>>;
struct Position {
    int row;
    int col;

    Position operator+(const Position& other) const {
        return {row + other.row, col + other.col};
    }
    bool operator==(const Position& other) const {
        return row == other.row && col == other.col;
    }
    bool operator !=(const Position& other) const {
        return !(*this == other);
    }
    bool valid(const Map& map) const {
        int rows = map.size();
        int cols = map[0].size();
        return 0 <= row && row < rows && 0 <= col && col < cols;
    }
};
std::ostream& operator<<(std::ostream& os, const Position& pos) {
    return os << "(" << pos.row << ", " << pos.col << ")";
}
struct PositionHash {
    std::size_t operator()(const Position& position) const {
        return std::hash<int>()(position.row) ^ (std::hash<int>()(position.col) << 1);
    }
};
struct DijkstraNode {
    Position position;
    float cost;
    
    bool operator>(const DijkstraNode& other) const {
        return cost > other.cost;
    }
};
struct PlanningProblem {
    Map map;
    Position start;
    Position end;
};
struct PlanningResult {
    float totalCost;
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

std::optional<PlanningResult> sequentialDijkstras(const PlanningProblem& problem) {
    Map map = problem.map;
    Position start = problem.start;
    Position end = problem.end;

    std::priority_queue<DijkstraNode, std::vector<DijkstraNode>, std::greater<DijkstraNode>> open;
    std::unordered_set<Position, PositionHash> closed;
    std::unordered_map<Position, Position, PositionHash> prev;
    std::unordered_map<Position, float, PositionHash> dist;

    open.push({start, 0});
    while (open.size() > 0) {
        DijkstraNode node = open.top();
        open.pop();

        if (closed.find(node.position) != closed.end()) {
            continue; 
        }
        closed.insert(node.position);

        if (node.position == end) {
            std::vector<Position> path;
            for (Position at = end; at != start; at = prev[at]) {
                path.push_back(at);
            }
            path.push_back(start);
            std::reverse(path.begin(), path.end());
            return PlanningResult{node.cost, path};
        }

        for (const Position& move : movements) {
            Position newPos = node.position + move;
            if (!newPos.valid(map) || closed.find(newPos) != closed.end()) { 
                continue; 
            }

            float newCost = node.cost + map[newPos.row][newPos.col] + 1;
            if (dist.find(newPos) == dist.end() || newCost < dist[newPos]) {
                dist[newPos] = newCost;
                prev[newPos] = node.position;
                open.push({newPos, newCost});
            }
        }
    }
    return std::nullopt;
}

struct Update {
    float newCost;
    Position position;
    int i;
};

std::optional<PlanningResult> deltaStepping(const PlanningProblem& problem, float delta) {
    Map map = problem.map;
    Position start = problem.start;
    Position end = problem.end;

    const int numBuckets = 10000;
    auto bucketIndex = [&](float dist) -> int {
        return static_cast<int>(std::floor(dist / delta));
    };

    std::vector<tbb::concurrent_vector<Position>> buckets(numBuckets);
    tbb::concurrent_unordered_map<Position, Position, PositionHash> prev;
    tbb::concurrent_unordered_map<Position, float, PositionHash> dist;

    buckets[0].push_back(start);
    dist[start] = 0;
    for (int bucketI = 0; bucketI < numBuckets; bucketI++) {

        tbb::concurrent_vector<Position> S;
        tbb::concurrent_vector<Position> R = buckets[bucketI];
        tbb::concurrent_vector<Position> nextR;

        while (R.size() > 0) {

            #pragma omp parallel for schedule(dynamic, 8)
            for (int i = 0; i < R.size(); i++) {
                const Position& position = R[i];
                S.push_back(position);

                float currCost = dist[position];
                for (const Position& movement : movements) {
                    const Position newPosition = position + movement;
                    if (!newPosition.valid(map)) continue;
    
                    float weight = map[newPosition.row][newPosition.col] + 1;
                    if (weight > delta) continue;
    
                    float newCost = currCost + weight;

                    if ((dist.count(newPosition) == 0 || newCost < dist[newPosition])) {
                        dist[newPosition] = newCost;
                        prev[newPosition] = position;
    
                        int nextBucketI = bucketIndex(newCost);
                        if (nextBucketI == bucketI) {
                            nextR.push_back(newPosition);
                        } else {
                            buckets[nextBucketI].push_back(newPosition);
                        }
                    }
                }
            }
            R = nextR;
            nextR.clear();
        }

        #pragma omp parallel for schedule(dynamic, 8)
        for (int i = 0; i < S.size(); i++) {
            const Position& position = S[i];
            float currCost = dist[position];

            for (const Position& movement : movements) {
                Position newPosition = position + movement;
                if (!newPosition.valid(map)) continue;

                float weight = map[newPosition.row][newPosition.col] + 1;
                if (weight <= delta) continue;

                float newCost = currCost + weight;
                if (dist.count(newPosition) == 0 || newCost < dist[newPosition]) {
                    dist[newPosition] = newCost;
                    prev[newPosition] = position;
                    buckets[bucketIndex(newCost)].push_back(newPosition);
                }
            }
        }
    }

    std::vector<Position> path;
    for (Position at = end; at != start; at = prev[at]) {
        path.push_back(at);
    }
    path.push_back(start);
    std::reverse(path.begin(), path.end());
    return PlanningResult{dist[end], path};
}

PlanningProblem loadProblem(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("Failed to open " + path);

    int32_t rows, cols;
    int32_t start_row, start_col, goal_row, goal_col;

    in.read(reinterpret_cast<char*>(&rows), sizeof(int32_t));
    in.read(reinterpret_cast<char*>(&cols), sizeof(int32_t));
    in.read(reinterpret_cast<char*>(&start_row), sizeof(int32_t));
    in.read(reinterpret_cast<char*>(&start_col), sizeof(int32_t));
    in.read(reinterpret_cast<char*>(&goal_row), sizeof(int32_t));
    in.read(reinterpret_cast<char*>(&goal_col), sizeof(int32_t));

    std::vector<float> flat(rows * cols);
    in.read(reinterpret_cast<char*>(flat.data()), flat.size() * sizeof(float));
    in.close();

    Map grid(rows, std::vector<float>(cols));
    for (int i = 0; i < rows; ++i) {
        std::copy(flat.begin() + i * cols, flat.begin() + (i + 1) * cols, grid[i].begin());
    }

    return PlanningProblem{
        .map = std::move(grid),
        .start = {start_row, start_col},
        .end = {goal_row, goal_col}
    };
}

void writePlan(const std::string& filename, const std::vector<Position>& plan) {
    std::ofstream out(filename, std::ios::binary);
    if (!out) throw std::runtime_error("Failed to open output file: " + filename);

    int32_t length = plan.size();
    out.write(reinterpret_cast<char*>(&length), sizeof(int32_t));

    for (const auto& pos : plan) {
        out.write(reinterpret_cast<const char*>(&pos.row), sizeof(int32_t));
        out.write(reinterpret_cast<const char*>(&pos.col), sizeof(int32_t));
    }

    out.close();
}

int main(int argc, char *argv[]) {
    omp_set_num_threads(8);
    PlanningProblem problem = loadProblem("problems/dc.bin");

    auto start = std::chrono::high_resolution_clock::now();
    std::optional<PlanningResult> result = deltaStepping(problem, 20.0f);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;

    if (result.has_value()) {
        writePlan("plans/dc.bin", result.value().positions);
        std::cout << "Found a path with " << result.value().totalCost << " total cost in " << diff.count() << " seconds." << std::endl;
    } else {
        std::cout << "No path found." << std::endl;
    }
    return 0;
}