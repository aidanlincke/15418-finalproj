/***
 *
 *
 *
 *
 * SOREN FILE
 *
 *
 *
 *
 *
 *
 */

#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
// #include <immintrin.h>
// #include <stm.h>
#include <tbb/concurrent_unordered_map.h>

struct Node {
    int value;
    Node* next;

    Node(int v, Node* n = nullptr)
        : value(v)
        , next(n)
    {
    }
};

struct LinkedList {
    Node* start;
    Node* end;

    LinkedList()
        : start(nullptr)
        , end(nullptr)
    {
    }

    void push_back(int value)
    {
        Node* new_node = new Node(value);
        if (start == nullptr) {
            start = new_node;
            end = new_node;
        } else {
            end->next = new_node;
            end = new_node;
        }
    }

    void append(LinkedList* other)
    {
        // for (Node* current = other->start; current != nullptr; current = current->next) {
        //     push_back(current->value);
        // }

        if (other->start == nullptr) {
            return;
        }

        if (start == nullptr) {
            start = other->start;
            end = other->end;
        } else {
            end->next = other->start;
            end = other->end;
        }
    }

    void print() const
    {
        Node* current = start;
        while (current != nullptr) {
            std::cout << current->value << " ";
            current = current->next;
        }
    }
};

// A custom structure to represent a weighted edge
struct Edge {
    float weight;

    // Constructor for convenience
    Edge(float w = 0.0f)
        : weight(w)
    {
    }
};

// A custom structure to represent a weighted edge with a path (linked list)
struct ShortcutEdge {
    float weight;
    LinkedList* path;

    // Constructor for convenience
    ShortcutEdge(float w = 0.0f)
        : weight(w)
        , path(new LinkedList())
    {
    }
    ShortcutEdge(float w, LinkedList* p)
        : weight(w)
        , path(p)
    {
    }
};

// Type 1: Simple Weighted Graph
// Each vertex is represented by an index, and its neighbors are stored in an unordered_map
// The unordered_map maps from neighbor vertex (int) to edge weight (float)
using SimpleGraph = std::vector<std::unordered_map<int, Edge>>;

// Type 2: Contracted Graph
// Contains both the graph data and a set of active nodes
struct ContractedGraph {
    std::vector<tbb::concurrent_unordered_map<int, ShortcutEdge>> edges;
    std::unordered_set<int> activeNodes;

    // Constructor to create from a SimpleGraph
    explicit ContractedGraph(const SimpleGraph& simple)
    {
        // Convert simple edges to shortcut edges
        edges.resize(simple.size());
        for (int v = 0; v < simple.size(); ++v) {
            for (const auto& [neighbor, edgeData] : simple[v]) {
                LinkedList* path = new LinkedList();
                edges[v][neighbor] = ShortcutEdge(edgeData.weight, path);
            }
            // Initially, all nodes are active
            activeNodes.insert(v);
        }
    }

    // Constructor with specified size
    explicit ContractedGraph(int size)
    {
        edges.resize(size);
        for (int i = 0; i < size; ++i) {
            activeNodes.insert(i);
        }
    }
};

// Function to load a SimpleGraph from a Pittsburgh graph file
SimpleGraph loadGraphFromFile(const std::string& filename)
{
    // Open the file
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return SimpleGraph();
    }

    // Read the number of nodes and edges
    int numNodes, numEdges;
    file >> numNodes >> numEdges;

    // Create graph with appropriate size
    SimpleGraph graph(numNodes);

    // Read each edge
    int numSelfLoops = 0;
    for (int i = 0; i < numEdges; ++i) {
        int from, to;
        float weight;

        file >> from >> to >> weight;

        if (from == to) {
            numSelfLoops++;
            continue;
        }

        // Add the edge to the graph (both directions for undirected graph)
        graph[from][to] = Edge(weight);
        graph[to][from] = Edge(weight);
    }
    numEdges -= numSelfLoops;

    std::cout << "Loaded graph with " << numNodes << " nodes and " << numEdges << " edges" << std::endl;
    return graph;
}

// Helper function to create a simple graph with a specified number of vertices
SimpleGraph createSimpleGraph(int numVertices)
{
    return SimpleGraph(numVertices);
}

// Helper function to add an edge to the simple graph
void addEdge(SimpleGraph& graph, int from, int to, float weight)
{
    graph[from][to] = Edge(weight);
    graph[to][from] = Edge(weight);
}

// Helper function to create a grid graph with rows*cols vertices
SimpleGraph createGridGraph(int rows, int cols)
{
    int numVertices = rows * cols;
    SimpleGraph graph = createSimpleGraph(numVertices);

    std::mt19937 rng(std::random_device {}());
    std::uniform_real_distribution<float> dist(1.0f, 10.0f);

    // Helper function to convert (row, col) to vertex index
    auto vertexIndex = [cols](int row, int col) {
        return row * cols + col;
    };

    // Add edges for the grid
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            int current = vertexIndex(row, col);

            // Add edge to right neighbor
            if (col + 1 < cols) {
                int right = vertexIndex(row, col + 1);
                float weight = dist(rng);
                addEdge(graph, current, right, weight);
            }

            // Add edge to bottom neighbor
            if (row + 1 < rows) {
                int bottom = vertexIndex(row + 1, col);
                float weight = dist(rng);
                addEdge(graph, current, bottom, weight);
            }

            // Add some diagonal edges for more interesting shortcuts
            if (row + 1 < rows && col + 1 < cols && dist(rng) > 7.0f) {
                int diagonal = vertexIndex(row + 1, col + 1);
                float weight = dist(rng);
                addEdge(graph, current, diagonal, weight);
            }
        }
    }

    return graph;
}

// Helper function to contract a node in the contracted graph
void contractNode(ContractedGraph& graph, int nodeToContract)
{
    std::vector<std::pair<int, ShortcutEdge>> edgeList(
        graph.edges[nodeToContract].begin(),
        graph.edges[nodeToContract].end());

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < edgeList.size(); ++i) {
        const auto& [start, edge_in] = edgeList[i];
        // for (const auto& [start, edge_in] : graph.edges[nodeToContract]) {

        for (const auto& [end, edge_out] : graph.edges[nodeToContract]) {

            if (start < end) {
                float new_weight = edge_in.weight + edge_out.weight;

                if (graph.edges[start].find(end) == graph.edges[start].end() || graph.edges[start][end].weight > new_weight) {
                    LinkedList* new_path = new LinkedList();

                    new_path->append(graph.edges[start][nodeToContract].path);
                    new_path->push_back(nodeToContract);
                    new_path->append(edge_out.path);

                    graph.edges[start][end] = ShortcutEdge(new_weight, new_path);

                    LinkedList* new_path_reverse = new LinkedList();
                    new_path_reverse->append(graph.edges[end][nodeToContract].path);
                    new_path_reverse->push_back(nodeToContract);
                    new_path_reverse->append(edge_in.path);

                    graph.edges[end][start] = ShortcutEdge(new_weight, new_path_reverse);
                }
            }
        }
    }

    for (const auto& [neighbor, edge] : graph.edges[nodeToContract]) {
        graph.edges[neighbor].unsafe_erase(nodeToContract);
    }

    graph.activeNodes.erase(nodeToContract);
}

ContractedGraph makeContraction(SimpleGraph& graph)
{
    ContractedGraph contractedGraph(graph);

    for (int i = 0; i < graph.size(); ++i) {
        if (i % 11 != 0) { //(contractedGraph.edges[i].size() < 4) {
            std::cout << "Contracting node " << i << std::endl;
            contractNode(contractedGraph, i);
        }
    }

    return contractedGraph;
}

// Helper function to print a simple graph
void printSimpleGraph(const SimpleGraph& graph)
{
    for (int v = 0; v < graph.size(); ++v) {
        std::cout << "Vertex " << v << ":\n";
        for (const auto& [neighbor, edge] : graph[v]) {
            std::cout << "  -> " << neighbor << " (weight: " << edge.weight << ")\n";
        }
    }
}

// Helper function to print a contracted graph
void printContractedGraph(const ContractedGraph& graph)
{
    std::cout << "Active nodes: ";
    for (int node : graph.activeNodes) {
        std::cout << node << " ";
    }
    std::cout << "\n";

    for (int node : graph.activeNodes) {
        std::cout << "Vertex " << node << ":\n";
        for (const auto& [neighbor, edge] : graph.edges[node]) {
            std::cout << "  -> " << neighbor << " (weight: " << edge.weight << ", path: ";
            edge.path->print();
            std::cout << ")\n";
        }
    }
    // for (int v = 0; v < graph.edges.size(); ++v) {
    //     std::cout << "Vertex " << v;
    //     if (graph.activeNodes.count(v) == 0) {
    //         std::cout << " (contracted)";
    //     }
    //     std::cout << ":\n";

    //     for (const auto& [neighbor, edge] : graph.edges[v]) {
    //         std::cout << "  -> " << neighbor << " (weight: " << edge.weight << ", path: ";
    //         edge.path->print();
    //         std::cout << ")\n";
    //     }
    // }
}

bool verifyGraph(SimpleGraph& graph)
{
    int numNodes = graph.size();
    for (int i = 0; i < numNodes; ++i) {
        for (const auto& [neighbor, edge] : graph[i]) {
            if (neighbor == i || neighbor < 0 || neighbor >= numNodes) {
                std::cout << "Invalid neighbor: " << neighbor << std::endl;
                return false;
            }
            if (graph[neighbor].find(i) == graph[neighbor].end() || graph[neighbor][i].weight != edge.weight) {
                std::cout << "Invalid edge weight: " << graph[neighbor][i].weight << " != " << edge.weight << std::endl;
                return false;
            }
        }
    }
    return true;
}

SimpleGraph simpleGraph1()
{
    SimpleGraph simpleGraph = createSimpleGraph(5);

    // Add some edges
    addEdge(simpleGraph, 0, 1, 4.5);
    addEdge(simpleGraph, 0, 2, 3.1);
    addEdge(simpleGraph, 1, 2, 2.0);
    addEdge(simpleGraph, 1, 3, 1.7);
    addEdge(simpleGraph, 2, 3, 2.5);
    addEdge(simpleGraph, 3, 4, 5.0);

    return simpleGraph;
}

/***
 *
 *
 *
 *
 * Aidan FILE
 *
 *
 *
 *
 *
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <omp.h>
#include <optional>
#include <queue>
#include <random>
#include <sstream>
#include <stdexcept>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_unordered_set.h>
#include <tbb/concurrent_vector.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using Map = std::vector<std::vector<float>>;
struct Position {
    int row;
    int col;

    Position operator+(const Position& other) const
    {
        return { row + other.row, col + other.col };
    }
    bool operator==(const Position& other) const
    {
        return row == other.row && col == other.col;
    }
    bool operator!=(const Position& other) const
    {
        return !(*this == other);
    }
    bool valid(const Map& map) const
    {
        int rows = map.size();
        int cols = map[0].size();
        return 0 <= row && row < rows && 0 <= col && col < cols;
    }
};
std::ostream& operator<<(std::ostream& os, const Position& pos)
{
    return os << "(" << pos.row << ", " << pos.col << ")";
}
struct PositionHash {
    std::size_t operator()(const Position& position) const
    {
        return std::hash<int>()(position.row) ^ (std::hash<int>()(position.col) << 1);
    }
};
struct DijkstraNode {
    Position position;
    float cost;

    bool operator>(const DijkstraNode& other) const
    {
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
    Position { 0, 0 },
    Position { 1, 0 },
    Position { -1, 0 },
    Position { 0, 1 },
    Position { 0, -1 },
    Position { 1, 1 },
    Position { -1, -1 },
    Position { 1, -1 },
    Position { -1, 1 }
};

std::optional<PlanningResult> sequentialDijkstras(const PlanningProblem& problem)
{
    Map map = problem.map;
    Position start = problem.start;
    Position end = problem.end;

    std::priority_queue<DijkstraNode, std::vector<DijkstraNode>, std::greater<DijkstraNode>> open;
    std::unordered_set<Position, PositionHash> closed;
    std::unordered_map<Position, Position, PositionHash> prev;
    std::unordered_map<Position, float, PositionHash> dist;

    open.push({ start, 0 });
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
            return PlanningResult { node.cost, path };
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
                open.push({ newPos, newCost });
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

std::pair<int, int> pickTwoRandom(const std::unordered_set<int>& activeNodes)
{
    if (activeNodes.size() < 2) {
        throw std::runtime_error("Not enough elements to pick two!");
    }

    // Move elements into a vector
    std::vector<int> nodes(activeNodes.begin(), activeNodes.end());

    // Random device and generator
    static std::random_device rd;
    static std::mt19937 gen(rd());

    std::uniform_int_distribution<> dis(0, nodes.size() - 1);

    int firstIndex = dis(gen);
    int secondIndex;
    do {
        secondIndex = dis(gen);
    } while (secondIndex == firstIndex);

    return { nodes[firstIndex], nodes[secondIndex] };
}

std::optional<PlanningResult> deltaSteppingStuff(const ContractedGraph& contractedGraph, float delta)
{
    std::pair<int, int> pair = pickTwoRandom(contractedGraph.activeNodes);
    int start = pair.first;
    int goal = pair.second;
    const int numBuckets = 10000;
    auto bucketIndex = [&](float dist) -> int {
        return static_cast<int>(std::floor(dist / delta));
    };

    std::vector<tbb::concurrent_vector<int>> buckets(numBuckets);
    tbb::concurrent_unordered_map<int, int> prev;
    tbb::concurrent_unordered_map<int, float> dist;

    buckets[0].push_back(start);
    dist[start] = 0;
    for (int bucketI = 0; bucketI < numBuckets; bucketI++) {

        tbb::concurrent_vector<int> S;
        tbb::concurrent_vector<int> R = buckets[bucketI];
        tbb::concurrent_vector<int> nextR;

        while (R.size() > 0) {

            #pragma omp parallel for schedule(dynamic, 8)
            for (int i = 0; i < R.size(); i++) {
                const int& position = R[i];
                S.push_back(position);

                float currCost = dist[position];
                tbb::concurrent_unordered_map<int, ShortcutEdge> edgesFromPosition = contractedGraph.edges[position];
                for (const auto& [newPosition, shortcutEdge] : edgesFromPosition) {
                    float weight = shortcutEdge.weight + 1;
                    if (weight > delta)
                        continue;

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
            const int& position = S[i];
            float currCost = dist[position];

            tbb::concurrent_unordered_map<int, ShortcutEdge> edgesFromPosition = contractedGraph.edges[position];
            for (const auto& [newPosition, shortcutEdge] : edgesFromPosition) {

                float weight = shortcutEdge.weight + 1;
                if (weight <= delta)
                    continue;

                float newCost = currCost + weight;
                if (dist.count(newPosition) == 0 || newCost < dist[newPosition]) {
                    dist[newPosition] = newCost;
                    prev[newPosition] = position;
                    buckets[bucketIndex(newCost)].push_back(newPosition);
                }
            }
        }
    }

    std::vector<int> path;
    for (int at = goal; at != start; at = prev[at]) {
        path.push_back(at);
    }
    path.push_back(start);
    std::reverse(path.begin(), path.end());
    return PlanningResult { dist[goal], path };
}

std::optional<PlanningResult> deltaStepping(const PlanningProblem& problem, float delta)
{
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
                    if (!newPosition.valid(map))
                        continue;

                    float weight = map[newPosition.row][newPosition.col] + 1;
                    if (weight > delta)
                        continue;

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
                if (!newPosition.valid(map))
                    continue;

                float weight = map[newPosition.row][newPosition.col] + 1;
                if (weight <= delta)
                    continue;

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
    return PlanningResult { dist[end], path };
}

PlanningProblem loadProblem(const std::string& path)
{
    std::ifstream in(path, std::ios::binary);
    if (!in)
        throw std::runtime_error("Failed to open " + path);

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

    return PlanningProblem {
        .map = std::move(grid),
        .start = { start_row, start_col },
        .end = { goal_row, goal_col }
    };
}

void writePlan(const std::string& filename, const std::vector<Position>& plan)
{
    std::ofstream out(filename, std::ios::binary);
    if (!out)
        throw std::runtime_error("Failed to open output file: " + filename);

    int32_t length = plan.size();
    out.write(reinterpret_cast<char*>(&length), sizeof(int32_t));

    for (const auto& pos : plan) {
        out.write(reinterpret_cast<const char*>(&pos.row), sizeof(int32_t));
        out.write(reinterpret_cast<const char*>(&pos.col), sizeof(int32_t));
    }

    out.close();
}

int main(int argc, char* argv[])
{
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