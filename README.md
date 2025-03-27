# 15418-finalproj

### Title
Parallelized Path Planning Algorithms, by Aidan Lincke (alincke) and Soren Dupont (sdupont)

### URL
https://aidanlincke.github.io/15418-finalproj/

### Summary
Our plan for this project is to implement a number of path planning algorithms on graphs that have some level of parallelism that can be exploited, but that is difficult to extract, and compare which algorithms best lend themselves to parallel speedup. We will use OpenMP to program familiar graph search algorithms such as Dijkstra's algorithm. Probabilistic roadmaps, which involves a higher degree of parallelism for generating random point, may be well suited for parallism via CUDA.

### Background


### The Challenge
This problem is challenging because while there is some parallelism in graph path-finding algorithms that can be exploited, there are many sensitive shared resources that processors must have access to in order to make progress, such as the graph itself. For Dijkstras algorithm, for example, there are also a shared set of visited nodes, a hashmap, and priority queue for which node should be accessed next. Allowing processors to work with slightly stale versions of these resources could be a solution, so long as each processors maintains consistent versions of those resources, as could changing how a resource like the priority queue is accessed (for example, perhaps elements of it could be assigned to certain nodes, and sent asynchronously to reduce contention over the priority queue). This algorithm introduces a complex set of tradeoffs that may need to be considered and tested for.

Algorithms like A* have parallelism much more inherently build in, but only offer heuristical solutions to pathfinding, and may present challenges for load balancing, as in A* static assignment is required. We are considering certain precomputations that might be done to help with load balancing in this case. Graph contraction techniques could play a role in this, and that itself is a difficult parallelizeable task, as it involves many nodes reading from and modifying a graph concurrently.

### Resources

### Goals and Deliverables

### Platform Choice

### Schedule
