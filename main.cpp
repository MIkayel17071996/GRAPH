#include <iostream>
#include "UnweightedGraph.hpp"
#include <utility>

int main() {
    Graph<int> graph(7);

    // Add edges to create a more complex undirected graph
    graph.addEdge(0, 1, true);
    graph.addEdge(0, 2, true);
    graph.addEdge(1, 3, true);
    graph.addEdge(1, 4, true);
    graph.addEdge(2, 4, true);
    graph.addEdge(3, 5, true);
    graph.addEdge(4, 5, true);
    graph.addEdge(5, 6, true);

    std::cout << std::endl;

    // Find all paths from node 0 to node 6
    std::vector<std::vector<int>> paths = graph.getAllPossiblePaths(0, 6);

    // Print the paths
    for (const auto& path : paths) {
        for (int node : path) {
            std::cout << node << " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;

    // Find the shortest path (optional)


    return 0;
}
