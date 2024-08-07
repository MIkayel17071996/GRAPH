#ifndef UNWEIGHTEDGRAPH_HPP
#define UNWEIGHTEDGRAPH_HPP

#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>

// Template class Graph to create a graph with vertices of type T
template <class T>
class Graph
{
public:
    // Type aliases for convenience
    using value_type = T;
    using pointer = T *;
    using size_type = std::size_t;
    using const_pointer = const T *;
    using reference = T &;
    using const_reference = const T &;

private:
    // Internal structure representing a vertex
    struct Vertex
    {
        value_type value;

        Vertex(const_reference val)
            : value(val)
        {
        }
    };

public:
    // Constructors and destructor
    Graph();                             // Default constructor
    explicit Graph(const size_type);     // Constructor with initial size
    Graph(const Graph &);                // Copy constructor
    Graph(Graph &&) noexcept;            // Move constructor
    ~Graph();                            // Destructor
    Graph &operator=(const Graph &);     // Copy assignment operator
    Graph &operator=(Graph &&) noexcept; // Move assignment operator

private:
    // Delete conversion operator
    operator T() = delete;

public:
    // Method to add a vertex without a value
    void addVertex();
    // Method to add a vertex with a value
    void addVertex(const value_type &);
    // Method to add an edge
    void addEdge(const size_type &, const size_type &, bool = true);
    // Transpose the graph
    void transpose();

    // Method to print adjacency matrix
    void printMatrix() const noexcept;
    // Method to print adjacency list
    void printList() const noexcept;
    // Method to print vertex values
    void printValues() const noexcept;

    // Get number of connected components
    const size_type getComponentCount() const;
    // Get number of nodes at a given level using BFS
    const size_type getNumberOfNodesAtAGivenLevelBFS(size_type, size_type = 0) const;
    // Get number of nodes at a given level using DFS
    const size_type getNumberOfNodesAtAGivenLevelDFS(size_type, size_type = 0) const;

    // Depth-first search method
    void dfs(size_type, bool = false) const;
    // Breadth-first search method
    void bfs(size_type) const;
    // Get shortest path
    const std::vector<int> getShortestPath(size_type, size_type) const;
    // Get all possible paths
    const std::vector<std::vector<int>> getAllPossiblePaths(size_type, size_type) const;
    // Check for cycle in an undirected graph
    bool isCycledUndirected() const;
    // Check for cycle in a directed graph
    bool isCycledDirected() const;

private:
    // Helper for component count
    void getComponentCountHelper(size_type, std::vector<bool> &) const;
    // Helper for nodes at level using DFS
    void getNumberOfNodesAtAGivenLevelDFSHelper(size_type, size_type, size_type, std::vector<bool> &, size_type &);
    // Helper for iterative DFS
    void dfsHelperIterative(size_type, std::vector<bool> &) const;
    // Helper for recursive DFS
    void dfsHelperRecursive(size_type, std::vector<bool> &) const;
    // Helper for all possible paths
    void getAllPossiblePathsHelper(size_type, size_type, std::vector<bool> &, std::vector<int> &, std::vector<std::vector<int>> &) const;
    // Helper for cycle detection in undirected graph
    bool isCycledUndirectedHelper(size_type, std::vector<bool> &, size_type) const;
    // Helper for cycle detection in directed graph
    bool isCycledDirectedHelper(size_type, std::vector<bool> &, std::vector<bool> &) const;
    // Swap method
    void swap(Graph &) noexcept;

private:
    size_type m_size;                                 // Number of vertices
    std::vector<std::vector<int>> m_adjacencyList;    // Adjacency list
    std::vector<std::vector<bool>> m_adjacencyMatrix; // Adjacency matrix
    std::vector<Vertex *> m_vertices;                 // List of vertices
};

#include "UnweightedGraph.cpp" // Include the implementation file
#endif                         // UNWEIGHTEDGRAPH_HPP
