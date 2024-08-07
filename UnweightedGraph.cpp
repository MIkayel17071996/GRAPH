#include <iostream>
#include "UnweightedGraph.hpp"
#include <exception>
#include <algorithm>
#include <stack>
#include <queue>

/********* Constructors and destructor *********/

// Graph constructor initializes the size to 0
template <typename T>
Graph<T>::Graph()
    : m_size(0)
{
}

// Graph constructor that initializes the graph with the given number of vertices, each without a concrete value
template <typename T>
Graph<T>::Graph(const size_type count)
    : m_size(0)
{
    for (int i = 0; i < count; ++i)
    {
        this->addVertex();
    }
}

// Graph copy constructor
template <typename T>
Graph<T>::Graph(const Graph &other)
    : m_size(other.m_size)
{
    m_vertices.resize(other.m_size);
    for (size_type i = 0; i < other.m_size; ++i)
    {
        m_vertices[i] = new Vertex(other.m_vertices[i]->value);
    }

    m_adjacencyList = other.m_adjacencyList;
    m_adjacencyMatrix = other.m_adjacencyMatrix;
}

// Graph move constructor
template <typename T>
Graph<T>::Graph(Graph &&other) noexcept
{
    this->m_size = other.m_size;
    this->m_adjacencyList = std::move(other.m_adjacencyList);
    this->m_adjacencyMatrix = std::move(other.m_adjacencyMatrix);
    this->m_vertices = std::move(other.m_vertices);
    other.m_size = 0;
}

// Graph destructor deletes all the vertex pointers
template <typename T>
Graph<T>::~Graph()
{
    for (auto &elem : m_vertices)
    {
        delete elem;
    }
}

// Graph copy assignment operator
template <typename T>
typename Graph<T>::Graph &Graph<T>::operator=(const Graph &other)
{
    if (this != &other)
    {

        for (auto &elem : m_vertices)
        {
            delete elem;
        }

        Graph<T> tmp(other);

        swap(tmp);
    }

    return *this;
}

// Graph move assignment operator
template <typename T>
typename Graph<T>::Graph &Graph<T>::operator=(Graph &&other) noexcept
{
    for (auto &elem : m_vertices)
    {
        delete elem;
    }

    m_size = other.m_size;
    m_adjacencyList = std::move(other.m_adjacencyList);
    m_adjacencyMatrix = std::move(other.m_adjacencyMatrix);
    m_vertices = std::move(other.m_vertices);

    other.m_size = 0;

    return *this;
}

/******************* Modifiers *******************/

// Method to add a vertex to the graph
template <typename T>
void Graph<T>::addVertex()
{
    m_vertices.push_back(new Vertex(value_type{}));

    m_adjacencyList.push_back(std::vector<int>());

    m_adjacencyMatrix.push_back(std::vector<bool>(m_size, 0));

    for (auto &row : m_adjacencyMatrix)
    {
        row.push_back(0);
    }

    ++m_size;
}

// Method to add a vertex to the graph
template <typename T>
void Graph<T>::addVertex(const value_type &value)
{
    m_vertices.push_back(new Vertex(value));

    m_adjacencyList.push_back(std::vector<int>());

    m_adjacencyMatrix.push_back(std::vector<bool>(m_size, 0));

    for (auto &row : m_adjacencyMatrix)
    {
        row.push_back(0);
    }

    ++m_size;
}

// Method to add an edge between two vertices
template <typename T>
void Graph<T>::addEdge(const size_type &from, const size_type &to, bool flag)
{
    try
    {
        if (from >= m_size)
        {
            throw std::out_of_range("'From' pointer is out of the range");
        }
        else if (to >= m_size)
        {
            throw std::out_of_range("'To' pointer is out of the range");
        }

        auto find = std::find(m_adjacencyList[from].begin(), m_adjacencyList[from].end(), to);

        if (find == m_adjacencyList[from].end())
        {
            if (flag) // For undirected graph
            {
                m_adjacencyList[from].push_back(to);
                m_adjacencyList[to].push_back(from);

                m_adjacencyMatrix[from][to] = true;
                m_adjacencyMatrix[to][from] = true;
            }
            else // For directed graph
            {
                m_adjacencyList[from].push_back(to);
                m_adjacencyMatrix[from][to] = true;
            }
        }
    }
    catch (const std::out_of_range &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}
// Method for transpose the
template <typename T>
void Graph<T>::transpose()
{
    Graph<T> transposedGraph(m_size);

    for (size_t from = 0; from < m_size; ++from)
    {
        for (const auto &to : m_adjacencyList[from])
        {
            transposedGraph.addEdge(to, from);
        }
    }

    swap(transposedGraph);
}

/************* Printing Operations *************/

// Method to print the adjacency matrix
template <typename T>
void Graph<T>::printMatrix() const noexcept
{
    for (int i = 0; i < m_adjacencyMatrix.size(); ++i)
    {
        for (int j = 0; j < m_adjacencyMatrix[i].size(); ++j)
        {
            std::cout << m_adjacencyMatrix[i][j] << "\t";
        }
        std::cout << std::endl;
    }
}

// Method to print the adjacency list
template <typename T>
void Graph<T>::printList() const noexcept
{
    for (int i = 0; i < m_size; ++i)
    {
        std::cout << i << "-> ";
        size_t columnSize = m_adjacencyList[i].size();

        for (int j = 0; j < columnSize; ++j)
        {
            if (j == columnSize - 1)
            {
                std::cout << m_adjacencyList[i][j];
                break;
            }

            std::cout << m_adjacencyList[i][j] << ", ";
        }

        std::cout << std::endl;
    }
}

// Method to print the values of the vertices
template <typename T>
void Graph<T>::printValues() const noexcept
{
    for (int i = 0; i < m_vertices.size(); ++i)
    {
        std::cout << m_vertices[i]->value << " ";
    }

    std::cout << std::endl;
}

/******************* Capacity *******************/
// PUBLIC METHODS

// Method to get the number of connected components in the graph
template <typename T>
const typename Graph<T>::size_type Graph<T>::getComponentCount() const
{
    if (m_adjacencyList.empty())
    {
        return 0;
    }

    std::vector<bool> visited(m_size, false);
    size_type componentCount = 0;

    for (size_type i = 0; i < m_size; ++i)
    {
        if (!visited[i])
        {
            ++componentCount;
            getComponentCountHelper(i, visited);
        }
    }

    return componentCount;
}

// Method for getting vertices count in particular level (using BFS)
template <typename T>
const typename Graph<T>::size_type
Graph<T>::getNumberOfNodesAtAGivenLevelBFS(size_type level, size_type source) const
{
    try
    {
        if (source >= m_size)
        {
            throw std::out_of_range("'Source' pointer is out of range");
        }
        if (level < 0)
        {
            throw std::out_of_range("Level cannot be negative");
        }

        std::vector<bool> visited(m_size, false);
        std::queue<std::pair<int, int>> queue;

        queue.push({source, 0});
        visited[source] = true;

        size_type nodesCount = 0;

        while (!queue.empty())
        {
            auto [currentNode, currentLevel] = queue.front();
            queue.pop();

            if (currentLevel == level)
            {
                ++nodesCount;
            }
            else if (currentLevel < level)
            {
                for (const auto &neighbor : m_adjacencyList[currentNode])
                {
                    if (!visited[neighbor])
                    {
                        visited[neighbor] = true;
                        queue.push({neighbor, currentLevel + 1});
                    }
                }
            }
        }
        return nodesCount;
    }
    catch (const std::out_of_range &e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
    }

    return 0;
}

// Method for getting vertices count in particular level (using DFS)
template <typename T>
const typename Graph<T>::size_type
Graph<T>::getNumberOfNodesAtAGivenLevelDFS(size_type level, size_type source) const
{
    try
    {
        if (source >= m_size)
        {
            throw std::out_of_range("'Source' pointer is out of the range");
        }
        else if (level < 0)
        {
            throw std::out_of_range("Level cannot be negative");
        }

        std::vector<bool> visited(m_size, false);
        size_type count = 0;

        getNumberOfNodesAtAGivenLevelDFSHelper(source, level, 0, visited, count);

        return count;
    }
    catch (const std::out_of_range &e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
    }

    return 0;
}

// PRIVATE (HELPER) METHODS FOR CAPACITY METHODS

//  Helper method to count components
template <typename T>
void Graph<T>::getComponentCountHelper(size_type index, std::vector<bool> &visited) const
{
    visited[index] = true;

    for (const auto &neighbor : m_adjacencyList[index])
    {
        if (!visited[neighbor])
        {
            getComponentCountHelper(neighbor, visited);
        }
    }
}

// Helper method for getNumberOfNodesAtAGivenLevelDFS
template <typename T>
void Graph<T>::getNumberOfNodesAtAGivenLevelDFSHelper(size_type current, size_type targetLevel, size_type currentLevel, std::vector<bool> &visited, size_type &count)
{
    visited[current] = true;

    if (currentLevel == targetLevel)
    {
        ++count;
    }
    else if (currentLevel < targetLevel)
    {
        for (const auto &neighbor : m_adjacencyList[current])
        {
            if (!visited[neighbor])
            {
                getNumberOfNodesAtAGivenLevelDFSHelper(neighbor, targetLevel, currentLevel + 1, visited, count);
            }
        }
    }
}

////////////////////////////////////////////////////////////////
/******************* Operations *******************/

// Depth-first search (DFS) method
template <typename T>
void Graph<T>::dfs(size_type source, bool flag) const
{
    try
    {
        if (source >= m_size)
        {
            throw std::out_of_range("The index is out of the range");
        }

        if (m_size > 100)
        {
            flag = true;
        }

        std::vector<bool> visited(m_size, 0);

        if (!flag)
        {
            dfsHelperRecursive(source, visited);

            for (int i = 0; i < m_size; ++i)
            {
                if (!visited[i])
                {
                    dfsHelperRecursive(i, visited);
                }
            }
        }
        else
        {
            dfsHelperIterative(source, visited);

            for (int i = 0; i < m_size; ++i)
            {
                if (!visited[i])
                {
                    dfsHelperIterative(i, visited);
                }
            }
        }
    }
    catch (const std::out_of_range &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

// Breadth-first search (BFS) method
template <typename T>
void Graph<T>::bfs(size_type source) const
{
    try
    {
        if (source >= m_size)
        {
            throw std::out_of_range("The Index out of the range");
        }

        std::vector<bool> visited(m_size, 0);
        std::queue<int> queue;
        visited[source] = true;

        queue.push(source);

        while (!queue.empty())
        {
            int vertexIndex = queue.front();
            queue.pop();

            std::cout << vertexIndex << " ";

            for (const auto &neighbor : m_adjacencyList[vertexIndex])
            {
                if (!visited[neighbor])
                {
                    visited[neighbor] = true;
                    queue.push(neighbor);
                }
            }
        }
    }
    catch (const std::out_of_range &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

// Method for getting shortes path from source to destination
template <typename T>
const std::vector<int> Graph<T>::getShortestPath(size_type from, size_type to) const
{
    try
    {
        if (from >= m_size)
        {
            throw std::out_of_range("'From' pointer is out of the range");
        }
        else if (to >= m_size)
        {
            throw std::out_of_range("'To' pointer is out of the range");
        }

        if (from == to)
        {
            return {static_cast<int>(from)};
        }

        std::vector<bool> visited(m_size, false);
        std::vector<int> parent(m_size, -1);
        std::queue<int> queue;

        queue.push(from);
        visited[from] = true;

        while (!queue.empty())
        {
            int current = queue.front();
            queue.pop();

            for (const auto &neigbhor : m_adjacencyList[current])
            {
                if (!visited[neigbhor])
                {
                    visited[neigbhor] = true;
                    parent[neigbhor] = current;
                    queue.push(neigbhor);

                    if (neigbhor == to)
                    {
                        std::vector<int> path;
                        for (int v = to; v != -1; v = parent[v])
                        {
                            path.push_back(v);
                        }
                        std::reverse(path.begin(), path.end());
                        return path;
                    }
                }
            }
        }
    }
    catch (const std::out_of_range &e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
    }
    return {};
}

// Method for getting all possible paths from sourceto destination
template <typename T>
const std::vector<std::vector<int>> Graph<T>::getAllPossiblePaths(size_type source, size_type destination) const
{
    std::vector<std::vector<int>> allPaths;
    if (source >= m_size || destination >= m_size)
    {
        std::cerr << "Source or destination out of range." << std::endl;
        return allPaths;
    }

    std::vector<bool> visited(m_size, false);
    std::vector<int> path;

    getAllPossiblePathsHelper(source, destination, visited, path, allPaths);

    return allPaths;
}

// Method to detect cycles in an undirected graph
template <typename T>
bool Graph<T>::isCycledUndirected() const
{
    std::vector<bool> visited(m_size, false);

    for (size_type i = 0; i < m_size; ++i)
    {
        if (!visited[i])
        {
            if (isCycledUndirectedHelper(i, visited, -1))
            {
                return true;
            }
        }
    }
    return false;
}

template <typename T>
bool Graph<T>::isCycledDirected() const
{
    std::vector<bool> visited(m_size, false);
    std::vector<bool> recStack(m_size, false);

    for (size_type i = 0; i < m_size; ++i)
    {
        if (isCycledDirectedHelper(i, visited, recStack))
        {
            return true;
        }
    }
    return false;
}

// PRIVATE (HELPER) METHODS FOR OPERATION METHODS

// Helper method for iterative DFS
template <typename T>
void Graph<T>::dfsHelperIterative(size_type source, std::vector<bool> &visited) const
{
    std::stack<int> stack;
    visited[source] = true;

    stack.push(source);

    while (!stack.empty())
    {
        int vertexIndex = stack.top();
        stack.pop();

        std::cout << vertexIndex << " ";

        for (const auto &neighbor : m_adjacencyList[vertexIndex])
        {
            if (!visited[neighbor])
            {
                visited[neighbor] = true;
                stack.push(neighbor);
            }
        }
    }
}

// Helper method for recursive DFS
template <typename T>
void Graph<T>::dfsHelperRecursive(size_type source, std::vector<bool> &visited) const
{
    std::cout << source << " ";
    visited[source] = true;

    for (const auto &neighbor : m_adjacencyList[source])
    {
        if (!visited[neighbor])
        {
            dfsHelperRecursive(neighbor, visited);
        }
    }
}

// Helper method for getAllPossiblePaths
template <typename T>
void Graph<T>::getAllPossiblePathsHelper(size_type current, size_type destination, std::vector<bool> &visited,
                                         std::vector<int> &path, std::vector<std::vector<int>> &allPaths) const
{
    visited[current] = true;
    path.push_back(current);

    if (current == destination)
    {
        allPaths.push_back(path);
    }
    else
    {
        for (const auto &neighbor : m_adjacencyList[current])
        {
            if (!visited[neighbor])
            {
                getAllPossiblePathsHelper(neighbor, destination, visited, path, allPaths);
            }
        }
    }

    path.pop_back();
    visited[current] = false;
}

// Helper method for isCycledUndirected
template <typename T>
bool Graph<T>::isCycledUndirectedHelper(size_type v, std::vector<bool> &visited, size_type parent) const
{
    visited[v] = true;

    for (const auto &neighbor : m_adjacencyList[v])
    {
        if (!visited[neighbor])
        {
            if (isCycledUndirectedHelper(neighbor, visited, v))
            {
                return true;
            }
        }
        else if (neighbor != parent)
        {
            return true;
        }
    }
    return false;
}

template <typename T>
bool Graph<T>::isCycledDirectedHelper(size_type v, std::vector<bool>& visited, std::vector<bool>& recStack) const
{
    if (!visited[v])
    {
        // Mark the current node as visited and add it to the recursion stack
        visited[v] = true;
        recStack[v] = true;

        // Recur for all the vertices adjacent to this vertex
        for (const auto& neighbor : m_adjacencyList[v])
        {
            if (!visited[neighbor] && isCycledDirectedHelper(neighbor, visited, recStack))
            {
                return true;
            }
            else if (recStack[neighbor])
            {
                return true;
            }
        }
    }

    // Remove the vertex from recursion stack
    recStack[v] = false;
    return false;
}

// Method to swap two graphs
template <typename T>
void Graph<T>::swap(Graph &second) noexcept
{
    using std::swap;
    swap(this->m_size, second.m_size);
    swap(this->m_adjacencyList, second.m_adjacencyList);
    swap(this->m_adjacencyMatrix, second.m_adjacencyMatrix);
    swap(this->m_vertices, second.m_vertices);
}
