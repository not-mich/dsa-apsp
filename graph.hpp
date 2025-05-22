#ifndef GRAPH_HPP_
#define GRAPH_HPP_

#include <iostream>
#include <fstream>
#include <utility>
#include <functional>
#include <vector>
#include <string>
#include <queue>
#include <unordered_map>
#include <limits>

template <typename T>
class Graph {
 private:
  std::vector<std::unordered_map<int, T> > adjList {};
  int numVertices {};

 public:
  // empty graph with N vertices
  explicit Graph(int N);

  // construct graph from edge list in filename
  explicit Graph(const std::string& filename);

  // add an edge directed from vertex i to vertex j with given weight
  void addEdge(int i, int j, T weight);

  // removes edge from vertex i to vertex j
  void removeEdge(int i, int j);

  // is there an edge from vertex i to vertex j?
  bool isEdge(int i, int j) const;

  // return weight of edge from i to j
  // will throw an exception if there is no edge from i to j
  T getEdgeWeight(int i, int j) const;

  // returns number of vertices in the graph
  int size() const;

  // return iterator to a particular vertex
  const std::unordered_map<int, T>& neighbours(int a) const {
    return adjList.at(a);
  }
};

template <typename T>
Graph<T>::Graph(int N) : adjList(N), numVertices {N} {}

template <typename T>
Graph<T>::Graph(const std::string& inputFile) {
  std::ifstream infile {inputFile};
  if (!infile) {
    std::cerr << inputFile << " could not be opened\n";
    return;
  }
  // first line has number of vertices
  infile >> numVertices;
  adjList.resize(numVertices);
  int i {};
  int j {};
  double weight {};
  // assume each remaining line is of form
  // origin dest weight
  while (infile >> i >> j >> weight) {
    addEdge(i, j, static_cast<T>(weight));
  }
}

template <typename T>
int Graph<T>::size() const {
  return numVertices;
}

template <typename T>
void Graph<T>::addEdge(int i, int j, T weight) {
  if (i < 0 or i >= numVertices or j < 0 or j >= numVertices) {
    throw std::out_of_range("invalid vertex number");
  }
  adjList[i].insert({j, weight});
}

template <typename T>
void Graph<T>::removeEdge(int i, int j) {
  // check if i and j are valid
  if (i >= 0 && i < numVertices && j >= 0 && j < numVertices) {
    adjList[i].erase(j);
  }
}

template <typename T>
bool Graph<T>::isEdge(int i, int j) const {
  if (i >= 0 && i < numVertices && j >= 0 && j < numVertices) {
    return adjList.at(i).contains(j);
  }
  return false;
}

template <typename T>
T Graph<T>::getEdgeWeight(int i, int j) const {
  return adjList.at(i).at(j);
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const Graph<T>& G) {
  for (int i = 0; i < G.size(); ++i) {
    out << i << ':';
    for (const auto& [neighbour, weight] : G.neighbours(i)) {
      out << " (" << i << ", " << neighbour << ")[" << weight << ']';
    }
    out << '\n';
  }
  return out;
}


// APSP functions
// Use this function to return an "infinity" value
// appropriate for the type T
template <typename T>
T infinity() {
  if (std::numeric_limits<T>::has_infinity) {
    return std::numeric_limits<T>::infinity();
  } else {
    return std::numeric_limits<T>::max();
  }
}

// implement an algorithm for determining if G
// has a negative weight cycle here
template <typename T>
bool existsNegativeCycle(const Graph<T>& G) {
  T inf = infinity<T>();
  const int numVertices = G.size();
  std::vector<T> shortestDist (numVertices, 0);

  // using Bellman-Ford's algorithm, repeat (Vertices - 1) times
  // relax the edges, look for the shortest path, dist[v] > dist[u] + weight
  for (int i = 0; i < numVertices; ++i) {
    for (int u = 0; u < numVertices; ++u) {
      for (const auto& [v, weight] : G.neighbours(u)) {
        if (shortestDist[u] < inf && shortestDist[v] > shortestDist[u] + weight) {
          shortestDist[v] = shortestDist[u] + weight;
        }
      }
    }
  }

  // if there is a negative weight cycle and return true, otherwise false
  for (int u = 0; u < numVertices; ++ u) {
    for (const auto& [v, weight] : G.neighbours(u)) {
      if (shortestDist[u] < inf && shortestDist[v] > shortestDist[u] + weight) {
        return true;
      }
    }
  }
  return false;
}

// implement Johnson's APSP algorithm here
template <typename T>
std::vector<std::vector<T> >
johnsonAPSP(const Graph<T>& G) {
  std::ignore = G;
  return {};
}

// implement the Floyd-Warshall APSP algorithm here
template <typename T>
std::vector<std::vector<T> >
floydWarshallAPSP(const Graph<T>& G) {
  T inf = infinity<T>();
  const int numVertices = G.size();
  std::vector<std::vector<T>> distanceMatrix(numVertices, std::vector<T>(numVertices, inf));

  // initialisation of distanceMatrix
  for (int i = 0; i < numVertices; i++) {
    distanceMatrix[i][i] = 0;
    for (const auto& [j, weight] : G.neighbours(i)) {
      distanceMatrix[i][j] = weight;
    }
  }

  // implementation of floyd-warshall
  // iterate through all of the pairs of vertices (i, j) and intermediate k vertex (outer-most loop) 
  for (int k = 0; k < numVertices; k++) {
    for (int i = 0; i < numVertices; i++) {
      for (int j = 0; j < numVertices; j++) {
        if (distanceMatrix[i][k] < inf && distanceMatrix[k][j] < inf) {
          // check path through vertices (i, j) and if path through k is shorter
          // if so, update the shortest distance
          if (distanceMatrix[i][k] + distanceMatrix[k][j] < distanceMatrix[i][j]) {
            distanceMatrix[i][j] = distanceMatrix[i][k] + distanceMatrix[k][j];
          }
        }
      }
    }
  }
  return distanceMatrix;
}

#endif      // GRAPH_HPP_
