#include "Salesman.h"

#include <algorithm>
#include <cmath>
#include <limits>

// Creates an empty object.
Salesman::Salesman() : vertices(0) {}

// Implements brute-force algorithm.
std::vector<std::vector<double>> Salesman::exact(std::vector<std::vector<double>> const & coords) {
    vertices = std::ssize(coords);
    std::vector<int> tour(vertices);
    for (int i = 0; i != vertices; ++i) {
        tour[i] = i;
    }
    double length = measure_length(tour, coords);
    std::vector<int> optimal_tour = tour;
    double optimal_length = length;
    while (std::next_permutation(tour.begin(), tour.end())) {
        length = measure_length(tour, coords);
        if (optimal_length > length) {
            optimal_tour = tour;
            optimal_length = length;
        }
    }
    return convert_tour(optimal_tour, coords);
}

// Implements 2-opt algorithm.
std::vector<std::vector<double>> Salesman::approximate(std::vector<std::vector<double>> const & coords) {
    vertices = std::ssize(coords);
    std::vector<std::vector<int>> mst;
    build_mst(mst, coords);
    std::vector<std::vector<std::vector<int>>> graph(vertices);
    build_graph(graph, mst);
    std::vector<int> path;
    build_path(0, path, graph);
    std::vector<int> tour;
    build_tour(tour, path);
    return convert_tour(tour, coords);
}

// Implements nearest neighbor algorithm.
std::vector<std::vector<double>> Salesman::heuristic_first(std::vector<std::vector<double>> const & coords) {
    vertices = std::ssize(coords);
    std::vector<std::vector<double>> graph(vertices, std::vector<double> (vertices, 0));
    for (int i = 0; i != vertices; ++i) {
        for (int j = i + 1; j != vertices; ++j) {
            graph[i][j] = graph[j][i] = measure_distance(i, j, coords);
        }
    }
    std::vector<bool> color(vertices, true);
    std::vector<int> tour(vertices);
    int vertex = 0;
    color[vertex] = false;
    tour[0] = vertex;
    for (int i = 1; i != vertices; ++i) {
        int neighbor = -1;
        for (int j = 0; j != vertices; ++j) {
            if (color[j] and (neighbor == -1 or graph[vertex][neighbor] > graph[vertex][j])) {
                neighbor = j;
            }
        }
        vertex = neighbor;
        color[vertex] = false;
        tour[i] = vertex;
    }
    return convert_tour(tour, coords);
}

// Implements lexicographic sorting algorithm.
std::vector<std::vector<double>> Salesman::heuristic_second(std::vector<std::vector<double>> const & coords) {
    vertices = std::ssize(coords);
    std::vector<int> tour(vertices);
    for (int i = 0; i != vertices; ++i) {
        tour[i] = i;
    }
    std::sort(tour.begin(), tour.end(), [&] (int lhs, int rhs) { return coords[lhs] < coords[rhs]; });
    return convert_tour(tour, coords);
}

// Measures the distance between two vertices.
double Salesman::measure_distance(int x, int y, std::vector<std::vector<double>> const & coords) const {
    return std::sqrt(std::pow(coords[x][0] - coords[y][0], 2) + std::pow(coords[x][1] - coords[y][1], 2));
}

// Converts the tour from one view to another.
std::vector<std::vector<double>> Salesman::convert_tour(std::vector<int> const & tour, std::vector<std::vector<double>> const & coords) const {
    std::vector<std::vector<double>> cycle(vertices);
    for (int i = 0; i != vertices; ++i) {
        cycle[i] = coords[tour[i]];
    }
    return cycle;
}

// Measures the length of the tour.
double Salesman::measure_length(std::vector<int> const & tour, std::vector<std::vector<double>> const & coords) const {
    double length = 0;
    for (int i = 0; i != vertices; ++i) {
        length += measure_distance(tour[i], tour[(i + 1) % vertices], coords);
    }
    return length;
}

// Builds the minimal spanning tree of the graph.
void Salesman::build_mst(std::vector<std::vector<int>> & mst, std::vector<std::vector<double>> const & coords) const {
    std::vector<std::vector<double>> graph(vertices, std::vector<double> (vertices, 0));
    for (int i = 0; i != vertices; ++i) {
        for (int j = i + 1; j != vertices; ++j) {
            graph[i][j] = graph[j][i] = measure_distance(i, j, coords);
        }
    }
    std::vector<bool> color(vertices, true);
    std::vector<int> ancestor(vertices, -1);
    std::vector<double> distance(vertices, std::numeric_limits<double>::max());
    for (int i = 0; i != vertices; ++i) {
        int vertex = -1;
        for (int j = 0; j != vertices; ++j) {
            if (color[j] and (vertex == -1 or distance[vertex] > distance[j])) {
                vertex = j;
            }
        }
        color[vertex] = false;
        if (ancestor[vertex] != -1) {
            mst.push_back({ancestor[vertex], vertex});
        }
        for (int j = 0; j != vertices; ++j) {
            if (distance[j] > graph[vertex][j]) {
                ancestor[j] = vertex;
                distance[j] = graph[vertex][j];
            }
        }
    }
}

// Builds the graph by doubling the minimal spanning tree's edges.
void Salesman::build_graph(std::vector<std::vector<std::vector<int>>> & graph, std::vector<std::vector<int>> const & mst) const {
    for (std::vector<int> const & edge : mst) {
        int first = std::ssize(graph[edge[0]]);
        int second = std::ssize(graph[edge[1]]);
        graph[edge[0]].push_back({edge[1], second});
        graph[edge[0]].push_back({edge[1], second + 1});
        graph[edge[1]].push_back({edge[0], first});
        graph[edge[1]].push_back({edge[0], first + 1});
    }
}

// Builds the eulerian path in the graph.
void Salesman::build_path(int vertex, std::vector<int> & path, std::vector<std::vector<std::vector<int>>> & graph) const {
    for (std::vector<int> & neighbor : graph[vertex]) {
        if (neighbor[0] != -1) {
            int adjacent = neighbor[0];
            graph[neighbor[0]][neighbor[1]][0] = -1;
            neighbor[0] = -1;
            build_path(adjacent, path, graph);
        }
    }
    path.push_back(vertex);
}

// Builds the hamiltonian path from the eulerian path.
void Salesman::build_tour(std::vector<int> & tour, std::vector<int> const & path) const {
    std::vector<bool> color(vertices, true);
    for (int vertex : path) {
        if (color[vertex]) {
            tour.push_back(vertex);
            color[vertex] = false;
        }
    }
}