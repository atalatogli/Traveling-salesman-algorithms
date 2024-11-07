#pragma once

#include <vector>

struct Salesman {
public:
    // Creates an empty object.
    Salesman();

    // Implements brute-force algorithm.
    std::vector<std::vector<double>> exact(std::vector<std::vector<double>> const & coords);

    // Implements 2-opt algorithm.
    std::vector<std::vector<double>> approximate(std::vector<std::vector<double>> const & coords);

    // Implements nearest neighbor algorithm.
    std::vector<std::vector<double>> heuristic_first(std::vector<std::vector<double>> const & coords);

    // Implements lexicographic sorting algorithm.
    std::vector<std::vector<double>> heuristic_second(std::vector<std::vector<double>> const & coords);
private:
    // Disables copy constructor for our class.
    Salesman(Salesman const & other);

    // Disables assignment operator for our class.
    Salesman & operator = (Salesman const & other);

    // Measures the distance between two vertices.
    double measure_distance(int x, int y, std::vector<std::vector<double>> const & coords) const;

    // Converts the tour from one view to another.
    std::vector<std::vector<double>> convert_tour(std::vector<int> const & tour, std::vector<std::vector<double>> const & coords) const;

    // Measures the length of the tour.
    double measure_length(std::vector<int> const & tour, std::vector<std::vector<double>> const & coords) const;
    
    // Builds the minimal spanning tree of the graph.
    void build_mst(std::vector<std::vector<int>> & mst, std::vector<std::vector<double>> const & coords) const;
    
    // Builds the graph by doubling the minimal spanning tree's edges.
    void build_graph(std::vector<std::vector<std::vector<int>>> & graph, std::vector<std::vector<int>> const & mst) const;
    
    // Builds the eulerian path in the graph.
    void build_path(int vertex, std::vector<int> & path, std::vector<std::vector<std::vector<int>>> & graph) const;

    // Builds the hamiltonian path from the eulerian path.
    void build_tour(std::vector<int> & tour, std::vector<int> const & path) const;
    
    int vertices;
};