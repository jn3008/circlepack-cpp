#pragma once
#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include "Vertex.h"

// -------------------------------------------------------------------
// Graph base class
// Though this class can be used to describe a general graph, the
// functionality assumes it is a triangulation of a closed topological
// disc.
//
// Member variables:
// - int n :
//      total number of vertices of the graph
//
// - std::vector<std::vector<int>> edges :
//      defines the adjacency of vertices in the class, can be
//      described as a triangular matrix since each std::vector<int>
//      in edges has increasing size, with incremements of 1,
//      beginning at 1. The reason to avoid a matrix is that the
//      graph is undirected so the matrix would be diagonally
//      symmetric.
//      The value of edges[i][j] denotes how many times i is
//      connected to j.
//
// - std::vector<Vertex *> vertices :
//      Pointers to Vertex objects that store information about each
//      vertex; its ordered and oriented petals (to avoid having to
//      calculate it more than once), its boundary status, and index.
//      We use pointers so that we can use polymorphism in the
//      derived graph class to pointer at derived vertex objects that
//      can hold more information such as radius and position of the
//      vertex (circle).
//
// Descriptions of each method are with their implementations.
class Graph
{
protected:
    int n;
    std::vector<Vertex *> vertices;

    virtual void add_adjacency(int vertex_index_1, int vertex_index_2);
    void remove_adjacency(int vertex_index_1, int vertex_index_2);

private:
    std::vector<std::vector<bool>> adjacencies;

    int dist(int vertex_index_1, int vertex_index_2) const;

    std::vector<int> find_unordered_petals(int vertex_index) const;
    std::vector<int> find_ordered_petals(int vertex_index) const;

    bool path_exists(std::vector<int> &vertex_indices, int current_idx, bool closed) const;
    void fix_orientation(int vertex_index, const std::array<int, 2> &dir, std::vector<bool> &fixed);

    void set_bound(int vertex_index, bool val);
    void set_petals(int vertex_index, const std::vector<int> &val);

public:
    // -------------------------------------------------------------------
    // Simple constructor for the class
    Graph()
    {
        n = 0;
    }

    //----------------------------------------------------------
    // Copy constructor for deep copy.
    Graph(const Graph &other) : n(other.n), adjacencies(other.adjacencies)
    {
        for (const Vertex *v : other.vertices)
            vertices.push_back(v->clone());
    }
    //----------------------------------------------------------
    // The pointers of the vertices vector point to heap objects
    // that need to be deleted manually.
    virtual ~Graph()
    {
        for (Vertex *v : vertices)
            delete v;
    }
    //----------------------------------------------------------
    // Assign operator
    Graph &operator=(const Graph &other)
    {
        if (this != &other)
        {
            // Assign dynamically allocated resources
            for (const Vertex *v : vertices)
                delete v;
            vertices.clear();
            for (const Vertex *v : other.vertices)
                vertices.push_back(v->clone());

            // Assign stack members
            n = other.n;
            adjacencies = other.adjacencies;
        }
        return *this;
    }
    //----------------------------------------------------------

    static Graph example_graph_AxB(int A, int B);

    const int &get_n() const { return n; }

    void auto_set_bounds();

    bool is_bound(int vertex_index) const;
    std::vector<int> get_petals(int vertex_index) const;
    Vertex *get_petal_pointer_of(int vertex_index, int petal_index) const;

    bool is_adjacent(int vertex_index_1, int vertex_index_2) const;
    bool exists(int vertex_index) const;

    void print_bounds() const;
    void print_adjacencies() const;
    void print_petals() const;

    void update_petals();
};