#pragma once
#include <iostream>
#include <vector>
#include <array>
#include <stack>
#include <algorithm>
#include <numeric>
#include "Graph.h"
#include "CirclePackVertex.h"

// -------------------------------------------------------------------
// Small struct to describe an operation to the graph, specifically,
// a connection added or removed between two vertices.
// 0 - Nothing (for completeness), 1 - add_adjacency, 2 - remove_adjacency
struct GraphOperation
{
    int operation_type, vertex_index_1, vertex_index_2;
};

// -------------------------------------------------------------------
// CirclePackGraph class
// Inherits from Graph class, adding specific functionality of the
// Graph in the context of computing a circle packing.
//
// Descriptions of each method are with the implementations
class CirclePackGraph : public Graph
{
private:
    // view_scale and view_translation allow the user to move around
    // and zoom into the packing.
    float view_scale = 1;
    ofVec2f view_translation;
    // history_stack allows the user to undo the last edit to the graph
    // Each edit is a vector of GraphOperation object which
    // represents the type of operation and two vertex indices.
    std::stack<std::vector<GraphOperation>> history_stack;

    static double angle(double v, double u, double w);
    double theta(int vertex_index) const;

    void compute_positions(int vertex_index, std::vector<bool> &fixed);

    void recenter_scale(float width, float height);

protected:
    void add_adjacency(int vertex_index_1, int vertex_index_2) override;

public:
    CirclePackGraph() : Graph(){};

    // ------------------------------------------------------------
    // Construct instance from an instance of the base class Graph.
    // Turn all Vertex* pointers which point to Vertex objects in
    // the original Graph instance to point to CirclePackVertex
    // objects which retain all vertex graph-related information
    // but now can hold information about position, label (radius)
    // and angle sum target.
    CirclePackGraph(const Graph &graph) : Graph(graph)
    {
        for (Vertex *&v : vertices)
        {
            Vertex *upgrade = new CirclePackVertex(*v);
            delete v;
            v = upgrade;
        }
    }

    double get_label(int vertex_index) const;
    void set_label(int vertex_index, double val);

    ofVec2f get_pos(int vertex_index) const;
    void set_pos(int vertex_index, const ofVec2f &val);

    double get_target_ang(int vertex_index) const;
    void set_target_ang(int vertex_index, double val);
    void adjust_target_ang(int boundary_vertex_index, double adjustment);
    void auto_set_targets();
    void set_ngon_target_ang(int ngon);

    void print_labels() const;

    void compute_labels_by_fixing_boundary_labels();
    void compute_labels_by_fixing_boundary_targets();
    void compute_positions();

    double get_error() const;

    void remove_and_fill(int vertex_index, std::vector<GraphOperation> &edit_to_graph);
    void add_vertex(int vertex_index_1, int vertex_index_2, int where);
    void flip_edge(int vertex_index_1, int vertex_index_2);

    

    Vertex *get_circle_at_position(ofVec2f pos) const;

    float get_view_scale() const;
    void set_view_scale(float val);
    ofVec2f get_view_translation() const;
    void set_view_translation(ofVec2f val);

    void recenter_scale();

    void undo();
};