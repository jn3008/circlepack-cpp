#include "CirclePackGraph.h"

// -------------------------------------------------------------------
// Getter and setter functions:
// Each has an assertion to ensure a valid vertex_index is given.
// These methods use methods of the CirclePackVertex class so
// the pointers of the base class Vertex need to be casted. We
// can use a static cast because we expect the pointers here to only
// point to instances of the derived class.
// -------------------------------------------------------------------
// Return copy of the label of a vertex with given index
double CirclePackGraph::get_label(int vertex_index) const
{
    assert(vertex_index >= 0 && vertex_index < get_n());
    return static_cast<CirclePackVertex *>(vertices[vertex_index])->get_label();
}
// -------------------------------------------------------------------
// Set the value of the label of a vertex with given index
void CirclePackGraph::set_label(int vertex_index, double val)
{
    assert(vertex_index >= 0 && vertex_index < get_n());
    static_cast<CirclePackVertex *>(vertices[vertex_index])->set_label(val);
}
// -------------------------------------------------------------------
// Return copy of the position of a vertex with given index
ofVec2f CirclePackGraph::get_pos(int vertex_index) const
{
    assert(vertex_index >= 0 && vertex_index < get_n());
    return static_cast<CirclePackVertex *>(vertices[vertex_index])->get_pos();
}
// -------------------------------------------------------------------
// Set the value of the position of a vertex with given index
void CirclePackGraph::set_pos(int vertex_index, const ofVec2f &val)
{
    assert(vertex_index >= 0 && vertex_index < get_n());
    static_cast<CirclePackVertex *>(vertices[vertex_index])->set_pos(val);
}
// -------------------------------------------------------------------
// Return copy of the target angle of a vertex with given index
double CirclePackGraph::get_target_ang(int vertex_index) const
{
    assert(vertex_index >= 0 && vertex_index < get_n());
    return static_cast<CirclePackVertex *>(vertices[vertex_index])->get_target_ang();
}
// -------------------------------------------------------------------
// Set the value of the target angle of a vertex with given index
void CirclePackGraph::set_target_ang(int vertex_index, double val)
{
    assert(vertex_index >= 0 && vertex_index < get_n());
    static_cast<CirclePackVertex *>(vertices[vertex_index])->set_target_ang(val);
}
// -------------------------------------------------------------------
// Adjust the angle sum target of a boundary vertex with given index and
// given adjustment. The total target angle sums should be constant
// (precisely (num_boundaries-2)*PI), so when adjusting this target angle
// sum, the difference should be recompensated across all other boundary vertices
void CirclePackGraph::adjust_target_ang(int boundary_vertex_index, double adjustment)
{
    assert(boundary_vertex_index >= 0 && boundary_vertex_index < get_n());
    set_target_ang(boundary_vertex_index,
                   get_target_ang(boundary_vertex_index) + adjustment);

    int n = get_n();
    // Count the number of boundary vertices
    int num_boundaries = 0;
    for (int i = 0; i < n; i++)
        num_boundaries += is_bound(i);

    for (int i = 0; i < n; i++)
        if (is_bound(i) && i != boundary_vertex_index)
            set_target_ang(i, get_target_ang(i) - adjustment / (num_boundaries - 1));
}
// -------------------------------------------------------------------
void CirclePackGraph::auto_set_targets()
{

    for (int i = 0; i < get_n(); i++)
    {
        if (!is_bound(i))
            set_target_ang(i, 2 * PI);
        else
            set_target_ang(i, theta(i));
    }
}

// -------------------------------------------------------------------
void CirclePackGraph::set_ngon_target_ang(int ngon)
{
    // Check there are vertices
    if (get_n() < 1)
        return;

    // Find all boundary vertices in order.
    // To do this, first find any one, then keeping adding
    // the latest one's first petal until we reach the first again
    int first_boundary_index = 0;
    while (!is_bound(first_boundary_index))
        first_boundary_index++;

    std::vector<int> boundary_indices{first_boundary_index};
    while (get_petals(boundary_indices.back())[0] != boundary_indices.front())
        boundary_indices.push_back(get_petals(boundary_indices.back())[0]);

    // Check there's enough boundary vertices
    if (boundary_indices.size() < ngon)
        return;

    // Now set the angle sum targets. We must choose 'ngon' vertices to be 
    // the corners, which will each have angle PI*(ngon-2) / ngon while
    // the rest will be PI (flat edge).

    float angle = PI*(ngon-2)/ngon;
    for (int i = 0; i < ngon; i++)
    {
        int rdm_index = rand()%boundary_indices.size();
        set_target_ang(boundary_indices[rdm_index], angle);
        boundary_indices.erase(boundary_indices.begin()+rdm_index);
    }
    while (!boundary_indices.empty())
    {
        set_target_ang(boundary_indices.back(), PI);
        boundary_indices.pop_back();
    }
}

// -------------------------------------------------------------------
// Print the label of each vertex
void CirclePackGraph::print_labels() const
{
    std::cout << "----------------vertex-labels-----------------" << std::endl;
    for (int i = 0; i < get_n(); i++)
        std::cout << i << ": " << get_label(i) << ", ";
    std::cout << std::endl;
}

// -------------------------------------------------------------------
// Override Graph::add_adjacency by calling it and turning any
// new pointers to Vertex object into CirclePackVertex objects.
void CirclePackGraph::add_adjacency(int vertex_index_1, int vertex_index_2)
{
    bool new_vertex_added = max(vertex_index_1, vertex_index_2) >= get_n();

    Graph::add_adjacency(vertex_index_1, vertex_index_2);

    if (!new_vertex_added)
        return;

    // If a new vertex is added, it will be a Vertex object so
    // turn it into a CirclePackVertex object
    Vertex *&v = vertices.back();
    Vertex *upgrade = new CirclePackVertex(*v);
    delete v;
    v = upgrade;
}

// -------------------------------------------------------------------
// Iterate through each interior vertex (not boundary vertex) to
// adjust its label, "relaxation" algorithm based on Uniform
// Neighbour Model
void CirclePackGraph::compute_labels_by_fixing_boundary_labels()
{
    for (int i = 0; i < get_n(); i++)
    {
        // Skip if the i'th vertex is a boundary or has no petals
        if (is_bound(i) || !exists(i))
            continue;

        int k = get_petals(i).size();
        double th = theta(i);
        double beta = sin(th * 0.5 / k);
        double delta = sin(PI / k);
        double r_hat = beta / (1 - beta) * get_label(i);
        double x = (1 - delta) / delta * r_hat;

        set_label(i, x);
    }
}

// -------------------------------------------------------------------
// Iterate through each vertex to adjust its label, "relaxation"
// algorithm based on Uniform Neighbour Model. Doesn't ignore boundary
// vertices, but compute the adjusted label based on its target angle
// which isn't 2*PI (unlike interior vertices which are 2*PI by default)
void CirclePackGraph::compute_labels_by_fixing_boundary_targets()
{
    for (int i = 0; i < get_n(); i++)
    {
        // Skip if the i'th vertex has no petals
        if (!exists(i))
            continue;

        int k = get_petals(i).size();
        double beta = sin(theta(i) * 0.5 / k);
        double target_ang = is_bound(i) ? get_target_ang(i) : PI * 2;
        double delta = sin(target_ang * 0.5 / k);
        double r_hat = beta / (1 - beta) * get_label(i);
        double x = (1 - delta) / delta * r_hat;
        set_label(i, x);
    }
}

// -------------------------------------------------------------------
// Calculate and update the positions of vertices (circles). This
// calculation is "trivial" given the labels are correct.
// It is done by first fixing the positions of two adjacent vertices,
// I picked the first vertex and its first neighbour.
// Then, recursively fixing the positions of the common vertices of
// each pair, by keeping track of indices of vertices whose positions
// have been fixed.
void CirclePackGraph::compute_positions()
{
    std::vector<bool> fixed(get_n(), false);

    // Fix the position of two vertices.
    int A = 0;
    while (!exists(A))
        A++;
    int B = get_petals(A)[0];
    set_pos(A, ofVec2f());
    set_pos(B, ofVec2f(get_label(A) + get_label(B), 0));
    fixed[A] = true;
    fixed[B] = true;

    // Call a recursive function to fix the positions of the
    // petals of vertex A
    compute_positions(A, fixed);

    // Translate and scale positions to fit the canvas nicely
    // (and scale labels accordingly).
    recenter_scale();
}

// -------------------------------------------------------------------
// Calculate and return the error of the packing. This is done by
// summing the difference of each interior vertex's angle sum and
// angle sum target (which is probably 2*PI).
double CirclePackGraph::get_error() const
{
    double error = 0;
    for (int i = 0; i < get_n(); i++)
        if (!is_bound(i) && exists(i))
            error += get_target_ang(i) - theta(i);
    return abs(error);
}

// -------------------------------------------------------------------
// Remove a vertex from the packing and "sew" the hole shut by
// connecting its petals together in a zig-zag fashion.
// The param 'edit_to_graph' is to information about the edit that will be made
// to graph, in order to push to the history stack, so that it can be undone.
// The reason it's passed as a references is that if this function needs to call itself
// again then we group all the operations to the graph as a single edit, so that they
// can all be undone at once.
void CirclePackGraph::remove_and_fill(int vertex_index, std::vector<GraphOperation> &edit_to_graph)
{
    assert(vertex_index >= 0 && vertex_index < get_n());
    std::vector<int> petals = get_petals(vertex_index);

    // Only push this edit to the history stack if the function
    // wasn't called by itself.
    bool push_edit_to_history_stack = edit_to_graph.empty();

    int len = petals.size();
    std::cout << "Remove" << (is_bound(vertex_index) ? " boundary " : " ")
              << "vertex " << vertex_index
              << ", #petals: " << len << std::endl;

    // check first whether adding these edges will cause a problem, namely
    // adding an adjacency that already exists, it should be impossible to
    // have a double connection on the plane for a circle packing (not for
    // planar graphs in general though). Keep trying with different offsets
    // for making the connections, until giving up.
    bool problem = false;
    if (len > 3 && !is_bound(vertex_index))
    {
        std::vector<std::array<int, 2>> pairs;
        for (int i = 0; i < len - 3; i++)
            pairs.push_back(std::array<int, 2>{1 + i % 2 + i / 2, len - i / 2 - 1});

        for (int offset = 0; offset < len; offset++)
        {
            problem = false;
            for (int i = 0; i < len - 3; i++)
            {
                if (is_adjacent(petals[(pairs[i][0] + offset) % len],
                                petals[(pairs[i][1] + offset) % len]))
                {
                    problem = true;
                    // std::cout << "problem, " << offset << std::endl;
                }
            }

            if (problem)
                continue;

            // Sew hole shut
            if (len > 3)
                for (int i = 0; i < len - 3; i++)
                {
                    int index_1 = petals[(pairs[i][0] + offset) % len];
                    int index_2 = petals[(pairs[i][1] + offset) % len];
                    add_adjacency(index_1, index_2);
                    edit_to_graph.push_back(GraphOperation{.operation_type = 1,
                                                           .vertex_index_1 = index_1,
                                                           .vertex_index_2 = index_2});
                }
            break;
        }
    }

    if (problem)
        return;

    // Remove the adjacencies of the vertex we're removing
    for (int i = 0; i < len; i++)
    {
        remove_adjacency(petals[i], vertex_index);
        edit_to_graph.push_back(GraphOperation{.operation_type = 2,
                                               .vertex_index_1 = petals[i],
                                               .vertex_index_2 = vertex_index});
    }

    // Update the graph structure
    auto_set_bounds();
    update_petals();

    // If removing the current vertex leaves another vertex with
    // two petals left, remove it also.
    for (int petal : petals)
        if (get_petals(petal).size() < 3)
        {
            std::cout << "Removing " << vertex_index << " also removes " << petal << std::endl;
            remove_and_fill(petal, edit_to_graph);
        }

    if (push_edit_to_history_stack)
        history_stack.push(edit_to_graph);
}

// -------------------------------------------------------------------
// Given an edge between vertices with specified indices, add a new vertex.
// If 'where == 0', add the vertex in the middle of the 'left' face (triangle)
// shared by both vertices. if 'where == 1', do the same on the 'right' side.
// Note, if both vertices are boundary vertices then only one face is available.
//      1
//    / │ \ 
//   L  │  R
//    \ │ /
//      2
// If 'where == 2', add a vertex between both these vertices by removing the edge
// and connected the new vertex to these two vertices and their two mutual petals
// (one petal if both vertices are boundary).
void CirclePackGraph::add_vertex(int vertex_index_1, int vertex_index_2, int where)
{
    if (!is_adjacent(vertex_index_1, vertex_index_2)) // Check function is well-defined.
        return;

    // Store information about the edit so that it can be undone later.
    std::vector<GraphOperation> edit_to_graph;

    // Make an array to access the indices more easily.
    // (indices[i] and indices[1-i])
    std::array<int, 2> indices{vertex_index_1, vertex_index_2};

    if (where < 2)
    {
        // If there is no face to add a vertex to.
        if (is_bound(indices[where]) &&
            is_bound(indices[1 - where]) &&
            get_petals(indices[1 - where])[0] == indices[where])
            return;

        std::vector<int> petals = get_petals(indices[1 - where]);
        std::vector<int>::iterator it = std::find(petals.begin(), petals.end(), indices[where]);
        int petal_index = std::distance(petals.begin(), it);
        int vertex_index_3 = petals[(petal_index + 1) % petals.size()];

        int new_vertex_index = get_n();
        add_adjacency(new_vertex_index, vertex_index_1);
        add_adjacency(new_vertex_index, vertex_index_2);
        add_adjacency(new_vertex_index, vertex_index_3);

        edit_to_graph.push_back(GraphOperation{.operation_type = 1,
                                               .vertex_index_1 = new_vertex_index,
                                               .vertex_index_2 = vertex_index_1});
        edit_to_graph.push_back(GraphOperation{.operation_type = 1,
                                               .vertex_index_1 = new_vertex_index,
                                               .vertex_index_2 = vertex_index_2});
        edit_to_graph.push_back(GraphOperation{.operation_type = 1,
                                               .vertex_index_1 = new_vertex_index,
                                               .vertex_index_2 = vertex_index_3});

        history_stack.push(edit_to_graph);
        update_petals();
        return;
    }
    else if (where == 2)
    {
        remove_adjacency(vertex_index_1, vertex_index_2);

        int new_vertex_index = get_n();
        add_adjacency(new_vertex_index, vertex_index_1);
        add_adjacency(new_vertex_index, vertex_index_2);

        edit_to_graph.push_back(GraphOperation{.operation_type = 2,
                                               .vertex_index_1 = vertex_index_1,
                                               .vertex_index_2 = vertex_index_2});
        edit_to_graph.push_back(GraphOperation{.operation_type = 1,
                                               .vertex_index_1 = new_vertex_index,
                                               .vertex_index_2 = vertex_index_1});
        edit_to_graph.push_back(GraphOperation{.operation_type = 1,
                                               .vertex_index_1 = new_vertex_index,
                                               .vertex_index_2 = vertex_index_2});

        for (int i = 0; i < 2; i++)
        {
            std::vector<int> petals = get_petals(indices[1 - i]);
            int len = petals.size();
            std::vector<int>::iterator it = std::find(petals.begin(), petals.end(), indices[i]);
            int petal_index = std::distance(petals.begin(), it);
            // There the vertex is a boundary and the other vertex is its last petal, skip
            if (is_bound(indices[i]) && is_bound(indices[1-i]) && petal_index == len - 1)
                continue;
            int vertex_index_3 = petals[(petal_index + 1) % len];
            add_adjacency(new_vertex_index, vertex_index_3);

            edit_to_graph.push_back(GraphOperation{.operation_type = 1,
                                                   .vertex_index_1 = new_vertex_index,
                                                   .vertex_index_2 = vertex_index_3});
        }

        // If a vertex is added between two boundary vertices, set the new
        // vertex's label to be similar to the other two and update the boundaries
        if (is_bound(vertex_index_1) && is_bound(vertex_index_2))
        {
            set_label(new_vertex_index, 0.5 * (get_label(vertex_index_1) + get_label(vertex_index_2)));
            auto_set_bounds();
        }

        history_stack.push(edit_to_graph);
        update_petals();
        return;
    }
}

// -------------------------------------------------------------------
// Given an edge between vertices with specified indices, remove it
// and add an edge between their mutual petals.
// Note, if both vertices are boundary vertices then this isn't possible and
// nothing will happen.
// Another note, need to be careful, need to check that we don't cause a
// vertex to be doubly connected to another vertex, so before adding an
// adjacency, check if they're already adjacent.
void CirclePackGraph::flip_edge(int vertex_index_1, int vertex_index_2)
{
    // Check function is well-defined.
    if (!is_adjacent(vertex_index_1, vertex_index_2) ||
        (is_bound(vertex_index_1) && is_bound(vertex_index_2)))
        return;

    std::vector<int> petals = get_petals(vertex_index_1);
    int len = petals.size();
    std::vector<int>::iterator it = std::find(petals.begin(), petals.end(), vertex_index_2);
    int petal_index = std::distance(petals.begin(), it);
    // There the vertex is a boundary and the other vertex is its last petal, skip

    int vertex_index_3 = petals[(petal_index + 1) % len];
    int vertex_index_4 = petals[(petal_index + len - 1) % len];

    // Check that we're not causing a double connection that will
    // break the graph (invalid triangulation)
    if (is_adjacent(vertex_index_3, vertex_index_4))
        return;

    remove_adjacency(vertex_index_1, vertex_index_2);
    add_adjacency(vertex_index_3, vertex_index_4);

    // Store information about the edit so that it can be undone later.
    std::vector<GraphOperation> edit_to_graph = {GraphOperation{.operation_type = 2,
                                                                .vertex_index_1 = vertex_index_1,
                                                                .vertex_index_2 = vertex_index_2},
                                                 GraphOperation{.operation_type = 1,
                                                                .vertex_index_1 = vertex_index_3,
                                                                .vertex_index_2 = vertex_index_4}};
    history_stack.push(edit_to_graph);
    update_petals();
}

// -------------------------------------------------------------------
// Static function: given 3 mutually tangent circles with radii
// u, v, w, calculate the angle of the triangle made by the centers
// of the circles, at the center of the circle with radius u.
double CirclePackGraph::angle(double v, double u, double w)
{
    return 2 * asin(sqrt(u * w / (u + v) / (w + v)));
}

// -------------------------------------------------------------------
// Calculate the 'angle sum' of a vertex, given its index. It is
// computed by summing the angle it makes with each pair of mutually
// adjacent petals
double CirclePackGraph::theta(int vertex_index) const
{
    if (!exists(vertex_index))
        return 0;
    std::vector<int> petals = get_petals(vertex_index);
    int len = petals.size();
    double ret = 0;

    for (int i = 0; i < len - 1; i++)
        ret += angle(get_label(vertex_index), get_label(petals[i]), get_label(petals[i + 1]));
    if (!is_bound(vertex_index))
        ret += angle(get_label(vertex_index), get_label(petals[len - 1]), get_label(petals[0]));

    return ret;
}

// -------------------------------------------------------------------
// A recursive function to calculate and update the position of the
// petals of the vertex whose index is given, by updating the boolean vector
// which keeps track of positions fixed.
void CirclePackGraph::compute_positions(const int vertex_index, std::vector<bool> &fixed)
{
    // The petals of vertex 'vertex_index'
    std::vector<int> petals = get_petals(vertex_index);
    int len = petals.size();

    // Find a petal which has had its position fixed (we expect at least one)
    std::vector<int>::iterator it = std::find_if(petals.begin(), petals.end(),
                                                 [&fixed](int index)
                                                 { return fixed[index]; });

    // This shouldn't happen but just in case we don't find one.
    if (it == petals.end())
        return;

    // Get the index from the iterator.
    int fixed_petal_index = std::distance(petals.begin(), it);

    if (fixed_petal_index > 0)
    {
        for (int i = fixed_petal_index - 1; i > 0; i--)
        {
            int unfixed_petal = petals[i];
            if (fixed[unfixed_petal])
                continue;

            int fixed_petal = petals[i + 1];

            // Calculate the next position to be fixed.
            ofVec2f pos = get_pos(fixed_petal) - get_pos(vertex_index);
            pos.normalize();
            pos *= (get_label(vertex_index) + get_label(unfixed_petal));
            pos.rotateRad(-angle(get_label(vertex_index), get_label(unfixed_petal), get_label(fixed_petal)));
            pos += get_pos(vertex_index);

            // Set the new position and mark it fixed.
            set_pos(unfixed_petal, pos);
            fixed[unfixed_petal] = true;

            compute_positions(unfixed_petal, fixed);
        }
    }
    if (fixed_petal_index < len - 1)
    {
        for (int i = fixed_petal_index + 1; i < len; i++)
        {
            int unfixed_petal = petals[i];
            if (fixed[unfixed_petal])
                continue;

            int fixed_petal = petals[i - 1];

            // Calculate the next position to be fixed.
            ofVec2f pos = get_pos(fixed_petal) - get_pos(vertex_index);
            pos.normalize();
            pos *= (get_label(vertex_index) + get_label(unfixed_petal));
            pos.rotateRad(angle(get_label(vertex_index), get_label(unfixed_petal), get_label(fixed_petal)));
            pos += get_pos(vertex_index);

            // Set the new position and mark it fixed.
            set_pos(unfixed_petal, pos);
            fixed[unfixed_petal] = true;

            compute_positions(unfixed_petal, fixed);
        }
    }
}

// -------------------------------------------------------------------
// Calls the recenter_scale method with width and height params
// chosen based on the canvas' size
void CirclePackGraph::recenter_scale()
{
    float margin = 0.05 * min(ofGetWidth(), ofGetHeight());
    recenter_scale(ofGetWidth() - margin, ofGetHeight() - margin);
}

// -------------------------------------------------------------------
// Re-center and scale the circle packing to fit the given dimensions
void CirclePackGraph::recenter_scale(float width, float height)
{
    // Calculate the bounding box for the circle packing with
    // the current positions and labels
    ofVec2f bound_max = ofVec2f(-9999, -9999);
    ofVec2f bound_min = ofVec2f(9999, 9999);
    for (int i = 0; i < get_n(); i++)
    {
        if (!exists(i))
            continue;
        ofVec2f pos = get_pos(i);
        float label = get_label(i);
        bound_max.x = max(bound_max.x, pos.x + label);
        bound_max.y = max(bound_max.y, pos.y + label);
        bound_min.x = min(bound_min.x, pos.x - label);
        bound_min.y = min(bound_min.y, pos.y - label);
    }

    // Scale factor
    float scale_factor = min(width / (bound_max.x - bound_min.x), height / (bound_max.y - bound_min.y));

    // Update the positions and labels.
    for (int i = 0; i < get_n(); i++)
    {
        ofVec2f pos = get_pos(i);
        pos -= bound_min * 0.5;
        pos -= bound_max * 0.5;
        pos *= scale_factor;
        pos += view_translation;
        pos *= view_scale;

        set_pos(i, pos);
        set_label(i, get_label(i) * scale_factor * view_scale);
    }
}

// -------------------------------------------------------------------
// Search for a vertex whose circle intersects the given position,
// if nothing was found return a nullptr.
Vertex *CirclePackGraph::get_circle_at_position(ofVec2f pos) const
{
    for (int i = 0; i < get_n(); i++)
    {
        if (!exists(i))
            continue;
        ofVec2f vertex_pos = get_pos(i);
        float vertex_label = get_label(i);
        if ((vertex_pos - pos).length() < vertex_label)
            return vertices[i];
    }
    return nullptr;
}

// -------------------------------------------------------------------
// Return the current view_scale, default is 1.
float CirclePackGraph::get_view_scale() const
{
    return view_scale;
}
// -------------------------------------------------------------------
// Set the view_scale.
void CirclePackGraph::set_view_scale(float val)
{
    view_scale = val;
}
// -------------------------------------------------------------------
// Return the current view_translation, default is (0,0).
ofVec2f CirclePackGraph::get_view_translation() const
{
    return view_translation;
}
// -------------------------------------------------------------------
// Set the view_translation.
void CirclePackGraph::set_view_translation(ofVec2f val)
{
    view_translation = val;
}
// -------------------------------------------------------------------

// -------------------------------------------------------------------
// Undo the latest edit to the graph
void CirclePackGraph::undo()
{
    // Check there is an edit to undo
    if (history_stack.empty())
        return;

    std::cout << "Undo" << std::endl;

    std::vector<GraphOperation> latest_edit = history_stack.top();
    history_stack.pop();

    for (GraphOperation operation : latest_edit)
    {
        if (operation.operation_type == 1)
        {
            remove_adjacency(operation.vertex_index_1, operation.vertex_index_2);
        }
        else if (operation.operation_type == 2)
        {
            add_adjacency(operation.vertex_index_1, operation.vertex_index_2);
        }
    }

    auto_set_bounds();
    update_petals();
}
// -------------------------------------------------------------------