#include "Graph.h"

// -------------------------------------------------------------------
// Static constructor, creating B rings of A vertices, each ring
// connected to the next, and the last ring "sewn shut" by connecting
// the vertices in the ring in a zig-zag fashion.
Graph Graph::example_graph_AxB(int A, int B)
{
    Graph graph;

    if (A < 3 || B < 1)
    {
        // if the parameters are invalid, return a W4 graph
        std::cout << "Please ensure A >= 3 and B >= 1" << std::endl;
        graph.add_adjacency(0, 1);
        graph.add_adjacency(1, 2);
        graph.add_adjacency(2, 3);
        graph.add_adjacency(3, 0);
        graph.add_adjacency(4, 0);
        graph.add_adjacency(4, 1);
        graph.add_adjacency(4, 2);
        graph.add_adjacency(4, 3);
    }
    else
    {
        int circum = A; // "circumference"
        int rings = B;

        for (int i = 0; i < circum; i++)
        {
            for (int j = 0; j < rings; j++)
            {
                // Connect vertices of a 'ring'.
                graph.add_adjacency((i + 0) % circum + circum * j, (i + 1) % circum + circum * j);
                if (j < rings - 1)
                {
                    // If it isn't the last ring, connect each vertex in the current
                    // ring to two vertices of the next ring.
                    graph.add_adjacency(i + circum * j, i + circum * (j + 1));
                    graph.add_adjacency(i + circum * j, (i + 1) % circum + circum * (j + 1));
                }
            }
        }

        // We're left with a hole if circum > 3.
        // "Sew" it shut by connecting the vertices of the last ring in a zig-zag fashion.
        if (circum > 3)
            for (int i = 0; i < circum - 3; i++)
                graph.add_adjacency(circum * (rings - 1) + 1 + i % 2 + i / 2,
                                    circum * (rings - 1) + circum - i / 2 - 1);
    }

    graph.auto_set_bounds();
    graph.print_bounds();
    graph.update_petals();

    return graph;
}

// -------------------------------------------------------------------
// Add a connection between a pair of vertices with given indices,
// that is, set their adjacency value to 'true'.
// If any of the indices are higher than the current total number of
// vertices of the graph (minus 1) then adjust the members 'n',
// 'adjacencies', 'vertices' accordingly.
// This is virtual so that the CirclePackGraph can turn any new pointers
// to Vertex objects into CirclePackVertex objects.
void Graph::add_adjacency(int vertex_index_1, int vertex_index_2)
{
    // std::cout << "Add adjacency : " << vertex_index_1 << "," << vertex_index_2 << std::endl;
    // Ensure no vertex index is negative
    assert(min(vertex_index_1, vertex_index_2) >= 0);

    // Handle whether we're adding an edge with a vertex index
    // higher than the current number of vertices.
    int old_n = n;
    int new_n = max(vertex_index_1, vertex_index_2) + 1;
    if (new_n > old_n)
    {
        for (int i = old_n; i < new_n; i++)
        {
            // Increase the size of the triangular matrix adjacencies
            std::vector<bool> new_col(i + 1, false);
            adjacencies.push_back(new_col);

            // Create new Vertex until the highest vertex index is reached
            vertices.push_back(new Vertex(i));
        }
        // Update the value of n accordingly.
        n = new_n;
    }

    adjacencies[max(vertex_index_1, vertex_index_2)][min(vertex_index_1, vertex_index_2)] = true;
}

// -------------------------------------------------------------------
// Remove a connection between a pair of vertices with given indices,
// that is, set their adjacency value to false.
void Graph::remove_adjacency(int vertex_index_1, int vertex_index_2)
{
    // std::cout << "Remove adjacency : " << vertex_index_1 << "," << vertex_index_2 << std::endl;
    assert(min(vertex_index_1, vertex_index_2) >= 0 && min(vertex_index_1, vertex_index_2) < n);
    adjacencies[max(vertex_index_1, vertex_index_2)][min(vertex_index_1, vertex_index_2)] = false;
}

// -------------------------------------------------------------------
// Automatically find which vertices are boundary vertices and set
// them as such. We decide that a vertex is boundary vertex if we
// can't find a closed path through its petals (neighbours).
// This only works if there are more than 3 boundary vertices.
void Graph::auto_set_bounds()
{
    update_petals();
    for (int i = 0; i < n; i++)
    {
        std::vector<int> petals = get_petals(i);
        set_bound(i, !path_exists(petals, 1, true) && exists(i));
    }
}

// -------------------------------------------------------------------
// Getter and setter functions:
// Each has an assertion that ensures a valid vertex index is given.
// -------------------------------------------------------------------
// Returns whether a vertex with given index is a boundary vertex
bool Graph::is_bound(int vertex_index) const
{
    assert(vertex_index >= 0 && vertex_index < n);
    return vertices[vertex_index]->is_bound();
}
// -------------------------------------------------------------------
// Sets the boundary status of a vertex with given index
void Graph::set_bound(int vertex_index, bool val)
{
    assert(vertex_index >= 0 && vertex_index < n);
    vertices[vertex_index]->set_bound(val);
    return;
}
// -------------------------------------------------------------------
// Gets the petals (neighbours) of a vertex with given index.
// They should be correctly ordered and oriented if 'update_petals'
// has been called after any modification to 'edges', assuming 'edges'
// describes a valid graph (that is a triangulation of a disc)
std::vector<int> Graph::get_petals(int vertex_index) const
{
    assert(vertex_index >= 0 && vertex_index < n);
    return vertices[vertex_index]->get_petals();
}
// -------------------------------------------------------------------
// Set the petals of a vertex with given index.
void Graph::set_petals(int vertex_index, const std::vector<int> &val)
{
    assert(vertex_index >= 0 && vertex_index < n);
    vertices[vertex_index]->set_petals(val);
    return;
}
// -------------------------------------------------------------------
// Get the petals of a vertex with given index.
Vertex *Graph::get_petal_pointer_of(int vertex_index, int petal_index) const
{
    assert(vertex_index >= 0 && vertex_index < n);
    std::vector<int> petals = get_petals(vertex_index);
    int len = petals.size();

    while (petal_index < 0)
        petal_index += len;
    petal_index %= len;

    return vertices[petals[petal_index]];
}
// // -------------------------------------------------------------------
// // Set the petals of a vertex with given index.
// Vertex *Graph::get_vertex_pointer(int vertex_index) const
// {
//     if (vertex_index < 0 || vertex_index >= n)
//         return nullptr;
//     return vertices[vertex_index];
// }
// -------------------------------------------------------------------

// -------------------------------------------------------------------
// Print indices of boundary vertices.
void Graph::print_bounds() const
{
    std::cout << "----------------boundary-vertices-----------------" << std::endl;
    for (int i = 0; i < n; i++)
        if (is_bound(i))
            std::cout << i << " ";
    std::cout << std::endl;
}

// -------------------------------------------------------------------
// Print the triangular matrix which described the adjacencies of
// vertices.
void Graph::print_adjacencies() const
{
    std::cout << "----------------edges-----------------" << std::endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < adjacencies[i].size(); j++)
            std::cout << (adjacencies[i][j] ? 1 : 0) << " ";
        std::cout << std::endl;
    }
}

// -------------------------------------------------------------------
// Print the petals of each vertex in the graph.
void Graph::print_petals() const
{
    std::cout << "----------------petals-----------------" << std::endl;
    for (const Vertex *v : vertices)
        std::cout << v->get_index() << ": " << v->get_petals_string() << std::endl;
}

// -------------------------------------------------------------------
// For each vertex in the graph, calculate and set its ordered petals.
// Then, make it so that orientation is consistent for the petals of
// each vertex. This is done by fixing the orientation of any one
// vertex, I chose vertex 0, and for each of its petals correct and
// fix their petals by calling the recursive function 'fix_orientation',
// keeping track of which vertices have already had their petals'
// orientation fixed with 'orientation_fixed' bool vector which is
// passed as a reference in the recursive function.
// Note; need to be extra careful when trying to find the ordered petals
// when a petal appears twice in the list. This means there's at least
// two ways to order the petals (up to reversing). To deal with this,
// we first calculate all petals without double connections, then
// do the rest, and check that no vertex is enclosed in the list of
// petals.
void Graph::update_petals()
{
    // Set each vertex's petals to be the computed ordered petals.
    for (int i = 0; i < n; i++)
    {
        set_petals(i, find_ordered_petals(i));
    }

    // Init a bool vector to keep track of indices of vertices whose
    // petals' orientation has been fixed in place.
    std::vector<bool> orientation_fixed(n, false);

    // Arbitrarily choose a first vertex.
    int first = 0;
    while (!exists(first))
        first++;
    orientation_fixed[first] = true;

    // Loop through each petal of the 'first' vertex (except for the
    // last petal if it is a boundary vertex).
    // Call fix_orientation for that petal, with 'dir' being the
    // 'first' vertex index and the index of the next consecutive petal
    std::vector<int> petals = get_petals(first);
    int len = petals.size();
    for (int i = 0; i < len - (is_bound(first) ? 1 : 0); i++)
    {
        int petal_idx = petals[i];
        std::array<int, 2> dir = {petals[(i + 1) % petals.size()], first};
        fix_orientation(petal_idx, dir, orientation_fixed);
    }

    print_petals();
}

// -------------------------------------------------------------------
// Returns whether the vertices with the two given indices are
// adjacent. There is an assertion to ensure that the indices are valid.
bool Graph::is_adjacent(int vertex_index_1, int vertex_index_2) const
{
    assert(min(vertex_index_1, vertex_index_2) >= 0 && max(vertex_index_1, vertex_index_2) > n);
    return adjacencies[max(vertex_index_1, vertex_index_2)][min(vertex_index_1, vertex_index_2)] > 0;
}

// -------------------------------------------------------------------
// Returns the graph distance between two vertices, given their
// indices.
int Graph::dist(int vertex_index_1, int vertex_index_2) const
{
    // Set distance to -1 to mark as unexplored
    std::vector<int> distance(n, -1);

    // distance from vertex vertex_index_1 to vertex vertex_index_1 is 0
    distance[vertex_index_1] = 0;
    // queue to mark the distance of a vertex's neighbours
    std::vector<int> q(1, vertex_index_1);

    while (q.size() > 0 && distance[vertex_index_2] < 0)
    {
        // pop front (get and remove first element of the vector)
        int current = q[0];
        q.erase(q.begin());
        for (int i : get_petals(current))
        {
            if (distance[i] < 0)
            {
                distance[i] = distance[current] + 1;
                q.push_back(i);
            }
        }
    }
    return distance[vertex_index_2];
}

// -------------------------------------------------------------------
// Find the unordered petal of a vertex with given index by using
// a simple loop that checks for adjacencies using 'is_edge'
std::vector<int> Graph::find_unordered_petals(int vertex_index) const
{
    std::vector<int> ret;
    for (int i = 0; i < n; i++)
        if (is_adjacent(vertex_index, i))
            ret.push_back(i);
    return ret;
}

// -------------------------------------------------------------------
// Find the ordered petals of a vertex with given index by first
// getting the unordered petals and then finding a path through them.
//
// Should the vertex be a boundary vertex, then its list of petals
// won't be closed (circular); it will have two petals that are boundary
// vertices, each of which will only have one neighbour in common with
// the vertex in question, hence, they should be the start and end of the
// list of ordered petals. So before finding a path through the list of
// unordered petals, we set the first petal to be one of the boundary
// vertex.
// Note: the returned petals don't guarantee orientation consistent
// with the rest of the vertices, which is why the orientation is
// fixed in another method.
std::vector<int> Graph::find_ordered_petals(int vertex_index) const
{
    std::vector<int> ret = find_unordered_petals(vertex_index);
    int n_petals = ret.size();

    if (n_petals > 0)
    {
        // If the vertex is boundary and its first (unordered) petal
        // isn't a boundary, loop and swap until it is.
        // std::cout <<vertex_index<<", "<< ret[0] << ", " << is_bound(vertex_index) << ", " << !is_bound(ret[0])  << std::endl;
        if (is_bound(vertex_index) && !is_bound(ret[0]))
        {
            for (int i = 1; i < n_petals; i++)
            {
                if (is_bound(ret[i]))
                {
                    std::swap(ret[0], ret[i]);
                    break;
                }
            }
        }

        // Search for a path through petals, search for a closed path
        // if it isn't a boundary vertex.
        path_exists(ret, 1, !is_bound(vertex_index));
    }

    return ret;
}

// -------------------------------------------------------------------
// Find a path between a given list of vertex indices
// This is done by performing recursive backtracking.
// The 'closed' parameter is for finding a circular path through these vertices.
bool Graph::path_exists(std::vector<int> &vertex_indices, int current_index, bool closed) const
{
    int len = vertex_indices.size();
    if (current_index == len)
    {
        if (!closed)
            return is_bound(vertex_indices.back());
        else
            return is_adjacent(vertex_indices.front(), vertex_indices.back());
    }

    for (int i = current_index; i < len; i++)
    {
        if (!is_adjacent(vertex_indices[current_index - 1], vertex_indices[i]))
            continue;

        std::swap(vertex_indices[current_index], vertex_indices[i]);
        if (path_exists(vertex_indices, current_index + 1, closed))
            return true;
        std::swap(vertex_indices[current_index], vertex_indices[i]);
    }
    return false;
}

// -------------------------------------------------------------------
// Fix the orientation of the petals of a vertex with given index.
// The orientation is determined by the 'dir' parameter which stands
// for 'direction' and represents two adjacent petals of the vertex in
// question. It signifies that the triangle made by the vertices with
// indices [vertex_index, dir[0], dir[1]] should be, say, clockwise.
// We check whether this is so in the current (ordered) petals,
// and if not, we simply reverse the petals. Then mark the petals'
// orientation as fixed and recursively call the function for those
// petals that haven't yet had their own petals' orientation fixed.
void Graph::fix_orientation(int vertex_index, const std::array<int, 2> &dir, std::vector<bool> &fixed)
{
    std::vector<int> petals = get_petals(vertex_index);
    int len = petals.size();
    for (int i = 0; i < len; i++)
        if (petals[i] == dir[0])
        {
            if (petals[(i + 1) % len] != dir[1])
            {
                std::reverse(petals.begin(), petals.end());
                set_petals(vertex_index, petals);
            }
            break;
        }

    fixed[vertex_index] = true;

    for (int i = 0; i < len - (is_bound(vertex_index) ? 1 : 0); i++)
        if (!fixed[petals[i]] && exists(petals[i]))
            fix_orientation(petals[i], std::array<int, 2>{petals[(i + 1) % len], vertex_index}, fixed);
}

// -------------------------------------------------------------------
// Check whether a vertex 'exists' in the graph. We say that it it's
// non-existant if it is isolated (has no neighbours), and use this to
// exlucde it from relevant calculations.
bool Graph::exists(int vertex_index) const
{
    if (vertex_index < 0 || vertex_index > n - 1)
        return false;
    return get_petals(vertex_index).size() > 0;
}
