#pragma once
#include "CirclePackGraph.h"

// CirclePacking class
// Used as an interface for CirclePackGraph class
// Handles creating the graph, updating the packing by computing the
// updated labels and positions, mouse and keyboard interactions.
class CirclePacking
{
private:
    CirclePackGraph graph;

    Vertex *hovered_vertex = nullptr,
           *selected_vertex = nullptr,
           *second_selected_vertex = nullptr;

    ofVec2f start_drag_position;

    bool theme = false, // false: dark, true: light
        mode = false;   // false: fix boundary labels, true: fix boundary angles

    std::array<std::vector<std::string>, 2> instructions{
        std::vector<std::string>{
            "BACKSPACE : undo",
            "RETURN : mode",
            "T : theme",
            "V : reset view",
            "R : remove vertex",
            "UP : edge to vertex",
            "DOWN : flip edge",
            "LEFT/RIGHT : face to vertex"},
        std::vector<std::string>{
            "BACKSPACE : undo",
            "RETURN : mode",
            "T : theme",
            "V : reset view",
            "UP/DOWN : increase/decrease target angle",
            "LEFT/RIGHT : navigate boundary vertices"}};

public:
    CirclePacking() {}

    CirclePacking(const Graph &graph_) : graph(graph_)
    {
        update_packing(300);

        graph.print_bounds();
        graph.print_adjacencies();
        graph.print_labels();
        graph.print_petals();
    }

    ofTrueTypeFont font_mono, font_mono_italic, font_mono_bold;

    void set_fonts(std::string address_mono,
                   std::string address_mono_italic,
                   std::string address_mono_bold, int size)
    {
        font_mono.load(address_mono, size);
        font_mono_italic.load(address_mono_italic, size);
        font_mono_bold.load(address_mono_bold, size);
    }

    void show(std::string last_key_pressed) const
    {
        ofTranslate(ofGetWidth() * 0.5, ofGetHeight() * 0.5);
        show_packing();

        ofTranslate(-ofGetWidth() * 0.5, -ofGetHeight() * 0.5);
        show_info(last_key_pressed);
    }

    void set_hovered(Vertex *v) { hovered_vertex = v; }
    void set_selected(Vertex *v) { selected_vertex = v; }
    void set_second_selected(Vertex *v) { second_selected_vertex = v; }

    bool get_theme() const { return theme; }

    // -------------------------------------------------------------------
    // Draw the circles of the packing and the edges of the underlying
    // graph
    void show_packing() const
    {

        for (int i = 0; i < graph.get_n(); i++)
        {
            if (!graph.exists(i))
                continue;
            ofVec2f pos = graph.get_pos(i);

            // Draw circles
            ofNoFill();
            ofSetColor(theme ? 0 : 255);
            ofSetLineWidth(theme ? 2 : 1);
            // Thicken the circle if it's selected
            if (selected_vertex)
                if (i == selected_vertex->get_index())
                    ofSetLineWidth(theme ? 5 : 4);
            if (second_selected_vertex)
                if (i == second_selected_vertex->get_index())
                    ofSetLineWidth(theme ? 5 : 4);

            // Darken the circle if it's hovered over
            if (hovered_vertex)
                if (i == hovered_vertex->get_index())
                    ofSetColor(125);

            ofDrawCircle(pos.x, pos.y, graph.get_label(i));
        }

        // Draw edges
        for (int i = 0; i < graph.get_n(); i++)
        {
            if (!graph.exists(i))
                continue;

            ofVec2f pos = graph.get_pos(i);
            ofSetLineWidth(theme ? 2 : 1);
            for (int j = i; j < graph.get_n(); j++)
            {
                if (graph.is_adjacent(i, j))
                {
                    ofSetColor(125);
                    ofDrawLine(pos, graph.get_pos(j));
                }
            }
        }
    }

    void show_info(std::string last_key_pressed) const
    {
        float height_offset = 0;
        float margin = 10;
        ofFill();

        // Draw text info for selected vertex
        if (selected_vertex)
        {
            ofTrueTypeFont font = hovered_vertex == selected_vertex ? font_mono_bold : font_mono;
            show_vertex_info(font, selected_vertex, height_offset, margin);
            height_offset += margin;
        }
        // Draw text info for second selected vertex
        if (second_selected_vertex)
        {
            ofTrueTypeFont font = hovered_vertex == second_selected_vertex ? font_mono_bold : font_mono;
            show_vertex_info(font, second_selected_vertex, height_offset, margin);
            height_offset += margin;
        }
        // Draw text info for hovered-over vertex, (if it's not already selected)
        if (hovered_vertex &&
            hovered_vertex != selected_vertex &&
            hovered_vertex != second_selected_vertex)
        {
            ofTrueTypeFont font = font_mono_italic;
            show_vertex_info(font, hovered_vertex, height_offset, margin);
        }

        // Draw text info for error of the packing
        std::string error = std::to_string(graph.get_error());
        // Add left-padding to the value
        error = "Error: " + std::string(max(0, 8 - static_cast<int>(error.length())), ' ') + error;
        font_mono.drawString(error, ofGetWidth() - margin - font_mono.stringWidth(error),
                             margin + font_mono.stringHeight(error));

        // Draw as text the last pressed key
        font_mono.drawString(last_key_pressed,
                             ofGetWidth() - margin - font_mono.stringWidth(last_key_pressed),
                             ofGetHeight() - margin);

        // Draw instructions; the set of instruction depends on the mode
        height_offset = 0;
        ofSetColor(128);
        for (std::string str : instructions[mode])
        {
            font_mono_italic.drawString(str, margin, ofGetHeight() - margin - height_offset);
            height_offset += font_mono_italic.getAscenderHeight();
        }
    }

    void show_vertex_info(const ofTrueTypeFont font, const Vertex *v, float &height_offset, const float margin) const
    {
        int vertex_index = v->get_index();

        std::string str = "Index:  " + std::to_string(vertex_index);
        height_offset += font.getAscenderHeight();
        font.drawString(str, margin, height_offset);

        str = "Petals: " + v->get_petals_string();
        height_offset += font.getAscenderHeight();
        font.drawString(str, margin, height_offset);

        if (!mode)
            return;
        str = "Target: " + std::to_string(graph.get_target_ang(vertex_index)) + "rad";
        height_offset += font.getAscenderHeight();
        font.drawString(str, margin, height_offset);
    }

    void update_hovered(ofVec2f mouse)
    {
        set_hovered(graph.get_circle_at_position(mouse));
    }
    void update_selected(int key)
    {
        // If we are in 'fix boundary angles' mode then we only want to be able to
        // select boundary vertices (or nothing, that is, deselect).
        if (mode && hovered_vertex && !graph.is_bound(hovered_vertex->get_index()))
            return;

        // Left-click, or nothing selected yet
        if (key == 0 || !selected_vertex)
        {
            // Select the first vertex, deselect the second
            set_selected(hovered_vertex);
            set_second_selected(nullptr);
        }
        // Right-click and (fix boundary labels mode)
        else if (key == 2 && !mode)
        {
            set_second_selected(hovered_vertex); // Select the second vertex
        }
    }

    void update_start_drag(ofVec2f mouse)
    {
        start_drag_position = mouse - graph.get_view_translation() * graph.get_view_scale();
    }
    void update_view_translation(ofVec2f mouse)
    {
        graph.set_view_translation((mouse - start_drag_position) / graph.get_view_scale());
        graph.recenter_scale();
    }
    void update_view_scale(float scrollY)
    {
        graph.set_view_scale(pow(2, log2(graph.get_view_scale()) + scrollY * 0.2));
        graph.recenter_scale();
    }

    void update_packing(int num_iterations)
    {
        // ------------------------------------------------------------
        for (int i = 0; i < num_iterations; i++)
        {
            if (mode)
                graph.compute_labels_by_fixing_boundary_targets();
            else
                graph.compute_labels_by_fixing_boundary_labels();
        }
        graph.compute_positions();
        // ------------------------------------------------------------
    }

    void key_pressed(int key)
    {
        // If 'r' is pressed and a vertex is selected, remove it.
        if (key == 'r' && selected_vertex && !mode)
        {
            std::vector<GraphOperation> edit_to_graph;
            graph.remove_and_fill(selected_vertex->get_index(), edit_to_graph);
            update_packing(200);
            set_selected(nullptr);
        }
        // If left or right arrow is pressed
        else if (key == OF_KEY_LEFT || key == OF_KEY_RIGHT)
        {
            if (mode) // 'Fix boundary target angles mode'
            {
                if (selected_vertex)
                    set_selected(graph.get_petal_pointer_of(
                        selected_vertex->get_index(),
                        key == OF_KEY_LEFT ? 0 : -1));
            }
            else if (second_selected_vertex) // 'Fix boundary target labels mode'
            {
                // If two vertices are selected, add a vertex in the face/triangle 'left'/'right'
                // of the edge made by those two vertices
                graph.add_vertex(selected_vertex->get_index(),
                                 second_selected_vertex->get_index(),
                                 key == OF_KEY_LEFT ? 0 : 1);
            }
            update_packing(200);
        }
        else if (key == OF_KEY_UP || key == OF_KEY_DOWN)
        {
            if (mode) // 'Fix boundary target angles mode'
            {
                if (selected_vertex)
                    graph.adjust_target_ang(selected_vertex->get_index(),
                                            (key == OF_KEY_UP ? 1 : -1) * 0.02);
            }
            else if (second_selected_vertex) // 'Fix boundary target labels mode'
            {
                if (key == OF_KEY_UP)
                {
                    // If up arrow pressed and two vertices are selected,
                    // add a vertex between those two vertices.
                    graph.add_vertex(selected_vertex->get_index(), second_selected_vertex->get_index(), 2);
                }
                else
                {
                    // If down arrow pressed and two vertices are selected,
                    // flip the edge to connect their two mutual petals instead.
                    // Note, at most one vertex can be boundary.
                    graph.flip_edge(selected_vertex->get_index(), second_selected_vertex->get_index());
                }
            }
            update_packing(200);
        }
        // If 'v' is pressed, reset the view.
        else if (key == 'v')
        {
            graph.set_view_scale(1);
            graph.set_view_translation(ofVec2f());
            graph.recenter_scale();
            update_start_drag(ofVec2f());
        }
        // If backspace is pressed, undo the latest edit to the graph
        else if (key == OF_KEY_BACKSPACE)
        {
            graph.undo();
            update_packing(200);
        }
        // If return/enter is pressed, toggle mode
        else if (key == OF_KEY_RETURN)
        {
            mode = !mode;
            // Set the angle sum targets of the boundary vertices
            // to be their current angle sums.
            graph.auto_set_targets();
            // Deselect both selections
            set_selected(nullptr);
            set_second_selected(nullptr);
        }
        // If backspace is pressed, toggle theme
        else if (key == 't')
        {
            theme = !theme;
        }
        // Otherwise, performs computations to reduce error of packing
        else
        {
            update_packing(50);
        }
    }
};