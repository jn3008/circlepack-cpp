#pragma once
// #include "ofMain.h"
#include "Vertex.h"

// -------------------------------------------------------------------
// CirclePackVertex class which inherits from Vertex.
// It serves to keep graph information about a circle in a 
// circle packing from the Vertex base class whilst also adding
// extra necessary information about the circle in the packing, namely,
// its label (radius), position (2D), and 'angle sum target' which is 
// 2*PI for interior vertices for our purposes. The sum of the 'angle
// sum target' for the m boundary vertices of a circle packing graph
// should be (m-2)*PI (like polygons: triangle:PI, square 2*PI). 
class CirclePackVertex : public Vertex
{
private:
    ofVec2f pos;
    double label, target;

public:

    CirclePackVertex(const Vertex &vertex) : Vertex(vertex)
    {
        pos = ofVec2f();
        label = 1;
        target = PI * 2;
    }

    // ------------------------------------------------------
    // Clone function for pointer, overrides the base class version
    CirclePackVertex* clone() const override
    {
        return new CirclePackVertex(*this);
    }

    // ------------------------------------------------------
    // Set the value of position
    void set_pos(const ofVec2f &val)
    {
        pos = val;
    }

    // ------------------------------------------------------
    // Return the vertex's (copied) position value.
    ofVec2f get_pos() const
    {
        return pos;
    }

    // ------------------------------------------------------
    // Set the value of label (radius).
    void set_label(const double &val)
    {
        label = val;
    }
    
    // ------------------------------------------------------
    // Get vertex's (copied) label value.
    double get_label() const
    {
        return label;
    }

    // ------------------------------------------------------
    // Set the value of target_ang (angle sum target)
    void set_target_ang(const double &val)
    {
        target = val;
    }
    
    // ------------------------------------------------------
    // Get vertex's (copied) target_ang value.
    double get_target_ang() const
    {
        return target;
    }
};