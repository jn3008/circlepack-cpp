#pragma once
#include "ofMain.h"

class Vertex
{
private:
    int index;
    bool is_bound_val;
    std::vector<int> petals;

public:
    Vertex(int index_) : index(index_)
    {
        is_bound_val = false;
    }

    // Clone function that needs to be overridden by the
    // derived CirclePackVertex class
    virtual Vertex* clone() const 
    {
        return new Vertex(*this);
    }

    // Destructor needs to be virtual to properly delete
    // CirlcePackVertex objects which have Vertex* pointers
    virtual ~Vertex() = default;

    int get_index() const
    {
        return index;
    }

    // ------------------------------------------------------
    // Set the status of a vertex as boundary or interior.
    void set_bound(const bool &val)
    {
        is_bound_val = val;
    }
    // ------------------------------------------------------
    // Return the vertex's (copied) boundary status.
    bool is_bound() const
    {
        return is_bound_val;
    }
    // ------------------------------------------------------
    // Set the (indices of) petals of the vertex.
    void set_petals(const std::vector<int> &val)
    {
        petals = val;
    }
    // ------------------------------------------------------
    // Return the (copied) vector of petal indices.
    std::vector<int> get_petals() const
    {
        return petals;
    }
    
    std::string get_petals_string() const {
        std::string ret = "";
        int len = petals.size();
        for (int i = 0; i < len; i++)
        {
            ret += std::to_string(petals[i]);
            ret += i < len-1 ? "-" : (is_bound_val ? " " : "-");
        }
        return ret;
    }
};