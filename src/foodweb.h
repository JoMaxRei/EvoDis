#ifndef FOODWEB_H_
#define FOODWEB_H_

#include <stdint.h>
#include <vector>
#include "species.h"

// AKA Patch
//
// Every foodweb has species  that exist in it
// Every foodweb calculates how species interact
class Foodweb
{
public:
    // Creates a new foodweb
    // @param resource adds a non-extinguishable species to the foodweb
    // @param id the id of this new foodweb
    Foodweb(Species *resource, size_t x, size_t y);

    // Adds a species to the foodweb
    //@return -1 if maximum species in this web was exceeded, position of the species otherwise
    int64_t add_species(Species *species);
    // Sum of all species dispersal rates on this foodweb
    int64_t m_local_dispersal_rate;

    // x position of the this foodweb
    size_t m_x;
    // y position of the this foodweb
    size_t m_y;

    // maximum species that can live in a foodweb
    static const size_t MAX_DIM = 100;

private:
    Species **m_species;
    size_t m_species_count;
    double *m_fitness;
};

#endif