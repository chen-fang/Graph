#include <map>
#include <iostream>
#include "Bipartite.hpp"

typedef vector< __vertex_type > graph_type;
typedef vector< vector<index_t> > seed_type;

int main()
{
    graph_type Graph;
    seed_type Seed;

    std::size_t Nx = 4;
    std::size_t Ny = 4;
    std::size_t Nz = 3;
    
    Bipartite_Cartesian( Graph, Nx, Ny, Nz );
    
    Seeding::Init(Graph, Seed);
    Seeding::Coloring( Graph );

}
