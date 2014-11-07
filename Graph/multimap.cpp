#include <vector>
#include <iostream>
#include "Bipartite.hpp"

#include "/home/blackicefc/ADETL_DISTRO/adetl/scalars/ADscalar.hpp"
#include "/home/blackicefc/ADETL_DISTRO/adetl/systems/ADvector.hpp"

typedef vector< __vertex_type > graph_type;
typedef vector< vector<index_t> > seed_type;

typedef adetl::ADscalar<> ADscalar;
typedef adetl::ADvector ADvector;

using std::vector;

int main()
{
   graph_type Graph;
   seed_type Seed;

   std::size_t Nx = 5;
   std::size_t Ny = 1;
   std::size_t Nz = 1;
    
   Bipartite_Cartesian( Graph, Nx, Ny, Nz );
   Sort( Graph );
    
   Seeding::Init(Graph, Seed);
   Seeding::Coloring( Graph );
   Print( Graph );



   std::size_t N = 5;
   ADvector vec1, vec2;
   vec1.resize( N );
   vec2.resize( N );

   for( std::size_t i = 0; i < N; ++i )
   {
      // vec1[i].make_independent( i );
      vec1[i].make_independent( Graph[i].color );
   }

   vec2[0] = vec1[0] + vec1[1];
   vec2[1] = vec1[0] + vec1[1] + vec1[2];
   vec2[2] = vec1[1] + vec1[2] + vec1[3];
   vec2[3] = vec1[2] + vec1[3] + vec1[4];
   vec2[4] = vec1[3] + vec1[4];

   for( std::size_t i = 0; i < N; ++i )
   {
      std::cout << vec2[i] << std::endl;
   }

   vector<double> resi;
   vector<int> row;
   vector<int> col;
   vector<double> data;

   vec2.extract_CSR( resi, row, col, data );

   for( std::size_t i = 0; i < resi.size(); ++i ) 
   {
      std::cout << resi[i] << std::endl;
   }
   std::cout << std::endl;

   for( std::size_t i = 0; i < row.size(); ++i ) 
   {
      std::cout << row[i] << std::endl;
   }
   std::cout << std::endl;

   for( std::size_t i = 0; i < col.size(); ++i ) 
   {
      std::cout << col[i] << std::endl;
   }
   std::cout << std::endl;

   for( std::size_t i = 0; i < data.size(); ++i ) 
   {
      std::cout << data[i] << std::endl;
   }



}
