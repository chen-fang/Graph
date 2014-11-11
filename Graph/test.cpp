#include <vector>
#include <iostream>
#include "Graph.hpp"

//#include "/home/blackicefc/ADETL_DISTRO/adetl/scalars/ADscalar.hpp"
//#include "/home/blackicefc/ADETL_DISTRO/adetl/systems/ADvector.hpp"

#include "/home/chenfang/Library/AD_Library/ADETL_DISTRO/adetl/scalars/ADscalar.hpp"
#include "/home/chenfang/Library/AD_Library/ADETL_DISTRO/adetl/systems/ADvector.hpp"


typedef adetl::ADscalar<> ADscalar;
typedef adetl::ADvector ADvector;
typedef GENSOL::CSR_Matrix<double,int> CSR_Matrix;

using std::vector;

int main()
{
   std::size_t Nx = 5;
   std::size_t Ny = 1;
   std::size_t Nz = 1;

   Seeding S;
    
   S.Build_BipartiteGraph( Nx, Ny, Nz );
   S.Coloring();
   S.Print();


   std::size_t N = Nx * Ny * Nz;
   ADvector vec1, vec2;
   vec1.resize( N );
   vec2.resize( N );

   for( std::size_t i = 0; i < N; ++i )
   {
      vec1[i].make_independent( S.Color(i) );
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

   Print_Vector( resi );
   Print_Vector( row );
   Print_Vector( col );
   Print_Vector( data );


   std::cout << "--------------- After Recovery -------------------" << std::endl;
   CSR_Matrix CSR;
   S.Recovery( vec2, CSR );
   
   Print_Vector( CSR.rowptr() );
   // Print_Vector( CSR.colind() );
   // Print_Vector( CSR.value() );
}
