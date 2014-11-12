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
   std::size_t Ny = 5;
   std::size_t Nz = 5;
   std::size_t Neqn = 2;

   Seeding S( Nx, Ny, Nz, Neqn );
   S.Print();


   std::size_t N = Nx * Ny * Nz * Neqn;
   ADvector vec1, vec2;
   vec1.resize( N );
   vec2.resize( N );

   for( std::size_t i = 0; i < N; ++i )
   {
      vec1[i].make_independent( S.Color(i) );
   }


   // [ 5 1 1 2 ]
   vec2[0] = vec1[0] + vec1[1] + vec1[2] + vec1[3];
   vec2[1] = vec1[0] + vec1[1] + vec1[2] + vec1[3];
   vec2[2] = vec1[0] + vec1[1] + vec1[2] + vec1[3] + vec1[4] + vec1[5];
   vec2[3] = vec1[0] + vec1[1] + vec1[2] + vec1[3] + vec1[4] + vec1[5];
   vec2[4] = vec1[2] + vec1[3] + vec1[4] + vec1[5] + vec1[6] + vec1[7];
   vec2[5] = vec1[2] + vec1[3] + vec1[4] + vec1[5] + vec1[6] + vec1[7];
   vec2[6] = vec1[4] + vec1[5] + vec1[6] + vec1[7] + vec1[8] + vec1[9];
   vec2[7] = vec1[4] + vec1[5] + vec1[6] + vec1[7] + vec1[8] + vec1[9];
   vec2[8] = vec1[6] + vec1[7] + vec1[8] + vec1[9];
   vec2[9] = vec1[6] + vec1[7] + vec1[8] + vec1[9];


   // // [ 5 1 1 ] Case #1
   // vec2[0] = vec1[0] + vec1[1];
   // vec2[1] = vec1[0] + vec1[1] + vec1[2];
   // vec2[2] = vec1[1] + vec1[2] + vec1[3];
   // vec2[3] = vec1[2] + vec1[3] + vec1[4];
   // vec2[4] = vec1[3] + vec1[4];

   // // [ 5 1 1 ] Case #2
   // vec2[0] = vec1[0];
   // vec2[1] = vec1[0] + vec1[1];
   // vec2[2] = vec1[1] + vec1[2];
   // vec2[3] = vec1[2] + vec1[3];
   // vec2[4] = vec1[3] + vec1[4];

   // [ 5 2 1 ]
   // vec2[0] = vec1[0] + vec1[1] + vec1[5];
   // vec2[1] = vec1[0] + vec1[1] + vec1[2] + vec1[6];
   // vec2[2] = vec1[1] + vec1[2] + vec1[3] + vec1[7];
   // vec2[3] = vec1[2] + vec1[3] + vec1[4] + vec1[8];
   // vec2[4] = vec1[3] + vec1[4] + vec1[9];
   // // continue
   // vec2[5] = vec1[0] + vec1[5] + vec1[6];
   // vec2[6] = vec1[1] + vec1[5] + vec1[6] + vec1[7];
   // vec2[7] = vec1[2] + vec1[6] + vec1[7] + vec1[8];
   // vec2[8] = vec1[3] + vec1[7] + vec1[8] + vec1[9];
   // vec2[9] = vec1[4] + vec1[8] + vec1[9];


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
   vector<double> residual;
   S.Recover_CSR( vec2, residual, CSR );
   
   Print_Vector( CSR.rowptr() );
   Print_Vector( CSR.colind() );
   Print_Vector( CSR.value() );
   Print_Vector( residual );

   S.Recover_CSR( vec2, residual, CSR );
}
