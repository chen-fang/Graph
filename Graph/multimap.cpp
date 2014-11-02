#include <map>
#include <iostream>
#include "Bipartite.hpp"


int main()
{
   multimap<std::size_t, std::size_t> Seed;
   multimap<std::size_t,std::size_t> Graph;

   std::size_t Nx = 2;
   std::size_t Ny = 2;
   std::size_t Nz = 2;
   std::size_t NumVar = Nx * Ny * Nz;
   std::size_t NumColor;
   Bipartite_Cartesian( Graph, Nx, Ny, Nz, NumColor );
   std::cout << NumColor << std::endl;


   Seeding::Get_Seed_Matrix( Graph, Seed, NumVar, NumColor );

}
