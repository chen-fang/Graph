#pragma once

#include <map>
#include <vector>

using std::size_t;
using std::multimap;
using std::make_pair;
using std::vector;

typedef int color_t;
typedef std::size_t index_t;
typedef multimap<index_t,index_t> Map;
typedef Map::iterator Iter;
typedef Map::const_iterator const_Iter;
typedef Map::size_type size_type;


void Bipartite_Cartesian ( Map& Graph, size_t _Nx, size_t _Ny, size_t _Nz, size_t& Num_Color)
{
   index_t Vertex1, Vertex2;

   // 1st direction in all layers
   for (index_t k=0; k < _Nz; ++k)
   {
      for (index_t j=0; j < _Ny; ++j)
      {
	 for (index_t i=0; i < _Nx-1; ++i)
	 {
	    Vertex1 = i + j*_Nx + k*_Nx*_Ny;
	    Vertex2 = Vertex1 + 1;
	    Graph.insert( make_pair(Vertex1,Vertex2) );
	    Graph.insert( make_pair(Vertex2,Vertex1) );
	 }
      }
   }
   // 2nd direction in all layers
   for (index_t k=0; k < _Nz; ++k)
   {
      for (index_t j=0; j < _Ny-1; ++j)
      {
	 for (index_t i=0; i < _Nx; ++i)
	 {
	    Vertex1 = i + j*_Nx + k*_Nx*_Ny;
	    Vertex2 = Vertex1 + _Nx;
	    Graph.insert( make_pair(Vertex1,Vertex2) );
	    Graph.insert( make_pair(Vertex2,Vertex1) );
	 }
      }
   }
   // 3rd direction
   for (index_t k=0; k < _Nz-1; ++k)
   {
      for (index_t j=0; j < _Ny; ++j)
      {
	 for (index_t i=0; i < _Nx; ++i)
	 {
	    Vertex1 = i + j*_Nx + k*_Nx*_Ny;
	    Vertex2 = Vertex1 + _Nx * _Ny;
	    Graph.insert( make_pair(Vertex1,Vertex2) );
	    Graph.insert( make_pair(Vertex2,Vertex1) );
	 }
      }
   }
   // Itself
   size_t N = _Nx * _Ny * _Nz;
   for (index_t i=0; i < N; ++i)
   {
      Graph.insert( make_pair(i,i) );
   }

   // Upper bound for the number of colors
   size_t Max_Degree = 0;
   for( index_t i=0; i < N; ++i )
   {
      size_t count = Graph.count(i);
      Max_Degree = Max_Degree >= count ? Max_Degree : count;
   }
   ++Max_Degree;
   Num_Color = Max_Degree <= N ? Max_Degree : N;
}




struct Seeding
{
   static void Init( const Map& Graph, Map& Seed, vector<color_t>& Color );
   static void Coloring( const Map& Graph, Map& Seed, vector<color_t>& Color, size_t NumColor );
   static void Get_Seed_Matrix( const Map& Graph, Map& Seed, size_t NumVar, size_t NumColor );
};


void Seeding::Init( const Map& Graph, Map& Seed, vector<color_t>& Color )
{
   color_t color(0);
   std::pair<const_Iter,const_Iter> p = Graph.equal_range( Graph.begin()->first );
   while( p.first != p.second )
   {
      index_t index = p.first->second;
      Seed.insert( make_pair(color,index) );
      Color[index] = color;
      ++p.first;
      ++color;
   }
   // verify
   for( Iter i=Seed.begin(); i != Seed.end(); ++i )
   {
      std::cout << i->first <<"\t"<< i->second << std::endl;
   }
}

void Seeding::Coloring( const Map& Graph, Map& Seed, vector<color_t>& Color, size_t NumColor )
{
   vector<bool> Palette;
   Palette.reserve( NumColor );
   Palette.assign( NumColor, true );

   index_t i;
   size_t Num_Var = Color.size();
   for( i=0; i<Num_Var; ++i )
   {
      std::cout << "====================================" << std::endl;
      if( Color[i] != -1 )
      {
	 std::pair<const_Iter,const_Iter> j = Graph.equal_range( i );
	 while( j.first != j.second )
	 {
	    index_t index = j.first->second;
	    color_t color = Color[index];
	    if( color != -1 )
	    {
	       Palette[color] = false;
	    }
	    ++j.first;
	 }
	 std::cout << "palette------" << std::endl;
	 for(color_t aa=0; aa<NumColor; ++aa)
	 {
	    std::cout << Palette[aa] << std::endl;
	 }
	 
	 // choose the minimum color
	 color_t color(0);
	 for(index_t k=0; k<NumColor; ++k)
	 {
	    if( Palette[k] != false )
	    {
	       color = k;
	       std::cout << "Selected:\t" << color << std::endl;
	       break;
	    }
	 }
	 std::cout << "end of palette------" << std::endl;
	 Seed.insert( make_pair(color,i) );
	 Color[i] = color;
	 // reset Palette
	 Color.assign( NumColor, true );
      }//end_of_if_statement
      std::cout << "====================================" << std::endl;
   }
}



void Seeding::Get_Seed_Matrix( const Map& Graph, Map& Seed, size_t NumVar, size_t NumColor)
{
   vector<color_t> Index_Color;
   Index_Color.reserve( NumVar );
   Index_Color.assign( NumVar, -1 );

   Init( Graph, Seed, Index_Color );
   Coloring( Graph, Seed, Index_Color, NumColor );
}
