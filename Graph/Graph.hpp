#pragma once

#include <vector>
#include <algorithm>
#include "CSR_Matrix.hpp"

using std::size_t;
using std::vector;

typedef int           color_t;
typedef size_t        index_t;

class Seeding
{
   struct __vertex_type
   {
      color_t color;
      vector<index_t> connection;
    
      void assign_color ( color_t _color )
      { color = _color; }

      void assign_vertex ( index_t _vertex_index )
      { connection.push_back( _vertex_index ); }

      void assign ( color_t _color, index_t _vertex_index )
      {
	 assign_color( _color );
	 assign_vertex( _vertex_index );
      }
   };

   // Data
   typedef vector< __vertex_type > graph_type;
   graph_type Graph;



public:
   Seeding ( size_t _Nx, size_t _Ny, size_t _Nz, size_t _Neqn = 1 );

   inline void Coloring ();

   template< typename T, typename V, typename M >
   inline void Recover_CSR ( const T& source, V& residual, M& CSR );

   // data access
   color_t Color ( index_t i ) const;
   size_t Max_Color () const;


   void Print ();

private:
   // policy can be added to tackle with various discretization domains
   inline void   __Build_BipartiteGraph ( size_t _Nx, size_t _Ny, size_t _Nz, size_t _Neqn );
   inline void   __Sort_Vertex_Connection ();
   inline void   __Initialize_Colors ();
   inline size_t __Theoretical_MaxColor () const;
};





// ========================================================= 
// ===================== public functions ==================
//
Seeding::Seeding ( size_t _Nx, size_t _Ny, size_t _Nz, size_t _Neqn )
{
   __Build_BipartiteGraph( _Nx, _Ny, _Nz, _Neqn );
   Coloring();
}


void Seeding::Coloring ()
{
   std::cout << "Coloring..." << std::endl;
   // initialize colors
   __Initialize_Colors();

   const size_t max_color = __Theoretical_MaxColor();
   vector<bool> Palette;
   Palette.resize( max_color );
   Palette.assign( max_color, true );

   const size_t N = Graph.size();
   index_t i;
   for( i = 0; i < N; ++i )
   {
//        std::cout << "Vertex --------------------------------------------" << i << std::endl;
      if( Graph[i].color == -1 )
      {
	 size_t n_distance1_index = Graph[i].connection.size();
	 index_t distance1_vertex;
	 size_t j;
	 for ( j = 0; j < n_distance1_index; ++j)
	 {
	    distance1_vertex = Graph[i].connection[j];
	    size_t n_distance2_vertex = Graph[distance1_vertex].connection.size();
                
	    index_t distance2_vertex;
	    size_t k;
	    color_t c_vertex;
	    for ( k = 0; k < n_distance2_vertex; ++k )
	    {
	       distance2_vertex = Graph[distance1_vertex].connection[k];
	       c_vertex = Graph[distance2_vertex].color;
	       if( c_vertex != -1 )
	       {
		  Palette[c_vertex] = false;
	       }
	    }
	 }
            
            
	 // choose the minimum color from Palatte
	 color_t min_color;
	 color_t k;
	 for( k = 0; k < max_color; ++k )
	 {
	    if( Palette[k] != false )
	    {
	       min_color = k;
	       break;
	    }
	 }

	 Graph[i].color = min_color;
//            Seed.insert( make_pair(color,i) );
//            Color[i] = color;
	 // reset Palette
	 Palette.assign( max_color, true );
      }//end_of_if_statement
      else
      {
	 // do nothing
      }
   }
}


/*
 * Recovery()
 *
 * Proposed:
 * Receive an ADvector parameter, and extract values & gradients
 * DIRECTLY from ADscalars and store them into "CSR_Matrix"
 *
 * Current implementation:
 * The default routine "extract_CSR" within the ADvector is called first
 * to obtain the "compressed" version of data in CSR format.
 * Then this "compressed" version of data is modified so that the values
 * sit into their real positions
 */
template< typename T, typename V, typename M >
void Seeding::Recover_CSR ( const T& source, V& residual, M& CSR )
{
   static bool is_first = true;
   std::cout << "is_first\t" << is_first << std::endl;

   // update CSR.rowptr() & CSR.colptr() ONLY for the first time
   if( is_first )
   {
      if( !CSR.rowptr().is_empty() )
	 CSR.rowptr().clear();
      CSR.rowptr().reserve( Graph.size() + 1 );

      if( !CSR.colind().is_empty() )
	 CSR.colind().clear();
      // CSR.colind().reserve( NNZ );

      // CSR.value() will be taken care of in "Recovery()"

      size_t counter_nnz = 0;
      CSR.rowptr().push_back( counter_nnz );

      for( size_t i = 0; i < Graph.size(); ++i )
      {
	 index_t col;
	 for( size_t j = 0; j < Graph[i].connection.size(); ++j )
	 {
	    col = Graph[i].connection[j];
	    CSR.colind().push_back( col );
	 }
	 counter_nnz += Graph[i].connection.size();
	 CSR.rowptr().push_back( counter_nnz );
      }

      is_first = false;
   }


   // update residual & CSR.value()
   if( !residual.empty() )
      residual.clear();
   residual.reserve( Graph.size() );

   if( !CSR.value().is_empty() )
      CSR.value().clear();
   CSR.value().reserve( CSR.colind().size() );


   for( size_t i = 0; i < Graph.size(); ++i )
   {
      residual.push_back( source[i].value() );

      index_t col;
      color_t color;
      for( size_t j = 0; j < Graph[i].connection.size(); ++j )
      {
	 col = Graph[i].connection[j];
	 // "color" tells which column it sits in the compressed array
	 color = Graph[col].color;
	 CSR.value().push_back( source[i].derivative(color) );
      }
   }

   CSR.M() = CSR.rowptr().size() - 1;
   CSR.N() = CSR.colind().size() - 1;
   CSR.NNZ() = CSR.value().size();
}



color_t Seeding::Color ( index_t i ) const
{
   return Graph[i].color;
}

size_t Seeding::Max_Color () const
{
   size_t counter = 0;
   for( index_t i = 0; i < Graph.size(); ++ i )
   {
      counter = ( counter >= Graph[i].color ) ? counter : Graph[i].color;
   }
   return counter;
}




// ============== end of public functions ==================




// ==========================================================
// ================= private functions ======================
/*
 * Build_BipartiteGraph
 * Current implementation:
 * Natural ordering / Cartesian
 */
void Seeding::__Build_BipartiteGraph ( size_t _Nx, size_t _Ny, size_t _Nz, size_t _Neqn )
{
   if ( !Graph.empty() )
   {
      Graph.clear();
   }
   const size_t n_cells = _Nx * _Ny * _Nz;
   Graph.resize( n_cells * _Neqn );

   index_t block_row, block_col;
   index_t row, col;
   size_t Nvar;

   // 1st direction in all layers
   for ( index_t k = 0; k < _Nz; ++k )
   {
      for ( index_t j = 0; j < _Ny; ++j )
      {
	 for ( index_t i = 0; i < _Nx-1; ++i )
	 {
	    block_row = i + j*_Nx + k*_Nx*_Ny;
	    block_col = block_row + 1;
                
	    for( index_t r = 0; r < _Neqn; ++r )
	    {
	       row = block_row * _Neqn + r;
	       Nvar = _Neqn;
	       for( index_t c = 0; c < Nvar; ++c )
	       {
		  col = block_col * Nvar + c;
		  Graph[row].assign(-1, col);
		  Graph[col].assign(-1, row);
	       }
	    }
	 }
      }
   }
   // 2nd direction in all layers
   for ( index_t k = 0; k < _Nz; ++k )
   {
      for ( index_t j = 0; j < _Ny-1; ++j )
      {
	 for ( index_t i = 0; i < _Nx; ++i )
	 {
	    block_row = i + j*_Nx + k*_Nx*_Ny;
	    block_col = block_row + _Nx;
                
	    for( index_t r = 0; r < _Neqn; ++r )
	    {
	       row = block_row * _Neqn + r;
	       Nvar = _Neqn;
	       for( index_t c = 0; c < Nvar; ++c )
	       {
		  col = block_col * Nvar + c;
		  Graph[row].assign(-1, col);
		  Graph[col].assign(-1, row);
	       }
	    }
	 }
      }
   }
   // 3rd direction
   for ( index_t k = 0; k < _Nz-1; ++k )
   {
      for ( index_t j = 0; j < _Ny; ++j )
      {
	 for ( index_t i = 0; i < _Nx; ++i )
	 {
	    block_row = i + j*_Nx + k*_Nx*_Ny;
	    block_col = block_row + _Nx * _Ny;
                
	    for( index_t r = 0; r < _Neqn; ++r )
	    {
	       row = block_row * _Neqn + r;
	       Nvar = _Neqn;
	       for( index_t c = 0; c < Nvar; ++c )
	       {
		  col = block_col * Nvar + c;
		  Graph[row].assign(-1, col);
		  Graph[col].assign(-1, row);
	       }
	    }

	 }
      }
   }
   // Itself
   for ( index_t i = 0; i < n_cells; ++i )
   {
      block_row = i;
      block_col = block_row;

      for( index_t r = 0; r < _Neqn; ++r )
      {
	 row = block_row * _Neqn + r;
	 Nvar = _Neqn;
	 for( index_t c = 0; c < Nvar; ++c )
	 {
	    col = block_col * Nvar + c;
	    Graph[row].assign( -1, col );
	 }
      }
   }
    
   __Sort_Vertex_Connection();
}




/*
 * Sort()
 * Sort the "connection" vector associated with each vertex
 * to provide convenience & efficient in following operations
 */
void Seeding::__Sort_Vertex_Connection ()
{
   vector<index_t>::iterator iter_begin, iter_end;
   for ( index_t i = 0; i < Graph.size(); ++i )
   {
      iter_begin = Graph[i].connection.begin();
      iter_end = Graph[i].connection.end();
      std::sort( iter_begin, iter_end );
   }
}

/*
 * Init()
 * Assign colors to vertices that are distance-1 to the very 1st vertex
 */
void Seeding::__Initialize_Colors ()
{
   size_t max_color = __Theoretical_MaxColor();
    
   color_t color(0);
   graph_type :: const_iterator iter = Graph.begin();
    
   size_t n_distance1_vertex = (*iter).connection.size();
   for ( index_t i = 0; i < n_distance1_vertex; ++i )
   {
      index_t distance1_vertex = (*iter).connection[i];
      Graph[ distance1_vertex ].assign_color(color);
      ++color;
   }
}


/*
 * Get_Max_Color()
 * Get the upper bound of colors to be assinged to vertices
 */
size_t Seeding::__Theoretical_MaxColor () const
{
   size_t max_degree = 0;
   size_t tmp;
   size_t N = Graph.size();

   for( index_t i = 0; i < N; ++i )
   {
      tmp = Graph[i].connection.size();
      max_degree = ( max_degree >= tmp ) ? max_degree : tmp;
   }
   size_t upper_bound = max_degree * max_degree + 1;
   //std::cout << "max_degree = " << max_degree << std::endl;
   return ( upper_bound <= N ? upper_bound : N );
}
// ============== end of private functions ====================
// ============================================================









// Print
template< typename T >
void Print_Vector ( const T& vec )
{
   for( size_t i = 0; i < vec.size(); ++i )
   {
      std::cout << vec[i] << std::endl;
   }
   std::cout << std::endl;
}

void Seeding::Print ()
{
   std::cout << "*********** Graph Summary ************"<< std::endl;
   std::cout << "Graph.size() = " << Graph.size() << std::endl;
   std::cout << "Theoretical Max Color = " << __Theoretical_MaxColor() << std::endl;
   std::cout << "Max Color = " << Max_Color() << std::endl;
   // for ( index_t i = 0; i < Graph.size(); ++i )
   // {
   //    std::cout << "-------- " << i << " ------> Color: " << Graph[i].color << std::endl;
   //    size_t sz = Graph[i].connection.size();
   //    for ( index_t j = 0; j < sz; ++j )
   //    {
   // 	 std::cout << Graph[i].connection[j] << std::endl;
   //    }
   // }
}





