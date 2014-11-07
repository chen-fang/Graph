#pragma once

#include <vector>
#include <algorithm>

using std::size_t;
using std::vector;

typedef int color_t;
typedef size_t index_t;

class Seeding
{
   struct __vertex_type
   {
      color_t color;
      vector<index_t> connection;
    
      void assign_color ( color_t _color )
      {
	 color = _color;
      }
    
      void assign_vertex ( index_t _vertex_index )
      {
	 connection.push_back(_vertex_index);
      }
    
      void assign ( color_t _color, index_t _vertex_index )
      {
	 assign_color(_color);
	 assign_vertex(_vertex_index);
      }
   };

   // Data
   typedef vector< __vertex_type > graph_type;
   graph_type Graph;

public:
   // policy can be added to tackle with various discretization domains
   inline void Build_BipartiteGraph ( size_t _Nx, size_t _Ny, size_t _Nz, size_t _n_eqn = 1 );
   inline void Coloring ();
   inline void Recovery ( const vector<double>& old_value,
			  vector<int>& new_rowptr,
			  vector<int>& new_colind, vector<double>& new_value );

   // data access
   color_t Color ( index_t i ) const;


   void Print ();

private:
   inline void Sort ();
   inline void Init ();
   inline size_t Get_Max_Color () const;
};





// ========================================================= 
// ===================== public functions ==================
//
void Seeding::Build_BipartiteGraph ( size_t _Nx, size_t _Ny, size_t _Nz, size_t _n_eqn )
{
   if ( !Graph.empty() )
   {
      Graph.clear();
   }
   size_t n_cells = _Nx * _Ny * _Nz;
   Graph.resize( n_cells * _n_eqn );
   index_t row, col;
   index_t i, j, k, m;
    
   // 1st direction in all layers
   for ( k = 0; k < _Nz; ++k )
   {
      for ( j = 0; j < _Ny; ++j )
      {
	 for ( i = 0; i < _Nx-1; ++i )
	 {
	    row = i + j*_Nx + k*_Nx*_Ny;
	    col = row + 1;
	    Graph[row].assign(-1, col);
	    Graph[col].assign(-1, row);
                
//                row = ( i + j*_Nx + k*_Nx*_Ny ) * _n_eqn;
//                col = ( row + 1 ) * _n_eqn;
//                for( m = 0; m < _n_eqn; ++m )
//                {
//                    Graph[row+m].assign(-1, col+m);
//                    Graph[col+m].assign(-1, row+m);
//                }
	 }
      }
   }
   // 2nd direction in all layers
   for ( k = 0; k < _Nz; ++k )
   {
      for ( j = 0; j < _Ny-1; ++j )
      {
	 for ( i = 0; i < _Nx; ++i )
	 {
	    row = i + j*_Nx + k*_Nx*_Ny;
	    col = row + _Nx;
	    Graph[row].assign(-1, col);
	    Graph[col].assign(-1, row);
                
//                row = ( i + j*_Nx + k*_Nx*_Ny ) * _n_eqn;
//                col = ( row + _Nx ) * _n_eqn;
//                for( m = 0; m < _n_eqn; ++m )
//                {
//                    Graph[row+m].assign(-1, col+m);
//                    Graph[col+m].assign(-1, row+m);
//                }
	 }
      }
   }
   // 3rd direction
   for ( k = 0; k < _Nz-1; ++k )
   {
      for ( j = 0; j < _Ny; ++j )
      {
	 for ( i = 0; i < _Nx; ++i )
	 {
	    row = i + j*_Nx + k*_Nx*_Ny;
	    col = row + _Nx * _Ny;
	    Graph[row].assign(-1, col);
	    Graph[col].assign(-1, row);
                
//                row = ( i + j*_Nx + k*_Nx*_Ny ) * _n_eqn;
//                col = ( row + _Nx * _Ny ) * _n_eqn;
//                for( m = 0; m < _n_eqn; ++m )
//                {
//                    Graph[row+m].assign(-1, col+m);
//                    Graph[col+m].assign(-1, row+m);
//                }

	 }
      }
   }
   // Itself
   for (i=0; i < n_cells; ++i)
   {
      row = i;
      col = row;
      Graph[row].assign(-1, col);
        
//        row = i * _n_eqn;
//        col = row;
//        for( m = 0; m < _n_eqn; ++m )
//        {
//            Graph[row+m].assign(-1, col+m);
//            Graph[col+m].assign(-1, row+m);
//        }
   }
    
   // Sort "connection"
   Sort();
}


void Seeding::Coloring ()
{
   // initialize colors
   Init();

   const size_t max_color = Get_Max_Color();
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
void Seeding::Recovery( const vector<double>& old_value,
			vector<int>& new_rowptr,
			vector<int>& new_colind, vector<double>& new_value )
{
   new_rowptr.clear();
   new_rowptr.reserve( Graph.size() );
   new_colind.clear();
   // ...
   new_value.clear();
   new_value.reserve( old_value.size() );

   size_t counter_nnz = 0;
   new_rowptr.push_back( counter_nnz );
   for( size_t i = 0; i < Graph.size(); ++i )
   {
      for( size_t j = 0; j < Graph[i].connection.size(); ++j )
      {
	 index_t col = Graph[i].connection[j];
	 new_colind.push_back( col ); // copy...faster?

	 // "color" tells which column it is in the compressed array
	 color_t color = Graph[col].color;
	 new_value.push_back( old_value[ counter_nnz + color ] );
      }

      counter_nnz += Graph[i].connection.size();
      new_rowptr.push_back( counter_nnz );

   }
}




color_t Seeding::Color ( index_t i ) const
{
   return Graph[i].color;
}




// ============== end of public functions ==================




// ==========================================================
// ================= private functions ======================
/*
 * Sort()
 * Sort the "connection" vector associated with each vertex
 * to provide convenience & efficient in following operations
 */
void Seeding::Sort ()
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
void Seeding::Init ()
{
   size_t max_color = Get_Max_Color();
    
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
size_t Seeding::Get_Max_Color () const
{
   size_t max_degree = 0;
   size_t tmp;
   size_t N = Graph.size();

   for( index_t i = 0; i < N; ++i )
   {
      tmp = Graph[i].connection.size();
      max_degree = ( max_degree >= tmp ) ? max_degree : tmp;
   }
   max_degree = max_degree * max_degree + 1;
   return ( max_degree <= N ? max_degree : N );
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
   for ( index_t i = 0; i < Graph.size(); ++i )
   {
      std::cout << "-------- " << i << " ------> Color: " << Graph[i].color << std::endl;
      size_t sz = Graph[i].connection.size();
      for ( index_t j = 0; j < sz; ++j )
      {
	 std::cout << Graph[i].connection[j] << std::endl;
      }
   }
}





