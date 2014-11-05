#pragma once

#include <vector>

using std::size_t;
using std::vector;

typedef int color_t;
typedef size_t index_t;

struct __vertex_type
{
    int color;
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

typedef vector< __vertex_type > graph_type;
typedef vector< vector<index_t> > seed_type;

void Bipartite_Cartesian ( graph_type& Graph, size_t _Nx, size_t _Ny, size_t _Nz )
{
    std::cout << "Bipartite_Cartesian" << std::endl;
    if ( !Graph.empty() )
    {
        Graph.clear();
    }
    size_t N = _Nx * _Ny * _Nz;
    Graph.resize( N );
    index_t _row, _col;
    index_t i, j, k;
    
    // 1st direction in all layers
    for (k=0; k < _Nz; ++k)
    {
        for (j=0; j < _Ny; ++j)
        {
            for (i=0; i < _Nx-1; ++i)
            {
                _row = i + j*_Nx + k*_Nx*_Ny;
                _col = _row + 1;
                Graph[_row].assign(-1, _col);
                Graph[_col].assign(-1, _row);
            }
        }
    }
    // 2nd direction in all layers
    for (k=0; k < _Nz; ++k)
    {
        for (j=0; j < _Ny-1; ++j)
        {
            for (i=0; i < _Nx; ++i)
            {
                _row = i + j*_Nx + k*_Nx*_Ny;
                _col = _row + _Nx;
                Graph[_row].assign(-1, _col);
                Graph[_col].assign(-1, _row);
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
                _row = i + j*_Nx + k*_Nx*_Ny;
                _col = _row + _Nx * _Ny;
                Graph[_row].assign(-1, _col);
                Graph[_col].assign(-1, _row);
            }
        }
    }
    // Itself

    for (i=0; i < N; ++i)
    {
        Graph[i].assign(-1, i);
    }
    
    // Print
    for ( i = 0; i < N; ++i )
    {
        size_t sz = Graph[i].connection.size();
        for ( j = 0; j < sz; ++j )
        {
            std::cout << i << "\t" << Graph[i].connection[j] << std::endl;
        }
    }
    std::cout << "SZ\t" << Graph.size() << std::endl;
}






struct Seeding
{
    static size_t Get_Max_Color( const graph_type& Graph );
    static void Init( graph_type& Graph, seed_type& Seed );
    static void Coloring( graph_type& Graph );
//    static void Get_Seed_Matrix( const Map& Graph, seed_type& Seed, size_t NumVar, size_t NumColor );
};

size_t Seeding::Get_Max_Color( const graph_type& Graph )
{
    size_t max_degree = 0;
    size_t tmp;
    index_t i = 0;
    size_t N = Graph.size();
    std::cout << "N:\t" << N << std::endl;
    for( i = 0; i < N; ++i )
    {
        tmp = Graph[i].connection.size();
        max_degree = ( max_degree >= tmp ) ? max_degree : tmp;
    }
    max_degree = max_degree * max_degree + 1;
    return ( max_degree <= N ? max_degree : N );
}

void Seeding::Init( graph_type& Graph, seed_type& Seed )
{
    std::cout << "Seeding::Init ----------" << std::endl;
    if( !Seed.empty() )
    {
        Seed.clear();
    }
    size_t max_color = Get_Max_Color( Graph );
    Seed.resize( max_color );
    
    std::cout << "max_color = " << max_color << std::endl;
    std::cout << "Seed.size = " << Seed.size() << std::endl;
    
    color_t color(0);
    graph_type :: const_iterator iter = Graph.begin();
    
    size_t n_vertex = (*iter).connection.size();
    size_t i;
    for ( i = 0; i < n_vertex; ++i )
    {
        index_t connected_vertex = (*iter).connection[i];
        Graph[ connected_vertex ].assign_color(color);
        Seed[color].push_back( connected_vertex );
        
        std::cout << "Seed:\t" << color << "\t" << connected_vertex << std::endl;
        
        ++color;
    }
}

void Seeding::Coloring( graph_type& Graph )
{
    const size_t max_color = Get_Max_Color( Graph );
    vector<bool> Palette;
    Palette.resize( max_color );
    Palette.assign( max_color, true );

    const size_t N = Graph.size();
    index_t i;
    for( i = 0; i < N; ++i )
    {
        std::cout << "Vertex --------------------------------------------" << i << std::endl;
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
                        std::cout << "Forbidden Color:\t" << c_vertex << std::endl;
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
                    std::cout << "selected color: " << k << std::endl;
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
            std::cout << "Already has color:\t" << Graph[i].color << std::endl;
        }
        std::cout << std::endl;
    }
}



//void Seeding::Get_Seed_Matrix( const Map& Graph, Map& Seed, size_t NumVar, size_t NumColor)
//{
//    vector<color_t> Index_Color;
//    Index_Color.reserve( NumVar );
//    Index_Color.assign( NumVar, -1 );
//    
//    Init( Graph, Seed, Index_Color );
//    Coloring( Graph, Seed, Index_Color, NumColor );
//}
