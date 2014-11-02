#pragma once

#include <vector>

using std::size_t;
using std::vector;

typedef size_t color_t;
typedef size_t index_t;

struct __vertex_type
{
    index_t color;
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

void Bipartite_Cartesian ( graph_type& Graph, size_t _Nx, size_t _Ny, size_t _Nz, size_t& Num_Color)
{
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
                Graph[_row].assign(0, _col);
                Graph[_col].assign(0, _row);
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
                Graph[_row].assign(0, _col);
                Graph[_col].assign(0, _row);
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
                Graph[_row].assign(0, _col);
                Graph[_col].assign(0, _row);
            }
        }
    }
    // Itself
    size_t N = _Nx * _Ny * _Nz;
    for (i=0; i < N; ++i)
    {
        Graph[i].assign(0, i);
    }
}

size_t Color_Boundary ( const graph_type& Graph )
{
    size_t max_degree = 0;
    size_t tmp;
    index_t i = 0;
    size_t N = Graph.size();
    for( i = 0; i < N; ++i )
    {
        tmp = Graph[i].connection.size();
        max_degree = max_degree >= tmp ? max_degree : tmp;
    }
    max_degree = max_degree * max_degree + 1;
    return ( max_degree <= N ? max_degree : N );
}





struct Seeding
{
    static void Init( graph_type& Graph, seed_type& Seed );
    static void Coloring( const graph_type& Graph, seed_type& Seed, vector<color_t>& Color, size_t NumColor );
    static void Get_Seed_Matrix( const Map& Graph, seed_type& Seed, size_t NumVar, size_t NumColor );
};


void Seeding::Init( graph_type& Graph, seed_type& Seed )
{
    color_t color(0);
    graph_type :: const_iterator iter = Graph.begin();
    
    size_t N = (*iter).connection.size();
    size_t i;
    for (i = 0; i < N; ++i)
    {
        index_t connected_vertex = (*iter).connection[i];
        Graph[ connected_vertex ].assign_color(color);
        
        std::cout << connected_vertex << "\t" << color << std::endl;
        
        ++color;
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
