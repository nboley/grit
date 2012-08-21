/* Copyright (c) 2011-2012 Nathan Boley */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

#define INITIAL_ARRAY_SIZE 0

#define MAX_CONTIG_NAME_LEN 255

struct contig_t {
    char* name;
    int size;
    double* values;
};

void
add_values_to_contig( struct contig_t* contig, int start, int stop, double value)
{
    if( stop > contig->size )
    {
        size_t curr_size = contig->size;
        contig->size = stop;
        contig->values = realloc( contig->values, sizeof( double )*stop );
        memset( contig->values + curr_size, 0, stop - curr_size);
    }
    
    int i;
    for( i = start; i < stop; i++ )
    {
        contig->values[ i ] += value;
    }
    
    return;
}

struct contigs_t {
    struct contig_t* contigs;
    int size;
};

void 
init_contigs( struct contigs_t** contigs )
{
    *contigs = malloc( sizeof(struct contigs_t) );
    (*contigs)->size = 0;
    (*contigs)->contigs = NULL;
}

int
add_new_contig( struct contigs_t* contigs, char* name )
{
    contigs->size += 1;
    contigs->contigs = realloc( 
        contigs->contigs, sizeof(struct contig_t)*contigs->size );
    
    struct contig_t* new_contig = &( contigs->contigs[ contigs->size-1 ] );
    /* copy the name over */
    new_contig->name = calloc( sizeof(char), (strlen(name)+1) );
    strncpy( new_contig->name, name, strlen(name) );
    
    new_contig->size = 0;
    new_contig->values = NULL;
    
    return contigs->size-1;
}


int 
find_contig_index( struct contigs_t* contigs, char* name ) 
{
    int i;
    for( i = 0; i < contigs->size; i++ )
    {
        if( 0 == strcmp( name, contigs->contigs[ i ].name )  )
            return i;
    }

    return -1;
};

int
load_bedgraph( char* fname, struct contigs_t** contigs )
{
    FILE* fp = fopen( fname, "r" );
    if( NULL == fp ){
        fprintf( stderr, "Couldnt open the bedgraph file for reading." );
        return 1;
    }
    
    init_contigs( contigs );
    
    while( !feof( fp ) )
    {
        // chr4 . . 4754 6133 . + . 1
        char chrm[ MAX_CONTIG_NAME_LEN ];
        int start = -1 ;
        int stop = -1;
        double value = -1;
        int rv1, rv2, rv3;
        
        fscanf( fp, "%s", chrm );
        if( 0 == strcmp( chrm, "track" ) )
        {
            while( !feof(fp) && fgetc( fp ) != '\n' );
            continue;
        }
        
        int contig_index = find_contig_index( *contigs, chrm );
        if( contig_index == -1 )
        {
            contig_index = add_new_contig( *contigs, chrm );
        }
        
        rv1 = fscanf( fp, "%i", &start );
        rv2 = fscanf( fp, "%i", &stop );
        rv3 = fscanf( fp, "%lf", &value );        
        if( isnan( value ) ) {
            fprintf( stderr, "WHAT???? %e\n", value );
        };
        
        if( rv1 == rv2 == rv3 == 1 )
        {
            assert( start >= 0 );
            assert( stop >= 0 );
            assert( !isnan(value) );
            add_values_to_contig(
                (*contigs)->contigs + contig_index, start,stop,value);
        } 
        
        while( !feof(fp) && fgetc( fp ) != '\n' );        
    }

    /*
    printf( "HERE: %e\n", (*contigs)->contigs[0].values[ 1717 ] );

    int i, j;
    for( i = 0; i < (*contigs)->size; i++ )
    {
        for( j = 0; j < (*contigs)->contigs[i].size; j++ )
        {
            assert( (*contigs)->contigs[i].values[ j ] >= 0 );
            printf( "%e\n", (*contigs)->contigs[i].values[ j ] );
        }
    }
    */
    
    return 0;
}

int main( int argc, char** argv )
{
    if( argc != 2 )
    {
        fprintf( stderr, "Need a bedgraph file.\n" );
        exit( -1 );
    }
    
    struct contigs_t* contigs;
    return load_bedgraph( argv[1], &contigs );

}
