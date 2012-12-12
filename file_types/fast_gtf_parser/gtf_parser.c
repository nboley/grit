/* Copyright (c) 2011-2012 Nathan Boley */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

enum ELEMENT_TYPE {
    UNKNOWN,
    TPUTR,
    FPUTR,
    UTR,
    CDS,
    GENE,
    TRANSCRIPT,
    START_CODON,
    STOP_CODON
};

struct gtf_line {
    char* chrm;
    int start;
    int stop;
    char strand;
    enum ELEMENT_TYPE element_type;
    int score;
    double rpkm;
    double rpk;
    char* gene_id;
    char* trans_id;
};

int 
cmp_gtf_lines( const void* i1, const void* i2 )
{
    int cmp = 0;
    struct gtf_line** line_1 = ( struct gtf_line** )i1;
    struct gtf_line** line_2 = ( struct gtf_line** )i2;
    
    cmp = strcmp( (*line_1)->gene_id, (*line_2)->gene_id );
    if( cmp != 0 )
        return cmp;
    
    cmp = strcmp( (*line_1)->trans_id, (*line_2)->trans_id );
    if( cmp != 0 )
        return cmp;
    
    return (*line_1)->start - (*line_2)->start;
}

void
get_id_ptr( 
    char* line_str, char* key_str, char** id )
{
    char* id_start;
    int id_len;
    // find the location of the id
    id_start = strstr( line_str, key_str );
    if( id_start == NULL ) {
        *id = NULL;
        return;
    }
    // Move past the gene id
    id_start += strlen( key_str );
    // skip to the next alpha numeric character
    while( !isalnum( *id_start ) )
        (id_start)++;
    
    // find the length of the gene identifier by moving
    // to the next ;, and then backing past non alpha 
    // numeric chars
    id_len = 0;
    while( ';' != *(id_start+id_len) )
        (id_len)++;
    
    while( !isalnum( *(id_start+id_len-1) ) && ')' != *(id_start+id_len-1) )
        (id_len)--;

    *id = calloc( sizeof(char), id_len + 1 );
    memcpy( *id, id_start, sizeof(char)*id_len );
    
    return;
}


enum ELEMENT_TYPE 
find_element_type_from_str( char* element_type_str ) {
    if( 0 == strcmp(element_type_str, ".") )
        return UNKNOWN;

    if( 0 == strcmp(element_type_str, "exon") )
        return UNKNOWN;

    if( 0 == strcmp(element_type_str, "gene") )
        return GENE;

    if( 0 == strcmp(element_type_str, "transcript") )
        return TRANSCRIPT;

    if( 0 == strcmp(element_type_str, "start_codon") )
        return START_CODON;

    if( 0 == strcmp(element_type_str, "stop_codon") )
        return STOP_CODON;

    if( 0 == strcmp(element_type_str, "UTR") )
        return UTR;
    
    if( 0 == strcmp(element_type_str, "3UTR") )
        return TPUTR;

    if( 0 == strcmp(element_type_str, "three_prime_UTR") )
        return TPUTR;
    
    if( 0 == strcmp(element_type_str, "5UTR") )
        return FPUTR;

    if( 0 == strcmp(element_type_str, "five_prime_UTR") )
        return TPUTR;
    
    if( 0 == strcmp(element_type_str, "CDS") )
        return CDS;
    
    fprintf(stderr, "WARNING: unrecognized element type: '");
    fprintf(stderr, "%s", element_type_str);
    fprintf(stderr, "' \n");
    
    return UNKNOWN;
}

struct gtf_line*
get_line_info( char* line ) {
    // Allocate space for the current line
    struct gtf_line* curr_line = 
        malloc( sizeof( struct gtf_line ) );
    
    char chrm[200];
    char element_type_str[200];
    char char_score[200];
    
    // get the exon location
    sscanf( line, "%s %*s %s %i %i %s %s", 
            chrm, element_type_str,
            &(curr_line->start), &(curr_line->stop), 
            char_score, &(curr_line->strand) );
    
    curr_line->chrm = calloc( sizeof(char), strlen( chrm ) + 1 );
    memcpy( curr_line->chrm, chrm, strlen( chrm ) );
    
    curr_line->element_type = find_element_type_from_str( element_type_str );
    
    if( 0 == strcmp( char_score, "." ) )
    {
        curr_line->score = -1;
    } else {
        curr_line->score = atoi( char_score );
    }
    
    // find the gene id
    get_id_ptr( line, "gene_id", &(curr_line->gene_id) );
    if( NULL == curr_line->gene_id ) {
        fprintf( stderr, "INVALID GTF LINE: no gene id\n" );
        free( curr_line );
        return NULL;
    }
    
    // find the transcript id
    get_id_ptr( line, "transcript_id", &(curr_line->trans_id) );
    if( NULL == curr_line->trans_id ) {
        fprintf( stderr, "INVALID GTF LINE: no transcript id\n" );
        free( curr_line );
        return NULL;
    }


    char* char_rpkm = NULL;
    get_id_ptr( line, "rpkm", &(char_rpkm) );
    if( NULL == char_rpkm ) {
        curr_line->rpkm = -1;
    } else {
        curr_line->rpkm = strtod( char_rpkm, NULL );
        free( char_rpkm );
        
        //fprintf( stderr, "INVALID GTF LINE: no transcript id\n" );
        //free( curr_line );
        //return NULL;
    }

    char* char_rpk = NULL;
    get_id_ptr( line, "rpk ", &(char_rpk) );
    if( NULL == char_rpk ) {
        curr_line->rpk = -1;
    } else {
        curr_line->rpk = strtod( char_rpk, NULL );
        free( char_rpk );
        
        //fprintf( stderr, "INVALID GTF LINE: no transcript id\n" );
        //free( curr_line );
        //return NULL;
    }

    
    return curr_line;
}

void
load_gtf_data( char* fname, struct gtf_line*** gtf_lines, int* num_gtf_lines )
{
    FILE* fp;
    char buffer[10000];
    
    *gtf_lines = NULL;
    *num_gtf_lines = 0;
    
    int num_lines_parsed = 0;
    int num_lines_alloced = 0;
    #define LINE_ALLOC_SIZE 10000;
    
    fp = fopen( fname, "r" );
    if( NULL == fp ){
        fprintf( stderr, "Couldnt open the gtf file for reading." );
        return;
    }

    while( !feof( fp ) )
    {
        // Read the next line
        char* res = fgets( buffer, 10000, fp );
        // if there was an error or eof
        if( res == NULL ) {
            // if eof, just break out of the loop
            if( feof( fp ) ) {
                break;
            } 
            // else, print the error
            else {
                perror ("Error opening file");
            }
        }
        
        // check to see if it's a track linme. If so, skip it
        char* track_substr = strstr( res, "track");
        if( track_substr == res )
            continue;
        
        struct gtf_line* curr_line = get_line_info( buffer );
        if( NULL == curr_line ) continue;
        
        if( num_lines_parsed >= num_lines_alloced )
        {
            num_lines_alloced += LINE_ALLOC_SIZE;
            *gtf_lines = realloc( 
                *gtf_lines, sizeof(struct gtf_line*)*num_lines_alloced );
        }
        
        (*gtf_lines)[ num_lines_parsed ] = curr_line;
        num_lines_parsed += 1;
    }

    // Realloc to save memory
    *gtf_lines = realloc( 
        *gtf_lines, sizeof(struct gtf_line*)*num_lines_parsed );
    
    *num_gtf_lines = num_lines_parsed;

    qsort( *gtf_lines, num_lines_parsed, sizeof(struct gtf_line*), cmp_gtf_lines );

    fclose( fp );
    
    return;
}

struct transcript {
    char* trans_id;
    int cds_start;
    int cds_stop;
    int num_exon_bnds;
    int* exon_bnds;
    int score;
    double rpkm;
    double rpk;
};

struct gene {
    char* gene_id;
    char* chrm;
    char strand;
    int min_loc;
    int max_loc;
    int num_transcripts;
    struct transcript** transcripts;
};

#define TRANS_ALLOC_SIZE 100;

void
add_transcript( struct transcript*** transcripts, 
                int* num_transcripts, int* num_allcd_trans,
                char* trans_id, 
                int num_exons, int* exon_bnds,
                int cds_start, int cds_stop,
                int score, double rpkm, double rpk )
{
    if( *num_transcripts >= *num_allcd_trans ) {
        *num_allcd_trans += TRANS_ALLOC_SIZE;
        *transcripts = realloc( 
            *transcripts, (*num_allcd_trans)*sizeof( struct transcript* ) );
    }
    
    (*transcripts)[ *num_transcripts ] = malloc( sizeof( struct transcript ) );
    
    (*transcripts)[ *num_transcripts ]->trans_id 
        = calloc(strlen(trans_id)+1, 1);
    strcpy( (*transcripts)[ *num_transcripts ]->trans_id, trans_id );
    
    (*transcripts)[ *num_transcripts ]->num_exon_bnds = num_exons;
    /* init memory for, and copy over the exon bounds */
    (*transcripts)[ *num_transcripts ]->exon_bnds
        = calloc( sizeof( int ), num_exons );
    
    memcpy( (*transcripts)[ *num_transcripts ]->exon_bnds, 
            exon_bnds, sizeof(int)*num_exons );
    
    (*transcripts)[ *num_transcripts ]->cds_start = cds_start;
    (*transcripts)[ *num_transcripts ]->cds_stop = cds_stop;
    
    (*transcripts)[ *num_transcripts ]->score = score;
    (*transcripts)[ *num_transcripts ]->rpkm = rpkm;
    (*transcripts)[ *num_transcripts ]->rpk = rpk;
    
    (*num_transcripts)++;
    
    return;
}

#define GENE_ALLOC_SIZE 2;

void
add_gene( struct gene*** genes, int* num_genes, int* num_allcd_genes,
          char* gene_id, char* chrm, char strand, 
          int min_loc, int max_loc,
          int num_trans, struct transcript** transcripts  )
{
    if( *num_genes >= *num_allcd_genes ) {
        *num_allcd_genes += GENE_ALLOC_SIZE;
        *genes = realloc( *genes, (*num_allcd_genes)*sizeof( struct gene* ) );
    }
    
    (*genes)[ *num_genes ] = malloc( sizeof( struct gene ) );
    
    (*genes)[ *num_genes ]->gene_id = calloc(strlen(gene_id)+1, 1);
    strcpy( (*genes)[ *num_genes ]->gene_id, gene_id );
    
    (*genes)[ *num_genes ]->chrm = calloc(strlen(chrm)+1, 1);
    strcpy( (*genes)[ *num_genes ]->chrm, chrm );
    
    (*genes)[ *num_genes ]->strand = strand;
    (*genes)[ *num_genes ]->num_transcripts = num_trans;
    (*genes)[ *num_genes ]->transcripts = malloc( 
        num_trans*sizeof(struct transcript*) );
    
    /* we need to copy the pointers because we may realloc transcripts */
    memcpy( (*genes)[ *num_genes ]->transcripts, 
            transcripts, num_trans*sizeof(struct transcript*) );
    
    (*genes)[ *num_genes ]->min_loc = min_loc;
    (*genes)[ *num_genes ]->max_loc = max_loc;
    
    (*num_genes)++;
    
    return;
}                

void
parse_gtf_data( struct gtf_line** gtf_lines, int num_lines, 
                struct gene*** genes_ptr, int* num_genes_ptr ) 
{
    char* prev_gene_id = gtf_lines[ 0 ]->gene_id;
    char* prev_trans_id = gtf_lines[ 0 ]->trans_id;
    char* prev_chrm = gtf_lines[ 0 ]->chrm;
    char prev_strand = gtf_lines[ 0 ]->strand;
    int curr_trans_exons[ 10000 ];
    memset( curr_trans_exons, 0, sizeof(int)*10000 );
    curr_trans_exons[0] = gtf_lines[ 0 ]->start;
    curr_trans_exons[1] = gtf_lines[ 0 ]->stop;
    int curr_trans_num_exons = 2;
    int curr_trans_CDS_start = -1;
    int curr_trans_CDS_stop = -1;
    if( gtf_lines[0]->element_type == CDS ) {
        curr_trans_CDS_start = gtf_lines[ 0 ]->start;
        curr_trans_CDS_stop = gtf_lines[ 0 ]->stop;
    }
    
    struct transcript** transcripts = NULL;
    int num_transcripts = 0;
    int num_trans_in_curr_gene = 0;
    int num_allcd_trans = 0;
    
    struct gene** genes = NULL;
    int num_genes = 0;
    
    int num_allcd_genes = 0;
    int curr_max_loc = -1;
    int curr_min_loc = (1 >> 30);
    
    int curr_score = -1;
    double curr_rpkm = -1;
    double curr_rpk = -1;
    
    int i;
    for( i = 1; i < num_lines; i++ )
    {
        struct gtf_line* line = gtf_lines[i];
        
        /* if we've moved to a new transcript, then add the 
           previous to the transcripts list, and reset the counters */
        if( 0 != strcmp( prev_trans_id, line->trans_id ) )
        {
            add_transcript( 
                &transcripts, &num_transcripts, &num_allcd_trans,
                prev_trans_id, curr_trans_num_exons, curr_trans_exons,
                curr_trans_CDS_start, curr_trans_CDS_stop, 
                curr_score, curr_rpkm, curr_rpk );
            
            // Move onto the next transcript
            prev_trans_id = line->trans_id;
            memset( curr_trans_exons, 0, sizeof(int)*10000 );
            curr_trans_num_exons = 0;
            curr_trans_CDS_start = -1;
            curr_trans_CDS_stop = -1;
            
            num_trans_in_curr_gene++;
        }
        
        /* if we've moved to a new gene */
        if( 0 != strcmp( prev_gene_id, line->gene_id ) )
        {
            struct transcript** curr_transcripts  = 
                transcripts + num_transcripts - num_trans_in_curr_gene;
            
            add_gene( &genes, &num_genes, &num_allcd_genes,
                      prev_gene_id, prev_chrm, prev_strand, 
                      curr_min_loc, curr_max_loc,
                      num_trans_in_curr_gene, curr_transcripts );
            
            // Move onto the next gene
            prev_gene_id = line->gene_id;
            prev_chrm = line->chrm;
            prev_strand = line->strand;

            num_trans_in_curr_gene = 0;
            curr_max_loc = -1;
            curr_min_loc = (1 << 30);            
        }
       
        /******** update the CDS *******/
        // if this is coding seqeunce
        if( line->element_type == CDS ) 
        {
            // and we havn't seen coding sequence yet, then update both the 
            // start and the stop
            if( -1 == curr_trans_CDS_start ) 
            {
                curr_trans_CDS_start = line->start;
                curr_trans_CDS_stop = line->stop;
            }
            // otherwise, just update the stop
            else {
                curr_trans_CDS_stop = line->stop;
            }
        }

        /* finally, add the exons and update the bounds */
        curr_score = line->score;
        curr_rpkm = line->rpkm;
        curr_rpk = line->rpk;
        curr_min_loc = MIN( curr_min_loc, line->start );
        curr_max_loc = MAX( curr_max_loc, line->stop );
        curr_trans_exons[ curr_trans_num_exons ] = line->start;
        curr_trans_exons[ curr_trans_num_exons+1 ] = line->stop;
        
        curr_trans_num_exons += 2;

    }

    /* Add the final transcript */
    add_transcript( &transcripts, &num_transcripts, &num_allcd_trans,
                    prev_trans_id, curr_trans_num_exons, curr_trans_exons,
                    curr_trans_CDS_start, curr_trans_CDS_stop, 
                    curr_score, curr_rpkm, curr_rpk );
    num_trans_in_curr_gene++;
    
    struct transcript** curr_transcripts  = 
        transcripts + num_transcripts - num_trans_in_curr_gene;

    add_gene( &genes, &num_genes, &num_allcd_genes,
              prev_gene_id, prev_chrm, prev_strand, 
              curr_min_loc, curr_max_loc,
              num_trans_in_curr_gene, curr_transcripts );
    
    free( transcripts );
    /* realloc to save memory from the block alloc */
    genes = realloc( genes, num_genes*sizeof( struct gene* ) );
    
    /* pass back the gene data, by ref. The transcript data is pointed 
       to by the genes. */
    *genes_ptr = genes;
    *num_genes_ptr = num_genes;
    
    return;
}

void
get_gtf_data( char* fname, struct gene*** arg_genes, int* arg_num_genes )
{
    struct gtf_line** gtf_lines = NULL;
    int num_lines = 0;
    load_gtf_data( fname, &gtf_lines, &num_lines );
    
    if( 0 == num_lines )
    {
        *arg_genes = NULL;
        *arg_num_genes = 0;
        return;
    }

    struct gene** genes;
    int num_genes;
    parse_gtf_data( gtf_lines, num_lines, &genes, &num_genes );
    
    int i;
    for( i = 0; i < num_lines; i++ )
    {
        free( gtf_lines[i]->chrm );
        free( gtf_lines[i]->gene_id );
        free( gtf_lines[i]->trans_id );
        free( gtf_lines[i] );
    }
    
    free( gtf_lines );
    
    *arg_genes = genes;
    *arg_num_genes = num_genes;
    
    return;
}


void
free_gtf_data( struct gene** genes, int num_genes )
{
    int i;
    for( i = 0; i < num_genes; i++ )
    {
        int j;
        for( j = 0; j < genes[i]->num_transcripts; j++ )
        {
            free( genes[i]->transcripts[j]->trans_id );
            free( genes[i]->transcripts[j]->exon_bnds );
            free( genes[i]->transcripts[j] );
        }
        
        free( genes[i]->chrm );
        free( genes[i]->gene_id );
        free( genes[i]->transcripts );
        free( genes[i] );
    }
    
    free( genes );
}

int main( int argc, char** argv )
{
    if( argc != 2 )
        exit( -1 );
    
    struct gene** genes;
    int num_genes;
    get_gtf_data( argv[1], &genes, &num_genes );

    
    int i;
    for( i = 0; i < num_genes; i++ ) 
    {
        printf( "NEW GENE: %s ( %i transcripts )\n", 
                genes[ i ]->gene_id,
                genes[ i ]->num_transcripts );
        
        int j;
        for( j = 0; j < genes[ i ]->num_transcripts; j++ )
        {
            printf( "\t%s %i %i", genes[ i ]->transcripts[j]->trans_id, 
                    genes[ i ]->transcripts[j]->cds_start, 
                    genes[ i ]->transcripts[j]->cds_stop );                
        }
        printf( "\n" );
    }
    
    
    free_gtf_data( genes, num_genes );
    
    return 0;
}
