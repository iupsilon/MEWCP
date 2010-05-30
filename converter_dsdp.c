/*****************************************************************************
 *  Project: Maximum Edge Weighted Clique Problem
 * 
 *  Authors:
 *  (c) 2009 Yari Melzani (yari.melzani@gmail.com) 
 *
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <stdbool.h>
#include "converter_dsdp.h"
#include "MEWCP_dsdp.h"
#include "MEWCP_tabu_definitions.h"


void MEWCP_generate_sdp_constraints(matrix_weights_t * matrix_weights, constraint_t ** sdp_constraints_matrix, double ** bi, constraint_t ** vect_mat_braching_constraints, double alpha,
                                    bool c_cardinality,
                                    bool c_simple_MC,
                                    bool c_improved_MC_A,
                                    bool c_improved_MC_B,
                                    bool c_improved_MC_C,
                                    bool c_4C2,
                                    bool c_4C3_A,
                                    bool c_4C3_B,
                                    unsigned int * out_num_contraints
                                   )
{



    unsigned int n,m,c,k;
    unsigned int number_constraints = 0;
    unsigned int num_blocks;
    unsigned int vector_length;






    /* Number of blocks is 1 */
    num_blocks = NUM_BLOCKS;

    /* Retreive n e m */
    /* loading di n,m */

    n = matrix_weights->n;
    m = matrix_weights->m;


    /* Compute C = n/m */
    if ((n % m) != 0) /* Checking that n is integer multiple of m */
    {
        fprintf(stderr,"** ERROR!, n must divide m **\n");
        exit(EXIT_FAILURE);
    }

    /* c is the cardinality of each partition Vk, k=1,..,m */
    c=n/m;

    /* The number of elements of the matrix is */
    vector_length = n*(n+1)/2;

    /* Scrivo le prime 5 righe di intestazione prima delle matrici */

    /* Calcolo il numero delle matrici dei vincoli */
    if (c_cardinality == true)
    {
        number_constraints += 1;
    }
    if ( c_simple_MC == true )
    {
        number_constraints += m;
    }
    if ( c_improved_MC_A == true )
    {
        number_constraints += n;
    }
    if ( c_improved_MC_B == true )
    {
        number_constraints += n*(m-1);
    }
    if ( c_improved_MC_C == true )
    {
        number_constraints += m;
    }
    if ( c_4C2 == true )
    {
        number_constraints += 1;
    }
    if ( c_4C3_A == true )
    {
        number_constraints += n;
    }
    if ( c_4C3_B == true )
    {
        number_constraints += n;
    }
    number_constraints += 1; /* I consider the branching constraints */



#if defined MEWCP_CONVERTER_DSDP_VERBOSE1

    printf("Constraints number: %d\n",number_constraints);
#endif

    /* End computing the number of contraints */

#if defined MEWCP_CONVERTER_DSDP_VERBOSE1

    printf("\n");
#endif



    /* I allocate the vect_mat_braching_constraints */
    *vect_mat_braching_constraints = MEWCP_allocate_vect_mat_branching_constraints(vector_length);

    /* I allocate sdp_constraints_matrix and vector bi */
    *sdp_constraints_matrix = MEWCP_allocate_sdp_constraints_matrix(n,number_constraints);
    *bi = MEWCP_allocate_bi(number_constraints);





    /******		BEGIN CONTRAINTS	********/
    /*	MATRIX related to variable i=0  is W  */

    k=0;

    MEWCP_write_constraints_weights_matrix(matrix_weights, *sdp_constraints_matrix,vector_length, n,k);
#if defined MEWCP_CONVERTER_DSDP_DEBUG

    MEWCP_print_contraints_matrix(*sdp_constraints_matrix,number_constraints,vector_length);
#endif

    /* I start writing constraints */
    k = 1;

    if (c_cardinality == true)
    {
        MEWCP_write_constraints_cardinality(*sdp_constraints_matrix,vector_length,n,k);
        MEWCP_write_b_cardinality(*bi,k-1,m);
        k += 1;
#if defined MEWCP_CONVERTER_DSDP_DEBUG

        MEWCP_print_contraints_matrix(*sdp_constraints_matrix,number_constraints,vector_length);
        MEWCP_print_vectorY(*bi,number_constraints);
#endif

    }
    if ( c_simple_MC == true )
    {
        MEWCP_write_constraints_simple_MC(*sdp_constraints_matrix,vector_length,n,c,k);

        MEWCP_write_b_simple_MC(*bi,k-1,n,m);
        k += m;


#if defined MEWCP_CONVERTER_DSDP_DEBUG

        MEWCP_print_contraints_matrix(*sdp_constraints_matrix,number_constraints,vector_length);
        MEWCP_print_vectorY(*bi,number_constraints);
#endif

    }
    if ( c_improved_MC_A == true )
    {
        MEWCP_write_constraints_improved_MC_A(*sdp_constraints_matrix,vector_length,n,c,k);
        MEWCP_write_b_improved_MC_A(*bi,k-1,n);

        k += n;
#if defined MEWCP_CONVERTER_DSDP_DEBUG

        MEWCP_print_contraints_matrix(*sdp_constraints_matrix,number_constraints,vector_length);
        MEWCP_print_vectorY(*bi,number_constraints);
#endif

    }
    if ( c_improved_MC_B == true )
    {
        MEWCP_write_constraints_improved_MC_B(*sdp_constraints_matrix,vector_length,n,m,c,k);
        MEWCP_write_b_improved_MC_B(*bi,k-1,n,m);

        k += n*(m-1);
#if defined MEWCP_CONVERTER_DSDP_DEBUG

        MEWCP_print_contraints_matrix(*sdp_constraints_matrix,number_constraints,vector_length);
        MEWCP_print_vectorY(*bi,number_constraints);
#endif

    }
    if ( c_improved_MC_C == true )
    {
        MEWCP_write_constraints_improved_MC_C(*sdp_constraints_matrix,vector_length,m,c,k);
        MEWCP_write_b_improved_MC_C(*bi,k-1,m);
        k += m;
#if defined MEWCP_CONVERTER_DSDP_DEBUG

        MEWCP_print_contraints_matrix(*sdp_constraints_matrix,number_constraints,vector_length);
        MEWCP_print_vectorY(*bi,number_constraints);
#endif

    }
    if ( c_4C2 == true )
    {
        MEWCP_write_constraints_4C2(*sdp_constraints_matrix,vector_length,n,m,c,k);
        MEWCP_write_b_4C2(*bi,k-1);
        k += 1;
#if defined MEWCP_CONVERTER_DSDP_DEBUG

        MEWCP_print_contraints_matrix(*sdp_constraints_matrix,number_constraints,vector_length);
        MEWCP_print_vectorY(*bi,number_constraints);
#endif

    }
    if ( c_4C3_A == true )
    {

        MEWCP_write_constraints_4C3_A(*sdp_constraints_matrix,vector_length,n,m,c,k);
        MEWCP_write_b_4C3_A(*bi,k-1,n);

        k += n;
#if defined MEWCP_CONVERTER_DSDP_DEBUG

        MEWCP_print_contraints_matrix(*sdp_constraints_matrix,number_constraints,vector_length);
        MEWCP_print_vectorY(*bi,number_constraints);
#endif

    }
    if ( c_4C3_B == true )
    {
        MEWCP_write_constraints_4C3_B( *sdp_constraints_matrix,vector_length,n,m,c,k)
        ;
        MEWCP_write_b_4C3_B(*bi,k-1,n,m);
        k += n;
#if defined MEWCP_CONVERTER_DSDP_DEBUG

        MEWCP_print_contraints_matrix(*sdp_constraints_matrix,number_constraints,vector_length);
        MEWCP_print_vectorY(*bi,number_constraints);
#endif

    }

    /* I write the bi of the branching constraints */
    MEWCP_write_b_branching(*bi,k-1);
    k += 1;

    /* 	Fine generazione matrici dei vincoli */

    /* I set the out variable */

    *out_num_contraints = number_constraints;   // branching constraints included


#if defined MEWCP_CONVERTER_DSDP_VERBOSE1

    printf("\n** End generating contraints **\n");
#endif
}


void MEWCP_generate_constraints_branch(list_blocked_nodes_t * list_blocked_nodes, constraint_t * vect_mat_contraints, const unsigned int vector_length, const unsigned int n, const unsigned int c)
{

    unsigned int i;
    int boundaries[4];
    unsigned int pos_vect_matrix;
    unsigned int num_blocked_nodes;
    unsigned int current_blocked_node;
    unsigned int pos_elem = 0;  /* The position in the constraint array */

    num_blocked_nodes = list_blocked_nodes->num_blocked_nodes;

    for (i=0; i<num_blocked_nodes;++i)
    {
        current_blocked_node = list_blocked_nodes->blocked_node[i];
        trova_boundaries_diagonale(c,current_blocked_node,boundaries);

        // I write  (cur_, curr_)= 1
        pos_vect_matrix = MEWCP_convert_coords_ij_to_vector_matrix(current_blocked_node+1,current_blocked_node+1);


        vect_mat_contraints->index[pos_elem] = pos_vect_matrix;
        vect_mat_contraints->weight[pos_elem] = (double) 1;
        ++pos_elem;

        /* Adding non diagonal elements */
        /*
        unsigned int j;
        for (j=i+1;j<n;++j)
        {
            if (  j > boundaries[1] ) // if (i,j) isn't in the diagonal box of (i,i)
            {

                pos_vect_matrix = MEWCP_convert_coords_ij_to_vector_matrix(j+1,current_blocked_node  +1);

                vect_mat_contraints->index[pos_elem] = pos_vect_matrix;
                vect_mat_contraints->weight[pos_elem] = (double) 1;
                ++ pos_elem;
            }
        }

        for(j=0;j<current_blocked_node;++j)
        {
            if (j<boundaries[0])
            {
                pos_vect_matrix = MEWCP_convert_coords_ij_to_vector_matrix(current_blocked_node +1,j+1);

                vect_mat_contraints->index[pos_elem] = pos_vect_matrix;
                vect_mat_contraints->weight[pos_elem] = (double) 1;
                ++ pos_elem;
            }
        }
        
        */
        /* End adding non diagonal elements */


    }

    vect_mat_contraints->num_nz = pos_elem;
}


void MEWCP_write_constraints_weights_matrix(matrix_weights_t * matrix_weights, constraint_t * sdp_constraints_matrix, const unsigned int vector_length, const unsigned int n, unsigned int k)
{

#if defined MEWCP_CONVERTER_DSDP_VERBOSE1
    printf("* MEWCP_write_constraints_weights_matrix *\n");
#endif

    unsigned int j;
    unsigned int i;
    weight_t w;
    unsigned int pos_vect_matrix;
    unsigned int pos_elem = 0;  /* The position in the constraint array */


    /* prendo i valori [i,j] valore */
    for(j=0;j<n;++j)
    {
        for(i=0;i<= j;++i)
        {
            w = matrix_weights->weight[j][i];
            /* printf("[%d,%d] %d\t",i+1,j+1,w); */

            pos_vect_matrix = MEWCP_convert_coords_ij_to_vector_matrix(j+1,i+1);
#if defined MEWCP_CONVERTER_DSDP_DEBUG

            printf("[%d][%d]-> %d\n",j,i,pos_vect_matrix);
#endif

            if (i != j)
            {

                sdp_constraints_matrix[k].index[pos_elem] = pos_vect_matrix;
                sdp_constraints_matrix[k].weight[pos_elem] = -(double) w/2;
                ++pos_elem;
            }
            else
            {
                sdp_constraints_matrix[k].index[pos_elem] = pos_vect_matrix;
                sdp_constraints_matrix[k].weight[pos_elem] = -(double) w;
                pos_elem++;
            }

#if defined MEWCP_CONVERTER_DSDP_VERBOSE1
            if (i != j)
            {
                printf("0 1 %d %d %7.1lf\n",i+1,j+1,(double) w/2);
            }
            else
            {
                printf("0 1 %d %d %7.1lf\n",i+1,j+1,(double) w);
            }
#endif

        }
    }
    sdp_constraints_matrix[k].num_nz = pos_elem;
}



void MEWCP_write_constraints_cardinality( constraint_t * sdp_constraints_matrix, const unsigned int vector_length, const unsigned int n,  unsigned int k)
{
#if defined MEWCP_CONVERTER_DSDP_VERBOSE1
    printf("\n* MEWCP_write_constraints_cardinality *\n");
#endif

    unsigned int i;
    unsigned int pos_vect_matrix;
    unsigned int pos_elem = 0;  /* The position in the constraint array */

    // Vincolo di cardinalità
    for (i=0; i<n ;++i)
    {
        pos_vect_matrix = MEWCP_convert_coords_ij_to_vector_matrix(i+1,i+1);

        sdp_constraints_matrix[k].index[pos_elem] = pos_vect_matrix;
        sdp_constraints_matrix[k].weight[pos_elem] = (double) 1;
        ++pos_elem;


#if defined MEWCP_CONVERTER_DSDP_VERBOSE1

        printf("%d %d 1\tcardinality %d\n",i+1, i +1,k);
#endif

    }
    sdp_constraints_matrix[k].num_nz = pos_elem;

    /* Fine vincolo di cardinalita */
}

void MEWCP_write_constraints_simple_MC( constraint_t * sdp_constraints_matrix, const unsigned int vector_length, const unsigned int n, const unsigned int c,  unsigned int k)
{
#if defined MEWCP_CONVERTER_DSDP_VERBOSE1
    printf("\n* MEWCP_write_constraints_simple_MC *\n");
#endif

    unsigned int i;
    unsigned int pos_vect_matrix;
    unsigned int pos_elem = 0;  /* The position in the constraint array */

    /* Simple multiple choice Sk*X = 1 */
    for (i=1; i<=n; ++i)
    {


        pos_vect_matrix = MEWCP_convert_coords_ij_to_vector_matrix(i,i);

        sdp_constraints_matrix[k].index[pos_elem] = pos_vect_matrix;
        sdp_constraints_matrix[k].weight[pos_elem] = (double) 1;
        ++ pos_elem;



#if defined MEWCP_CONVERTER_DSDP_VERBOSE1

        printf("%d 1 %d %d 1\t\tSimple MC: %d\n",k,i,i,k);
#endif

        if ((i % c) == 0)
        {
            sdp_constraints_matrix[k].num_nz = pos_elem;
            pos_elem = 0;
            ++k;

        }

    }
    /* Fine vincoli di multiple choice */

}

void MEWCP_write_constraints_improved_MC_A(constraint_t * sdp_constraints_matrix, const unsigned int vector_length, const unsigned int n, const unsigned int c,  unsigned int k)
{
#if defined MEWCP_CONVERTER_DSDP_VERBOSE1
    printf("\n* MEWCP_write_constraints_improved_MC_A *\n");
#endif

    unsigned int i;
    unsigned int j;
    unsigned int pos_vect_matrix;
    int boundaries[4];
    unsigned int pos_elem = 0;  /* The position in the constraint array */


    // For each partition Sk i generate for each vetex belonging to it
    // M_ik*X = 0   IMC A
    /* i belongs to partition k */

    for (i=0;i<n;++i)  // I consider Xii=1 belonging to partition Sk
    {
        trova_boundaries_diagonale(c,i,boundaries);

        for (j=boundaries[0];j<=boundaries[1];++j) // metto gli 1 sulla riga 1 nella partizione Sk
        {
            if (i<j)
            {
                pos_vect_matrix = MEWCP_convert_coords_ij_to_vector_matrix(j+1,i+1);

                sdp_constraints_matrix[k].index[pos_elem] = pos_vect_matrix;
                sdp_constraints_matrix[k].weight[pos_elem] = (double) 1;
                ++ pos_elem;

#if defined MEWCP_CONVERTER_DSDP_VERBOSE1

                printf("%d 1 %d %d %d\t\tIMC A: %d\n",k,i+1,j+1,1,k);
#endif

            }
        }
        for (j=boundaries[2];j<=boundaries[3];++j)  // I put 1 on the column 1 of partition Sk
        {
            if (i>j)
            {

                pos_vect_matrix = MEWCP_convert_coords_ij_to_vector_matrix(i+1,j+1);

                sdp_constraints_matrix[k].index[pos_elem] = pos_vect_matrix;
                sdp_constraints_matrix[k].weight[pos_elem] = (double) 1;
                ++ pos_elem;

#if defined MEWCP_CONVERTER_DSDP_VERBOSE1

                printf("%d 1 %d %d %d\t\tIMC A: %d\n",k,j+1,i+1,1,k);
#endif

            }
        }
        sdp_constraints_matrix[k].num_nz = pos_elem;
        pos_elem = 0;
        ++k;
    }

}


void MEWCP_write_constraints_improved_MC_B( constraint_t * sdp_constraints_matrix, const unsigned int vector_length, const unsigned int n, const unsigned int m, const unsigned int c,  unsigned int k)
{
#if defined MEWCP_CONVERTER_DSDP_VERBOSE1
    printf("\n* MEWCP_write_constraints_improved_MC_B *\n");
#endif

    /* Vincolo di Multiple Choice IMC B */
    //  caso in cui i non appartiene alla partizione k
    // MI creo (c*m)*(m-1) vincoli

    unsigned int i;
    unsigned int u,v;
    int boundaries_u[4];
    int boundaries_v[4];
    unsigned int pos_vect_matrix;
    unsigned int pos_elem = 0;  /* The position in the constraint array */

    for (v=0; v<n; ++v)
    {
        trova_boundaries_diagonale(c,v,boundaries_v);

        for (u=0;u<m;++u)
        {
            trova_boundaries_diagonale(c,c*u,boundaries_u);
            if( (boundaries_u[0] == boundaries_v[0]) &&  (boundaries_u[1] == boundaries_v[1]) ) // considero stessa partizione
            {
                continue;
            }

            // I set a -1 in poszition (v,v)

            pos_vect_matrix = MEWCP_convert_coords_ij_to_vector_matrix(v+1,v+1);

            sdp_constraints_matrix[k].index[pos_elem] = pos_vect_matrix;
            sdp_constraints_matrix[k].weight[pos_elem] = (double) -1;
            ++ pos_elem;


#if defined MEWCP_CONVERTER_DSDP_VERBOSE1

            printf("%d 1 %d %d %d\t\tIMC B: %d\n",k,v+1,v+1,-1,k);
#endif


            for(i=boundaries_u[0];i<=boundaries_u[1];++i)
            {
                if (v<i)
                {

                    pos_vect_matrix = MEWCP_convert_coords_ij_to_vector_matrix(i+1,v+1);

                    sdp_constraints_matrix[k].index[pos_elem] = pos_vect_matrix;
                    sdp_constraints_matrix[k].weight[pos_elem] = (double) 0.5;
                    ++ pos_elem;


#if defined MEWCP_CONVERTER_DSDP_VERBOSE1

                    printf("%d 1 %d %d %f\t\tIMC B: %d\n",k,v+1,i+1,0.5,k);
#endif

                }
            }
            for(i=boundaries_u[2];i<=boundaries_u[3];++i)
            {
                if (i<v)
                {
                    pos_vect_matrix = MEWCP_convert_coords_ij_to_vector_matrix(v+1,i+1);
                    sdp_constraints_matrix[k].index[pos_elem] = pos_vect_matrix;
                    sdp_constraints_matrix[k].weight[pos_elem] = (double) 0.5;
                    ++ pos_elem;


#if defined MEWCP_CONVERTER_DSDP_VERBOSE1

                    printf("%d 1 %d %d %f\t\tIMC B: %d\n",k,i+1,v+1,0.5,k);
#endif

                }
            }
            sdp_constraints_matrix[k].num_nz = pos_elem;
            pos_elem = 0;
            ++k;  // cambio matrice
        }
    }

}

void MEWCP_write_constraints_improved_MC_C( constraint_t * sdp_constraints_matrix, const unsigned int vector_length, const unsigned int m, const unsigned int c,  unsigned int k)
{
#if defined MEWCP_CONVERTER_DSDP_VERBOSE1
    printf("\n* MEWCP_write_constraints_improved_MC_C *\n");
#endif

    /* Vincolo di Multiple choice rafforzati IMC C */
    // Scorro le partizioni con la variabile v e considero gli elementi i,j

    unsigned int i,j,z;
    int boundaries[4];
    unsigned int pos_vect_matrix;
    unsigned int pos_elem = 0;  /* The position in the constraint array */

    for(i=0;i<m;++i)
    {
        trova_boundaries_diagonale(c,c*i,boundaries);

        for(j=boundaries[0];j<=boundaries[1];++j)
        {
            for (z=boundaries[2];z<=boundaries[3];++z)
            {
                if (z<j)
                {
                    pos_vect_matrix = MEWCP_convert_coords_ij_to_vector_matrix(j+1,z+1);
                    sdp_constraints_matrix[k].index[pos_elem] = pos_vect_matrix;
                    sdp_constraints_matrix[k].weight[pos_elem] = (double) 1;
                    ++ pos_elem;


#if defined MEWCP_CONVERTER_DSDP_VERBOSE1

                    printf("%d 1 %d %d %d\t\tIMC C: %d\n",k,z+1,j+1,1,k);
#endif

                }
            }
        }
        sdp_constraints_matrix[k].num_nz = pos_elem;
        pos_elem = 0;
        ++k;  // This way I create m contraints instead of 1, should be the same
    }

}


void MEWCP_write_constraints_4C2( constraint_t * sdp_constraints_matrix, const unsigned int vector_length, const unsigned int n, const unsigned int m, const unsigned int c,  unsigned int k)
{
#if defined MEWCP_CONVERTER_DSDP_VERBOSE1
    printf("\n* MEWCP_write_constraints_4C2 *\n");
#endif

    /* Scrivo i vincoli 4C''  (J - MI) = 0 */
    unsigned int i,j;
    int boundaries[4];
    unsigned int pos_vect_matrix;
    unsigned int pos_elem = 0;  /* The position in the constraint array */


    for (i=0;i<n;++i)
    {
        trova_boundaries_diagonale(c,i,boundaries);

        for (j=0;j<n;++j)
        {

            if (i == j)
            {

                pos_vect_matrix = MEWCP_convert_coords_ij_to_vector_matrix(i+1,j+1);
                sdp_constraints_matrix[k].index[pos_elem] = pos_vect_matrix;
                sdp_constraints_matrix[k].weight[pos_elem] = 1-(double) m;
                ++ pos_elem;




#if defined MEWCP_CONVERTER_DSDP_VERBOSE1

                printf("%d 1 %d %d %d\t\t4C'' (J-Im)=0: %d\n",k,j+1,i+1,1-m,k);
#endif

            }
            else if (  j > boundaries[1] ) // if (i,j) isn't in the diagonal box of (i,i)
            {

                pos_vect_matrix = MEWCP_convert_coords_ij_to_vector_matrix(j+1,i+1);

                sdp_constraints_matrix[k].index[pos_elem] = pos_vect_matrix;
                sdp_constraints_matrix[k].weight[pos_elem] = (double) 1;
                ++ pos_elem;


#if defined MEWCP_CONVERTER_DSDP_VERBOSE1

                printf("%d 1 %d %d %f\t\t4C'' (J-Im)=0: %d\n",k,i+1,j+1,1.0,k);
#endif

            }
        }
    }
    sdp_constraints_matrix[k].num_nz = pos_elem;
}

void MEWCP_write_constraints_4C3_A( constraint_t * sdp_constraints_matrix, const unsigned int vector_length, const unsigned int n, const unsigned int m, const unsigned int c,  unsigned int k)
{
#if defined MEWCP_CONVERTER_DSDP_VERBOSE1
    printf("\n* MEWCP_write_constraints_4C3_A *\n");
#endif

    unsigned int i,j;
    int boundaries[4];
    unsigned int pos_vect_matrix;
    unsigned int pos_elem = 0;  /* The position in the constraint array */

    for (i=0; i<n;++i)
    {
        trova_boundaries_diagonale(c,i,boundaries);

        // I write  (i,i)= 1-m
        pos_vect_matrix = MEWCP_convert_coords_ij_to_vector_matrix(i+1,i+1);

        sdp_constraints_matrix[k].index[pos_elem] = pos_vect_matrix;
        sdp_constraints_matrix[k].weight[pos_elem] = 1- (double ) m;
        ++ pos_elem;





#if defined MEWCP_CONVERTER_DSDP_VERBOSE1

        printf("%d 1 %d %d %d\t\t4C'''A\n",k,i+1,i+1,1-m);
#endif

        // I write in the same line in the others boxes  1/2
        for(j=0;j<n;++j)
        {
            if (((j<boundaries[0]) || (j>boundaries[1])) && i<=j)
            {
                pos_vect_matrix = MEWCP_convert_coords_ij_to_vector_matrix(j+1,i+1);
                sdp_constraints_matrix[k].index[pos_elem] = pos_vect_matrix;
                sdp_constraints_matrix[k].weight[pos_elem] = (double) 0.5;
                ++ pos_elem;




#if defined MEWCP_CONVERTER_DSDP_VERBOSE1

                printf("%d 1 %d %d %f\n",k,i+1,j+1,0.5);
#endif

            }
        }


        // I write on the rows of other boxes 1/2
        for(j=0;j<n;++j)
        {
            if (((j<boundaries[2]) || (j>boundaries[3])) && j<=i)
            {
                pos_vect_matrix = MEWCP_convert_coords_ij_to_vector_matrix(i+1,j+1);

                sdp_constraints_matrix[k].index[pos_elem] = pos_vect_matrix;
                sdp_constraints_matrix[k].weight[pos_elem] = (double) 0.5;
                ++ pos_elem;



#if defined MEWCP_CONVERTER_DSDP_VERBOSE1

                printf("%d 1 %d %d %f\n",k,j+1,i+1,0.5);
#endif

            }
        }



        sdp_constraints_matrix[k].num_nz = pos_elem;
        pos_elem = 0;
        ++k;
    }
    /* End constraints  4C'''A */

}

void MEWCP_write_constraints_4C3_B(constraint_t* sdp_constraints_matrix, const unsigned int vector_length, const unsigned int n, const unsigned int m, const unsigned int c,  unsigned int k)
{

#if defined MEWCP_CONVERTER_DSDP_VERBOSE1
    printf("\n* MEWCP_write_constraints_4C3_B *\n");
#endif

    unsigned int i,j;
    int boundaries[4];
    unsigned int pos_vect_matrix;
    unsigned int pos_elem = 0;  /* The position in the constraint array */


    for (i=0; i<n;++i)
    {
        trova_boundaries_diagonale(c,i,boundaries);

        // I write (i,i)= m
        pos_vect_matrix = MEWCP_convert_coords_ij_to_vector_matrix(i+1,i+1);

        sdp_constraints_matrix[k].index[pos_elem] = pos_vect_matrix;
        sdp_constraints_matrix[k].weight[pos_elem] = (double) m;
        ++ pos_elem;



#if defined MEWCP_CONVERTER_DSDP_VERBOSE1

        printf("%d 1 %d %d %d\t\t4C'''B\n",k,i,i,m);
#endif

        // I write 1 on the main diagonal save in pos (i,i)
        for (j=0;j<n;++j)
        {
            if ( j != i )
                //if ( (j<boundaries[0]) || (j>boundaries[1]))
            {
                pos_vect_matrix = MEWCP_convert_coords_ij_to_vector_matrix(j+1,j+1);

                sdp_constraints_matrix[k].index[pos_elem] = pos_vect_matrix;
                sdp_constraints_matrix[k].weight[pos_elem] = (double) 1;
                ++ pos_elem;




#if defined MEWCP_CONVERTER_DSDP_VERBOSE1

                printf("%d 1 %d %d %d\n",k,j+1,j+1,1);
#endif

            }
        }

        // I write 1/2 on the same row of other boxes
        for(j=0;j<n;++j)
        {
            if (((j<boundaries[0]) || (j>boundaries[1])) && i<=j)
            {
                pos_vect_matrix = MEWCP_convert_coords_ij_to_vector_matrix(j+1,i+1);

                sdp_constraints_matrix[k].index[pos_elem] = pos_vect_matrix;
                sdp_constraints_matrix[k].weight[pos_elem] = (double) -0.5;
                ++ pos_elem;



#if defined MEWCP_CONVERTER_DSDP_VERBOSE1

                printf("%d 1 %d %d %f\n",k,i+1,j+1,-0.5);
#endif

            }
        }

        // Scrivo sulla riga verticale delle altre sk -1/2
        for(j=0;j<n;++j)
        {
            if (((j<boundaries[2]) || (j>boundaries[3])) && j<=i)
            {
                pos_vect_matrix = MEWCP_convert_coords_ij_to_vector_matrix(i+1,j+1);

                sdp_constraints_matrix[k].index[pos_elem] = pos_vect_matrix;
                sdp_constraints_matrix[k].weight[pos_elem] = (double) -0.5;
                ++ pos_elem;




#if defined MEWCP_CONVERTER_DSDP_VERBOSE1

                printf("%d 1 %d %d %f\n",k,j+1,i+1,-0.5);
#endif

            }
        }
        sdp_constraints_matrix[k].num_nz = pos_elem;
        pos_elem = 0;
        ++k;
    }

    /* End of constraints 4c''' (B) */

}


/*************************************************************************
 * 
 * filling vector bi 
 * 
 *************************************************************************/

void MEWCP_write_b_cardinality(double * bi, const int ishift, const unsigned m)
{

    bi[ishift] = (double) m;

#if defined MEWCP_CONVERTER_DSDP_VERBOSE1

    printf("%d ",m);
#endif

}

void MEWCP_write_b_simple_MC(double * bi, const int ishift, const unsigned n, const unsigned m)
{
    unsigned int i;

    /* b simple multiple choice  */
    for (i=0; i < m ; ++i)
    {
        bi[ishift + i] = 1;

#if defined MEWCP_CONVERTER_DSDP_VERBOSE1

        printf("1");
        printf(" ");
#endif

    }
}


void MEWCP_write_b_improved_MC_A(double * bi, const int ishift, const unsigned n)
{
    /* Improved Multiple choice A */

    unsigned int i;

    for (i=0; i < n ; ++i)
    {
        bi[ishift + i] = 0;

#if defined MEWCP_CONVERTER_DSDP_VERBOSE1

        printf("0");
        printf(" ");
#endif

    }
}


void MEWCP_write_b_improved_MC_B(double * bi, const int ishift, const unsigned n, const unsigned m)
{
    /* Improved Multiple choice B */

    unsigned int i;

    for (i=0; i < (n*(m-1)) ; ++i)
    {
        bi[ishift + i] = 0;


#if defined MEWCP_CONVERTER_DSDP_VERBOSE1

        printf("0");
        printf(" ");
#endif

    }
}

void MEWCP_write_b_improved_MC_C(double * bi, const int ishift,const unsigned int m)
{
    unsigned int i;

    /* Improved Multiple choice C */
    for (i=0; i < m ; ++i)
    {
        bi[ishift + i] = 0;

#if defined MEWCP_CONVERTER_DSDP_VERBOSE1

        printf("0");
        printf(" ");
#endif

    }
}

void MEWCP_write_b_4C2(double * bi, const int ishift)
{
    /* b constraints 4C''  (J-mI)=0 */

    bi[ishift] = 0;



#if defined MEWCP_CONVERTER_DSDP_VERBOSE1

    printf("0");
    printf(" ");
#endif
}

void MEWCP_write_b_4C3_A(double * bi, const int ishift, const unsigned n)
{
    /* b 4c'''A */

    unsigned int i;

    for (i=0; i < n ; ++i)
    {
        bi[ishift +i] = 0;

#if defined MEWCP_CONVERTER_DSDP_VERBOSE1

        printf("0");
        printf(" ");
#endif

    }
}

void MEWCP_write_b_4C3_B(double * bi, const int ishift,const unsigned n, const unsigned m)
{
    /* b 4c'''B */

    unsigned int i;

    for (i=0; i < n ; ++i)
    {

        bi[ishift +i] = m;

#if defined MEWCP_CONVERTER_DSDP_VERBOSE1

        printf("%d",m);
        printf(" ");
#endif

    }
}

void MEWCP_write_b_branching(double * bi, const int ishift)
{
    /* b branching constraints */
    bi[ishift] = (double) 0;

}


