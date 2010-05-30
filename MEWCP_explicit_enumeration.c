/*****************************************************************************
 *  Project: Maximum Edge Weighted Clique Problem with multiple choice contraints
 * 
 *  Authors:
 *  (c) 2009 Yari Melzani (yari.melzani@gmail.com) 
 *
 ****************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <stdbool.h>
#include <memory.h>
#include <time.h>
#include <sys/times.h>

#include "MEWCP_explicit_enumeration.h"
#include "MEWCP_dsdp.h"




/*
 * Prototypes 
 */


void MEWCP_bound_explicit(open_node_t * open_node, matrix_weights_t * matrix_weights,
                          const unsigned num_partitions,
                          const unsigned int cardinality_partitions)

{

#if defined MEWCP_EXPLICIT_ENUMERATION_DEBUG
    printf("* MEWCP_bound_explicit: %d *\n",open_node->serial_node);
#endif

    double z_best;
    open_node->list_nodes_solution = MEWCP_allocate_list_nodes_solution(num_partitions);

    MEWCP_compute_explicit_enumeration(open_node->list_blocked_nodes,
                                       matrix_weights,
                                       num_partitions,
                                       cardinality_partitions,
                                       open_node->list_nodes_solution,
                                       &z_best);

    open_node->PB = z_best;
    open_node->DB = z_best;
#if defined MEWCP_BOUNDING_VERBOSE1

    printf("<< EXPLICIT Enumeration: (%d) (P) Z: %.2lf \t level: %u\n",open_node->serial_node, z_best, open_node->depth_level);
#endif

}


void MEWCP_compute_explicit_enumeration(list_blocked_nodes_t * list_blocked_nodes,
                                        matrix_weights_t * matrix_weights,
                                        const unsigned num_partitions,
                                        const unsigned int cardinality_partitions,
                                        int * out_best_list_nodes_solution,
                                        double * out_z_best_solution)
{
#if defined MEWCP_EXPLICIT_ENUMERATION_DEBUG
    printf("* MEWCP_compute_explicit_enumeration *\n");
#endif

    int cin;
    unsigned int i,k;
    double z_best;
    double z_ith;
    unsigned long long int num_solutions;

    unsigned int pos_k[num_partitions];
    int ithsolution[num_partitions];

#if defined MEWCP_EXPLICIT_ENUMERATION_DEBUG

    printf("Cardinality: %d\n",cardinality_partitions);
#endif

    list_free_variables_partitions_t * list_free_variables_partitions;


    list_free_variables_partitions = MEWCP_allocate_list_free_varibles_partitions(num_partitions,cardinality_partitions);
    MEWCP_generate_list_free_varibles_partitions(list_blocked_nodes ,list_free_variables_partitions,num_partitions,cardinality_partitions);

    /* Let's count the solutions */
    num_solutions = list_free_variables_partitions[0].number_free_nodes;
    for(i=1;i<num_partitions;++i)
    {
        num_solutions *= list_free_variables_partitions[i].number_free_nodes;
    }

#if defined MEWCP_EXPLICIT_ENUMERATION_DEBUG
    printf("Number of soluzions: %lld\n",num_solutions);
#endif

    z_best = MEWCP_MIN_DOUBLE;

    // Initializes the position
    for (i =0; i< num_partitions;++i)
    {

        pos_k[i] = 0;
    }


    do
    {

        for(i=0;i<num_partitions;++i)
        {
            ithsolution[i] = list_free_variables_partitions[i].free_nodes[pos_k[i] ];
            //printf("%d ",ithsolution[i]);

        }
        z_ith = MEWCP_evaluate_list_nodes_solution(ithsolution,matrix_weights,num_partitions);


        //printf(" z: %.2lf\n",z_ith);



        if  (( z_ith - z_best) > MEWCP_EPSILON)
        {
            z_best = z_ith;
            MEWCP_clone_list_nodes_solution(ithsolution,out_best_list_nodes_solution,num_partitions);
        }


        cin = -1;
        do
        {
            cin += 1;
            //printf("cin: %d ",cin);
            pos_k[cin] = (pos_k[cin] + 1) % list_free_variables_partitions[cin].number_free_nodes;

        }
        while ((pos_k[cin ] == 0) && (cin<num_partitions-1) );


    }
    while( !( (cin == num_partitions - 1)  && (pos_k[num_partitions-1] == 0 ) ));




#if defined MEWCP_EXPLICIT_ENUMERATION_DEBUG
    /* Now let's see what is in the list of free varoiables */
    for(k=0; k<num_partitions;++k)
    {
        list_free_variables_partitions[k].position_free_nodes = 0;
        printf("Partition %d: ",k);
        for(i=0; i<list_free_variables_partitions[k].number_free_nodes; ++i)
        {
            printf("%d ",list_free_variables_partitions[k].free_nodes[i]);
        }
        printf("\n");

    }
#endif

    *out_z_best_solution = z_best;
    MEWCP_free_list_free_varibles_partitions(list_free_variables_partitions,num_partitions);


}


void MEWCP_generate_list_free_varibles_partitions( list_blocked_nodes_t * list_blocked_nodes,
        list_free_variables_partitions_t *  list_free_variables_partitions,
        const unsigned int num_partitions,
        const unsigned int cardinality_partitions)
{

#if defined MEWCP_EXPLICIT_ENUMERATION_DEBUG
    printf("* MEWCP_generate_list_free_varibles_partitions *\n");
#endif

    unsigned int i,k;
    unsigned int pos_k;

    /* starting from the first partitions  k=0,.., m-1 */
    for(k=0;k<num_partitions;++k)
    {
        list_free_variables_partitions[k].position_free_nodes = 0;
        list_free_variables_partitions[k].number_free_nodes = 0;

        for(i=k*cardinality_partitions ; i<k*cardinality_partitions + cardinality_partitions; ++i)
        {
            /* I' working on the k-th element of the list_free_variables */

            if(list_blocked_nodes->bool_list[i] == false )
            {
                pos_k = list_free_variables_partitions[k].position_free_nodes;
                list_free_variables_partitions[k].free_nodes[pos_k] = i;
                list_free_variables_partitions[k].position_free_nodes +=1;
                list_free_variables_partitions[k].number_free_nodes +=1;
            }
        }
    }


}


bool MEWCP_is_node_little_enough(list_blocked_nodes_t * list_blocked_nodes,
                                 const unsigned int num_partitions,
                                 const unsigned int cardinality_partitions,
                                 const unsigned int max_solutions)
{


    unsigned int dimension[num_partitions];
    unsigned int k,i;
    long long unsigned int number_solutions = 0;

    for(k=0;k<num_partitions;++k)
    {
        dimension[k] = 0; // initialize

        for(i=k*cardinality_partitions ; i<k*cardinality_partitions + cardinality_partitions; ++i)
        {
            /* I' working on the k-th element of the list_free_variables */

            if(list_blocked_nodes->bool_list[i] == false )
            {

                dimension[k] += 1;

            }
        }
    }

    number_solutions = dimension[0];
    for(i=1;i<num_partitions;++i)
    {
        number_solutions *= dimension[i];
        // if block that avoid overflow
        if (number_solutions > max_solutions )
        {
#if defined MEWCP_EXPLICIT_ENUMERATION_DEBUG
            printf("number_solutions: %llu\n",number_solutions);
#endif

            return false;
        }
        // end block if
    }

#if defined MEWCP_EXPLICIT_ENUMERATION_VERBOSE1
    printf("\t\t\t# Number_solutions: %llu\n",number_solutions);
#endif

    if (number_solutions > max_solutions )
    {
        return false;
    }
    else
    {
        return true;
    }

}




/*
 * ALLOCATION FUNCTIONS 
 */

list_free_variables_partitions_t * MEWCP_allocate_list_free_varibles_partitions(const unsigned num_partitions, const unsigned int cardinality_partitions )
{
#if defined MEWCP_EXPLICIT_ENUMERATION_DEBUG
    printf("* MEWCP_allocate_list_free_varibles_partitions *\n");
#endif

    unsigned int i;
    list_free_variables_partitions_t * list_free_variables_partitions;

    list_free_variables_partitions = (list_free_variables_partitions_t *) calloc(num_partitions, sizeof(list_free_variables_partitions_t));

    if ( list_free_variables_partitions == NULL)
    {
        printf("*******  Allocation MEWCP_allocate_list_free_varibles_partitions FAILED!\n");
        exit(EXIT_FAILURE);
    }

    /* I also can allocate for each partition the maximum number of free variables */
    for(i=0; i < num_partitions; ++i )
    {
        list_free_variables_partitions[i].free_nodes = (int *) calloc(cardinality_partitions, sizeof(int));
        list_free_variables_partitions[i].number_free_nodes = 0;
        list_free_variables_partitions[i].position_free_nodes = 0;
    }
    return list_free_variables_partitions;
}

/*
 * FREE FUNCTIONS 
 */

void MEWCP_free_list_free_varibles_partitions(list_free_variables_partitions_t *  list_free_variables_partitions,
        const unsigned int num_partitions)
{
#if defined MEWCP_EXPLICIT_ENUMERATION_DEBUG
    printf("* MEWCP_free_list_free_varibles_partitions *\n");
#endif

    unsigned int i;

    for(i=0; i < num_partitions; ++i)
    {
        free(list_free_variables_partitions[i].free_nodes );
    }
    free(list_free_variables_partitions);

}



