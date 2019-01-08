#ifndef MEWCP_EXPLICIT_ENUMERATION_H_
#define MEWCP_EXPLICIT_ENUMERATION_H_
/*****************************************************************************
 *  Project: Maximum Edge Weighted Clique Problem with multiple choice contraints
 * 
 *  Authors:
 *  (c) 2009 Yari Melzani (yari.melzani@gmail.com) 
 *
 ****************************************************************************/

#include "MEWCP_dsdp.h"

/*
* Definitions and Data structures 
*/

/*
 * This struct stands for explicit enumeration of solutions
 * For each partition i keep:
 * - the nodes not blocked, thus free
 * - the number of free nodes
 * - the position cursor that help the convolution.
 */
typedef struct element_free_variables_partitions_s
{
    int * free_nodes;
    int  number_free_nodes;
    int  position_free_nodes;

}
list_free_variables_partitions_t;


/*
 * Prototypes 
 */
void MEWCP_bound_explicit(open_node_t * open_node, matrix_weights_t * matrix_weights, 
							const unsigned num_partitions, 
							const unsigned int cardinality_partitions);

void MEWCP_generate_list_free_varibles_partitions( list_blocked_nodes_t * list_blocked_nodes,
        list_free_variables_partitions_t *  list_free_variables_partitions,
        const unsigned int num_partitions,
        const unsigned int cardinality_partitions);


void MEWCP_compute_explicit_enumeration(list_blocked_nodes_t * list_blocked_nodes,
                                        matrix_weights_t * matrix_weights,
                                        const unsigned num_partitions,
                                        const unsigned int cardinality_partitions,
                                        int * out_best_list_nodes_solution,
                                        double * out_z_best_solution);

bool MEWCP_is_node_little_enough(list_blocked_nodes_t * list_blocked_nodes, 
								const unsigned int num_partitions,
								const unsigned int cardinality_partitions,
								const unsigned int max_solutions);
/*
 * ALLOCATION FUNCTIONS 
 */

list_free_variables_partitions_t * MEWCP_allocate_list_free_varibles_partitions(const unsigned num_partitions, 
																				const unsigned int cardinality_partitions );

/*
 * FREE FUNCTIONS 
 */

void MEWCP_free_list_free_varibles_partitions(list_free_variables_partitions_t *  list_free_variables_partitions,
        const unsigned int num_partitions);





#endif /*MEWCP_EXPLICIT_ENUMERATION_H_*/
