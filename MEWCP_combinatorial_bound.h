#ifndef MEWCP_COMBINATORIAL_BOUND_H_
#define MEWCP_COMBINATORIAL_BOUND_H_

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

#include "MEWCP_combinatorial_bound.h"
#include "MEWCP_dsdp.h"

 
 /*
 * Prototypes 
 */
 
 void MEWCP_compute_combinatorial_preprocessing(open_node_t * open_node, matrix_weights_t * matrix_weights,
 												const unsigned int num_partitions,
 												const unsigned int cardinality_partitions,
 												const double initial_best_primal);
 
 /* given a node position and its partition generates the list of blocked nodes consisting of other nodes of the same partition*/
list_blocked_nodes_t * MEWCP_blocked_other_nodes_of_partition(const unsigned int pos, 
																const unsigned int partition, 
																const unsigned int num_partitions, 
																const unsigned int cardinality_partitions);

 
 void MEWCP_bound_combinatorial(open_node_t * open_node,
 								matrix_weights_t * matrix_weights,
 								const unsigned int num_partitions,
 								const unsigned int cardinality_partitions,
 								const double best_primal);
 								
 
 
 /*
 * ALLOCATION FUNCTIONS 
 */

#endif /*MEWCP_COMBINATORIAL_BOUND_H_*/
