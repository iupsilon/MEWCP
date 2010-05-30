/*****************************************************************************
 *  Project: Maximum Edge Weighted Clique Problem
 * 
 *  Authors:
 *  (c) 2009 Yari Melzani (yari.melzani@gmail.com) 
 *
 ****************************************************************************/

#ifndef MEWCP_TABU_DEFINITIONS_H_
#define MEWCP_TABU_DEFINITIONS_H_


#include "stdbool.h"
/*******************************************************************
 * 		DEFINITIONS
 * *****************************************************************/


#define MEWCP_DBL_EPSILON 1E-6
#define MEWCP_MAX_NEG_WEIGHT -1E6

typedef int weight_t;   /* Edges and vertex costs are integer */


#define NULL_POINTER -1
typedef int pointer_node_t;
typedef unsigned int iteration_t;

/*******************************************************************
 * 		STRUCTURES
 ******************************************************************/

typedef struct matrix_weights_s
{
    unsigned int n;  // Number of vertex
    unsigned int m;  // Number of partitions
    unsigned int c;  // Cardinality of each partition
    weight_t ** weight;
}
matrix_weights_t;

typedef struct node_s
{
    bool belongsM; // true if node belongs to set of solution
    pointer_node_t next; // next element
    pointer_node_t prev; // prev element
    weight_t sum_din; // somma distanze del nodo con tutti gli altri M
    weight_t sum_dout; // somma distanze del nodo con gli altri N non M
}
node_t;

typedef struct node_list_s
{
    node_t * node;  // Vector of nodes
	
	pointer_node_t * selected_node_partition;  /* Tells for each partition which is the node taken */
	
    // NULL cursor is equal to -1
    pointer_node_t N_head; // points to the head of N
    pointer_node_t M_head; // points to the head of M (nodes belonging to solution)
    pointer_node_t N_tail; // points to the tail N
    pointer_node_t M_tail; // points to the tail fo M (nodes belonging to the solution)

    weight_t Z; // objective function corresponding to the current solution (set M)

    unsigned int card_N; // cardinality set N
    unsigned int card_M; // cardinality set M
}
node_list_t;


typedef struct solution_s
{
	weight_t Z;	/* The value of OBJ */
	pointer_node_t * node_solution;   /* Keeps the list of solution id nodes */
} solution_t;



#endif /*MEWCP_DEFINITIONS_H_*/
