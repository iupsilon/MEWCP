/*****************************************************************************
 *  Project: Maximum Edge Weighted Clique Problem
 * 
 *  Authors:
 *  (c) 2009 Yari Melzani (yari.melzani@gmail.com) 
 *
 ****************************************************************************/

#ifndef MEWCP_TABU_H_
#define MEWCP_TABU_H_

#include "MEWCP_tabu_definitions.h"

/*******************************************************************
 * 		DEFINITIONS
 ******************************************************************/

//#define MEWCP_TABU_VERBOSE1
//#define MEWCP_TABU_VERBOSE2
//#define MEWCP_TABU_DEBUG1
//#define ASSERT


#define MAX_WORSENING_ITERATIONS 1000
#define TABU_IN_ITERATIONS 8	
#define TABU_OUT_ITERATIONS 1


/*******************************************************************
 * 		STRUCTURES
 ******************************************************************/

typedef struct list_iterations_s
{
    iteration_t current_iteration;
    unsigned int number_iterations_limit;         // Tells the iterations limit
    weight_t * Z_iteration;      // Per ogni iterazione mi dice quale è stato il valore Z trovato
}
list_iterations_t;

typedef struct tabu_node_state_s
{
    /* Tells the last iteration in whis the node has been involved into a move */
    iteration_t iteration_in;
    iteration_t iteration_out;
    unsigned int num_times_tabu; /* Number of times that node has not been taken because of tabu state */
}
tabu_node_state_t;


/* It is the tabu structure containing the overall tabu state */
typedef struct tabu_node_list_s
{
    tabu_node_state_t * tabu_node_state;
}
tabu_node_list_t;



typedef struct tabu_result_s
{
    iteration_t last_improvement_iteration;
    double last_improvement_time;
    solution_t solution;
}
tabu_result_t;


/*******************************************************************
 * 		PROTOTYPES
 ******************************************************************/



tabu_result_t MEWCP_compute_tabu_search(const unsigned int num_iterations, matrix_weights_t * matrix_weights, node_list_t * node_list);



weight_t MEWCP_compute_iteration_tabu_search(	matrix_weights_t * matrix_weights,
												node_list_t * node_list , 
												solution_t * best_solution, 
												list_iterations_t * list_iterations,
												tabu_node_list_t * tabu_node_list );
												
/* It gives the Z pretending the swap of [n1,n2] */
weight_t MEWCP_evaluate_swap_nodes(const pointer_node_t n1, const pointer_node_t n2, matrix_weights_t * matrix_weights, node_list_t * node_list);

void MEWCP_create_matrix_weights(const unsigned int n, const unsigned int m, matrix_weights_t * matrix_weights );
void MEWCP_load_AMPL_instance(char * filename_AMPL, matrix_weights_t * matrix_weights); /* loads an AMPL instance into a given empty matrix_weights */

void MEWCP_create_node_list(matrix_weights_t * matrix_weights, node_list_t * node_list); /* Initizlize the node_list */
void MEWCP_initialize_node_list(matrix_weights_t * matrix_weights, node_list_t * node_list);

void MEWCP_create_solution(matrix_weights_t * matrix_weights, solution_t * solution);

/* compute a starting solution taking the first element of each partition */
void MEWCP_compute_starting_solution(matrix_weights_t * matrix_weights, node_list_t * node_list);

/* Add node n to the set solution M in the node_list */
void MEWCP_add_node_to_solution_list(pointer_node_t n, matrix_weights_t * matrix_weights, node_list_t * node_list);
/* Remove node n from the set solution M in the node_list */
void MEWCP_remove_node_from_solution_list(pointer_node_t n, matrix_weights_t * matrix_weights, node_list_t * node_list);
/* Recalculate Z according to node_list */
weight_t MEWCP_update_Z_node_list(node_list_t * node_list, matrix_weights_t * matrix_weights);

/* Returns the partition number of node i */
pointer_node_t MEWCP_get_node_partition(pointer_node_t i, matrix_weights_t * matrix_weights);

/* Save the current solution from the node_list to solution */
void MEWCP_dump_tabu_solution(node_list_t * node_list, solution_t * solution);

/* clone a solution to an othere */
void MEWCP_clone_solution(solution_t * src_solution, solution_t * dest_solution, const unsigned int m);

/* Criteria to dermine if a node can go in or out the solution */
bool MEWCP_is_tabu_in(const pointer_node_t node, tabu_node_list_t * tabu_node_list, list_iterations_t * list_iterations);
bool MEWCP_is_tabu_out(const pointer_node_t node, tabu_node_list_t * tabu_node_list, list_iterations_t * list_iterations);



/******
 * Functions for the initializzation of tabu structures *
 ******/
/* Receives structs and allocate and initilize structures*/
void MEWCP_initialize_tabu_structures(	matrix_weights_t * matrix_weights,
                                       const unsigned int num_iterations_limit,
                                       list_iterations_t * list_iterations,
                                       solution_t * best_solution,
                                       tabu_node_list_t * tabu_node_list);

void MEWCP_create_iterations_list(const unsigned int num_iterations_limit, list_iterations_t * list_iterations);
void MEWCP_create_tabu_node_list(matrix_weights_t * matrix_weights, tabu_node_list_t * tabu_node_list);


/* Free functions */
void MEWCP_free_matrix_weights(matrix_weights_t *);
void MEWCP_free_node_list(node_list_t * node_list);
void MEWCP_free_solution(solution_t * solution);

void EWCP_free_iterations_list(list_iterations_t * list_iterations );
void EWCP_free_tabu_node_list(tabu_node_list_t * tabu_node_list);
void MEWCP_free_tabu_node_list(tabu_node_list_t * tabu_node_list);


/* Print functions */

void MEWCP_print_matrix_weights(matrix_weights_t *  matrix_weights);
void MEWCP_print_node_list(node_list_t * node_list);
void MEWCP_print_solution(matrix_weights_t * matrix_weights, solution_t * solution);


#endif /*MEWCP_TABU_H_*/
