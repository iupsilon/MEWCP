#ifndef MEWCP_DSDP_H_
#define MEWCP_DSDP_H_

/*****************************************************************************
 *  Project: Maximum Edge Weighted Clique Problem
 * 
 *  Authors:
 *  (c) 2009 Yari Melzani (yari.melzani@gmail.com) 
 *
 ****************************************************************************/

#include "dsdp/dsdp5.h"
#include "MEWCP_tabu_definitions.h"



/* LOG DEFINITION */
//#define ASSERT

//#define MEWCP_CONVERTER_DSDP_VERBOSE1
//#define MEWCP_CONVERTER_DSDP_VERBOSE2
//#define MEWCP_CONVERTER_DSDP_DEBUG



//#define MEWCP_DSDP_DEBUG
#define MEWCP_DSDP_VERBOSE1
#define MEWCP_DSDP_VERBOSE2


/* Verbosing for Branching */
//#define MEWCP_BRANCHING_DEBUG


/* Verbosing for Bounding */
//#define MEWCP_BOUNDING_DEBUG
//#define MEWCP_BOUNDING_VERBOSE2
#define MEWCP_BOUNDING_VERBOSE1

/* DSDP PARAMETERS */

#define MEWCP_GAP_TOLERANCE 0.001
#define MEWCP_POTENTIAL_PARAMETER 5
#define MEWCP_REUSE_MATRIX 1
#define MEWCP_SET_PNORM_TOLERANCE 1.0
#define MEWCP_ALPHA 1.0

#define MEWCP_EPSILON 10E-4
#define MEWCP_MIN_DOUBLE -10E7 /* Is the (double) -infinity */


/*
 * DSDP Data structures 
 */

typedef struct solution_bb_s
{
    double z_opt;
    unsigned int number_explored_nodes;
    unsigned int id_node_best_primal;
    int * list_nodes_best_solution;
    
    /* informations at the root node */
    double DB_root;
    double PB_root_bestK;
    double gap_root;
    double timestamp_user_time_root;
    double timestamp_system_time_root;
}
solution_bb_t;

#define NUM_BLOCKS 1


typedef struct constraint_s
{
    int * index;
    double * weight;
    int num_nz;
}
constraint_t;

typedef struct list_blocked_node_s
{
    unsigned int * blocked_node;
    bool * bool_list;	/* Is the same as the blocked_node but is rapresented as a vector of bools */
    unsigned int  num_blocked_nodes;
}
list_blocked_nodes_t;


typedef struct open_node_s
{

    int id_node;	/* if i is father's id, id_node is left:(2*i +1) right:(2*i +2) */
    double DB;		/* my dual bound */
    double PB;		/* my primal bound */
    list_blocked_nodes_t * list_blocked_nodes;
    constraint_t * vect_mat_branching_contraint;
    double * vect_y;

    /* For rounding i need the fractional value of sdp,
     * I save the solution nodes of rounding into the list_nodes_solution 
     */
    double * diagX;
    int * list_nodes_solution;
}
open_node_t;


typedef struct branching_open_node_s
{
    open_node_t * open_node;

    struct branching_open_node_s * next_branching_node;
    struct branching_open_node_s * prev_branching_node;

}
branching_open_node_t;

typedef struct list_branching_s
{
    branching_open_node_t * branching_open_node_t;

    double best_primal;
    unsigned int number_open_nodes;
    unsigned int number_explored_nodes;
    unsigned int id_node_best_primal;
    int * list_nodes_best_solution;

    branching_open_node_t * head;
    branching_open_node_t * tail;

}
list_branching_t;



/*
 * BRANCHING FUNCTIONS
 */

/* returns false if all partitions have only one fractional value */
bool MEWCP_generate_equi_branch_node(double * diag_X, const unsigned int n, const unsigned int m,  int * out_num_part,  int * out_id_node);

/*
 * Generate a new list of blocked nodes based on an old blocked list end equibranch decision 
 * branch_num_part is the number of partition 0,...,m-1
 * brach_id_node is the id of the node 0,...,n-1
 * father_blocked_nodes is the given list of blocked nodes (associated to the father)
 * left_son_blocked_nodes is the list of blocked nodes for the left son in the branching tree
 * right_son_blocked_nodes is the list of blocked nodes for the left son in the branching tree
 */
void MEWCP_generate_list_blocked_nodes_branching_sons(list_blocked_nodes_t * father_blocked_nodes,
        list_blocked_nodes_t * left_son_blocked_nodes,
        list_blocked_nodes_t * right_son_blocked_nodes,
        const int branch_num_part,
        const int branch_id_node,
        const unsigned int num_nodes,
        const unsigned int cardinality_partition  );

/* Execute sd to the open_node */
void MEWCP_bound(open_node_t * open_node, constraint_t * constraints_matrix,matrix_weights_t * matrix_weigths,double * bi,
                 const unsigned int num_constraints,
                 const unsigned int dim_matrix,
                 const unsigned int num_nodes,
                 const unsigned int num_partitions);


bool MEWCP_branch( open_node_t * open_node,
                   const unsigned int dim_matrix,
                   const unsigned int num_nodes,
                   const unsigned int num_partitions,
                   const unsigned int num_constraints,
                   open_node_t ** out_left_son,
                   open_node_t ** out_right_son);



solution_bb_t * MEWCP_branch_and_bound(open_node_t * open_root_node, constraint_t * constraints_matrix,  matrix_weights_t * matrix_weigths ,double * bi,
                            const unsigned int num_constraints,
                            const unsigned int dim_matrix,
                            const unsigned int num_nodes,
                            const unsigned int num_partitions,
                            double best_primal_obj,
                            int * list_node_best_solution);


void MEWCP_close_open_node(list_branching_t * list_branching, open_node_t * open_node);

/*
 * BRANCHING LIST FUNCIONS 
 */

/* It just allocates a struct and set pointers to NULL */
branching_open_node_t * MEWCP_allocate_branching_open_node(void);
list_branching_t * MEWCP_allocate_list_branching(void);
void MEWCP_free_branching_open_node(branching_open_node_t * branching_open_node);
void MEWCP_free_list_branching(list_branching_t * list_branching);

/* PUSH and POP: the list is a queue FIFO, push to the head, pop from the tail */
void MEWCP_push_open_node(open_node_t * open_node, list_branching_t * list_branching);
open_node_t * MEWCP_pop_open_node( list_branching_t * list_branching);
open_node_t * MEWCP_pop_specific_open_node(branching_open_node_t * branching_open_node , list_branching_t * list_branching);

bool MEWCP_is_list_branching_empty(list_branching_t * list_branching);
branching_open_node_t * MEWCP_find_worst_bound_element(list_branching_t * list_branching );




/*
 * UTILS
 */
// Trova i limiti che determinano il quadrato sulla diagonale del box in cui si trova il punto (i,i)
// Scrive in OUT in boundaries
/* It simply add a node to a list of blocked nodes */
void MEWCP_add_blocked_node(const unsigned int id_node, list_blocked_nodes_t * list_blocked_nodes);
void trova_boundaries_diagonale(const unsigned int c, const unsigned int i, int * boundaries);
unsigned int MEWCP_convert_coords_ij_to_vector_matrix(const unsigned int i, const unsigned j);
void MEWCP_clone_list_blocked_modes(list_blocked_nodes_t * list_to_be_cloned, list_blocked_nodes_t * list_cloned, const unsigned int num_nodes);
void MEWCP_clone_vect_y(double * vect_y_to_be_cloned, double * vect_y_cloned, const unsigned int num_constraints);
void MEWCP_clone_list_nodes_solution( int * list_to_be_cloned,  int * list_cloned, const unsigned int num_nodes);
void MEWCP_dump_diag_X(SDPCone * sdpcone, double * dst_diag_X, const unsigned int num_nodes);
void MEWCP_dump_vect_y(DSDP * dsdp, double * dst_vect_y, const unsigned int num_constraints);
void MEWCP_compute_sdp_rounding(double * diag_X,  int * list_nodes_rounded_solution, const unsigned int num_nodes, const unsigned int c);

/* Calculate the Objective function related to a list_node_soluztion, the weights matrix is given into a vector of n*(n+1) elements */
double MEWCP_evaluate_list_nodes_solution( int * list_node_solution, matrix_weights_t * matrix_weights ,const unsigned int m);
int sort_compare (const void * a, const void * b);
void Take_Time(double * user_time, double * system_time);


/* ALLOCATION FUNCTIONS */
double * MEWCP_allocate_vect_y(const unsigned int num_constraints);
constraint_t * MEWCP_allocate_sdp_constraints_matrix(const unsigned int n, const unsigned int num_constraints );
double * MEWCP_allocate_bi(const unsigned int num_contraints );
list_blocked_nodes_t *  MEWCP_allocate_list_blocked_nodes(const unsigned int num_nodes );
double * MEWCP_allocate_diag_X(const unsigned int length);
int * MEWCP_allocate_list_nodes_solution( const unsigned int m);
constraint_t * MEWCP_allocate_vect_mat_branching_constraints(const unsigned int dim_matrix);
open_node_t * MEWCP_allocate_open_node(void);
solution_bb_t * MEWCP_allocate_solution_bb(unsigned int num_partitions);

// Print functions
void MEWCP_print_contraints_matrix(double ** matrix, const unsigned int length_i, const unsigned int length_j);
void MEWCP_print_vectorY(double * vector, const unsigned int dim);
void MEWCP_print_diag_X(double * diag_X, const unsigned int num_nodes);
void MEWCP_print_list_blocked_nodes(list_blocked_nodes_t * list_blocked_modes);
void MEWCP_print_list_nodes_solution( int * list_nodes_solution, const unsigned int m);
void MEWCP_print_list_nodes_solution_cplex( int * list_nodes_solution, const unsigned int m);

/* FREE FUNCTIONS */
void MEWCP_free_list_blocked_nodes(list_blocked_nodes_t * list_blocked_nodes );
void MEWCP_free_sdp_constraints_matrix(constraint_t * sdp_constraints_matrix, const unsigned int num_constraints );
void MEWCP_free_bi(double * bi);
void MEWCP_free_diag_X(double * diag_X);
void MEWCP_free_vect_y(double * vect_y);
void MEWCP_free_list_nodes_solution( int * list_nodes_solution);
void MEWCP_free_vect_mat_branching_constraints(constraint_t * vect_mat_branching_contraints);
void MEWCP_free_open_node(open_node_t * open_node);
void MEWCP_free_solution_bb(solution_bb_t * solution_bb);


#endif /*MEWCP_DSDP_H_*/
