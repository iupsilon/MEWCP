/*****************************************************************************
 *  Project: Maximum Edge Weighted Clique Problem
 * 
 *  Authors:
 *  (c) 2009 Yari Melzani (yari.melzani@gmail.com) 
 *
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <sys/times.h>


#include "dsdp/dsdp5.h"
#include "converter_dsdp.h"
#include "MEWCP_dsdp.h"
#include "MEWCP_tabu.h"
#include "MEWCP_tabu_definitions.h"


/*******************************************
 * PROTOTYPES
 ******************************************/

void show_usage(void);

int main (int argc, char * argv[])
{
    /* Variables to decide which constrains are enabled */
    bool c_cardinality;
    bool c_simple_MC;
    bool c_improved_MC_A;
    bool c_improved_MC_B;
    bool c_improved_MC_C;
    bool c_4C2;
    bool c_4C3_A;
    bool c_4C3_B;

    c_cardinality =         false;
    c_simple_MC =           true;
    c_improved_MC_A =       true;
    c_improved_MC_B =       false;
    c_improved_MC_C =       false;
    c_4C2 =                 true;
    c_4C3_A	 =              false;
    c_4C3_B =               false;


    char * filename_in;
    
    
    unsigned int dim_matrix;
    FILE * file_out;
    char * filename_out = NULL;

    unsigned int num_constraints;
    unsigned int num_blocks;
    unsigned int num_nodes;
    unsigned int num_partitions, cardinality_partition;

    double  * bi;
    constraint_t * vect_mat_branching_constraints;
    constraint_t * constraints_matrix;
    solution_bb_t * solution_bb;

    /* variables to take time */
    double t_user;
    double t_system;
    double t_user_start;
    double t_system_start;



    open_node_t * open_node;
    num_blocks = NUM_BLOCKS;


    /* Now I load the instance */
    iteration_t iterations;

    matrix_weights_t matrix_weights;
    node_list_t node_list;

    tabu_result_t tabu_result;

    /* Checking parameters number */
    if(argc != 4 && argc != 3)
    {
        show_usage();
        return EXIT_SUCCESS;
    }


    iterations = atoi(argv[1]);
    filename_in = argv[2];
    if (argc == 4)
    {
        
        filename_out = argv[3];
    }

    /* I take starting time */
    Take_Time(&t_user_start,&t_system_start);

    MEWCP_load_AMPL_instance(filename_in,&matrix_weights);

    num_nodes = matrix_weights.n;
    num_partitions  = matrix_weights.m;
    cardinality_partition = matrix_weights.c;

    MEWCP_create_node_list(&matrix_weights,&node_list);
    MEWCP_initialize_node_list(&matrix_weights,&node_list);
    MEWCP_compute_starting_solution(&matrix_weights,&node_list);

    tabu_result = MEWCP_compute_tabu_search(iterations,&matrix_weights,&node_list);


    //MEWCP_print_solution(&matrix_weights,&tabu_result.solution);
    //printf("%d %d\n",tabu_result.solution.Z,tabu_result.last_improvement_iteration);



    MEWCP_generate_sdp_constraints(&matrix_weights, &constraints_matrix,&bi, &vect_mat_branching_constraints,MEWCP_ALPHA,
                                   c_cardinality,
                                   c_simple_MC,c_improved_MC_A,
                                   c_improved_MC_B,c_improved_MC_C,
                                   c_4C2,
                                   c_4C3_A,
                                   c_4C3_B,
                                   &num_constraints);


    dim_matrix = num_nodes*(num_nodes +1)/2;



    open_node = MEWCP_allocate_open_node();

    open_node->id_node = 0;

    open_node->list_blocked_nodes = MEWCP_allocate_list_blocked_nodes(num_nodes);
    open_node->vect_mat_branching_contraint = vect_mat_branching_constraints;
    MEWCP_generate_constraints_branch(open_node->list_blocked_nodes, vect_mat_branching_constraints,dim_matrix,num_nodes, cardinality_partition);


    /* I simply allocate a 0 vect_y */
    open_node->vect_y = MEWCP_allocate_vect_y(num_constraints);
    MEWCP_clone_vect_y(bi,open_node->vect_y,num_constraints);


    solution_bb = MEWCP_branch_and_bound(open_node,constraints_matrix, &matrix_weights,bi,num_constraints,dim_matrix,num_nodes,num_partitions,tabu_result.solution.Z,tabu_result.solution.node_solution);




    /* I take final time */
    Take_Time(&t_user,&t_system);

#if defined MEWCP_DSDP_VERBOSE1


    printf("Z_opt: %.2lf\tr_best_PB: %.2lf\tr_DB: %.2lf\tr_gap: %.2lf %%\ttime_root: %.2lf\t Best_node: %u\t Exp_nodes: %u \t Time: %.2lf\n",
           solution_bb->z_opt,
           solution_bb->PB_root_bestK,
           solution_bb->DB_root,
           solution_bb->gap_root,
           (solution_bb->timestamp_user_time_root - t_user_start ) + ( solution_bb->timestamp_system_time_root - t_system_start),
           solution_bb->id_node_best_primal,
           solution_bb->number_explored_nodes,
           (t_user-t_user_start) + ( t_system-t_system_start) );


#endif

    if (argc == 4)
    {
        if ((file_out = fopen(filename_out,"ab")) == NULL)
        {
            show_usage();
            exit(EXIT_FAILURE);
        }
        fprintf(file_out,"Z_opt: %.2lf\tr_best_PB: %.2lf\tr_DB: %.2lf\tr_gap: %.2lf %%\ttime_root: %.2lf\t Best_node: %u\t Exp_nodes: %u \t Time: %.2lf\n",
                solution_bb->z_opt,
                solution_bb->PB_root_bestK,
                solution_bb->DB_root,
                solution_bb->gap_root,
                (solution_bb->timestamp_user_time_root - t_user_start ) + ( solution_bb->timestamp_system_time_root - t_system_start),
                solution_bb->id_node_best_primal,
                solution_bb->number_explored_nodes,
                (t_user-t_user_start) + ( t_system-t_system_start) );

        fclose(file_out);
    }

    MEWCP_free_matrix_weights(&matrix_weights);
    MEWCP_free_node_list(&node_list);
    MEWCP_free_solution_bb(solution_bb);

    return EXIT_SUCCESS;
}


void show_usage(void)
{
    printf(" ***************	*****************************************\n");
    printf(" *\t\t\t\t\t\t\t*\n");
    printf(" *  Project: Maximum Edge Weighted Clique Problem\t*\n");		
    printf(" *  Authors:\t\t\t\t\t\t*\n");
    printf(" *  (c) 2009 Yari Melzani (yari.melzani@gmail.com)\t*\n");
    printf(" *\t\t\t\t\t\t\t*\n");
    printf(" ********************************************************\n\n");
    printf("MEWCP Parameters:\n");
    printf("Parameters:\n");
    printf("\t\t1) number tabu search iterations\n");
    printf("\t\t2) instance in format .dat\n");
    printf("\t\t3) <output file (optional)>\n");
    printf("\n");


}
