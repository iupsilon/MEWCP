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

#include "converter_dsdp.h"
#include "MEWCP_explicit_enumeration.h"
#include "MEWCP_combinatorial_bound.h"
#include "MEWCP_tabu.h"
#include "MEWCP_dsdp.h"
#include "dsdp/dsdp5.h"


solution_bb_t * MEWCP_branch_and_bound(open_node_t * open_root_node, constraint_t * constraints_matrix, matrix_weights_t * matrix_weigths, double * bi,
                                       const unsigned int num_constraints,
                                       const unsigned int dim_matrix,
                                       const unsigned int num_nodes,
                                       const unsigned int num_partitions,
                                       double best_primal_obj,
                                       int * list_node_best_solution, double time_limit)
{

#if defined MEWCP_CONVERTER_DSDP_VERBOSE1
    printf("* MEWCP_branch_and_bound\n\tHere we goooooooooooooo!!!! *\n");
#endif

    unsigned int cardinality_partitions;
    bool new_best_PB_found;
    bool active_combinatorial_bound = false;

    list_branching_t * list_branching;

    /* Branching variables */
    open_node_t * son_left;
    open_node_t * son_right;
    branching_open_node_t * branching_open_worst_bound;
    open_node_t * open_node_worst_bound;
    bool possible_branch;

    float gap; /* is the current % gap */
    bool left_to_be_closed;
    bool right_to_be_closed;
    
    double time_tmp = 0;  //useful to take time

    cardinality_partitions = num_nodes/num_partitions;

    /* Root node bounded, ready to gather some result informations */
    solution_bb_t * solution_bb;
    solution_bb = MEWCP_allocate_solution_bb(num_partitions);
    
    list_branching = MEWCP_allocate_list_branching();
    list_branching->list_nodes_best_solution = MEWCP_allocate_list_nodes_solution(num_partitions);


    /* First I compute Tabu search in order to get a good Primal Bound */
    list_branching->best_primal = best_primal_obj;
    MEWCP_clone_list_nodes_solution(list_node_best_solution,list_branching->list_nodes_best_solution, num_partitions);


    /* Let's consider root node */

#if defined PREPROCESSING_ACTIVE
    /* COMBINATORIAL PREPROCESSING */
    MEWCP_compute_combinatorial_preprocessing(open_root_node,matrix_weigths,num_partitions,cardinality_partitions,list_branching->best_primal);

    /* Now the root node has been modified with a new list of blocked nodes and a related
     * brancing constraint 
     */
#endif // combinatorial preprocessing

    /* Explicit Enumeration */
    if ( MEWCP_is_node_little_enough(open_root_node->list_blocked_nodes,num_partitions,num_nodes/num_partitions,MEWCP_MAX_EXPLICIT_SOLUTIONS) == true)
    {
        MEWCP_bound_explicit(open_root_node,matrix_weigths,num_partitions,num_nodes/num_partitions);

        // I update the best valueif needed
        new_best_PB_found = MEWCP_is_new_best_PB_and_update(open_root_node,list_branching,num_partitions);

#if defined MEWCP_DSDP_VERBOSE1

        if( new_best_PB_found  == true)
        {
            printf("\t*****(Explicit enumeration)  Node: %d\tNew best PB: %.2lf\n",open_root_node->serial_node,open_root_node->PB );
        }
#endif

    }  // end explicit enumeration
    else
    {
#if defined COMBINATORIAL_BOUND_ACTIVE
        // Combinatorial BOUND
        
        time_tmp = get_cpu_time();  // I take the time in order to determine how long the combinatorial bounding takes
        MEWCP_bound_combinatorial(open_root_node,matrix_weigths,num_partitions,num_nodes/num_partitions,list_branching->best_primal);
        new_best_PB_found = MEWCP_is_new_best_PB_and_update(open_root_node,list_branching,num_partitions);
		
		// I set the root DB due to combinatorial
		solution_bb->DB_root_combinatorial = open_root_node->DB_comb;
		solution_bb->PB_root_combinatorial = open_root_node->PB_comb;
		solution_bb->time_DB_root_combinatorial = get_cpu_time() - time_tmp;
		

#if defined MEWCP_DSDP_VERBOSE1

        if( new_best_PB_found  == true)
        {
            printf("\t*****(Comb Bound)  Node: %d\tNew best PB: %.2lf\n",open_root_node->serial_node,open_root_node->PB );
        }
#endif
        //end combinatorial Bound
#endif // end combinatorial bound condition

        /* Semidefinite BOUND */
        MEWCP_bound(open_root_node,constraints_matrix,matrix_weigths, bi,num_constraints, dim_matrix, num_nodes,num_partitions, list_branching->best_primal);
        new_best_PB_found = MEWCP_is_new_best_PB_and_update(open_root_node,list_branching,num_partitions);
		
		// I set the root DB due to semidefinite
		solution_bb->DB_root_semidefinite = open_root_node->DB_SDP;
		solution_bb->PB_root_semidefinite = open_root_node->PB_SDP;
		solution_bb->time_DB_root_semidefinite = get_cpu_time() - solution_bb->time_DB_root_combinatorial;
		
#if defined MEWCP_DSDP_VERBOSE1

        if( new_best_PB_found  == true)
        {
            printf("\t*****(SDP Bound)  Node: %d\tNew best PB: %.2lf\n",open_root_node->serial_node,open_root_node->PB );
        }
#endif

        /* end semidefinite bound */
    }



	/* gathering some root informations */
    solution_bb->DB_root = open_root_node->DB;
    solution_bb->PB_root_bestK = list_branching->best_primal;
    solution_bb->gap_root = (open_root_node->DB - list_branching->best_primal)/list_branching->best_primal*100;
    Take_Time(&solution_bb->timestamp_user_time_root, &solution_bb->timestamp_system_time_root);


    /* Let's start! */

    MEWCP_push_open_node(open_root_node,list_branching);



    /* All done, I'm ready to start with branching procedure */


    while(MEWCP_is_list_branching_empty(list_branching) == false)
    {
    	#if defined ROOT_NODE_SIMULATION_ONLY
    	// I only want the root simulation... i break!
    	break;
    	#endif
    	
    	if ((get_cpu_time() - time_limit) > MEWCP_EPSILON)
    	{
    		//Time limit exceeded
    		#if defined MEWCP_DSDP_VERBOSE1
    		printf("\n**************\tTime limit reached!\t****************\n\n");
    		#endif
    		solution_bb->best_bound_left = open_root_node->DB;
    		break;
    	}


        branching_open_worst_bound = MEWCP_find_worst_bound_element(list_branching);
        open_node_worst_bound = MEWCP_pop_specific_open_node(branching_open_worst_bound,list_branching);


        if ( (open_node_worst_bound->DB - list_branching->best_primal) > MEWCP_EPSILON )
        {

            gap = (open_node_worst_bound->DB - list_branching->best_primal)/list_branching->best_primal * 100;
#if defined MEWCP_DSDP_VERBOSE1

            printf("\n++ Current DB(%d): %.2lf \t level: %u \tBest P(%d): %.2lf \t open_nodes: %d\t explored: %u \t gap: %.3f %%\n",open_node_worst_bound->serial_node,open_node_worst_bound->DB, open_node_worst_bound->depth_level, list_branching->serial_node_best_primal,list_branching->best_primal, list_branching->number_open_nodes, list_branching->number_explored_nodes,gap);
#endif

            possible_branch = MEWCP_branch(open_node_worst_bound,dim_matrix,num_nodes,num_partitions,num_constraints,open_node_worst_bound->depth_level, &list_branching->current_serial_number, &son_left,&son_right);



            if (possible_branch == true)
            {

                /* I decide wether the node has to be closed */
                left_to_be_closed = false;
                right_to_be_closed = false;

#if defined COMBINATORIAL_BOUND_ACTIVE
                /* I want to know if combinatorial bound was convinient on the father */
                if( ((open_node_worst_bound->DB_SDP - open_node_worst_bound->DB_comb) > MEWCP_EPSILON) && open_node_worst_bound->DB_comb != 0)
                {
                    active_combinatorial_bound = true;
#if defined  MEWCP_DSDP_VERBOSE2

                    printf("combinatorial ACTIVE\t DB_SDP= %.2lf\t DB_comb= %.2lf\n", open_node_worst_bound->DB_SDP,open_node_worst_bound->DB_comb);
#endif

                }
                else
                {
                    active_combinatorial_bound = false;
#if defined  MEWCP_DSDP_VERBOSE2

                    printf("combinatorial NOT active\t  DB_SDP= %.2lf\t DB_comb= %.2lf\n",open_node_worst_bound->DB_SDP,open_node_worst_bound->DB_comb);
#endif

                }

#endif // if combinatorial bound is active

                /* I check the depth */
                if ((open_node_worst_bound->depth_level +1) > list_branching->max_exploration_level )
                {
                    list_branching->max_exploration_level =  (open_node_worst_bound->depth_level +1);
                }

                /* I decide what type of bound use */
                if ( MEWCP_is_node_little_enough(son_left->list_blocked_nodes,num_partitions,cardinality_partitions,MEWCP_MAX_EXPLICIT_SOLUTIONS) == true)
                {
                    MEWCP_bound_explicit(son_left,matrix_weigths,num_partitions,cardinality_partitions);


                    // I update the best value if needed
                    new_best_PB_found = MEWCP_is_new_best_PB_and_update(son_left,list_branching,num_partitions);

#if defined MEWCP_DSDP_VERBOSE1

                    if( new_best_PB_found  == true)
                    {
                        printf("\t*****(Explicit enumeration)  Node: %d\tNew best PB: %.2lf\n",son_left->serial_node,son_left->PB );
                    }
#endif
                    left_to_be_closed = true;
                }
                else
                {
#if defined COMBINATORIAL_BOUND_ACTIVE
                    // Combinatorial Bound
                    if (active_combinatorial_bound == true)
                    {
                        MEWCP_bound_combinatorial(son_left,matrix_weigths,num_partitions,cardinality_partitions,list_branching->best_primal);
                        new_best_PB_found = MEWCP_is_new_best_PB_and_update(son_left,list_branching,num_partitions);

#if defined MEWCP_DSDP_VERBOSE1

                        if( new_best_PB_found  == true)
                        {
                            printf("\t*****(Comb Bound)  Node: %d\tNew best PB: %.2lf\n",son_left->serial_node,son_left->PB );
                        }
#endif



                        // end Combinatorial Bound


                        // check if left son has to be closed
                        if ( (list_branching->best_primal - son_left->DB ) > MEWCP_EPSILON)
                        {
                            left_to_be_closed = true;
                        }
                    }
#endif /* end combinatorial */

                    /* Semidefinite Bound */
                    if (left_to_be_closed == false)
                    {
                        MEWCP_bound(son_left,constraints_matrix,matrix_weigths, bi,num_constraints,dim_matrix,num_nodes,num_partitions,  list_branching->best_primal );
                        /* I check if PB is improved */
                        new_best_PB_found = MEWCP_is_new_best_PB_and_update(son_left,list_branching,num_partitions);

#if defined MEWCP_DSDP_VERBOSE1

                        if( new_best_PB_found  == true)
                        {
                            printf("\t*****(SDP Bound)  Node: %d\tNew best PB: %.2lf\n",son_left->serial_node,son_left->PB );
                        }
#endif

                    }
                    /* end semidefinite bound */

                    /* check if left son has to be closed */
                    if ( (list_branching->best_primal - son_left->DB ) > MEWCP_EPSILON)
                    {
                        left_to_be_closed = true;
                    }

                }

                /****** BOUND RIGHT NODE ********/

                if ( MEWCP_is_node_little_enough(son_right->list_blocked_nodes,num_partitions,num_nodes/num_partitions,MEWCP_MAX_EXPLICIT_SOLUTIONS) == true)
                {
                    MEWCP_bound_explicit(son_right,matrix_weigths,num_partitions,num_nodes/num_partitions);

                    // I update the best value if needed
                    new_best_PB_found = MEWCP_is_new_best_PB_and_update(son_right,list_branching,num_partitions);

#if defined MEWCP_DSDP_VERBOSE1

                    if( new_best_PB_found  == true)
                    {
                        printf("\t*****(Explicit enumeration)  Node: %d\tNew best PB: %.2lf\n",son_right->serial_node,son_right->PB );
                    }

#endif


                    right_to_be_closed = true;
                }
                else
                {
#if defined COMBINATORIAL_BOUND_ACTIVE
                    if (active_combinatorial_bound == true)
                    {
                        // Combinatorial Bound
                        MEWCP_bound_combinatorial(son_right,matrix_weigths,num_partitions,num_nodes/num_partitions,list_branching->best_primal);

                        new_best_PB_found = MEWCP_is_new_best_PB_and_update(son_right,list_branching,num_partitions);

#if defined MEWCP_DSDP_VERBOSE1

                        if( new_best_PB_found  == true)
                        {
                            printf("\t*****(Comb Bound)  Node: %d\tNew best PB: %.2lf\n",son_right->serial_node,son_right->PB );
                        }
#endif
                        // end Combinatorial Bound

                        // check if left son has to be closed
                        if ( (list_branching->best_primal - son_right->DB ) > MEWCP_EPSILON)
                        {
                            right_to_be_closed = true;
                        }
                    }
#endif /* end combinatorial */

                    /* Semidefinite Bound */

                    if (right_to_be_closed == false)
                    {
                        MEWCP_bound(son_right,constraints_matrix,matrix_weigths, bi,num_constraints,dim_matrix,num_nodes,num_partitions, list_branching->best_primal );

                        /* I check if PB is improved */
                        new_best_PB_found = MEWCP_is_new_best_PB_and_update(son_right,list_branching,num_partitions);

#if defined MEWCP_DSDP_VERBOSE1

                        if( new_best_PB_found  == true)
                        {
                            printf("\t*****(SDP Bound)  Node: %d\tNew best PB: %.2lf\n",son_right->serial_node,son_right->PB );
                        }
#endif

                    }
                    /* end semidefinite bound */

                    /* check if left son has to be closed */
                    if ( (list_branching->best_primal - son_right->DB ) > MEWCP_EPSILON)
                    {
                        right_to_be_closed = true;
                    }
                }

                /* let's close or add nodes to the open nodes */
                if ( left_to_be_closed == true)
                {
                    MEWCP_close_open_node(list_branching,son_left);
                }
                else
                {
                    MEWCP_push_open_node(son_left,list_branching);
                }


                if ( right_to_be_closed == true)
                {
                    MEWCP_close_open_node(list_branching,son_right);
                }
                else
                {
                    MEWCP_push_open_node(son_right,list_branching);
                }

            }
        }

        /* I can close the node! */

        MEWCP_close_open_node(list_branching,open_node_worst_bound);
    }


    /* Before I leave I'd let you know the best solution */
#if defined MEWCP_DSDP_VERBOSE2


    printf("\nBranch & Bound Z*: %.2lf\tExplored_nodes: %d\t Best node: %d\t\n",list_branching->best_primal,list_branching->number_explored_nodes, list_branching->serial_node_best_primal);
#endif
#if defined MEWCP_DSDP_VERBOSE2

    MEWCP_print_list_nodes_solution_cplex(list_branching->list_nodes_best_solution, num_partitions);
#endif

    /* Now let's fill the solution! */
    solution_bb->node_best_primal = list_branching->serial_node_best_primal;
    solution_bb->depth_best_primal = list_branching->depth_node_best_primal;
    solution_bb->max_exploration_depth = list_branching->max_exploration_level;
    solution_bb->number_explored_nodes = list_branching->number_explored_nodes;
    solution_bb->z_opt = list_branching->best_primal;
    MEWCP_clone_list_nodes_solution(list_branching->list_nodes_best_solution, solution_bb->list_nodes_best_solution, num_partitions);


    /* Freeing structures */

    MEWCP_free_list_branching(list_branching);

    return solution_bb;
}

void MEWCP_close_open_node(list_branching_t * list_branching, open_node_t * open_node)
{
#if defined MEWCP_DSDP_VERBOSE1
    printf("-- Close node: %d\n",open_node->serial_node);
#endif

    list_branching->number_explored_nodes += 1;
    MEWCP_free_open_node(open_node);
}

void MEWCP_bound(open_node_t * open_node, constraint_t * constraints_matrix,matrix_weights_t * matrix_weigths, double * bi,
                 const unsigned int num_constraints,
                 const unsigned int dim_matrix,
                 const unsigned int num_nodes,
                 const unsigned int num_partitions,
                 const double best_PB)
{

#if defined MEWCP_CONVERTER_DSDP_VERBOSE1
    printf("* MEWCP_bound *\n");
#endif


    unsigned int i;
    unsigned cardinality_partition;

    /* SDP structures */
    SDPCone sdpcone;
    DSDP dsdp;
    DSDPTerminationReason reason;
    int info;
    unsigned int num_blocks;
    double sol_traceX;	/* Is the trace of X */
    double pobj;   /* Is the value of sd relax */


    /* varibles for rounding */
    double z_rouded;

    num_blocks = NUM_BLOCKS;
    cardinality_partition = num_nodes/num_partitions;

    info = DSDPCreate(num_constraints,&dsdp);
    info = DSDPCreateSDPCone(dsdp,num_blocks,&sdpcone);
    info = SDPConeSetBlockSize(sdpcone, 0, num_nodes);  /* dimension of block is n */



    /* I feel the DSP structures */
    for (i=0;i<num_constraints;++i)
    {
        info = DSDPSetDualObjective(dsdp,i+1,bi[i]);
    }

    for (i=0; i < num_constraints +1; ++i) /* Matrix are number of constraints +1 because of matrix W */
    {
        if (i != num_constraints ) /* It's not the last constraint matrix */
        {
            SDPConeSetASparseVecMat(sdpcone, 0, i, num_nodes, MEWCP_ALPHA, 0, constraints_matrix[i].index, constraints_matrix[i].weight, constraints_matrix[i].num_nz);

        }
        else /* It's the last matrix containing branching contraints */
        {
            SDPConeSetASparseVecMat(sdpcone, 0, i, num_nodes, MEWCP_ALPHA, 0, open_node->vect_mat_branching_contraint->index, open_node->vect_mat_branching_contraint->weight, open_node->vect_mat_branching_contraint->num_nz);

        }

#if defined MEWCP_BOUNDING_DEBUG
        printf("%d Another contrain matrix:\n",i);
        SDPConeViewDataMatrix(sdpcone, 0, i);
#endif

    }

    /* set DSDP parameters */

    info=DSDPSetGapTolerance(dsdp,MEWCP_GAP_TOLERANCE);
    info=DSDPSetPotentialParameter(dsdp,MEWCP_POTENTIAL_PARAMETER);
    info=DSDPReuseMatrix(dsdp,MEWCP_REUSE_MATRIX);
    info=DSDPSetPNormTolerance(dsdp,MEWCP_SET_PNORM_TOLERANCE);
    //info = DSDPSetR0(dsdp,);


    /* I stop the computation when DD is greater to -best_PB */
    info = DSDPSetDualBound( dsdp, -best_PB);



#if defined MEWCP_BOUNDING_VERBOSE2

    DSDPSetStandardMonitor(dsdp, 1); /* verbose each iteration */
    DSDPLogInfoAllow(1,0);
#endif




    /* Now I set the initial values of the variables y in (D) */
    for (i=0; i< num_constraints; ++i)
    {
        DSDPSetY0(dsdp, i+1, open_node->vect_y[i]);

    }

    DSDPSetup(dsdp);
    DSDPSolve(dsdp);

    DSDPComputeX(dsdp);

    DSDPStopReason(dsdp, &reason);
    DSDPGetPObjective(dsdp, &pobj);
    DSDPGetTraceX(dsdp, &sol_traceX);

#if defined MEWCP_BOUNDING_DEBUG

    if ((reason != 1) && (reason != 5))
    {
        printf("DSDP Error!: %d\n",reason);
        getchar();
    }
#endif
    /* I take the negative pobj */
    pobj = -pobj;

#if defined MEWCP_BOUNDING_DEBUG

    double * sol_vect_X;
    int sol_dim_vect_X;


    SDPConeGetXArray(sdpcone, 0, &sol_vect_X, &sol_dim_vect_X);

    printf("\n");
    SDPConeViewX(sdpcone, 0, num_nodes, sol_vect_X, sol_dim_vect_X);
    printf("\n");

#endif


    /* Now I get the diagonal of solution X */

    open_node->diagX = MEWCP_allocate_diag_X(num_nodes);
    MEWCP_dump_diag_X(&sdpcone,open_node->diagX,num_nodes);

    /* Now I get the value of Y variables */
    open_node->vect_y = MEWCP_allocate_vect_y(num_constraints);
    MEWCP_dump_vect_y(&dsdp, open_node->vect_y, num_constraints);


#if defined MEWCP_BOUNDING_DEBUG

    MEWCP_print_diag_X(open_node->diagX,num_nodes);
    MEWCP_print_vectorY(open_node->vect_y,num_constraints);

#endif


    /* I compute the rounding of diagonan_X */

    int * list_nodes_rounded = MEWCP_allocate_list_nodes_solution(num_partitions);


    MEWCP_compute_sdp_rounding(open_node->diagX,list_nodes_rounded,num_nodes,cardinality_partition);
#if defined MEWCP_BOUNDING_VERBOSE2

    MEWCP_print_list_nodes_solution(list_nodes_rounded,num_partitions);
#endif

    z_rouded = MEWCP_evaluate_list_nodes_solution(list_nodes_rounded,matrix_weigths,num_partitions);



#if defined MEWCP_BOUNDING_DEBUG

    printf("Trace X: %.2lf\n",sol_traceX);
#endif

#if defined MEWCP_BOUNDING_VERBOSE1

    printf("(BB)  Bound SDP: (%d) DB: %.2lf \t Rounded: %.2lf \t level: %u \t Status: %d\n",open_node->serial_node, pobj,z_rouded, open_node->depth_level,reason);
#endif




    /* Update open_node */

    open_node->DB_SDP = pobj;
    open_node->PB_SDP = z_rouded;

    if ((open_node->DB - pobj) > MEWCP_EPSILON )
    {
        open_node->DB = pobj;
    }

    if ( (z_rouded - open_node->PB) > MEWCP_EPSILON )
    {
        open_node->PB = z_rouded;
        open_node->list_nodes_solution = MEWCP_allocate_list_nodes_solution(num_partitions);
        MEWCP_clone_list_nodes_solution(list_nodes_rounded,open_node->list_nodes_solution,num_partitions);
    }


    MEWCP_free_list_nodes_solution(list_nodes_rounded);
    DSDPDestroy ( dsdp);


}

bool MEWCP_branch( open_node_t * open_node,
                   const unsigned int dim_matrix,
                   const unsigned int num_nodes,
                   const unsigned int num_partitions,
                   const unsigned int num_constraints,
                   const unsigned int father_depth_level,
                   unsigned int * serial_number_node,
                   open_node_t ** out_left_son,
                   open_node_t ** out_right_son)

{
#if defined MEWCP_DSDP_DEBUG
    printf("* MEWCP_branching *\n");
#endif

    unsigned cardinality_partition;


    /* varibles for branching */
    int out_num_part;
    int out_id_node;
    bool possible_branch;
    open_node_t * son_left;  /* In case of branching is possible */
    open_node_t * son_right;


    cardinality_partition = num_nodes/num_partitions;



    /* WELL Now let's generate branching sons */
    possible_branch = MEWCP_generate_perfect_equi_branch(open_node->diagX,num_nodes,num_partitions,&out_num_part,&out_id_node);
    //possible_branch = MEWCP_generate_max_fractional_branch(open_node->diagX,num_nodes,num_partitions,&out_num_part,&out_id_node);


#if defined MEWCP_DSDP_DEBUG

    printf("Possible branch: %d \t partition: %d \t id_node: %d \n",possible_branch, out_num_part, out_id_node);
#endif

    if (possible_branch == true) /* Branching is possible */
    {

        son_left = MEWCP_allocate_open_node();
        son_right = MEWCP_allocate_open_node();

        son_left->id_node = 2*open_node->id_node +1;
        son_right->id_node = 2*open_node->id_node +2;

        *serial_number_node +=1;
        son_left->serial_node = *serial_number_node;
        *serial_number_node +=1;
        son_right->serial_node = *serial_number_node;


        son_left->list_blocked_nodes = MEWCP_allocate_list_blocked_nodes(num_nodes);
        son_right->list_blocked_nodes = MEWCP_allocate_list_blocked_nodes( num_nodes);

        /* I generate blocked nodes */
        MEWCP_generate_list_blocked_nodes_branching_sons(open_node->list_blocked_nodes, son_left->list_blocked_nodes, son_right->list_blocked_nodes,out_num_part, out_id_node,num_nodes,cardinality_partition);

#if defined MEWCP_DSDP_DEBUG

        MEWCP_print_list_blocked_nodes(son_left->list_blocked_nodes);
        MEWCP_print_list_blocked_nodes(son_right->list_blocked_nodes);
#endif

        /* generate branching constraints based on the blocked nodes */
        son_left->vect_mat_branching_contraint = MEWCP_allocate_vect_mat_branching_constraints(dim_matrix);
        son_right->vect_mat_branching_contraint = MEWCP_allocate_vect_mat_branching_constraints( dim_matrix);

        MEWCP_generate_constraints_branch(son_left->list_blocked_nodes, son_left->vect_mat_branching_contraint,dim_matrix,num_nodes,cardinality_partition);
        MEWCP_generate_constraints_branch(son_right->list_blocked_nodes, son_right->vect_mat_branching_contraint, dim_matrix,num_nodes,cardinality_partition);


        /* I clone the vect_y */
        son_left->vect_y = MEWCP_allocate_vect_y(num_constraints);
        son_right->vect_y = MEWCP_allocate_vect_y(num_constraints);

        MEWCP_clone_vect_y(open_node->vect_y, son_left->vect_y,num_constraints);
        MEWCP_clone_vect_y(open_node->vect_y, son_right->vect_y,num_constraints);

        /* I set the level */
        son_left->depth_level = father_depth_level +1;
        son_right->depth_level = father_depth_level +1;

        *out_left_son = son_left;
        *out_right_son = son_right;

        return true;
    }
    else
    {
        return false;
    }

}

branching_open_node_t * MEWCP_find_worst_bound_element(list_branching_t * list_branching )
{
#if defined MEWCP_DSDP_DEBUG
    printf("* MEWCP_find_worst_bound_element *\n");
#endif

    branching_open_node_t * worst_element = NULL;
    branching_open_node_t * element;
    double worst_bound;

    worst_bound = MEWCP_MIN_DOUBLE;

#if defined ASSERT

    assert(MEWCP_is_list_branching_empty(list_branching) == false);
#endif

    /* I search worse element from the tail */
    for(element = list_branching->tail; element != NULL; element = element->prev_branching_node)
    {
        if( (element->open_node->DB -  worst_bound) >= MEWCP_EPSILON )
        {
            worst_bound = element->open_node->DB;
            worst_element = element;
        }
    }
    return worst_element;
}

bool MEWCP_is_list_branching_empty(list_branching_t * list_branching)
{
    bool is_empty;

    /* I check that pointer to head and tail are NULL */

    if ( (list_branching->head == NULL) &&  ( list_branching->tail == NULL  ) )
    {
        is_empty = true;
    }
    else
    {
        is_empty = false;
    }
    return is_empty;
}

void MEWCP_push_open_node(open_node_t * open_node, list_branching_t * list_branching)
{
    /* I insert the element to the head */
    branching_open_node_t * element;
    branching_open_node_t * element_tmp;

    /* I create a list element */
    element = MEWCP_allocate_branching_open_node();

    element->open_node = open_node;

    /* I have to append the element to the head */
    /* The list is empty */
    if( MEWCP_is_list_branching_empty(list_branching) == true)
    {
        list_branching->head = element;
        list_branching->tail = element;

        element->prev_branching_node = NULL;
        element->next_branching_node = NULL;
    }
    else /* the list isn't empty */
    {

        element_tmp = list_branching->head;

        list_branching->head = element;

        element->prev_branching_node = NULL;
        element->next_branching_node = element_tmp;

        element_tmp->prev_branching_node = element;
    }

    list_branching->number_open_nodes += 1;


}

open_node_t * MEWCP_pop_open_node( list_branching_t * list_branching)
{
    open_node_t * open_node;
    branching_open_node_t * element;

#if defined ASSERT

    assert( MEWCP_is_list_branching_empty(list_branching) == false);
#endif

    /* I extract the element from the tail */

    /* If the element is the first */
    if ( list_branching->tail->prev_branching_node == NULL)
    {
#if defined ASSERT
        /* The head and the tail point to the same first element */
        assert(  list_branching->head == list_branching->tail );
#endif

        element = list_branching->tail;
        list_branching->head = NULL;
        list_branching->tail = NULL;
    }
    else /* element isn't the first */
    {
#if defined ASSERT
        /* I have at least one previous element  */
        assert (  list_branching->tail->prev_branching_node != NULL );
#endif

        element = list_branching->tail;

        list_branching->tail = element->prev_branching_node;
        list_branching->tail->next_branching_node = NULL;

    }

    open_node = element->open_node;
    list_branching->number_open_nodes -= 1;
    free(element); /* I only free the element with the pointers but not the open_node structure */
    return open_node;
}

open_node_t * MEWCP_pop_specific_open_node(branching_open_node_t * branching_open_node , list_branching_t * list_branching)
{
#if defined MEWCP_DSDP_DEBUG
    printf("* MEWCP_pop_specific_open_node: %d *\n",branching_open_node->open_node->id_node);
#endif

    open_node_t * open_node;

#if defined ASSERT

    assert( MEWCP_is_list_branching_empty(list_branching) == false);
#endif

    /* I extract the given  branching_open_node element from the list */

    /* If the element is the only element in the list */

    if ( (branching_open_node->prev_branching_node == NULL) && ( branching_open_node->next_branching_node == NULL) )
    {
#if defined ASSERT
        /* The head and the tail point to the same first element */
        assert(  list_branching->head == list_branching->tail );
        assert(  list_branching->head ==  branching_open_node);

#endif

        list_branching->head = NULL;
        list_branching->tail = NULL;
    }
    else if (branching_open_node->prev_branching_node == NULL) /* element is the first but not the only*/
    {
#if defined ASSERT
        /* I have at least one follow element  */
        assert (  list_branching->head == branching_open_node);
        assert ( branching_open_node->next_branching_node != NULL );
#endif

        list_branching->head = branching_open_node->next_branching_node;
        list_branching->head->prev_branching_node = NULL;


    }
    else if (branching_open_node->next_branching_node == NULL) /* element is the last but not the only*/
    {
#if defined ASSERT
        /* I have at least one follow element  */
        assert (  list_branching->tail == branching_open_node);
        assert ( branching_open_node->prev_branching_node != NULL );
#endif

        list_branching->tail = branching_open_node->prev_branching_node;
        list_branching->tail->next_branching_node = NULL;
    }
    else /* Element is not the tail nor the head of the list */
    {
#if defined ASSERT

        assert ( branching_open_node->prev_branching_node != NULL );
        assert ( branching_open_node->next_branching_node != NULL );
#endif

        branching_open_node->prev_branching_node->next_branching_node = branching_open_node->next_branching_node;
        branching_open_node->next_branching_node->prev_branching_node = branching_open_node->prev_branching_node;
    }

    open_node = branching_open_node->open_node;
    list_branching->number_open_nodes -= 1;
    free(branching_open_node); /* I only free the element with the pointers but not the open_node structure */
    return open_node;
}


void MEWCP_generate_list_blocked_nodes_branching_sons(list_blocked_nodes_t * father_blocked_nodes,
        list_blocked_nodes_t * left_son_blocked_nodes,
        list_blocked_nodes_t * right_son_blocked_nodes,
        const int branch_num_part,
        const int branch_id_node,
        const unsigned int num_nodes,
        const unsigned int cardinality_partition )
{
    int  boundaries[4];
    int i;


    /* As first thing I clone the father blocked list */
    MEWCP_clone_list_blocked_modes(father_blocked_nodes, left_son_blocked_nodes,num_nodes);
    MEWCP_clone_list_blocked_modes(father_blocked_nodes, right_son_blocked_nodes,num_nodes);

    /* I consider left son the one with nodes blocked on the right of branch_id_node in the same partition */
    /* I consider right son the one with nodes blocked on the left of branch_id_node (also blocked) in the same partition */

    /* [ left-son +| branch_id_node || 0 | 0 | 0 ]
     * [ 0 | 0 | 0 | branch_id_node || right-son ]
     */

    /* I have in boundaries[0]=a and boundaries[1]=b   [a  X  X  X  b] */
    trova_boundaries_diagonale(cardinality_partition,branch_id_node,boundaries);

    /* I fill in the right son first */
    for(i=boundaries[0]; i<= branch_id_node; ++i )
    {
        MEWCP_add_blocked_node(i,right_son_blocked_nodes);
    }
    /* Now, I fill in the left son */
    for(i=branch_id_node+1; i<= boundaries[1]; ++i )
    {
        MEWCP_add_blocked_node(i,left_son_blocked_nodes);
    }

}

bool MEWCP_generate_equi_branch_node(double * diag_X, const unsigned int n, const unsigned int m,  int * out_num_part,  int * out_id_node)
{
#if defined MEWCP_CONVERTER_DSDP_VERBOSE1
    printf("* MEWCP_generate_equi_branch_node *\n");
#endif

    unsigned int i,j,c;
    double sum_partition[m];
    int  id_node_part[m];
    int num_fract_var[m];
    double min_diff;
    double cur_min_diff;

    double  sum_cur;
    int max_num_fract_var =0;

    c = n/m;

    for(i=0;i<m;++i)
    {

        cur_min_diff = 1;
        sum_cur = 0;

        num_fract_var[i] = 0;


        for(j=i*c;j<(i*c +c); ++j )
        {
            /* I find the number of fractional variables */
            sum_cur += diag_X[j];


            if ( (fabs(sum_cur - 0.5) < cur_min_diff) && (diag_X[j] > MEWCP_EPSILON )  && ( (1.0 - diag_X[j])  > MEWCP_EPSILON)  )
            {
                cur_min_diff = fabs(sum_cur - 0.5);

                id_node_part[i] = j;
                sum_partition[i] = sum_cur;

                num_fract_var[i] += 1;

                if (num_fract_var[i] > max_num_fract_var)
                {
                    /* I set the maximum fractional number of variables */
                    max_num_fract_var = num_fract_var[i];
                }

            }
#if defined MEWCP_BRANCHING_DEBUG

            printf("Part: %d \t diag[%d]: %.5lf \t Sum: %.5lf num_fract: %d\n",i,j,diag_X[j],sum_cur,num_fract_var[i]);
#endif

        }
#if defined MEWCP_BRANCHING_DEBUG
        printf("Ho preso: %d \t somma: %.2lf\n",id_node_part[i],sum_partition[i]);
#endif

    }
    min_diff = 1;
    for (i=0;i<m;++i)
    {
        if ( num_fract_var[i] == max_num_fract_var )
        {



            /* I take the partition node nearest to sum 0.5 */

            if (  ( fabs(sum_partition[i] - 0.5) < min_diff)   )
            {


                min_diff = fabs(sum_partition[i] - 0.5);
#if defined MEWCP_BRANCHING_DEBUG

                printf("sum_part[%d]: %.2lf \t mindiff: %.2lf\n",i,sum_partition[i], min_diff);
#endif

                *out_id_node = id_node_part[i];
                *out_num_part = i;
            }

        }
    }



    if ( max_num_fract_var < 1 )
    {

        return false;
    }
    else
    {
#if defined MEWCP_BRANCHING_DEBUG
        printf("taken\t part: %d \t node: %d\n",*out_num_part, *out_id_node);
#endif

        return true;
    }

}

bool MEWCP_generate_perfect_equi_branch(double * diag_X, const unsigned int n, const unsigned int m,  int * out_num_part,  int * out_id_node)
{
#if defined MEWCP_CONVERTER_DSDP_VERBOSE1
    printf("* MEWCP_generate_perfect_equi_branch_node *\n");
#endif


    unsigned int i,j,c;
    double sum_partition[m];
    int  id_node_part[m];
    int num_fract_var[m];
    double min_diff;
    double cur_min_diff;

    double  sum_cur;
    int max_num_fract_var =0;

    c = n/m;

    for(i=0;i<m;++i)
    {

        cur_min_diff = 1;
        sum_cur = 0;

        num_fract_var[i] = 0;


        for(j=i*c;j<(i*c +c); ++j )
        {
            /* I find the number of fractional variables */
            sum_cur += diag_X[j];


            if( (diag_X[j] - 0 ) >= MEWCP_EPSILON )
            {
                num_fract_var[i] += 1;

                if (num_fract_var[i] > max_num_fract_var)
                {
                    /* I set the maximum fractional number of variables */
                    max_num_fract_var = num_fract_var[i];
                }
            }





            if ( fabs(sum_cur - 0.5) < cur_min_diff )
            {
                cur_min_diff = fabs(sum_cur - 0.5);

                id_node_part[i] = j;
                sum_partition[i] = sum_cur;

            }
#if defined MEWCP_BRANCHING_DEBUG

            printf("Part: %d \t diag[%d]: %.4lf \t Sum: %.2lf num_fract: %d\n",i,j,diag_X[j],sum_cur,num_fract_var[i]);
#endif

        }
#if defined MEWCP_BRANCHING_DEBUG
        printf("Ho preso: %d \t somma: %.2lf\n",id_node_part[i],sum_partition[i]);
#endif

    }
    min_diff = 1;
    for (i=0;i<m;++i)
    {
        if ( num_fract_var[i] == max_num_fract_var )
        {



            /* I take the partition node nearest to sum 0.5 */

            if (  ( fabs(sum_partition[i] - 0.5) < min_diff)   )
            {


                min_diff = fabs(sum_partition[i] - 0.5);
#if defined MEWCP_BRANCHING_DEBUG

                printf("sum_part[%d]: %.2lf \t mindiff: %.2lf\n",i,sum_partition[i], min_diff);
#endif

                *out_id_node = id_node_part[i];
                *out_num_part = i;
            }

        }
    }



    if ( max_num_fract_var < 2 )
    {

        return false;
    }
    else
    {
#if defined MEWCP_BRANCHING_DEBUG
        printf("taken\t part: %d \t node: %d\n",*out_num_part, *out_id_node);
#endif

        return true;
    }
}

/*
 * I select the partition with the most number of fractional variables
 * I take the element that divides the partition in 2 equi number of fractional
 * If I have more than one partition with the same max number of fract variables I take the one with pivot element having nearest sum to 1/2
 */
bool MEWCP_generate_max_fractional_branch(double * diag_X, const unsigned int n, const unsigned int m,  int * out_num_part,  int * out_id_node)
{
#if defined MEWCP_CONVERTER_DSDP_VERBOSE1
    printf("* MEWCP_generate_max_fractional_branch *\n");
#endif


    unsigned int i,j,c;


    int num_fract_var[m];
    double min_diff;


    double  sum_cur;
    int max_num_fract_var =0;

    c = n/m;
    int fract[m][c];  // number of fractional varibles until that node
    double sum[m][c];  //density until that element

    /* I have to clean the vectors */
    for(i=0;i<m;++i)
    {
        for (j = 0; j < c; ++j)
        {
            fract[i][j] = 0;
            sum[i][j] = 0.0;
        }
    }


    for(i=0;i<m;++i)
    {


        sum_cur = 0;

        num_fract_var[i] = 0;


        for(j=i*c;j<(i*c +c); ++j )
        {
            /* I find the number of fractional variables */
            sum_cur += diag_X[j];
            sum[i][j%c] = sum_cur;


            if( (diag_X[j] - 0 ) >= MEWCP_EPSILON )
            {
                num_fract_var[i] += 1;
                fract[i][j%c] = num_fract_var[i];

                if (num_fract_var[i] > max_num_fract_var)
                {
                    /* I set the maximum fractional number of variables */
                    max_num_fract_var = num_fract_var[i];
                }
            }
#if defined MEWCP_BRANCHING_DEBUG
            printf("part: %d id: %d\t diag: %.2lf  sum: %.2lf  fract: %d\n",i,j,diag_X[j],sum[i][j%c], fract[i][j%c]);
#endif

        }
#if defined MEWCP_BRANCHING_DEBUG
        printf("\n");
#endif

    }
    int half_pivot;
    int candidate_node = -1;
    int candidate_partition = -1;
    min_diff = 1;

    half_pivot = (int) (floor ((float) max_num_fract_var / 2));
    for (i=0;i<m;++i)
    {
        if ( num_fract_var[i] == max_num_fract_var )
        {

            for(j=0;j<c;++j)
            {
                if(fract[i][j] == half_pivot)
                {

                    if( (fabs(sum[i][j] - 0.5) < min_diff) )
                    {
                        min_diff = fabs((sum[i][j] - 0.5));
                        candidate_node = j;
                        candidate_partition = i;

                    }
                    break;
                }
            }


            *out_id_node = (candidate_partition*c) + candidate_node;
            *out_num_part = candidate_partition;
        }

    }




    if ( max_num_fract_var < 2 )
    {

        return false;
    }
    else
    {
#if defined MEWCP_BRANCHING_DEBUG
        printf("taken\t part: %d \t node: %d\n",*out_num_part, *out_id_node);
#endif

        return true;
    }

}




void MEWCP_add_blocked_node(const unsigned int id_node, list_blocked_nodes_t * list_blocked_nodes)
{
    unsigned int pos_insertion;

    pos_insertion = list_blocked_nodes->num_blocked_nodes;

    if (list_blocked_nodes->bool_list[id_node] == false)
    {
        list_blocked_nodes->bool_list[id_node] = true;
        list_blocked_nodes->blocked_node[pos_insertion] = id_node;
        list_blocked_nodes->num_blocked_nodes += 1;
    }
#if defined MEWCP_BRANCHING_DEBUG
    else
    {
        printf("Node: %d is already blocked!\n",id_node);


    }
#endif

    /* All structures are updated */
}

double MEWCP_evaluate_list_nodes_solution(int * list_node_solution, matrix_weights_t * matrix_weights ,const unsigned int m)
{
#if defined MEWCP_CONVERTER_DSDP_VERBOSE1
    printf("* MEWCP_evaluate_list_nodes_solution *\n");
#endif

    unsigned int i,j;

    double z_tmp;
    int pos_i;
    int pos_j;

    z_tmp = 0;



    for (i=0;i<m;++i)
    {
        for(j=i; j<m;++j)
        {
            pos_i = list_node_solution[i];
            pos_j = list_node_solution[j];
            z_tmp += (double) matrix_weights->weight[pos_i][pos_j];
        }
    }
#if defined MEWCP_CONVERTER_DSDP_DEBUG
    printf("Objective value of node_list: %.2lf\n",z_tmp);
#endif

    return z_tmp;
}

int sort_compare (const void * a, const void * b)
{
    return ( *(int*)a - *(int*)b );

}

void Take_Time(double * user_time, double * system_time)
{
    struct tms buff;
    times(&buff);
    *system_time = buff.tms_stime;
    *user_time = buff.tms_utime;
    *system_time /= 100;
    *user_time /= 100;
    
    
}

double get_cpu_time(void)
{
	double user_time;
	double  system_time;
	struct tms buff;
    
    times(&buff);
    system_time = buff.tms_stime;
    user_time = buff.tms_utime;
    system_time /= 100;
    user_time /= 100;
    
    return user_time + system_time; 
	
}

open_node_t * MEWCP_allocate_open_node(void)
{
    open_node_t * open_node;

    open_node = (open_node_t *) calloc(1,sizeof (open_node_t));
    if ( open_node == NULL)
    {
        printf("!!! ERROR allocation open_node! \n");
        exit(EXIT_FAILURE);
    }

    open_node->id_node = -1;
    open_node->DB = (double) MEWCP_MAX_DOUBLE;
    open_node->PB = (double) MEWCP_MIN_DOUBLE;



    open_node->list_nodes_solution = NULL;
    open_node->diagX = NULL;
    open_node->list_blocked_nodes = NULL;
    open_node->vect_mat_branching_contraint = NULL;

    return open_node;

}

solution_bb_t * MEWCP_allocate_solution_bb(unsigned int num_partitions)
{
#if defined MEWCP_DSDP_DEBUG
    printf("* MEWCP_allocate_solution_bb *\n");
#endif

    solution_bb_t * solution_bb;

    solution_bb = ( solution_bb_t *) calloc(1,sizeof(solution_bb_t) );
    if ( solution_bb == NULL)
    {
        printf("!!! ERROR allocation MEWCP_allocate_solution_bb! \n");
        exit(EXIT_FAILURE);
    }

    solution_bb->list_nodes_best_solution= MEWCP_allocate_list_nodes_solution(num_partitions);

    return solution_bb;

}


list_blocked_nodes_t *  MEWCP_allocate_list_blocked_nodes(const unsigned int num_nodes )
{
#if defined MEWCP_CONVERTER_DSDP_VERBOSE1
    printf("* MEWCP_allocate_list_blocked_nodes *\n");
#endif

    list_blocked_nodes_t * list_blocked_nodes;

    list_blocked_nodes = (list_blocked_nodes_t *) calloc(1, sizeof( list_blocked_nodes_t));
    if (list_blocked_nodes == NULL)
    {
        printf("!!! ERROR allocation list_blocked_nodes! \n");
        exit(EXIT_FAILURE);
    }

    list_blocked_nodes -> blocked_node = (unsigned int *) calloc(num_nodes, sizeof(unsigned int));
    if (list_blocked_nodes -> blocked_node  == NULL)
    {
        printf("!!! ERROR allocation list_blocked_nodes -> blocked_node! \n");
        exit(EXIT_FAILURE);
    }
    list_blocked_nodes->bool_list = (bool *) calloc(num_nodes, sizeof(bool));

    if ( list_blocked_nodes->bool_list == NULL)
    {
        printf("!!! ERROR allocation list_blocked_nodes->bool_list! \n");
        exit(EXIT_FAILURE);
    }

    list_blocked_nodes->num_blocked_nodes = 0;

    return list_blocked_nodes;

}

branching_open_node_t * MEWCP_allocate_branching_open_node(void)
{
#if defined MEWCP_DSDP_DEBUG
    printf("* MEWCP_allocate_branching_open_node *\n");
#endif

    branching_open_node_t * branching_open_node;

    branching_open_node = (branching_open_node_t *) calloc(1,sizeof(branching_open_node_t) );

    branching_open_node->next_branching_node = NULL;
    branching_open_node->prev_branching_node = NULL;
    branching_open_node->open_node = NULL;

    return branching_open_node;
}


list_branching_t * MEWCP_allocate_list_branching(void)
{

#if defined MEWCP_DSDP_DEBUG
    printf("* MEWCP_allocate_list_branching *\n");
#endif

    list_branching_t * list_branching;

    list_branching = (list_branching_t *) calloc(1,sizeof(list_branching_t));
    list_branching->list_nodes_best_solution = NULL;
    list_branching->branching_open_node_t = NULL;
    list_branching->head = NULL;
    list_branching->tail = NULL;

    list_branching->number_open_nodes = 0;

    return list_branching;


}


constraint_t * MEWCP_allocate_sdp_constraints_matrix(const unsigned int n, const unsigned int num_constraints )
{
#if defined MEWCP_CONVERTER_DSDP_VERBOSE1
    printf("* MEWCP_allocate_sdp_constraints_matrix *\n");
#endif

    unsigned int i;
    unsigned int dim_vect_matrix;

    dim_vect_matrix = n*(n+1)/2;

    constraint_t * sdp_constraints_matrix;

    /* I allocate the constraints matrix that is card(b) -1
    * I don't consider the last matrix of branching constraints */
    sdp_constraints_matrix = (constraint_t *) calloc(num_constraints , sizeof(constraint_t));
    if (sdp_constraints_matrix == NULL)
    {
        printf("Allocation ERROR: sdp_constraints_matrix\n");
        exit (EXIT_FAILURE);
    }

    for (i=0; i < num_constraints ; ++i )
    {
        sdp_constraints_matrix[i].index = (int *) calloc(dim_vect_matrix, sizeof(int));
        sdp_constraints_matrix[i].weight = (double *) calloc(dim_vect_matrix, sizeof(double));

        if ( (sdp_constraints_matrix[i].index == NULL) || (sdp_constraints_matrix[i].weight == NULL))
        {
            printf("Allocation ERROR: sdp_constraints_matrix\n");
            exit (EXIT_FAILURE);
        }
    }

    return sdp_constraints_matrix;
}



double * MEWCP_allocate_bi(const unsigned int num_contraints )
{
#if defined MEWCP_CONVERTER_DSDP_VERBOSE1
    printf("* MEWCP_allocate_bi *\n");
#endif

    double * bi;

    bi = (double *) calloc(num_contraints, sizeof(double));
    return bi;
}

unsigned int MEWCP_convert_coords_ij_to_vector_matrix(const unsigned int i, const unsigned j)
{
    /* i >= j
     * i,j = 1,...,n
     * */
#if defined ASSERT
    assert(i>=j);
#endif

    return i*(i-1)/2 +j -1;
}

constraint_t * MEWCP_allocate_vect_mat_branching_constraints(const unsigned int dim_matrix)
{
    constraint_t * vect_mat_branching_contraints;

    vect_mat_branching_contraints = (constraint_t *) calloc(1, sizeof(constraint_t));
    if (vect_mat_branching_contraints  == NULL)
    {
        printf("*******  Allocation MEWCP_allocate_vect_mat_branching_constraints FAILED!\n");
        exit(EXIT_FAILURE);
    }
    vect_mat_branching_contraints->index = (int *) calloc(dim_matrix,sizeof(int));
    vect_mat_branching_contraints->weight = (double *) calloc(dim_matrix,sizeof(double));
    if ( vect_mat_branching_contraints->index  == NULL)
    {
        printf("*******  Allocation MEWCP_allocate_vect_mat_branching_constraints FAILED!\n");
        exit(EXIT_FAILURE);
    }

    if ( vect_mat_branching_contraints->weight  == NULL)
    {
        printf("*******  Allocation MEWCP_allocate_vect_mat_branching_constraints FAILED!\n");
        exit(EXIT_FAILURE);
    }



    return vect_mat_branching_contraints;
}

int * MEWCP_allocate_list_nodes_solution( const unsigned int m)
{
    int * list_nodes_solution;
    list_nodes_solution = ( int *) calloc (m, sizeof(int));
    if (list_nodes_solution  == NULL)
    {
        printf("*******  Allocation MEWCP_allocate_list_nodes_solution FAILED!\n");
        exit(EXIT_FAILURE);
    }
    return list_nodes_solution;
}

double * MEWCP_allocate_diag_X(const unsigned int length)
{
    double * diagX;

    diagX = (double *)  calloc( length, sizeof(double));
    if (diagX  == NULL)
    {
        printf("*******  Allocation MEWCP_allocate_diag_X FAILED!\n");
        exit(EXIT_FAILURE);
    }

    return diagX;
}

double * MEWCP_allocate_vect_y(const unsigned int num_constraints)
{
    double * vect_y;

    vect_y = (double *) calloc(num_constraints, sizeof(double));
    if (vect_y  == NULL)
    {
        printf("*******  Allocation MEWCP_allocate_vect_y FAILED!\n");
        exit(EXIT_FAILURE);
    }


    return vect_y;
}

void trova_boundaries_diagonale(const unsigned int c, const unsigned int i, int * boundaries)
{
    // Devo stabilire in che partizione Sk si trova il punto i,j
    // mi aspetto i,j i=j
    // trovo le coordinate interne che mi delimitano la partizione
    unsigned int k; // mi dice in quale partizione sono (0,...,m-1)

    k=i/c;  /* mi aspetto k=0,...,(m-1) */

#if defined LOG2

    printf("Sono nella aprtizione: %d\n",k);
#endif

    // mi trovo xa,xb
    boundaries[0] = k *c;  //xa
    boundaries[1] = c*(k+1) -1;  //xb

    boundaries[2] = boundaries[0]; //ya
    boundaries[3] = boundaries[1]; //yb

#if defined LOG2

    printf("I'm in the partition: %d\n",k);
    printf("Boundaries: xa:%d\txb:%d\tya:%d\tyb:%d\n",boundaries[0],boundaries[1],boundaries[2],boundaries[3]);
#endif

}

/*
 * Cloned the blocked nodes list.
 */
void MEWCP_clone_list_blocked_modes(list_blocked_nodes_t * list_to_be_cloned, list_blocked_nodes_t * list_cloned, const unsigned int num_nodes)
{

    /* I copy the list */
    memcpy(list_cloned->blocked_node, list_to_be_cloned->blocked_node, num_nodes * sizeof(unsigned int) );
    memcpy(list_cloned->bool_list, list_to_be_cloned->bool_list, num_nodes * sizeof(bool));

    list_cloned->num_blocked_nodes = list_to_be_cloned->num_blocked_nodes;
}


void MEWCP_clone_list_nodes_solution( int * list_to_be_cloned,  int * list_cloned, const unsigned int num_nodes)
{
    /* I copy the list */
    memcpy(list_cloned,list_to_be_cloned,num_nodes*sizeof( int));
}

void MEWCP_clone_vect_y(double * vect_y_to_be_cloned, double * vect_y_cloned, const unsigned int num_constraints)
{
    memcpy(vect_y_cloned, vect_y_to_be_cloned, num_constraints*sizeof(double));

}

void MEWCP_dump_diag_X(SDPCone * sdpcone, double * dst_diag_X, const unsigned int num_nodes)
{
    unsigned int i;
    int dim_matrix;
    double * sol_X;
    unsigned int pos;

    (void) SDPConeGetXArray(*sdpcone, 0, &sol_X, &dim_matrix);



    for(i=0; i< num_nodes; ++i)
    {
        pos = MEWCP_convert_coords_ij_to_vector_matrix(i+1,i+1);
        dst_diag_X[i] = sol_X[pos];
    }

}


void MEWCP_dump_vect_y(DSDP  * dsdp, double * dst_vect_y, const unsigned int num_constraints)
{

    DSDPGetY (*dsdp, dst_vect_y , num_constraints);
}


bool MEWCP_is_new_best_PB_and_update(open_node_t * open_node, list_branching_t * list_branching, const unsigned int num_partitions)
{
    bool has_been_updated = false;

    /* if PB improved */
    if( (open_node->PB - list_branching->best_primal) >= MEWCP_EPSILON)
    {
        list_branching->best_primal = open_node->PB;
        list_branching->id_node_best_primal = open_node->id_node;
        list_branching->serial_node_best_primal = open_node->serial_node;
        list_branching->depth_node_best_primal = open_node->depth_level;
        MEWCP_clone_list_nodes_solution(open_node->list_nodes_solution,list_branching->list_nodes_best_solution, num_partitions);

        /* Ok update completed */
        has_been_updated = true;
    }
    return has_been_updated;
}


void MEWCP_compute_sdp_rounding(double * diag_X,  int * list_nodes_rounded_solution, const unsigned int num_nodes, const unsigned int c)
{
#if defined MEWCP_CONVERTER_DSDP_VERBOSE1
    printf("* MEWCP_compute_sdp_rounding *\n");
#endif

    unsigned int i;
    unsigned int j;
    unsigned int m;
    double val_max;
    unsigned pos_max;

    m = num_nodes / c;

    for (i=0;i<m;++i)
    {
        val_max = 0;
        pos_max = i*c;

        for (j=i*c;j<(i*c+c); ++j)
        {
            if( (diag_X[j] - val_max) > MEWCP_EPSILON)
            {
                val_max = diag_X[j] ;
                pos_max = j;
            }
        }

        list_nodes_rounded_solution[i] = pos_max;
    }

}

/* ******************************************
 * PRINTING FUNCTIONS 
 * ******************************************/


void MEWCP_print_contraints_matrix(double ** matrix, const unsigned int length_i, const unsigned int length_j)
{
    printf("[Printing matrix of contraints ]\n");
    unsigned int i,j;

    for (i=0;i<length_i;++i)
    {
        for(j=0;j<length_j;++j)
        {
            printf("[%d][%d] %.2lf\t",i,j,matrix[i][j]);

        }
        printf("\n");
    }
}

void MEWCP_print_vectorY(double * vector, const unsigned int dim)
{
    unsigned int i;
    printf("Varibili Y: ");

    for (i=0;i<dim;++i)
    {
        printf("[%d] %.2lf\t",i,vector[i]);
    }
    printf("\n");

}

void MEWCP_print_diag_X(double * diag_X, const unsigned int num_nodes)
{
    unsigned int i;

    printf("Diag X: ");
    for (i=0;i<num_nodes; ++i)
    {
        printf("%.2lf ",diag_X[i]);
    }
    printf("\n");
}


void MEWCP_print_list_blocked_nodes(list_blocked_nodes_t * list_blocked_modes)
{
    unsigned int i;
    unsigned int n;

    n=list_blocked_modes->num_blocked_nodes;
    printf("Blocked nodes: ");

    for (i=0;i<n;++i)
    {
        printf("%u ",list_blocked_modes->blocked_node[i]);
    }
    printf("\n");
}

void MEWCP_print_list_nodes_solution( int * list_nodes_solution, const unsigned int m)
{
    unsigned int i;
    printf("Nodes list solution: : ");

    for (i=0;i<m;++i)
    {
        printf("%d ",list_nodes_solution[i]);
    }
    printf("\n");
}


void MEWCP_print_list_nodes_solution_cplex( int * list_nodes_solution, const unsigned int m)
{
    unsigned int i;
    printf("Nodes list solution (cplex): ");

    qsort (list_nodes_solution, m, sizeof(int), sort_compare);


    for (i=0;i<m;++i)
    {
        printf("%d ",list_nodes_solution[i]+1);
    }
    printf("\n");
}

/* ******************************************
 * FREE FUNCTIONS 
 * ******************************************/
void MEWCP_free_list_blocked_nodes(list_blocked_nodes_t * list_blocked_nodes)
{
#if defined MEWCP_DSDP_DEBUG
    printf("* MEWCP_free_list_blocked_nodes *\n");
#endif

    free( list_blocked_nodes->blocked_node);
    free (list_blocked_nodes->bool_list);
}

void MEWCP_free_sdp_constraints_matrix(constraint_t * sdp_constraints_matrix, const unsigned int num_constraints )
{
#if defined MEWCP_DSDP_DEBUG
    printf("* MEWCP_free_sdp_constraints_matrix *\n");
#endif

    unsigned int i;
    for(i=0;i<num_constraints;++i) /* matrix are num_contraints +1 but we do not consider the branching constraint matrix */
    {
        free(sdp_constraints_matrix[i].index);
        free(sdp_constraints_matrix[i].weight);
    }
}
void MEWCP_free_bi(double * bi)
{
#if defined MEWCP_DSDP_DEBUG
    printf("* MEWCP_free_bi *\n");
#endif

    free(bi);
}

void MEWCP_free_diag_X(double * diag_X)
{
#if defined MEWCP_DSDP_DEBUG
    printf("* MEWCP_free_diag_X *\n");
#endif

    free(diag_X);
}

void MEWCP_free_vect_y(double * vect_y)
{
#if defined MEWCP_DSDP_DEBUG
    printf("* MEWCP_free_vect_y *\n");
#endif

    free(vect_y);
}

void MEWCP_free_list_nodes_solution( int * list_nodes_solution)
{
#if defined MEWCP_DSDP_DEBUG
    printf("* MEWCP_free_list_nodes_solution *\n");
#endif

    free(list_nodes_solution);
}

void MEWCP_free_vect_mat_branching_constraints(constraint_t * vect_mat_branching_contraints)
{
#if defined MEWCP_DSDP_DEBUG
    printf("* MEWCP_free_vect_mat_branching_constraints *\n");
#endif

    free(vect_mat_branching_contraints->index);
    free(vect_mat_branching_contraints->weight);
    free(vect_mat_branching_contraints);
}

void MEWCP_free_open_node(open_node_t * open_node)
{
#if defined MEWCP_DSDP_DEBUG
    printf("* MEWCP_free_open_node *\n");
#endif


    MEWCP_free_list_blocked_nodes (open_node->list_blocked_nodes);
    MEWCP_free_vect_mat_branching_constraints(open_node->vect_mat_branching_contraint);
    MEWCP_free_diag_X(open_node->diagX);
    MEWCP_free_list_nodes_solution(open_node->list_nodes_solution);
    MEWCP_free_vect_y(open_node->vect_y);

    free(open_node);
}


void MEWCP_free_branching_open_node(branching_open_node_t * branching_open_node)
{
#if defined MEWCP_DSDP_DEBUG
    printf("* MEWCP_free_branching_open_node *\n");
#endif

    MEWCP_free_open_node(branching_open_node->open_node);
    free(branching_open_node);

}
void MEWCP_free_list_branching(list_branching_t * list_branching)
{
#if defined MEWCP_DSDP_DEBUG
    printf("* MEWCP_free_branching *\n");
#endif

    if (MEWCP_is_list_branching_empty(list_branching) == false)
    {
        printf("**\t ERROR! Cannot free list_branching because isn't empty! Check it out!\n");
        //exit(EXIT_FAILURE);
    }

    MEWCP_free_list_nodes_solution(list_branching->list_nodes_best_solution);
    /* As the list is empty I simply free the list structure */
    free(list_branching);

}

void MEWCP_free_solution_bb(solution_bb_t * solution_bb)
{
#if defined MEWCP_DSDP_DEBUG
    printf("* MEWCP_free_solution_bb *\n");
#endif

    free(solution_bb->list_nodes_best_solution);
    free(solution_bb );
}
