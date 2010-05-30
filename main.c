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

#include "dsdp/dsdp5.h"
#include "converter_dsdp.h"
#include "MEWCP_dsdp.h"



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

    /*
    unsigned  int i;
    int info;
    */

    unsigned int dim_matrix;
    /*
        double pobj;
        double sol_traceX;
        double * sol_vect_X;
        int sol_dim_vect_X;
        int tmp_num_variables;
        */


    unsigned int num_constraints;
    unsigned int num_blocks;
    unsigned int num_nodes;
    unsigned int num_partitions, cardinality_partition;


    double  * bi;
    double * vect_mat_branching_constraints;
    double ** constraints_matrix;

    /*
    SDPCone sdpcone;
    DSDP dsdp;
    DSDPTerminationReason reason;
    */

    open_node_t * open_node;

    filename_in = argv[1];


    num_blocks = NUM_BLOCKS;

    MEWCP_generate_sdp_constraints(filename_in,&constraints_matrix,&bi, &vect_mat_branching_constraints,MEWCP_ALPHA,
                                   c_cardinality,
                                   c_simple_MC,c_improved_MC_A,
                                   c_improved_MC_B,c_improved_MC_C,
                                   c_4C2,
                                   c_4C3_A,
                                   c_4C3_B,
                                   &num_nodes,
                                   &num_partitions,
                                   &cardinality_partition,
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
	

    MEWCP_branch_and_bound(open_node,constraints_matrix,bi,num_constraints,dim_matrix,num_nodes,num_partitions);

	#if defined MEWCP_DSDP_DEBUG
	printf("\tEnd Branch and bound... \n");	
	#endif

    //MEWCP_close_open_node(open_node);
    
    /*
    MEWCP_bound( open_node,  constraints_matrix,  bi,
                                   num_constraints,
                                   dim_matrix,
                                   num_nodes,
                                   num_partitions);

    */

    //MEWCP_print_vectorY(bi,num_constraints);

    /*
        info = DSDPCreate(num_constraints,&dsdp);
        info = DSDPCreateSDPCone(dsdp,num_blocks,&sdpcone);
        info = SDPConeSetBlockSize(sdpcone, 0, num_nodes);  // dimension of block is n 
    */

    /*
     
        for (i=0;i<num_constraints;++i)
        {
            info = DSDPSetDualObjective(dsdp,i+1,bi[i]);
     
        }
     
        for (i=0; i < num_constraints +1; ++i) // Matrix are number of constraints +1 because of matrix W 
        {
            if (i != num_constraints ) // It's not the last constraint matrix 
            {
                SDPConeSetADenseVecMat(sdpcone, 0, i, num_nodes, MEWCP_ALPHA , constraints_matrix[i], dim_matrix);
            }
            else // It's the last matrix containing branching contraints 
            {
                SDPConeSetADenseVecMat(sdpcone, 0, i, num_nodes, 1 , vect_mat_branching_constraints, dim_matrix);
            }
     
    #if defined MEWCP_CONVERTER_DSDP_DEBUG
     
            SDPConeViewDataMatrix(sdpcone, 0, i);
    #endif
     
        }
     
    */

    /* Get read to go */

    /*
        info=DSDPSetGapTolerance(dsdp,MEWCP_GAP_TOLERANCE);
        info=DSDPSetPotentialParameter(dsdp,MEWCP_POTENTIAL_PARAMETER);
        info=DSDPReuseMatrix(dsdp,MEWCP_REUSE_MATRIX);
        info=DSDPSetPNormTolerance(dsdp,MEWCP_SET_PNORM_TOLERANCE);
     
     
    #if defined MEWCP_CONVERTER_DSDP_VERBOSE2
     
        DSDPSetStandardMonitor(dsdp, 1); // verbose each iteration 
        DSDPLogInfoAllow(1,0);
    #endif
     
    #if defined MEWCP_CONVERTER_DSDP_DEBUG
     
        double vect_y[num_constraints ];
        DSDPGetY (dsdp,  vect_y, num_constraints);
        MEWCP_print_vectorY(vect_y,num_constraints);
    #endif
     
     
        DSDPSetup(dsdp);
        DSDPSolve(dsdp);
     
        DSDPComputeX(dsdp);
     
     
     
     
     
     
    #if defined MEWCP_CONVERTER_DSDP_DEBUG
     
        SDPConeGetXArray(sdpcone, 0, &sol_vect_X, &sol_dim_vect_X);
     
        printf("\n");
        SDPConeViewX(sdpcone, 0, num_nodes, sol_vect_X, sol_dim_vect_X);
        printf("\n");
     
        DSDPGetY (dsdp,  vect_y, num_constraints);
        MEWCP_print_vectorY(vect_y,num_constraints);
    #endif
     
     
        double * diagX;
        MEWCP_allocate_diag_X(&diagX, num_nodes);
        MEWCP_dump_diag_X(&sdpcone,diagX,num_nodes);
        //MEWCP_print_diag_X(diagX,num_nodes);
     
        unsigned int * list_nodes_solution;
     
        MEWCP_allocate_list_nodes_solution(&list_nodes_solution,num_partitions);
        MEWCP_compute_sdp_rounding(diagX,list_nodes_solution,num_nodes,cardinality_partition);
        MEWCP_print_list_nodes_solution(list_nodes_solution,num_partitions);
     
        double z_rouded = MEWCP_evaluate_list_nodes_solution(list_nodes_solution,constraints_matrix[0],num_partitions);
    */

    /* generating branching */
    /*
        int out_num_part;
        int out_id_node;
        bool possible_branch;
        possible_branch = MEWCP_generate_equi_branch_node(diagX,num_nodes,num_partitions,&out_num_part,&out_id_node);
     
        printf("Possible branch: %d \t partition: %d \t id_node: %d \n",possible_branch, out_num_part, out_id_node);
     
        list_blocked_nodes_t left_son_blocked_nodes;
        list_blocked_nodes_t right_son_blocked_nodes;
     
        MEWCP_allocate_list_blocked_nodes(&left_son_blocked_nodes,num_nodes);
        MEWCP_allocate_list_blocked_nodes(&right_son_blocked_nodes,num_nodes);
     
        MEWCP_generate_list_blocked_nodes_branching_sons(&list_blocked_nodes,&left_son_blocked_nodes,&right_son_blocked_nodes,out_num_part, out_id_node,num_nodes,cardinality_partition);
     
        MEWCP_print_list_blocked_nodes(&left_son_blocked_nodes);
        MEWCP_print_list_blocked_nodes(&right_son_blocked_nodes);
    */
    /* generate branching constraints based on the blocked nodes */

    /*
        double * left_son_vect_mat_branching_contsraints;
        double * right_son_vect_mat_branching_contsraints;
        MEWCP_allocate_vect_mat_branching_constraints(&left_son_vect_mat_branching_contsraints,dim_matrix);
        MEWCP_allocate_vect_mat_branching_constraints(&right_son_vect_mat_branching_contsraints,dim_matrix);
        MEWCP_generate_constraints_branch(&left_son_blocked_nodes,left_son_vect_mat_branching_contsraints,dim_matrix,num_nodes,cardinality_partition);
        MEWCP_generate_constraints_branch(&right_son_blocked_nodes,right_son_vect_mat_branching_contsraints,dim_matrix,num_nodes,cardinality_partition);
     
        printf("Left: \n");
        SDPConeViewX(sdpcone, 0, num_nodes, left_son_vect_mat_branching_contsraints, dim_matrix);
        printf("Right: \n");
        SDPConeViewX(sdpcone, 0, num_nodes, right_son_vect_mat_branching_contsraints, dim_matrix);
    */

    /* Ok */

    /*
        // (void)	SDPConeRemoveDataMatrix(sdpcone, 0, num_constraints);
     
        //(void) 	SDPConeAddDataMatrix ( sdpcone, 0, num_constraints, num_nodes, char format, struct DSDPDataMat_Ops *dsdpdataops, void *data)
        printf("\n");
     
        SDPConeRemoveDataMatrix (sdpcone, 0, num_constraints);
        SDPConeAddADenseVecMat ( sdpcone, 0, num_constraints, num_nodes,MEWCP_ALPHA , left_son_vect_mat_branching_contsraints, dim_matrix);
     
     
        SDPConeViewDataMatrix(sdpcone, 0, num_constraints);
     
     
     
        DSDPSetup(dsdp);
        DSDPSolve(dsdp);
        DSDPComputeX(dsdp);
     
     
     
     
        SDPConeRemoveDataMatrix (sdpcone, 0, num_constraints);
        SDPConeSetADenseVecMat(sdpcone, 0, num_constraints, num_nodes, 1 , right_son_vect_mat_branching_contsraints, dim_matrix);
        SDPConeViewDataMatrix(sdpcone, 0, num_constraints);
     
        DSDPSetup(dsdp);
        DSDPSolve(dsdp);
     
        DSDPComputeX(dsdp);
     
     
     
        printf("\n");
     
     
     
    */


    /*
    #if defined MEWCP_CONVERTER_DSDP_DEBUG
     
        SDPConeView3(sdpcone);
    #endif
     
        DSDPStopReason(dsdp, &reason);
        DSDPGetPObjective(dsdp,&pobj);
     
     
        DSDPGetTraceX(dsdp, &sol_traceX);
        printf("(P) Obj_value: %.2lf \t Rounded: %.2lf \t Trace(X): %.2lf \t Termination value: %d\n",-pobj,z_rouded, sol_traceX,reason);
        //printf("Trace(X): %.2lf\n",sol_traceX);
        //printf("\n");
     
    #if defined MEWCP_CONVERTER_DSDP_DEBUG
     
        printf("\n\n");
        DSDPView(dsdp);
    #endif
     
        // ************************** END printing solution
     
        DSDPDestroy ( dsdp);
     */

    return EXIT_SUCCESS;
}
