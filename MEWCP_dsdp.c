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
#include <memory.h>

#include "converter_dsdp.h"
#include "MEWCP_dsdp.h"
#include "dsdp/dsdp5.h"


void MEWCP_branch_and_bound(open_node_t * open_root_node, double ** constraints_matrix, double * bi,
                            const unsigned int num_constraints,
                            const unsigned int dim_matrix,
                            const unsigned int num_nodes,
                            const unsigned int num_partitions)
{

#if defined MEWCP_CONVERTER_DSDP_VERBOSE1
    printf("* MEWCP_branch_and_bound\n\tHere we goooooooooooooo!!!! *\n");
#endif

    list_branching_t * list_branching;

    /* Branching variables */
    open_node_t * son_left;
    open_node_t * son_right;
    branching_open_node_t * branching_open_worst_bound;
    open_node_t * open_node_worst_bound;
    bool possible_branch;


    /* First I compute Tabu search in order to get a good Primal Bound */
    // ******************* TABU SEARCH!


    list_branching = MEWCP_allocate_list_branching();

    /* Let's consider root node */
    MEWCP_bound(open_root_node,constraints_matrix,bi,num_constraints,dim_matrix,num_nodes,num_partitions);

    MEWCP_push_open_node(open_root_node,list_branching);
    /* All done, I'm ready to start with branching procedure */


    while(MEWCP_is_list_branching_empty(list_branching) == false)
    {


        branching_open_worst_bound = MEWCP_find_worst_bound_element(list_branching);
        open_node_worst_bound = MEWCP_pop_specific_open_node(branching_open_worst_bound,list_branching);


        if ( (open_node_worst_bound->DB - list_branching->best_primal) > MEWCP_EPSILON )
        {
            possible_branch = MEWCP_branch(open_node_worst_bound,dim_matrix,num_nodes,num_partitions,&son_left,&son_right);

            if (possible_branch == true)
            {
                MEWCP_bound(son_left,constraints_matrix,bi,num_constraints,dim_matrix,num_nodes,num_partitions);
                MEWCP_bound(son_right,constraints_matrix,bi,num_constraints,dim_matrix,num_nodes,num_partitions);

                /* I check if PB is improved */
                if( (son_left->PB - list_branching->best_primal) > MEWCP_EPSILON)
                {
                    list_branching->best_primal = son_left->PB;
                }
                if( (son_right->PB - list_branching->best_primal) > MEWCP_EPSILON)
                {
                    list_branching->best_primal = son_right->PB;
                }

                /* I insert new nodes into the branching list if possible */
                if ( (list_branching->best_primal - son_left->PB ) > MEWCP_EPSILON)
                {
                    MEWCP_close_open_node(son_left);
                }
                else
                {
                    MEWCP_push_open_node(son_left,list_branching);
                }

                if ( (list_branching->best_primal - son_right->PB ) > MEWCP_EPSILON)
                {
                    MEWCP_close_open_node(son_right);
                }
                else
                {
                    MEWCP_push_open_node(son_right,list_branching);
                }

            }
            else
            {
                MEWCP_close_open_node(open_node_worst_bound);
            }

        }

        /* I can close the node! */

        MEWCP_close_open_node(open_node_worst_bound);



    }














    MEWCP_free_list_branching(list_branching);


}

void MEWCP_close_open_node(open_node_t * open_node)
{
#if defined MEWCP_DSDP_VERBOSE1
    printf("-- Close node: %d\n",open_node->id_node);
#endif

    MEWCP_free_open_node(open_node);
}


void MEWCP_bound(open_node_t * open_node, double ** constraints_matrix, double * bi,
                 const unsigned int num_constraints,
                 const unsigned int dim_matrix,
                 const unsigned int num_nodes,
                 const unsigned int num_partitions)
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
            SDPConeSetADenseVecMat(sdpcone, 0, i, num_nodes, MEWCP_ALPHA , constraints_matrix[i], dim_matrix);
        }
        else /* It's the last matrix containing branching contraints */
        {
            SDPConeSetADenseVecMat(sdpcone, 0, i, num_nodes, 1 , open_node->vect_mat_branching_contraint, dim_matrix);
        }


#if defined MEWCP_DSDP_DEBUG

        SDPConeViewDataMatrix(sdpcone, 0, i);
#endif

    }

    /* set DSDP parameters */
    info=DSDPSetGapTolerance(dsdp,MEWCP_GAP_TOLERANCE);
    info=DSDPSetPotentialParameter(dsdp,MEWCP_POTENTIAL_PARAMETER);
    info=DSDPReuseMatrix(dsdp,MEWCP_REUSE_MATRIX);
    info=DSDPSetPNormTolerance(dsdp,MEWCP_SET_PNORM_TOLERANCE);

#if defined MEWCP_CONVERTER_DSDP_VERBOSE2

    DSDPSetStandardMonitor(dsdp, 1); /* verbose each iteration */
    DSDPLogInfoAllow(1,0);
#endif


    DSDPSetup(dsdp);
    DSDPSolve(dsdp);
    DSDPComputeX(dsdp);

    DSDPStopReason(dsdp, &reason);
    DSDPGetPObjective(dsdp, &pobj);
    DSDPGetTraceX(dsdp, &sol_traceX);

    /* I take the negative pobj */
    pobj = -pobj;

#if defined MEWCP_DSDP_DEBUG

    double * sol_vect_X;
    int sol_dim_vect_X;
    double vect_y[num_constraints ];

    SDPConeGetXArray(sdpcone, 0, &sol_vect_X, &sol_dim_vect_X);

    printf("\n");
    SDPConeViewX(sdpcone, 0, num_nodes, sol_vect_X, sol_dim_vect_X);
    printf("\n");

    DSDPGetY (dsdp,  vect_y, num_constraints);
    MEWCP_print_vectorY(vect_y,num_constraints);
#endif


    /* Now I get the diagonal of solution X */

    open_node->diagX = MEWCP_allocate_diag_X(num_nodes);
    MEWCP_dump_diag_X(&sdpcone,open_node->diagX,num_nodes);

#if defined MEWCP_DSDP_DEBUG

    MEWCP_print_diag_X(open_node->diagX,num_nodes);
#endif


    /* I compute the rounding of diagonan_X */
    open_node->list_nodes_solution = MEWCP_allocate_list_nodes_solution(num_partitions);
    MEWCP_compute_sdp_rounding(open_node->diagX,open_node->list_nodes_solution,num_nodes,cardinality_partition);
    MEWCP_print_list_nodes_solution(open_node->list_nodes_solution,num_partitions);

    z_rouded = MEWCP_evaluate_list_nodes_solution(open_node->list_nodes_solution,constraints_matrix[0],num_partitions);

#if defined MEWCP_DSDP_VERBOSE1

    printf("(P) Obj_value: %.2lf \t Rounded: %.2lf \t Trace(X): %.2lf \t Termination value: %d\n",pobj,z_rouded, sol_traceX,reason);
#endif




    /* Update open_node */
    open_node->DB = pobj;
    open_node->PB = z_rouded;


    DSDPDestroy ( dsdp);


}

bool MEWCP_branch( open_node_t * open_node,
                   const unsigned int dim_matrix,
                   const unsigned int num_nodes,
                   const unsigned int num_partitions,
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
    possible_branch = MEWCP_generate_equi_branch_node(open_node->diagX,num_nodes,num_partitions,&out_num_part,&out_id_node);

#if defined MEWCP_DSDP_DEBUG

    printf("Possible branch: %d \t partition: %d \t id_node: %d \n",possible_branch, out_num_part, out_id_node);
#endif

    if (possible_branch == true) /* Branching is possible */
    {

        son_left = MEWCP_allocate_open_node();
        son_right = MEWCP_allocate_open_node();

        son_left->id_node = 2*open_node->id_node +1;
        son_right->id_node = 2*open_node->id_node +2;

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

        /* sdpcone not needed by the function, waiting for a proper printing funcion
        #if defined  MEWCP_DSDP_DEBUG
         
                printf("Left: \n");
                SDPConeViewX(sdpcone, 0, num_nodes, son_left->vect_mat_branching_contraint,dim_matrix );
                printf("Right: \n");
                SDPConeViewX(sdpcone, 0, num_nodes,son_right->vect_mat_branching_contraint, dim_matrix);
        #endif
        */
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

    branching_open_node_t * worst_element;
    branching_open_node_t * element;
    double worst_bound;

    worst_bound = MEWCP_MIN_DOUBLE;

#if defined ASSERT

    assert(MEWCP_is_list_branching_empty(list_branching) == false);
#endif

    for(element = list_branching->head; element != NULL; element = element->next_branching_node)
    {
        if( (element->open_node->DB -  worst_bound) > MEWCP_MIN_DOUBLE )
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
    printf("* MEWCP_pop_specific_open_node *\n");
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
    double max_sum_tmp;


    double sum_prev, sum_cur;

    c = n/m;

    for(i=0;i<m;++i)
    {
        sum_cur = 0;
        sum_prev =0;

        for(j=i*c;j<(i*c +c); ++j )
        {
            sum_cur = sum_prev + diag_X[j];

            if ( (1- diag_X[j] <=  MEWCP_EPSILON) ) /* only one element not zero I cannot branch! */
            {
#if defined MEWCP_CONVERTER_DSDP_DEBUG
                printf("Partition: %d only one node is 1\n",i);
#endif

                id_node_part[i] = -1;
                sum_partition[i] = -1;
                break;

            }
            else if ( (sum_cur - 0.5 > MEWCP_EPSILON) && j != (i*c +c -1) ) /* OK I have 1/2 and more considering not all the elements */
            {
                id_node_part[i] = j;
                sum_partition[i] = sum_cur;
            }
            else if ( (sum_cur - 0.5 > MEWCP_EPSILON) && j == (i*c +c -1) ) /* BRRR I have 1/2 and more considering not the elements */
            {
                id_node_part[i] = j-1;
                sum_partition[i] = sum_prev;
            }


            sum_prev = sum_cur;
        }
    }

    max_sum_tmp = 0;
    for (i=0;i<m;++i)
    {
        if (  (sum_partition[i] > max_sum_tmp) && ( (sum_partition[i] - 1.0)<MEWCP_EPSILON) )
        {
            max_sum_tmp = sum_partition[i];


            *out_id_node = id_node_part[i];
            *out_num_part = i;
        }
    }
#if defined MEWCP_CONVERTER_DSDP_DEBUG

    for(i=0;i<m;++i)
    {
        printf("Partition: %d \t Sum: %.3lf \t id_node: %d\n",i, sum_partition[i], id_node_part[i]);
    }

#endif



    if ( max_sum_tmp < MEWCP_EPSILON)
    {
        return false;
    }
    else
    {
        return true;
    }
}


void MEWCP_add_blocked_node(const unsigned int id_node, list_blocked_nodes_t * list_blocked_nodes)
{
    unsigned int pos_insertion;

    pos_insertion = list_blocked_nodes->num_blocked_nodes;

#if defined DEBUG

    assert(list_blocked_nodes->bool_list[id_node] == false);
#endif

    list_blocked_nodes->bool_list[id_node] = true;
    list_blocked_nodes->blocked_node[pos_insertion] = id_node;
    list_blocked_nodes->num_blocked_nodes += 1;

    /* All structures are updated */
}

double MEWCP_evaluate_list_nodes_solution(unsigned int * list_node_solution, double * vect_matrix_weights,const unsigned int m)
{
#if defined MEWCP_CONVERTER_DSDP_VERBOSE1
    printf("* MEWCP_evaluate_list_nodes_solution *\n");
#endif

    unsigned int i,j;
    unsigned int pos_vect;
    double z_tmp;

    z_tmp = 0;

    for (i=0;i<m;++i)
    {
        for(j=i; j<m;++j)
        {
            if(list_node_solution[i] >= list_node_solution[j])
            {
                pos_vect = MEWCP_convert_coords_ij_to_vector_matrix(list_node_solution[i]+1, list_node_solution[j]+1);
            }
            else
            {
                pos_vect = MEWCP_convert_coords_ij_to_vector_matrix(list_node_solution[j]+1, list_node_solution[i]+1);
            }

            if (i==j)
            {
                z_tmp += vect_matrix_weights[pos_vect];
            }
            else /* I have to consider 2* because of in the W matris non diagonal element are considered /2 */
            {
                z_tmp += 2*vect_matrix_weights[pos_vect];
            }

        }
    }
#if defined MEWCP_CONVERTER_DSDP_DEBUG
    printf("Objective value of node_list: %.2lf\n",-z_tmp);
#endif

    return -z_tmp;
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
    open_node->DB = (double) 0;
    open_node->PB = (double) 0;

    open_node->list_nodes_solution = NULL;
    open_node->diagX = NULL;
    open_node->list_blocked_nodes = NULL;
    open_node->vect_mat_branching_contraint = NULL;

    return open_node;

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

    list_branching->branching_open_node_t = NULL;
    list_branching->head = NULL;
    list_branching->tail = NULL;

    list_branching->number_open_nodes = 0;

    return list_branching;


}


double ** MEWCP_allocate_sdp_constraints_matrix(const unsigned int n, const unsigned int num_constraints )
{
#if defined MEWCP_CONVERTER_DSDP_VERBOSE1
    printf("* MEWCP_allocate_sdp_constraints_matrix *\n");
#endif

    unsigned int i;
    unsigned int dim_vect_matrix;

    dim_vect_matrix = n*(n+1)/2;

    double ** sdp_constraints_matrix;

    /* I allocate the constraints matrix that is card(b) -1
    * I don't consider the last matrix of branching constraints */
    sdp_constraints_matrix = (double **) calloc(num_constraints , sizeof(double *));
    if (sdp_constraints_matrix == NULL)
    {
        printf("Allocation ERROR: sdp_constraints_matrix\n");
        exit (EXIT_FAILURE);
    }

    for (i=0; i < num_constraints ; ++i )
    {
        sdp_constraints_matrix[i] = (double *) calloc(dim_vect_matrix, sizeof(double));
        if (sdp_constraints_matrix[i] == NULL)
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

double * MEWCP_allocate_vect_mat_branching_constraints(const unsigned int dim_matrix)
{
    double * vect_mat_branching_contraints;

    vect_mat_branching_contraints = (double *) calloc(dim_matrix, sizeof(double));
    if (vect_mat_branching_contraints  == NULL)
    {
        printf("*******  Allocation MEWCP_allocate_vect_mat_branching_constraints FAILED!\n");
        exit(EXIT_FAILURE);
    }
    return vect_mat_branching_contraints;
}

unsigned int * MEWCP_allocate_list_nodes_solution( const unsigned int m)
{
    unsigned int * list_nodes_solution;
    list_nodes_solution = (unsigned int *) calloc (m, sizeof(unsigned int));
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

    /* I copu the list */
    memcpy(list_cloned->blocked_node, list_to_be_cloned->blocked_node, num_nodes * sizeof(unsigned int) );
    memcpy(list_cloned->bool_list, list_to_be_cloned->bool_list, num_nodes * sizeof(bool));

    list_cloned->num_blocked_nodes = list_to_be_cloned->num_blocked_nodes;
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

void MEWCP_compute_sdp_rounding(double * diag_X, unsigned int * list_nodes_rounded_solution, const unsigned int num_nodes, const unsigned int c)
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
        printf("%lf ",diag_X[i]);
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

void MEWCP_print_list_nodes_solution(unsigned int * list_nodes_solution, const unsigned int m)
{
    unsigned int i;
    printf("Nodes list solution: : ");

    for (i=0;i<m;++i)
    {
        printf("%u ",list_nodes_solution[i]);
    }
    printf("\n");
}

/* ******************************************
 * FREE FUNCTIONS 
 * ******************************************/
void MEWCP_free_list_blocked_nodes(list_blocked_nodes_t * list_blocked_nodes, const unsigned int num_nodes )
{
#if defined MEWCP_DSDP_DEBUG
    printf("* MEWCP_free_list_blocked_nodes *\n");
#endif

    free( list_blocked_nodes->blocked_node);
}

void MEWCP_free_sdp_constraints_matrix(double ** sdp_constraints_matrix, const unsigned int num_constraints )
{
#if defined MEWCP_DSDP_DEBUG
    printf("* MEWCP_free_sdp_constraints_matrix *\n");
#endif

    unsigned int i;
    for(i=0;i<num_constraints;++i) /* matrix are num_contraints +1 but we do not consider the branching constraint matrix */
    {
        free(sdp_constraints_matrix[i]);
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

void MEWCP_free_list_nodes_solution(unsigned int * list_nodes_solution)
{
#if defined MEWCP_DSDP_DEBUG
    printf("* MEWCP_free_list_nodes_solution *\n");
#endif

    free(list_nodes_solution);
}

void MEWCP_free_vect_mat_branching_constraints(double * vect_mat_branching_contraints)
{
#if defined MEWCP_DSDP_DEBUG
    printf("* MEWCP_free_vect_mat_branching_constraints *\n");
#endif

    free(vect_mat_branching_contraints);
}

void MEWCP_free_open_node(open_node_t * open_node)
{
#if defined MEWCP_DSDP_DEBUG
    printf("* MEWCP_free_open_node *\n");
#endif

    free(open_node->list_blocked_nodes);
    free(open_node->vect_mat_branching_contraint);

    free(open_node->diagX);
    free(open_node->list_nodes_solution);

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
        exit(EXIT_FAILURE);
    }

    /* As the list is empty I simply free the list structure */
    free(list_branching);

}
