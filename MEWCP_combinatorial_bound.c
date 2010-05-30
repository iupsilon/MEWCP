/*****************************************************************************
 *  Project: Maximum Edge Weighted Clique Problem with multiple choice contraints
 * 
 *  Authors:
 *  (c) 2009 Yari Melzani (yari.melzani@gmail.com) 
 *
 ****************************************************************************/

#include "MEWCP_combinatorial_bound.h"
#include "MEWCP_tabu.h"
#include "converter_dsdp.h"

void MEWCP_compute_combinatorial_preprocessing(open_node_t * root_node, matrix_weights_t * matrix_weights,
        const unsigned int num_partitions,
        const unsigned int cardinality_partitions,
        const double initial_best_primal)
{
#if defined MEWCP_COMBINATORIAL_PREPROCESSING_DEBUG
    printf("\n* MEWCP_compute_combinatorial_preprocessing *\n");
#endif

    unsigned int i,k,n,dim_vect;
    n = num_partitions * cardinality_partitions;
    double vect_dual[n];  // I save the DB for each iteration

    list_blocked_nodes_t * current_list_blocked_nodes;
    list_blocked_nodes_t * final_list_blocked_nodes;
    open_node_t * current_open_node;
    double best_primal;

    best_primal = initial_best_primal;
    final_list_blocked_nodes = MEWCP_allocate_list_blocked_nodes(num_partitions*cardinality_partitions);

    for (k = 0; k<num_partitions; ++k)
    {
        for (i = k*cardinality_partitions; i <  (k*cardinality_partitions + cardinality_partitions); ++i)
        {
            current_list_blocked_nodes = MEWCP_blocked_other_nodes_of_partition(i,k,num_partitions,cardinality_partitions);
            current_open_node = MEWCP_allocate_open_node();

            current_open_node->list_blocked_nodes = current_list_blocked_nodes;

            MEWCP_bound_combinatorial(current_open_node,matrix_weights,num_partitions,cardinality_partitions,best_primal);

            vect_dual[i] = current_open_node->DB;
#if defined MEWCP_COMBINATORIAL_PREPROCESSING_DEBUG

            printf("---> node: %d \t DB: %.2lf\n",i,vect_dual[i]);
#endif

            if ( (current_open_node->PB - best_primal) >  MEWCP_EPSILON)
            {
                best_primal = current_open_node->PB;
#if defined MEWCP_COMBINATORIAL_PREPROCESSING_DEBUG

                printf("Preprocessing New best PB k: %d i: %d\t PB: %.2lf\n",k,i,best_primal);
#endif

            }

            /* all done, I can deallocate the current node */
            MEWCP_free_list_blocked_nodes(current_list_blocked_nodes);
            free(current_open_node); /* It contains no vectors */


        }
    }
#if defined MEWCP_COMBINATORIAL_PREPROCESSING_VERBOSE1

    printf("\n Best_primal: %.2lf\n",best_primal);
#endif

    /* now let's check if the node is about to be close
                 * in this case I add the only free element of the partition to the final_blocked_list */
    for (i=0;i<n;++i)
    {
        if ((best_primal - vect_dual[i] ) > MEWCP_EPSILON  )
        {
#if defined MEWCP_COMBINATORIAL_PREPROCESSING_DEBUG
            printf("XXX Preprocessing Node to be close! i: %d\t DB: %.2lf\n",i,vect_dual[i]);
#endif
            /* I add the node i to the list of final blocked nodes */
            MEWCP_add_blocked_node(i,final_list_blocked_nodes);

        }
        else
        {
#if defined MEWCP_COMBINATORIAL_PREPROCESSING_DEBUG
            printf("OOOOO Preprocessing Node survives:  i: %d\t DB: %.2lf\n",i,vect_dual[i]);
#endif

        }
    }


    /* Before I go I modify the root node given with the final blocked list nodes and
     * I generate the correspondant branching constraint according to the blocked list
     */

#if defined MEWCP_COMBINATORIAL_PREPROCESSING_VERBOSE1

    printf("Preprocessing has blocked the following nodes:\n");
    MEWCP_print_list_blocked_nodes(final_list_blocked_nodes);
#endif

    /* I free the root node fields i have to replace */
    MEWCP_free_list_blocked_nodes(root_node->list_blocked_nodes);
    MEWCP_free_vect_mat_branching_constraints(root_node->vect_mat_branching_contraint);


    root_node->list_blocked_nodes = final_list_blocked_nodes;


    dim_vect = n*(n+1)/2;
    MEWCP_allocate_vect_mat_branching_constraints(dim_vect);
    MEWCP_generate_constraints_branch(final_list_blocked_nodes,root_node->vect_mat_branching_contraint,dim_vect,n,cardinality_partitions);


}



list_blocked_nodes_t * MEWCP_blocked_other_nodes_of_partition(const unsigned int pos,
        const unsigned int partition,
        const unsigned int num_partitions,
        const unsigned int cardinality_partitions)
{
    unsigned int i;

    list_blocked_nodes_t * list_blocked_nodes;
    list_blocked_nodes = MEWCP_allocate_list_blocked_nodes(num_partitions*cardinality_partitions);

    for (i = partition*cardinality_partitions; i < ( partition*cardinality_partitions + cardinality_partitions); ++i)
    {
        if (i != pos)
        {

            MEWCP_add_blocked_node(i,list_blocked_nodes);
#if defined  MEWCP_COMBINATORIAL_PREPROCESSING_DEBUG

            //printf("other node part: %d\tblocked: %d\n",partition,i);
#endif

        }
    }



    return list_blocked_nodes;

}


void MEWCP_bound_combinatorial(open_node_t * open_node,
                               matrix_weights_t * matrix_weights,
                               const unsigned int num_partitions,
                               const unsigned int cardinality_partitions,
                               const double best_primal)
{
#if defined MEWCP_BOUNDING_COMBINATORIAL_DEBUG
    printf("\n* MEWCP_bound_combinatorial *\n");
#endif

    unsigned int i,j,k,t,z,r;
    double z_primal;
    double sum_i,sum_j;
    double total_sum = 0;
    double sum_tmp;

    /* Vector of selected nodes for each partition Vk, k=1,..,m */
    int selected_node[num_partitions];
    double value_k[num_partitions];

    int tmp_selected_edges[num_partitions];

    /* for each partition Vk I say which are the nodes connected to selected node in Vk */
    int selected_edges[num_partitions][num_partitions];

    for(i=0;i<num_partitions;++i)
    {
        for(j=0;j<num_partitions;++j)
        {
            selected_edges[i][j] = -1;
        }
    }



#if defined MEWCP_BOUNDING_COMBINATORIAL_DEBUG
    printf("m: %d\tc: %d\n",num_partitions, cardinality_partitions);
    MEWCP_print_list_blocked_nodes(open_node->list_blocked_nodes);
#endif

    for ( k = 0; k<num_partitions; ++k)
    {

        sum_i = MEWCP_MIN_DOUBLE;
#if defined MEWCP_BOUNDING_COMBINATORIAL_DEBUG

        printf("Azzero sum_i\n");
#endif

        for(i = cardinality_partitions*k; i< (cardinality_partitions*k +cardinality_partitions); ++i)
        {
            if (open_node->list_blocked_nodes->bool_list[i] == true)
            {
                //#if defined MEWCP_BOUNDING_COMBINATORIAL_DEBUG
                //              printf("Nodo i: %d bloccato\n",i);
                //#endif

                continue;
            }
#if defined MEWCP_BOUNDING_COMBINATORIAL_DEBUG
            printf("i: %d\t k: %d\n",i,k);
#endif

            for(z=0;z<num_partitions;++z)
            {
                tmp_selected_edges[z] = -1;
            }

            for(t=((k+1)%num_partitions), r=0;  r<num_partitions-1; t = ((t+1) % num_partitions), ++r)
            {

#if defined ASSERT
                assert(t != k);
#endif

                sum_j = MEWCP_MIN_DOUBLE;

                for(j = cardinality_partitions*t ; j<(cardinality_partitions*t + cardinality_partitions); ++j)
                {

                    if (open_node->list_blocked_nodes->bool_list[j] == true)
                    {
                        //#if defined MEWCP_BOUNDING_COMBINATORIAL_DEBUG
                        //                        printf("Nodo j: %d bloccato\n",j);
                        //#endif

                        continue;
                    }
                    /* I keep the maximum weight edge*/

                    sum_tmp = matrix_weights->weight[i][j];
                    sum_tmp += (matrix_weights->weight[j][j])/(num_partitions -1);


                    if ( (sum_tmp - sum_j) > MEWCP_EPSILON)
                    {
                        sum_j = sum_tmp;
                        tmp_selected_edges[t] = j;
                    }
            
                }
#if defined MEWCP_BOUNDING_COMBINATORIAL_DEBUG
                printf("Partizione k: %d t: %d \tnode: %d\t sum_tmp: %.2lf\n",k,t, tmp_selected_edges[t],sum_tmp);
#endif


            }

            /* Before I change i I want to know the value of this sum */
            sum_tmp = 0;
            for(z=0;z<num_partitions;++z)
            {
                if(tmp_selected_edges[z] != -1)
                {
                    sum_tmp += matrix_weights->weight[i][tmp_selected_edges[z]];
                    sum_tmp += matrix_weights->weight[tmp_selected_edges[z]][tmp_selected_edges[z]]/(num_partitions -1);
                }
            }
            sum_tmp += matrix_weights->weight[i][i];

            if ((sum_tmp - sum_i) > MEWCP_EPSILON)
            {
                selected_node[k] = i;
                for(z=0;z<num_partitions;++z)
                {
                    selected_edges[k][z] = tmp_selected_edges[z];
                }
                sum_i = sum_tmp;
                value_k[k] = sum_i;
#if defined MEWCP_BOUNDING_COMBINATORIAL_DEBUG

                printf("Update sum_i: %.2lf\n",sum_i);
#endif

            }

#if defined MEWCP_BOUNDING_COMBINATORIAL_DEBUG
            printf("Partition:  k: %d \t taken nodo: %d\t value: %.lf \t best_va: %.2lf\n Edges:",k, selected_node[k],sum_tmp, sum_i);
            for(z=0;z<num_partitions;++z)
            {
                printf("%d ",selected_edges[k][z]);
            }
            printf("\n\n");
#endif

        }

    }


    /* Now let's evaluate the value */
    total_sum = 0;
    for (i=0;i<num_partitions;++i)
    {
        total_sum += value_k[i];
    }
    total_sum = total_sum / 2;

#if defined MEWCP_BOUNDING_COMBINATORIAL_DEBUG
    /* Let's see nodes and edges taken */


    printf("Taken nodes: ");
    for (i=0;i<num_partitions;++i)
    {
        printf("%d ",selected_node[i]);
    }
    printf("\n");

    for(i=0; i<num_partitions;++i)
    {
        printf("Partition: %d \t value: %.2lf\nEdges:",i,value_k[i]);
        for(j=0;j<num_partitions;++j)
        {
            printf(" %d",selected_edges[i][j]);
        }
        printf("\n");
    }

#endif


    z_primal = MEWCP_evaluate_list_nodes_solution(selected_node,matrix_weights,num_partitions);


#if defined MEWCP_BOUNDING_COMBINATORIAL_DEBUG

    printf("Combinatorial: %.2lf \t Value Z primal: %.2lf\n",total_sum,z_primal);
#endif

#if defined MEWCP_BOUNDING_VERBOSE1

    printf("(BB)  Bound Combinatorial: (%d) DB: %.2lf \t PB: %.2lf \t level: %u\n",open_node->serial_node, total_sum, z_primal, open_node->depth_level);
#endif

#if defined ASSERT
    assert((total_sum - z_primal) >= 0);
#endif
    /* Update values DB PB */

    open_node->DB_comb = total_sum;
    open_node->PB_comb = z_primal;

    if ((open_node->DB - total_sum) > MEWCP_EPSILON )
    {
        open_node->DB = total_sum;
    }

    if ( (z_primal - open_node->PB) > MEWCP_EPSILON )
    {
        open_node->PB = z_primal;
        open_node->list_nodes_solution = MEWCP_allocate_list_nodes_solution(num_partitions);
        MEWCP_clone_list_nodes_solution(selected_node,open_node->list_nodes_solution,num_partitions);
    }


}


