/*****************************************************************************
 *  Project: Maximum Edge Weighted Clique Problem
 * 
 *  Authors:
 *  (c) 2009 Yari Melzani (yari.melzani@gmail.com) 
 *
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <string.h>


#include "MEWCP_tabu_definitions.h"
#include "MEWCP_tabu.h"

tabu_result_t MEWCP_compute_tabu_search(const unsigned int num_iterations, matrix_weights_t * matrix_weights, node_list_t * node_list)
{
#if defined MEWCP_TABU_VERBOSE1
    printf(".: Compute Tabu Search :.\n");
#endif

    unsigned int i;
    unsigned int counter_last_improvement;
    iteration_t last_improvement_iteration;

    solution_t best_solution;
    list_iterations_t  list_iterations;
    tabu_node_list_t  tabu_node_list;

    tabu_result_t tabu_result;
    weight_t Z_iteration_tabu;
    weight_t Z_best;


    tabu_result.last_improvement_iteration = 0;
    tabu_result.last_improvement_time = (double) 0;
    MEWCP_create_solution(matrix_weights,&tabu_result.solution);

    MEWCP_initialize_tabu_structures(matrix_weights, num_iterations,&list_iterations,&best_solution,&tabu_node_list);


    MEWCP_dump_tabu_solution(node_list,&best_solution);


    Z_best = best_solution.Z;
    list_iterations.Z_iteration[0] = node_list->Z;
    list_iterations.current_iteration = 0;
    counter_last_improvement = 0;

    for (i = 1; i < num_iterations; ++i)
    {
        list_iterations.current_iteration = i;
#if defined MEWCP_TABU_DEBUG1

        MEWCP_print_node_list(node_list);
#endif

        if (counter_last_improvement < MAX_WORSENING_ITERATIONS )
        {
            Z_iteration_tabu = MEWCP_compute_iteration_tabu_search (matrix_weights,node_list,&best_solution,&list_iterations,&tabu_node_list);

            if (Z_iteration_tabu - Z_best > MEWCP_DBL_EPSILON) /* Zbest Improved! */
            {
                counter_last_improvement = 0;
                Z_best = Z_iteration_tabu;
                last_improvement_iteration = i;
                MEWCP_dump_tabu_solution(node_list,&best_solution);

#if defined MEWCP_TABU_VERBOSE1

                printf("\t*** Z_best Improved!\tIter: %d, Z: %d *\n",i,Z_best);
#endif

            }
            else
            {
                counter_last_improvement +=1;  /* Not improved */
            }
        }
        else
        {
            break; /* No improvents for MAX_WORSENING_ITERATIONS */
        }
    }

    /* OK, all done, now I have to return */

    tabu_result.last_improvement_iteration = last_improvement_iteration;
    MEWCP_clone_solution(&best_solution, &tabu_result.solution, matrix_weights->m);
    MEWCP_free_solution(&best_solution);
    MEWCP_free_tabu_node_list(&tabu_node_list);

    return tabu_result;


}

weight_t MEWCP_compute_iteration_tabu_search(	matrix_weights_t * matrix_weights,
        node_list_t * node_list ,
        solution_t * best_solution,
        list_iterations_t * list_iterations,
        tabu_node_list_t * tabu_node_list )
{

#if defined MEWCP_TABU_DEBUG1
    printf("\n.: Iteration: %d :.\n", list_iterations->current_iteration);
#endif

    weight_t Z_prev;
    weight_t Z_current;
    weight_t Z_best;
    weight_t Z_swap;


    int m,c;
    int i,node_i,node_j,partition;
    pointer_node_t n1,n2;

    m = matrix_weights->m;
    c = matrix_weights->c;

    n1=NULL_POINTER;
    n2=NULL_POINTER;

    Z_current = (weight_t) MEWCP_MAX_NEG_WEIGHT;
    Z_prev = node_list->Z;
    Z_best = best_solution->Z;

    /* For each partition I evaluate the swap of element in solution with others of the same partition */
    for ( i = 0; i< matrix_weights->m; ++i)
    {

        node_i = node_list->selected_node_partition[i];
        partition = MEWCP_get_node_partition(node_i,matrix_weights); /* I get the node selected in the partition i-th */

        for(node_j = c*partition; node_j< (c*partition +c); ++node_j) /* I try swap with nodes of same partitions */
        {
            if (node_j != node_i )
            {
                Z_swap = MEWCP_evaluate_swap_nodes(node_i,node_j,matrix_weights,node_list);
#if defined MEWCP_TABU_DEBUG1

                printf("+++ Swap gain: [%d,%d]-> %d\n",node_i,node_j,Z_swap);
#endif


                if ( ((Z_swap - Z_current) > MEWCP_DBL_EPSILON ) && ( MEWCP_is_tabu_in(node_j,tabu_node_list,list_iterations) == false )  &&  ( MEWCP_is_tabu_out(node_i,tabu_node_list,list_iterations) == false ) )
                {
                    Z_current = Z_swap;
                    n1 = node_i;
                    n2 = node_j;

                    if ((Z_swap - Z_best) > MEWCP_DBL_EPSILON ) // New Z_best!

                    {

                        Z_best = Z_swap;


#if defined MEWCP_TABU_VERBOSE2

                        (void) printf("Z_current: %d\tNew Z_best: %d\n",Z_current,Z_best );
#endif

                    }

                }
                else if ( ( (Z_swap - Z_current) > MEWCP_DBL_EPSILON ) &&
                          ( ( MEWCP_is_tabu_in(node_j,tabu_node_list,list_iterations) == true )  ||  ( MEWCP_is_tabu_out(node_i,tabu_node_list,list_iterations) == true ) ) )
                {
#if defined MEWCP_TABU_DEBUG1
                    printf("--- Nodes: [%d,%d] not taken because of tabu: IN: %d\t OUT: %d\n",node_i,node_j,
                           (list_iterations->current_iteration - tabu_node_list->tabu_node_state[node_i].iteration_in),
                           (list_iterations->current_iteration- tabu_node_list->tabu_node_state[node_j].iteration_out) );
#endif

                    continue;
                }
                else  /* The swap is not good */
                {
                    continue;
                }
            }
        }
    }

    /* I have not found a suitable SWAP because all swaps are TABU */
    if ( (n1 == NULL_POINTER) || (n2 == NULL_POINTER))
    {

        //#if defined MEWCP_TABU_VERBOSE1
        printf("\n\t\t### ALL swaps are tabu!!!  ### \n");
        getchar();
        //#endif

        list_iterations->Z_iteration[list_iterations->current_iteration] = Z_current;
        return Z_current;
    }


    /* I handle the condition in which all swaps are tabu */
    /*
    #if defined ASSERT
        assert(n1 != NULL_POINTER);
        assert(n2 != NULL_POINTER);
    #endif
    */


    /* Now I perform the swap: n1 exits and n2 enters the solution */
    MEWCP_remove_node_from_solution_list(n1,matrix_weights,node_list);
    MEWCP_add_node_to_solution_list(n2,matrix_weights,node_list);

    tabu_node_list->tabu_node_state[n1].iteration_out = list_iterations->current_iteration;
    tabu_node_list->tabu_node_state[n2].iteration_in = list_iterations->current_iteration;
    list_iterations->Z_iteration[list_iterations->current_iteration] = Z_current;


#if defined MEWCP_TABU_VERBOSE1

    printf("%u Swap out: %u \tin: %u \tZ_prev: %d \tZ_current= %d  DeltaZ= %d\tDelta_Zbest= %d\n",
           list_iterations->current_iteration,
           n1,n2,Z_prev,Z_current,
           Z_current-Z_prev,
           Z_best-best_solution->Z);
#endif

    return Z_current;

}


bool MEWCP_is_tabu_in(const pointer_node_t node, tabu_node_list_t * tabu_node_list, list_iterations_t * list_iterations)
{
    if (tabu_node_list->tabu_node_state[node].iteration_out == 0)
    {
        return false;
    }
    else if( (list_iterations->current_iteration - tabu_node_list->tabu_node_state[node].iteration_out) > TABU_IN_ITERATIONS)
    {
        return false;
    }
    else
    {
        return true;
    }
}


bool MEWCP_is_tabu_out(const pointer_node_t node, tabu_node_list_t * tabu_node_list, list_iterations_t * list_iterations)
{
    if (tabu_node_list->tabu_node_state[node].iteration_in == 0)
    {
        return false;
    }
    else if( (list_iterations->current_iteration - tabu_node_list->tabu_node_state[node].iteration_in) > TABU_OUT_ITERATIONS)
    {
        return false;
    }
    else
    {
        return true;
    }
}

weight_t MEWCP_evaluate_swap_nodes(const pointer_node_t n1, const pointer_node_t n2, matrix_weights_t * matrix_weights, node_list_t * node_list)
{
    weight_t Z_swap;

#if defined ASSERT

    assert(node_list->node[n1].belongsM == true );
#endif

    Z_swap = node_list->Z;

    Z_swap  -= node_list->node[n1].sum_din;
    Z_swap -= matrix_weights->weight[n1][n1];
    Z_swap  += node_list->node[n2].sum_din;
    Z_swap += matrix_weights->weight[n2][n2];

    Z_swap -= matrix_weights->weight[n1][n2];


    return Z_swap;
}

void MEWCP_clone_solution(solution_t * src_solution, solution_t * dest_solution, const unsigned int m)
{
    dest_solution->Z = src_solution->Z;
    memcpy(dest_solution->node_solution, src_solution->node_solution, m * sizeof(pointer_node_t));
}

void MEWCP_dump_tabu_solution(node_list_t * node_list, solution_t * solution)
{
	#if defined MEWCP_TABU_DEBUG1
	printf ("\t\tDDDDD DUMP best solution\n");
	#endif   
    unsigned int i;
    pointer_node_t pos=0;

    solution->Z = node_list->Z;

    for(i = node_list->M_head ; i != (unsigned) NULL_POINTER ; i = node_list->node[i].next)
    {
        solution->node_solution[pos] = i;
        ++pos;
    }

}

void MEWCP_create_matrix_weights(const unsigned int n, const unsigned int m, matrix_weights_t * matrix_weights )
{
    unsigned int i;

    matrix_weights->n = n;
    matrix_weights->m = m;
    matrix_weights->c = n/m;

    matrix_weights->weight = (weight_t **) malloc(sizeof(weight_t *) * n);
    for(i = 0; i< n; ++ i)
    {
        matrix_weights -> weight[i] = (weight_t *) malloc(sizeof(weight_t) * n);
    }
}

void MEWCP_free_matrix_weights(matrix_weights_t * matrix_weights)
{
    unsigned int i;

    for(i = 0; i < matrix_weights->n; ++ i)
    {
        free( matrix_weights -> weight[i]);
    }
    free(matrix_weights->weight);
    matrix_weights -> n=0;
    matrix_weights -> m=0;
}

void MEWCP_free_node_list(node_list_t * node_list)
{
    free(node_list->node);
    free(node_list->selected_node_partition);
    node_list->node = NULL;
    node_list->selected_node_partition = NULL;

}

void MEWCP_load_AMPL_instance(char * filename_AMPL, matrix_weights_t * matrix_weights)
{
    unsigned int n,m,c;
    unsigned int i,j;
    weight_t weight_ij;
    char  notcare[6];
    char  notcare1[2];
    char  notcare2[3];
    char  notcare3[2];
    char  notcare4[12]; // Reads until [10000,10000]
    FILE * file_AMPL;


    if ((file_AMPL = fopen(filename_AMPL,"r")) == NULL)
    {
        printf("Open AMPL file %s failed!\n",filename_AMPL);
        exit(EXIT_FAILURE);
    }


    // reading n,m
    (void) fscanf(file_AMPL,"%s %s %s %u %s\n",notcare,notcare1,notcare2,&n,notcare3);
    (void) fscanf(file_AMPL,"%s %s %s %u %s\n",notcare,notcare1,notcare2,&m,notcare3);

    // I compute c from m,n
    c= n/m;

    // I throw away lines I don't care
    for (i=0; i< (2*c*m +7); ++i)
    {
        fscanf(file_AMPL,"%s\n",notcare4);
#if defined MEWCP_TABU_DEBUG2

        printf("%s\n",notcare4);
#endif

    }


    MEWCP_create_matrix_weights(n,m,matrix_weights);


    // Now I take matrix weights [i,j]
    for(i=0;i<n;++i)
    {
        for(j=0;j<=i;++j)
        {
            (void) fscanf(file_AMPL,"%s %lf\t",notcare4,&weight_ij);
            matrix_weights->weight[i][j] = weight_ij;

#if defined MEWCP_TABU_DEBUG1

            printf("%s ",notcare4);
            printf("%d ",weight_ij);
#endif

            if (i!=j)
            {
                matrix_weights->weight[j][i] = weight_ij;
            }
        }
#if defined MEWCP_TABU_DEBUG1

        printf("\n");
#endif

    }
}


void MEWCP_create_node_list(matrix_weights_t * matrix_weights, node_list_t * node_list)
{
#if defined MEWCP_TABU_VERBOSE1
    printf(".: Create Node List :.\n");
#endif

    unsigned int N = matrix_weights->n;
    unsigned int M = matrix_weights->m;

    node_list->node = (node_t *) malloc (sizeof (node_t) * N);
    node_list->selected_node_partition = (pointer_node_t *) malloc(sizeof(pointer_node_t) * M);

    MEWCP_initialize_node_list(matrix_weights,node_list);
}

void MEWCP_initialize_node_list(matrix_weights_t * matrix_weights, node_list_t * node_list)
{

#if defined MEWCP_TABU_VERBOSE1
    printf(".: Initialize Node List :.\n");
#endif

    unsigned int i,j;
    weight_t sum_temp;
    unsigned int N = matrix_weights->n;
    unsigned int M = matrix_weights->m;

    node_list->Z = (weight_t) 0; // Z is 0
    node_list->N_head = (pointer_node_t) 0;
    node_list->M_head = NULL_POINTER;  // M is empty so the head points to NULL
    node_list->N_tail = (pointer_node_t) (N -1); // punto all'ultimo elemento
    node_list->M_tail = NULL_POINTER; // M is empty so the tail points to NULL

    node_list->card_N = N; // cardinality N
    node_list->card_M = 0; // cardinality M = 0 not yet nodes in M

    /*	I initialize the cursor of selected element for each partition	*/
    for (i = 0; i < M; ++i)
    {
        node_list->selected_node_partition[i]= NULL_POINTER;
    }

    for (i = 0; i < N; ++i)
    {
        node_list->node[i].belongsM = false;
        node_list->node[i].next = i+1;
        node_list->node[i].prev = i-1;

        sum_temp = 0;
        for (j = 0; j < N; ++j)
        {
            /* I consider edge weights and vertex weights*/
            sum_temp += matrix_weights->weight[i][j];

        }
        node_list->node[i].sum_din = 0;
        node_list->node[i].sum_dout = sum_temp;
    }
    // Last element of N point to NULL
    node_list->node[N-1].next = NULL_POINTER;
}

void MEWCP_create_solution(matrix_weights_t * matrix_weights, solution_t * solution)
{
    unsigned int m;

    m= matrix_weights->m;
    solution->Z = (weight_t) 0;
    solution->node_solution = (pointer_node_t *) calloc(m, sizeof(pointer_node_t));

    if (solution-> node_solution == NULL)
    {
        fprintf(stderr,"ERROR: allocation Solution FAILED!\n");
    }

}

void MEWCP_add_node_to_solution_list(pointer_node_t n, matrix_weights_t * matrix_weights, node_list_t * node_list)
{
    unsigned int i,partition;
    weight_t weight_ni;

#if defined ASSERT

    assert(node_list->node[n].belongsM != true);
#endif

    // I remove node n from set N
    if (node_list->node[n].prev == NULL_POINTER) // il node n è il primo della node_list
    {
        // I patch the cursor to the head of N
        node_list->N_head = node_list->node[n].next;
    }
    else // n is not the firs element
    {
        node_list->node[node_list->node[n].prev].next = node_list->node[n].next;
    }
    if (node_list->node[n].next == NULL_POINTER) // se il n è l'ultimo elemento
    {
        // I have to path tail of prev element
        node_list->N_tail = node_list->node[n].prev;
    }
    else // n is not the last element
    {
        node_list->node[node_list->node[n].next].prev = node_list->node[n].prev;
    }

    // I add n to the tail of list M

    // check if M is empty
    if (node_list->M_tail == NULL_POINTER)
    {
#if defined ASSERT
        // if tail points to NULL then head also has to
        assert(node_list->M_head == NULL_POINTER);
#endif

        node_list->M_head = n; // the head of M now points to element
        node_list->M_tail = n; // M tail is now n
        node_list->node[n].prev = NULL_POINTER;
        node_list->node[n].next = NULL_POINTER;

    }
    else // at least one element still exists
    {
        node_list->node[n].prev = node_list->M_tail;
        node_list->node[node_list->M_tail].next = n;
        node_list->M_tail = n;
        node_list->node[n].next = NULL_POINTER;
    }

    // I update cardinality of M,N lists
    node_list->card_N -= 1;
    node_list->card_M += 1;
    // Now I set n belonging to M
    node_list->node[n].belongsM = true;


    // Now I update the sum of weights from each element of N and M
    weight_ni = 0;
    for ( i=0 ; i<matrix_weights->n ; ++ i)
    {

        // Update sum of weights  IN  increase of d(n,i)
        // Update sum of weights  IN  decrease of d(n,i)
        /*if (n == i)
        {
            continue ;
        } */
        weight_ni = matrix_weights->weight[n][i];
        node_list->node[i].sum_din += weight_ni;
        node_list->node[i].sum_dout -= weight_ni;
    }

    /* node n is selected for its partition */
    partition = MEWCP_get_node_partition(n,matrix_weights);
    node_list->selected_node_partition[partition] = n;

    // aggiorno il nuovo valore della funzione_obiettivo e la scrivo in node_list
    (void) MEWCP_update_Z_node_list(node_list,matrix_weights);

}

void MEWCP_remove_node_from_solution_list(pointer_node_t n, matrix_weights_t * matrix_weights, node_list_t * node_list)
{
    unsigned int i,partition;
    weight_t weight_ni;

#if defined ASSERT
    // Voglio essere sicuro che il node che tolgo appartenga alla soluzione
    assert(node_list->node[n].belongsM == true);
#endif

    // Sgancio il node n dall'insieme M dei nodi soluzione

    if (node_list->node[n].prev == NULL_POINTER) // il node n è il primo della node_list
    {
        // Sistemo il cursore alla testa di M
        node_list->M_head = node_list->node[n].next;
    }
    else // n non è il primo elemento e quindi posso considerare i rif. dell'elem precedente
    {
        node_list->node[node_list->node[n].prev].next = node_list->node[n].next;
    }
    if (node_list->node[n].next == NULL_POINTER) // se il n è l'ultimo elemento
    {
        // la coda deve puntare all'elemento precedente
        node_list->M_tail = node_list->node[n].prev;
    }
    else // non è l'ultimo elemento della lista
    {
        node_list->node[node_list->node[n].next].prev = node_list->node[n].prev;
    }

    // Aggancio l'elemento n all'insieme N dei nodi non soluzione in coda

    // controllo se la lista N è vuota
    if (node_list->N_tail == NULL_POINTER)
    {
#if defined ASSERT
        // se la coda punta a null, allora anche la testa punta a null
        assert(node_list->N_head == NULL_POINTER);
#endif

        node_list->N_head = n; // aggancio l'elemento alla testa
        node_list->N_tail = n; // aggancio l'elemento in coda
        node_list->node[n].prev = NULL_POINTER;
        node_list->node[n].next = NULL_POINTER;

    }
    else // esiste almeno un elemento
    {
        node_list->node[n].prev = node_list->N_tail;
        node_list->node[node_list->N_tail].next = n;
        node_list->N_tail = n;
        node_list->node[n].next = NULL_POINTER;
    }

    // Sistemo le cardinalità degli insiemi:  incremento N decremento M
    node_list->card_N += 1;
    node_list->card_M -= 1;

    // Ora posso dire che l'elemento n non appartiene più alla soluzione M
    node_list->node[n].belongsM = false;


    // Aggiorno la somma delle distanze verso N e M dei nodi
    weight_ni = 0;
    for ( i=0 ; i<matrix_weights->n ; ++ i)
    {

        // La somma delle distanze OUT crescono di d(n,i)
        // La somma delle distanze IN  decrescono di d(n,i)
        /*
        if (n == i)
        {
            continue;
        }
        */
        weight_ni = matrix_weights->weight[n][i];
        node_list->node[i].sum_dout += weight_ni;
        node_list->node[i].sum_din -= weight_ni;

    }

    /* node n is selected for its partition */
    partition = MEWCP_get_node_partition(n,matrix_weights);
#if defined ASSERT

    assert(node_list->selected_node_partition[partition] != NULL_POINTER);
#endif

    node_list->selected_node_partition[partition] = NULL_POINTER;

    // aggiorno il nuovo valore della funzione_obiettivo e la scrivo in node_list; IN REALTA' CAMBIA DI DIN[n]
    (void) MEWCP_update_Z_node_list(node_list,matrix_weights);
}




weight_t MEWCP_update_Z_node_list(node_list_t * node_list, matrix_weights_t * matrix_weights)
{
    int i;

    weight_t objective_function = 0;
    weight_t overall_weights = 0;
    weight_t vertex_weights = 0;

    for (i = node_list->M_head; i != NULL_POINTER ; i = node_list->node[i].next)
    {
        /* I consider 2*(edge weights) + vertex weights ...*/
        overall_weights += node_list->node[i].sum_din;
        vertex_weights += matrix_weights->weight[i][i];
    }

    // I have to get edge weights + vertex weights
#if defined ASSERT
    assert(  (overall_weights - vertex_weights)  % 2 == 0);
#endif

    // I update the value of Z related to the current nodes position
    objective_function = overall_weights - vertex_weights;
    objective_function = objective_function /2 ;
    objective_function = objective_function + vertex_weights ;

    node_list->Z = objective_function;

    return objective_function;
}

pointer_node_t MEWCP_get_node_partition(pointer_node_t i, matrix_weights_t * matrix_weights)
{
    unsigned int n,m,c;

    n = matrix_weights->n;
    m = matrix_weights->m;
    c = n/m;

    // Returns the integer division from i and c
    return (pointer_node_t) i/c;
}

void MEWCP_compute_starting_solution(matrix_weights_t * matrix_weights, node_list_t * node_list)
{
#if defined MEWCP_TABU_VERBOSE1
    printf(".: Compute starting solution :.\n");
#endif

    unsigned int i;
    unsigned int m,c;

    m = matrix_weights->m;
    c = matrix_weights->c;

    for (i=0; i<m ; ++i )
    {
        MEWCP_add_node_to_solution_list(i*c,matrix_weights,node_list);
    }
}


void MEWCP_create_iterations_list(const unsigned int num_iterations_limit, list_iterations_t * list_iterations)
{
    list_iterations->current_iteration = 0;
    list_iterations->number_iterations_limit = num_iterations_limit;
    list_iterations->Z_iteration = (weight_t *) calloc(num_iterations_limit, sizeof(weight_t));

    if ( list_iterations->Z_iteration == NULL)
    {
        fprintf(stderr,"ERROR: allocation Z_iteration failed!\n");
        exit(EXIT_FAILURE);
    }

}

void MEWCP_create_tabu_node_list(matrix_weights_t * matrix_weights, tabu_node_list_t * tabu_node_list)
{

    unsigned int n;

    n= matrix_weights->n;

    tabu_node_list->tabu_node_state = (tabu_node_state_t *) calloc(n, sizeof(tabu_node_state_t));
    if ( tabu_node_list->tabu_node_state  == NULL )
    {
        (void) fprintf(stderr,"ERROR: allocation tabu_node_list FAILED\n");
        exit(EXIT_FAILURE);
    }



}

void MEWCP_initialize_tabu_structures(	matrix_weights_t * matrix_weights,
                                       const unsigned int num_iterations_limit,
                                       list_iterations_t * list_iterations,
                                       solution_t * best_solution,
                                       tabu_node_list_t * tabu_node_list)
{
    MEWCP_create_iterations_list(num_iterations_limit,list_iterations);
    MEWCP_create_tabu_node_list(matrix_weights,tabu_node_list);
    MEWCP_create_solution(matrix_weights,best_solution);
}

void EWCP_free_iterations_list(list_iterations_t * list_iterations )
{
    list_iterations->current_iteration = 0;
    list_iterations->number_iterations_limit = 0;
    free(list_iterations->Z_iteration);
}

void EWCP_free_tabu_node_list(tabu_node_list_t * tabu_node_list)
{
    free(tabu_node_list->tabu_node_state);
}

void MEWCP_free_solution(solution_t * solution)
{
    solution->Z = (weight_t) 0;
    free(solution->node_solution);
    solution->node_solution = NULL;
}

void MEWCP_free_tabu_node_list(tabu_node_list_t * tabu_node_list)
{
    free(tabu_node_list->tabu_node_state);
    tabu_node_list->tabu_node_state = NULL;
}

/********************************************************************************
 * PRINT FUNCTIONS
 * *****************************************************************************/


void MEWCP_print_matrix_weights(matrix_weights_t *  matrix_weights)
{
    unsigned int i,j;

    (void) printf("n: %u\nm: %u\nc: %u\n",matrix_weights->n,matrix_weights->m,matrix_weights->c);

    for (i = 0; i < matrix_weights->n ; ++i)
    {
        printf("\n");
        for (j = 0; j < matrix_weights->n ; ++j)
        {
            printf("%.2lf ", matrix_weights->weight[i][j]);
        }
    }
    printf("\n");

}

void MEWCP_print_node_list(node_list_t * node_list)
{
    int i;

    printf("\nNodes in list N\n");
    printf("|N|: %u\n",node_list->card_N);

    printf("Number  Is_M  Next   Sum_Din  Sum_Dout\n");
    for(i = node_list->N_head ; i != NULL_POINTER ; i = node_list->node[i].next)
    {
        printf("%u\t%d\t%d\t%.2lf\t%.2lf\n",i,node_list->node[i].belongsM,
               node_list->node[i].next,
               node_list->node[i].sum_din,
               node_list->node[i].sum_dout);
    }

    printf("\nNodes in list M\n");
    printf("|M|: %u\n",node_list->card_M);

    printf("Number  Is_M  Next   Sum_Din  Sum_Dout\n");
    for(i = node_list->M_head ; i != NULL_POINTER ; i = node_list->node[i].next)
    {
        printf("%u\t%d\t%d\t%.2lf\t%.2lf\n",i,node_list->node[i].belongsM,
               node_list->node[i].next,
               node_list->node[i].sum_din,
               node_list->node[i].sum_dout);
    }

    printf("\nOBJ Function:  %.2lf\n\n",node_list->Z);
}

void MEWCP_print_solution(matrix_weights_t * matrix_weights, solution_t * solution)
{
    unsigned int i,m;


    m= matrix_weights->m;

    printf("X: {");
    for (i=0;i<m;++i)
    {
        printf(" %d",solution->node_solution[i]);
    }
    printf(" }\n");

    printf("Z*: %.2lf\n",solution->Z);

}

