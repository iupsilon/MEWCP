#if !defined CONVERTER_DSDP_H
#define CONVERTER_DSDP_H

/*****************************************************************************
 *  Project: Maximum Edge Weighted Clique Problem
 * 
 *  Authors:
 *  (c) 2009 Yari Melzani (yari.melzani@gmail.com) 
 *
 ****************************************************************************/


#include "MEWCP_dsdp.h"
#include "MEWCP_tabu_definitions.h"



void MEWCP_generate_sdp_constraints(matrix_weights_t * matrix_weights, constraint_t ** sdp_constraints_matrix, double ** bi, constraint_t ** vect_mat_braching_constraints, double alpha,
                                    bool c_cardinality,
                                    bool c_simple_MC,
                                    bool c_improved_MC_A,
                                    bool c_improved_MC_B,
                                    bool c_improved_MC_C,
                                    bool c_4C2,
                                    bool c_4C3_A,
                                    bool c_4C3_B,
                                    unsigned int * out_num_contraints
                                   );





/**********************************************************************/
/*		Functions for contraints
 * ********************************************************************/

void MEWCP_write_constraints_weights_matrix(matrix_weights_t * matrix_weights, constraint_t * sdp_constraints_matrix, const unsigned int vector_length, const unsigned int n, unsigned int k);
void MEWCP_write_constraints_cardinality(constraint_t * sdp_constraints_matrix, const unsigned int vector_length, const unsigned int n,  unsigned int k);
void MEWCP_write_constraints_simple_MC(constraint_t * sdp_constraints_matrix, const unsigned int vector_length, const unsigned int n, const unsigned int c,  unsigned int k);
void MEWCP_write_constraints_improved_MC_A(constraint_t * sdp_constraints_matrix, const unsigned int vector_length, const unsigned int n, const unsigned int c,  unsigned int k);
void MEWCP_write_constraints_improved_MC_B(constraint_t * sdp_constraints_matrix, const unsigned int vector_length, const unsigned int n, const unsigned int m, const unsigned int c,  unsigned int k);
void MEWCP_write_constraints_improved_MC_C(constraint_t * sdp_constraints_matrix, const unsigned int vector_length, const unsigned int m, const unsigned int c,  unsigned int k);
void MEWCP_write_constraints_4C2(constraint_t * sdp_constraints_matrix, const unsigned int vector_length, const unsigned int n, const unsigned int m, const unsigned int c,  unsigned int k);
void MEWCP_write_constraints_4C3_A(constraint_t * sdp_constraints_matrix, const unsigned int vector_length, const unsigned int n, const unsigned int m, const unsigned int c,  unsigned int k);
void MEWCP_write_constraints_4C3_B(constraint_t * sdp_constraints_matrix, const unsigned int vector_length, const unsigned int n, const unsigned int m, const unsigned int c,  unsigned int k);
void MEWCP_generate_constraints_branch(list_blocked_nodes_t * list_blocked_nodes, constraint_t * vect_mat_contraints, const unsigned int vector_length, const unsigned int n, const unsigned int c);

/*	FUNZIONI PER SCRIVERE I TERMINI NOTI DEI VINCOLI */
void MEWCP_write_b_cardinality(double * bi, const int ishift, const unsigned m);
void MEWCP_write_b_simple_MC(double * bi, const int ishift, const unsigned n, const unsigned m);
void MEWCP_write_b_improved_MC_A(double * bi, const int ishift, const unsigned n);
void MEWCP_write_b_improved_MC_B(double * bi, const int ishift, const unsigned n, const unsigned m);
void MEWCP_write_b_improved_MC_C(double * bi, const int ishift,const unsigned int m);
void MEWCP_write_b_4C2(double * bi, const int ishift);
void MEWCP_write_b_4C3_A(double * bi, const int ishift, const unsigned n);
void MEWCP_write_b_4C3_B(double * bi, const int ishift,const unsigned n, const unsigned m);
void MEWCP_write_b_branching(double * bi, const int ishift);


#endif /*CONVERTER_DSDP_H*/

