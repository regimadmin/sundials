/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2026, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * Copyright (c) 2013-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * Copyright (c) 2002-2013, Lawrence Livermore National Security.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * This is the (private) header file for the ARKODE rooted-tree
 * module used to express the elementary differentials of an
 * autonomous ODE y' = f(y).  The elementary differentials of
 * order n (f, f'f, f''(f,f), f'f'f, ...) are in one-to-one
 * correspondence with the unlabeled rooted trees having n nodes,
 * whose counts per order form the OEIS sequence A000081
 * (1, 1, 2, 4, 9, 20, 48, 115, 286, 719, ...).  Each tree t
 * gives one Runge--Kutta order condition Phi(t) = 1/gamma(t),
 * where Phi is the elementary weight of a Butcher table and
 * gamma is the tree density.
 *
 * Trees are represented canonically by their *level sequence*:
 * the array of node depths (root depth 1) in a preorder walk of
 * the tree, with the subtrees of every node listed in
 * non-increasing lexicographic order.
 *--------------------------------------------------------------*/

#ifndef _ARKODE_BUTCHER_TREES_IMPL_H
#define _ARKODE_BUTCHER_TREES_IMPL_H

#include <arkode/arkode_butcher.h>
#include <stdio.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*---------------------------------------------------------------
  Type : ARKodeButcherTreeIter
  ---------------------------------------------------------------
  Iterator over all unlabeled rooted trees with a fixed number of
  nodes.  Trees are visited in decreasing lexicographic order of
  their canonical level sequences, beginning with the "path" tree
  {1,2,...,n} and ending with the "bush" tree {1,2,...,2}, using
  the constant-amortized-time successor rule of Beyer and
  Hedetniemi, "Constant time generation of rooted trees", SIAM J.
  Comput. 9 (1980).  The number of trees visited for each order n
  equals A000081(n).
  ---------------------------------------------------------------*/
typedef struct
{
  int order;                /* number of nodes in each tree                */
  int* levels;              /* level sequence of the current tree          */
  sunbooleantype have_tree; /* SUNTRUE if levels holds a tree  */
} ARKodeButcherTreeIter;

/* Initialize an iterator over all rooted trees with the given number
   of nodes (order >= 1).  Returns ARK_SUCCESS, ARK_ILL_INPUT for an
   invalid order, or ARK_MEM_FAIL on an allocation failure. */
int arkButcherTreeIterInit(int order, ARKodeButcherTreeIter* iter);

/* Advance the iterator to the next tree, placing its level sequence
   in iter->levels.  Returns SUNTRUE if a tree was produced, and
   SUNFALSE once all A000081(order) trees have been visited. */
sunbooleantype arkButcherTreeIterNext(ARKodeButcherTreeIter* iter);

/* Free all memory allocated within the iterator */
void arkButcherTreeIterFree(ARKodeButcherTreeIter* iter);

/* Return the number of rooted trees with the given number of nodes,
   A000081(order), computed with the standard divisor-sum recurrence
     a(n+1) = (1/n) * sum_{k=1}^{n} ( sum_{d|k} d*a(d) ) * a(n-k+1),
   or 0 if order < 1 or the count would overflow a long int. */
long int arkButcherTreeCount(int order);

/* Return the order r(t) (number of nodes) of the tree with the given
   level sequence */
int arkButcherTreeOrder(const int* levels, int order);

/* Return the symmetry sigma(t) of the tree, i.e., the size of its
   automorphism group, computed recursively over the multiset of
   subtrees of each node */
long int arkButcherTreeSigma(const int* levels, int order);

/* Return the density gamma(t) of the tree, i.e., the product over
   all nodes of the number of nodes in the subtree rooted there.
   Each tree yields the order condition Phi(t) = 1/gamma(t). */
long int arkButcherTreeGamma(const int* levels, int order);

/* Return the elementary differential F(t) of the tree rendered as a
   string, e.g., "f", "f'f", "f''(f,f)", "f'f'f", or NULL on an
   allocation failure.  The caller is responsible for freeing the
   returned string. */
char* arkButcherTreeDiffString(const int* levels, int order);

/* Return the elementary weight Phi(t) of the tree rendered as a
   string in MATLAB-like notation, e.g., "b'*e", "b'*c",
   "b'*((A*c.^2).*(A*c))", or NULL on an allocation failure ("e"
   denotes the vector of all ones).  The caller is responsible for
   freeing the returned string. */
char* arkButcherTreeWeightString(const int* levels, int order);

/* Evaluate the elementary weight Phi(t) of a Butcher table for the
   tree with the given level sequence, i.e., the value of b'*(...)
   obtained by mapping each leaf to c and each internal node to
   A*(elementwise product of its children).  If phi_hat is non-NULL
   and B has an embedding, also evaluate the embedded weight d'*(...).
   Returns ARK_SUCCESS, ARK_ILL_INPUT for invalid inputs, or
   ARK_MEM_FAIL on an allocation failure. */
int arkButcherTreePhi(const int* levels, int order, ARKodeButcherTable B,
                      sunrealtype* phi, sunrealtype* phi_hat);

#ifdef __cplusplus
}
#endif

#endif
