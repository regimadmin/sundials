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
 * Private header for the ARKODE rooted-tree utilities used to
 * generate Runge--Kutta order conditions.
 *
 * The elementary differentials of an autonomous ODE y' = f(y) are
 * in one-to-one correspondence with unlabeled rooted trees (Butcher
 * theory), whose count per order is OEIS sequence A000081 (1, 1, 2,
 * 4, 9, 20, 48, 115, 286, 719, ...).  Each rooted tree t of order
 * r(t) = n yields one Runge--Kutta order condition
 *    Phi(t) = 1/gamma(t),
 * where Phi(t) is the elementary weight of the Butcher table and
 * gamma(t) is the density of the tree.
 *
 * Each tree is stored canonically as a "level sequence": the depth
 * of each node in a preorder (depth-first) traversal, with the root
 * at level 1.  For example, [1,2,3,2] is the order-4 tree whose
 * root has two children, the first of which has one child.
 *--------------------------------------------------------------*/

#ifndef _ARKODE_BUTCHER_TREES_IMPL_H
#define _ARKODE_BUTCHER_TREES_IMPL_H

#include <stddef.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*---------------------------------------------------------------
  Iterator over all rooted trees of a given order (number of
  nodes).  Trees are produced as canonical level sequences using
  the Beyer--Hedetniemi successor algorithm.
  ---------------------------------------------------------------*/
typedef struct
{
  int order;   /* number of nodes in each tree            */
  int* levels; /* current level sequence, levels[0] = 1   */
} ARKodeButcherTreeIter;

/* Initialize an iterator over all rooted trees with the given number
   of nodes (order >= 1); returns 0 on success, -1 on failure.  Upon
   successful return the iterator holds the first tree (the "path"
   tree [1,2,...,order]). */
int arkodeButcherTrees_IterInit(ARKodeButcherTreeIter* iter, int order);

/* Advance the iterator to the next canonical tree; returns SUNTRUE if
   a new tree is available, SUNFALSE if the enumeration is complete. */
sunbooleantype arkodeButcherTrees_IterNext(ARKodeButcherTreeIter* iter);

/* Release iterator resources */
void arkodeButcherTrees_IterFree(ARKodeButcherTreeIter* iter);

/*---------------------------------------------------------------
  Tree counting and per-tree functionals
  ---------------------------------------------------------------*/

/* Number of rooted trees with the given number of nodes (OEIS
   A000081), computed from the standard convolution recurrence;
   returns -1 for non-positive order */
long int arkodeButcherTrees_Count(int order);

/* Density gamma(t) of the tree: gamma(leaf) = 1 and gamma(t) =
   r(t) * prod(gamma(subtrees of the root's children)) */
long int arkodeButcherTrees_Density(const int* levels, int order);

/* Symmetry sigma(t) of the tree: the order of the automorphism group,
   sigma(t) = prod over classes of equal child subtrees of
   (multiplicity! * sigma(subtree)^multiplicity) */
long int arkodeButcherTrees_Symmetry(const int* levels, int order);

/*---------------------------------------------------------------
  Elementary weights
  ---------------------------------------------------------------*/

/* Evaluate the elementary weight Phi(t) of the tree for the Butcher
   table data (A, b, c) with s stages.  The arrays b, A and c are
   indexed by "color" so that the same routine serves both single
   tables and additive (ARK) pairs:

   * colors == NULL: all nodes use color 0 (single-table weight)
   * colors != NULL: colors[i] in {0,1,...} gives the color of node i
     (in level-sequence order); the root color selects the b (or d)
     vector, every other node v selects the A matrix applied at v,
     and leaves select the corresponding c vector.

   The work array must hold at least 2*(order+1)*s entries.  Returns
   0 on success and -1 on failure (NULL pointers). */
int arkodeButcherTrees_Weight(const int* levels, const int* colors, int order,
                              sunrealtype* const* b, sunrealtype** const* A,
                              sunrealtype* const* c, int s, sunrealtype* work,
                              sunrealtype* phi);

/*---------------------------------------------------------------
  Expression rendering (for diagnostic output and documentation)
  ---------------------------------------------------------------*/

/* Write the elementary differential of the tree, e.g. "f", "f'f",
   "f''(f,f)", "f'f'f", into str (at most len characters, always
   NUL-terminated); returns 0 on success, -1 on failure */
int arkodeButcherTrees_ElementaryDiff(const int* levels, int order, char* str,
                                      size_t len);

/* Write the elementary weight of the tree in vector notation, e.g.
   "b'*e", "b'*c", "b'*(c.*c)", "b'*A*c", "b'*((A*(c.*c)).*(A*c))",
   into str (at most len characters, always NUL-terminated); returns
   0 on success, -1 on failure */
int arkodeButcherTrees_WeightString(const int* levels, int order, char* str,
                                    size_t len);

#ifdef __cplusplus
}
#endif

#endif
