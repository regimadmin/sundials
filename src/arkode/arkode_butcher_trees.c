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
 * Rooted-tree utilities used to generate Runge--Kutta order
 * conditions from the elementary differentials of an ODE (see
 * arkode_butcher_trees_impl.h for an overview).
 *
 * References:
 *
 * J. C. Butcher, "Numerical Methods for Ordinary Differential
 * Equations," 2nd ed., Wiley, 2008 (chapter 3).
 *
 * E. Hairer, S. P. Norsett and G. Wanner, "Solving Ordinary
 * Differential Equations I: Nonstiff Problems," 2nd ed., Springer,
 * 1993 (section II.2).
 *
 * T. Beyer and S. M. Hedetniemi, "Constant time generation of
 * rooted trees," SIAM J. Comput., 9(4):706-712, 1980.
 *
 * OEIS Foundation Inc., "The On-Line Encyclopedia of Integer
 * Sequences," entry A000081, https://oeis.org/A000081.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_math.h>

#include "arkode_butcher_trees_impl.h"

/*---------------------------------------------------------------
  Private helper routines
  ---------------------------------------------------------------*/

/* Number of nodes in the subtree rooted at position pos: the
   contiguous run of following entries with a deeper level */
static int treeSubtreeSize(const int* levels, int order, int pos)
{
  int end = pos + 1;
  while ((end < order) && (levels[end] > levels[pos])) { end++; }
  return (end - pos);
}

/* Compare the subtrees (as level-sequence slices) rooted at two
   children of a common parent for equality */
static sunbooleantype treeSubtreesEqual(const int* levels, int pos1, int size1,
                                        int pos2, int size2)
{
  int i;
  if (size1 != size2) { return (SUNFALSE); }
  for (i = 0; i < size1; i++)
  {
    if (levels[pos1 + i] != levels[pos2 + i]) { return (SUNFALSE); }
  }
  return (SUNTRUE);
}

/* Density gamma of the subtree rooted at position pos */
static long int treeSubtreeDensity(const int* levels, int order, int pos)
{
  int size, child;
  long int gamma;
  size  = treeSubtreeSize(levels, order, pos);
  gamma = size;
  child = pos + 1;
  while (child < pos + size)
  {
    gamma *= treeSubtreeDensity(levels, order, child);
    child += treeSubtreeSize(levels, order, child);
  }
  return (gamma);
}

/* Symmetry sigma of the subtree rooted at position pos */
static long int treeSubtreeSymmetry(const int* levels, int order, int pos)
{
  int size, child, child2, csize, c2size;
  long int sigma, mult;
  size  = treeSubtreeSize(levels, order, pos);
  sigma = 1;
  child = pos + 1;
  while (child < pos + size)
  {
    csize = treeSubtreeSize(levels, order, child);
    sigma *= treeSubtreeSymmetry(levels, order, child);

    /* multiplicity factor: children equal to this one seen so far
       (including itself), so that a class of m equal subtrees
       contributes 1*2*...*m = m! overall */
    mult   = 1;
    child2 = pos + 1;
    while (child2 < child)
    {
      c2size = treeSubtreeSize(levels, order, child2);
      if (treeSubtreesEqual(levels, child, csize, child2, c2size)) { mult++; }
      child2 += c2size;
    }
    sigma *= mult;

    child += csize;
  }
  return (sigma);
}

/* Compute the elementwise product over the children ch of the node at
   position pos of the stage vectors psi(ch), where psi(leaf) = c and
   psi(internal node v) = A * (product of psi over children of v); an
   empty product yields the vector of ones.  The work array must
   provide 2*s entries per tree level below pos. */
static int treeWeightChildren(const int* levels, const int* colors, int order,
                              int pos, sunrealtype** const* A,
                              sunrealtype* const* c, int s, sunrealtype* prod,
                              sunrealtype* work)
{
  int i, j, child, csize, size, color, retval;
  sunrealtype* psi = work;
  sunrealtype* sub = work + s;

  for (i = 0; i < s; i++) { prod[i] = SUN_RCONST(1.0); }
  size  = treeSubtreeSize(levels, order, pos);
  child = pos + 1;
  while (child < pos + size)
  {
    csize = treeSubtreeSize(levels, order, child);
    color = (colors == NULL) ? 0 : colors[child];
    if (csize == 1) /* leaf */
    {
      for (i = 0; i < s; i++) { psi[i] = c[color][i]; }
    }
    else /* internal node: psi = A * (product over its children) */
    {
      retval = treeWeightChildren(levels, colors, order, child, A, c, s, sub,
                                  work + 2 * s);
      if (retval) { return (retval); }
      for (i = 0; i < s; i++)
      {
        psi[i] = SUN_RCONST(0.0);
        for (j = 0; j < s; j++) { psi[i] += A[color][i][j] * sub[j]; }
      }
    }
    for (i = 0; i < s; i++) { prod[i] *= psi[i]; }
    child += csize;
  }
  return (0);
}

/* Append text to the output cursor, returning -1 if it does not fit */
static int strAppend(char** pos, size_t* rem, const char* text)
{
  size_t n = strlen(text);
  if (n + 1 > *rem) { return (-1); }
  memcpy(*pos, text, n);
  *pos += n;
  **pos = '\0';
  *rem -= n;
  return (0);
}

/* Render the elementary differential of the subtree rooted at pos */
static int treeDiffNode(const int* levels, int order, int pos, char** out,
                        size_t* rem)
{
  int size, child, csize, nchildren, retval;
  char deriv[16];

  size      = treeSubtreeSize(levels, order, pos);
  nchildren = 0;
  child     = pos + 1;
  while (child < pos + size)
  {
    nchildren++;
    child += treeSubtreeSize(levels, order, child);
  }

  if (nchildren == 0) { return (strAppend(out, rem, "f")); }

  if (nchildren == 1)
  {
    retval = strAppend(out, rem, "f'");
    if (retval) { return (retval); }
    return (treeDiffNode(levels, order, pos + 1, out, rem));
  }

  /* nchildren >= 2 */
  if (nchildren == 2) { strcpy(deriv, "f''"); }
  else if (nchildren == 3) { strcpy(deriv, "f'''"); }
  else { snprintf(deriv, sizeof(deriv), "f^(%d)", nchildren); }
  retval = strAppend(out, rem, deriv);
  if (retval) { return (retval); }
  retval = strAppend(out, rem, "(");
  if (retval) { return (retval); }
  child = pos + 1;
  while (child < pos + size)
  {
    csize  = treeSubtreeSize(levels, order, child);
    retval = treeDiffNode(levels, order, child, out, rem);
    if (retval) { return (retval); }
    if (child + csize < pos + size)
    {
      retval = strAppend(out, rem, ",");
      if (retval) { return (retval); }
    }
    child += csize;
  }
  return (strAppend(out, rem, ")"));
}

/* Render the stage-vector expression psi for the subtree rooted at a
   non-root node: "c" for a leaf, otherwise "A*" applied to the
   (possibly parenthesized) product over children */
static int treeWeightNode(const int* levels, int order, int pos, char** out,
                          size_t* rem);

/* Render the elementwise product of psi over the children of the node
   at pos, joining factors with ".*" and parenthesizing non-leaf
   factors; nchildren must be >= 1 */
static int treeWeightProduct(const int* levels, int order, int pos, char** out,
                             size_t* rem)
{
  int size, child, csize, retval;

  size  = treeSubtreeSize(levels, order, pos);
  child = pos + 1;
  while (child < pos + size)
  {
    csize = treeSubtreeSize(levels, order, child);
    if (csize == 1) { retval = strAppend(out, rem, "c"); }
    else
    {
      retval = strAppend(out, rem, "(");
      if (retval) { return (retval); }
      retval = treeWeightNode(levels, order, child, out, rem);
      if (retval) { return (retval); }
      retval = strAppend(out, rem, ")");
    }
    if (retval) { return (retval); }
    if (child + csize < pos + size)
    {
      retval = strAppend(out, rem, ".*");
      if (retval) { return (retval); }
    }
    child += csize;
  }
  return (0);
}

static int treeWeightNode(const int* levels, int order, int pos, char** out,
                          size_t* rem)
{
  int size, child, nchildren, retval;

  size = treeSubtreeSize(levels, order, pos);
  if (size == 1) { return (strAppend(out, rem, "c")); }

  nchildren = 0;
  child     = pos + 1;
  while (child < pos + size)
  {
    nchildren++;
    child += treeSubtreeSize(levels, order, child);
  }

  retval = strAppend(out, rem, "A*");
  if (retval) { return (retval); }
  if (nchildren == 1)
  {
    return (treeWeightNode(levels, order, pos + 1, out, rem));
  }
  retval = strAppend(out, rem, "(");
  if (retval) { return (retval); }
  retval = treeWeightProduct(levels, order, pos, out, rem);
  if (retval) { return (retval); }
  return (strAppend(out, rem, ")"));
}

/*---------------------------------------------------------------
  Exported (library-private) routines
  ---------------------------------------------------------------*/

int arkodeButcherTrees_IterInit(ARKodeButcherTreeIter* iter, int order)
{
  int i;
  if ((iter == NULL) || (order < 1)) { return (-1); }
  iter->order  = order;
  iter->levels = (int*)malloc(order * sizeof(int));
  if (iter->levels == NULL) { return (-1); }
  /* the first canonical tree is the path [1,2,...,order] */
  for (i = 0; i < order; i++) { iter->levels[i] = i + 1; }
  return (0);
}

sunbooleantype arkodeButcherTrees_IterNext(ARKodeButcherTreeIter* iter)
{
  int i, p, q;
  int* levels;
  if ((iter == NULL) || (iter->levels == NULL)) { return (SUNFALSE); }
  levels = iter->levels;

  /* Beyer--Hedetniemi successor: find the last position p with level
     greater than 2 (if none, the enumeration ended at the "star"
     tree), locate the parent q of node p, then repeat the pattern
     starting at q through the end of the sequence */
  p = -1;
  for (i = iter->order - 1; i > 0; i--)
  {
    if (levels[i] > 2)
    {
      p = i;
      break;
    }
  }
  if (p < 0) { return (SUNFALSE); }

  q = -1;
  for (i = p - 1; i >= 0; i--)
  {
    if (levels[i] == levels[p] - 1)
    {
      q = i;
      break;
    }
  }
  if (q < 0) { return (SUNFALSE); } /* unreachable for valid sequences */

  for (i = p; i < iter->order; i++) { levels[i] = levels[i - (p - q)]; }
  return (SUNTRUE);
}

void arkodeButcherTrees_IterFree(ARKodeButcherTreeIter* iter)
{
  if (iter == NULL) { return; }
  free(iter->levels);
  iter->levels = NULL;
  iter->order  = 0;
}

long int arkodeButcherTrees_Count(int order)
{
  int n, k, d;
  long int inner, total, count;
  long int* a;

  if (order < 1) { return (-1); }
  if (order == 1) { return (1); }

  /* a(1) = 1 and a(n+1) = (1/n) * sum_{k=1}^{n} [ sum_{d|k} d*a(d) ] * a(n-k+1) */
  a = (long int*)calloc(order + 1, sizeof(long int));
  if (a == NULL) { return (-1); }
  a[1] = 1;
  for (n = 1; n < order; n++)
  {
    total = 0;
    for (k = 1; k <= n; k++)
    {
      inner = 0;
      for (d = 1; d <= k; d++)
      {
        if (k % d == 0) { inner += d * a[d]; }
      }
      total += inner * a[n - k + 1];
    }
    a[n + 1] = total / n;
  }
  count = a[order];
  free(a);
  return (count);
}

long int arkodeButcherTrees_Density(const int* levels, int order)
{
  if ((levels == NULL) || (order < 1)) { return (-1); }
  return (treeSubtreeDensity(levels, order, 0));
}

long int arkodeButcherTrees_Symmetry(const int* levels, int order)
{
  if ((levels == NULL) || (order < 1)) { return (-1); }
  return (treeSubtreeSymmetry(levels, order, 0));
}

int arkodeButcherTrees_Weight(const int* levels, const int* colors, int order,
                              sunrealtype* const* b, sunrealtype** const* A,
                              sunrealtype* const* c, int s, sunrealtype* work,
                              sunrealtype* phi)
{
  int i, color, retval;
  sunrealtype* prod;

  if ((levels == NULL) || (order < 1) || (b == NULL) || (A == NULL) ||
      (c == NULL) || (s < 1) || (work == NULL) || (phi == NULL))
  {
    return (-1);
  }

  prod   = work;
  retval = treeWeightChildren(levels, colors, order, 0, A, c, s, prod, work + s);
  if (retval) { return (retval); }

  color = (colors == NULL) ? 0 : colors[0];
  *phi  = SUN_RCONST(0.0);
  for (i = 0; i < s; i++) { *phi += b[color][i] * prod[i]; }
  return (0);
}

int arkodeButcherTrees_ElementaryDiff(const int* levels, int order, char* str,
                                      size_t len)
{
  char* pos = str;
  size_t rem;
  if ((levels == NULL) || (order < 1) || (str == NULL) || (len < 1))
  {
    return (-1);
  }
  rem  = len;
  *pos = '\0';
  return (treeDiffNode(levels, order, 0, &pos, &rem));
}

int arkodeButcherTrees_WeightString(const int* levels, int order, char* str,
                                    size_t len)
{
  int nchildren, child, retval;
  char* pos = str;
  size_t rem;

  if ((levels == NULL) || (order < 1) || (str == NULL) || (len < 1))
  {
    return (-1);
  }
  rem  = len;
  *pos = '\0';

  nchildren = 0;
  child     = 1;
  while (child < order)
  {
    nchildren++;
    child += treeSubtreeSize(levels, order, child);
  }

  if (nchildren == 0) { return (strAppend(&pos, &rem, "b'*e")); }
  retval = strAppend(&pos, &rem, "b'*");
  if (retval) { return (retval); }
  if (nchildren == 1) { return (treeWeightNode(levels, order, 1, &pos, &rem)); }
  retval = strAppend(&pos, &rem, "(");
  if (retval) { return (retval); }
  retval = treeWeightProduct(levels, order, 0, &pos, &rem);
  if (retval) { return (retval); }
  return (strAppend(&pos, &rem, ")"));
}

/*---------------------------------------------------------------
  EOF
  ---------------------------------------------------------------*/
