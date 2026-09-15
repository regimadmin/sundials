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
 * This is the implementation file for the ARKODE rooted-tree
 * module expressing the elementary differentials of an autonomous
 * ODE y' = f(y) (see arkode_butcher_trees_impl.h for details).
 *--------------------------------------------------------------*/

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arkode_butcher_trees_impl.h"
#include "arkode_impl.h"
#include "sundials_utils.h"

/* Multiply two positive long ints, returning 0 on overflow */
static long int mul_guard(long int x, long int y)
{
  if (x <= 0 || y <= 0 || x > LONG_MAX / y) { return 0; }
  return x * y;
}

/* Add two nonnegative long ints, returning -1 on overflow */
static long int add_guard(long int x, long int y)
{
  if (x < 0 || y < 0 || x > LONG_MAX - y) { return -1; }
  return x + y;
}

/* Checks that levels is a well-formed level sequence of a rooted tree in
 * preorder: every non-root node is deeper than the root and at most one
 * level deeper than its predecessor */
static sunbooleantype tree_valid(const int* levels, int order)
{
  if (levels == NULL || order < 1) { return SUNFALSE; }
  for (int i = 1; i < order; i++)
  {
    if (levels[i] <= levels[0] || levels[i] > levels[i - 1] + 1)
    {
      return SUNFALSE;
    }
  }
  return SUNTRUE;
}

/* Returns the index one past the end of the subtree rooted at index i */
static int subtree_end(const int* levels, int order, int i)
{
  int end = i + 1;
  while (end < order && levels[end] > levels[i]) { end++; }
  return end;
}

/*---------------------------------------------------------------
  Tree enumeration (Beyer--Hedetniemi successor rule)
  ---------------------------------------------------------------*/

int arkButcherTreeIterInit(int order, ARKodeButcherTreeIter* iter)
{
  if (iter == NULL) { return ARK_ILL_INPUT; }
  iter->order     = order;
  iter->levels    = NULL;
  iter->have_tree = SUNFALSE;
  if (order < 1) { return ARK_ILL_INPUT; }
  iter->levels = (int*)malloc(order * sizeof(*iter->levels));
  if (iter->levels == NULL) { return ARK_MEM_FAIL; }
  return ARK_SUCCESS;
}

sunbooleantype arkButcherTreeIterNext(ARKodeButcherTreeIter* iter)
{
  if (iter == NULL || iter->levels == NULL || iter->order < 1)
  {
    return SUNFALSE;
  }

  /* The first tree is the "path" tree {1,2,...,n} */
  if (!iter->have_tree)
  {
    for (int i = 0; i < iter->order; i++) { iter->levels[i] = i + 1; }
    iter->have_tree = SUNTRUE;
    return SUNTRUE;
  }

  /* Find the last node deeper than level 2; if none remain then the
   * "bush" tree {1,2,...,2} was the final tree of this order */
  int p = -1;
  for (int i = iter->order - 1; i > 0; i--)
  {
    if (iter->levels[i] > 2)
    {
      p = i;
      break;
    }
  }
  if (p < 0)
  {
    iter->have_tree = SUNFALSE;
    return SUNFALSE;
  }

  /* Find the parent q of node p, then replace the remainder of the
   * sequence with copies of the subsequence starting at q */
  int q = p - 1;
  while (iter->levels[q] != iter->levels[p] - 1) { q--; }
  for (int i = p; i < iter->order; i++)
  {
    iter->levels[i] = iter->levels[i - (p - q)];
  }
  return SUNTRUE;
}

void arkButcherTreeIterFree(ARKodeButcherTreeIter* iter)
{
  if (iter == NULL) { return; }
  free(iter->levels);
  iter->levels    = NULL;
  iter->have_tree = SUNFALSE;
}

/*---------------------------------------------------------------
  Tree counting (A000081 divisor-sum recurrence)
  ---------------------------------------------------------------*/

long int arkButcherTreeCount(int order)
{
  if (order < 1) { return 0; }

  long int* a = (long int*)calloc(order + 1, sizeof(*a));
  if (a == NULL) { return 0; }

  a[1] = 1;
  for (int n = 1; n < order; n++)
  {
    /* a(n+1) = (1/n) * sum_{k=1}^{n} ( sum_{d|k} d*a(d) ) * a(n-k+1) */
    long int total = 0;
    for (int k = 1; k <= n; k++)
    {
      long int inner = 0;
      for (int d = 1; d <= k; d++)
      {
        if (k % d == 0) { inner = add_guard(inner, mul_guard(d, a[d])); }
      }
      long int term = (inner < 0) ? 0 : mul_guard(inner, a[n - k + 1]);
      total         = (term == 0) ? -1 : add_guard(total, term);
      if (total < 0) { break; }
    }
    if (total < 0)
    {
      free(a);
      return 0;
    }
    a[n + 1] = total / n;
  }

  long int count = a[order];
  free(a);
  return count;
}

/*---------------------------------------------------------------
  Per-tree functionals: order, symmetry, and density
  ---------------------------------------------------------------*/

int arkButcherTreeOrder(const int* levels, int order)
{
  return tree_valid(levels, order) ? order : 0;
}

/* Recursion for sigma(t): the product over the children of each node of
 * sigma(child), times m! for every group of m identical children.  In a
 * canonical level sequence identical subtrees are adjacent. */
static long int tree_sigma(const int* levels, int order)
{
  long int sigma = 1;
  long int mult  = 1;
  int prev       = -1;
  int child      = 1;
  while (child < order && sigma > 0)
  {
    int end = subtree_end(levels, order, child);
    sigma   = mul_guard(sigma, tree_sigma(&levels[child], end - child));
    if (prev >= 0 && end - child == child - prev &&
        memcmp(&levels[child], &levels[prev], (end - child) * sizeof(*levels)) ==
          0)
    {
      mult++;
      sigma = mul_guard(sigma, mult);
    }
    else { mult = 1; }
    prev  = child;
    child = end;
  }
  return sigma;
}

long int arkButcherTreeSigma(const int* levels, int order)
{
  if (!tree_valid(levels, order)) { return 0; }
  return tree_sigma(levels, order);
}

/* Recursion for gamma(t): the number of nodes in the tree times the
 * product of gamma over the subtrees of the root's children */
static long int tree_gamma(const int* levels, int order)
{
  long int gamma = order;
  int child      = 1;
  while (child < order && gamma > 0)
  {
    int end = subtree_end(levels, order, child);
    gamma   = mul_guard(gamma, tree_gamma(&levels[child], end - child));
    child   = end;
  }
  return gamma;
}

long int arkButcherTreeGamma(const int* levels, int order)
{
  if (!tree_valid(levels, order)) { return 0; }
  return tree_gamma(levels, order);
}

/*---------------------------------------------------------------
  Expression rendering
  ---------------------------------------------------------------*/

/* A minimal dynamically-sized string builder */
typedef struct
{
  char* str;
  size_t len;
  size_t cap;
} tree_strbuf;

static sunbooleantype strbuf_append(tree_strbuf* buf, const char* s)
{
  size_t slen = strlen(s);
  if (buf->len + slen + 1 > buf->cap)
  {
    size_t new_cap = 2 * (buf->len + slen + 1);
    char* new_str  = (char*)realloc(buf->str, new_cap);
    if (new_str == NULL) { return SUNFALSE; }
    buf->str = new_str;
    buf->cap = new_cap;
  }
  memcpy(&buf->str[buf->len], s, slen + 1);
  buf->len += slen;
  return SUNTRUE;
}

static sunbooleantype strbuf_append_long(tree_strbuf* buf, long int value)
{
  char digits[32];
  snprintf(digits, sizeof(digits), "%ld", value);
  return strbuf_append(buf, digits);
}

/* Renders the elementary differential of a subtree: a leaf is "f", a node
 * with a single child is "f'" followed by the child, and a node with m > 1
 * children is the m-th derivative of f applied to the children */
static sunbooleantype diff_render(const int* levels, int order, tree_strbuf* buf)
{
  int children = 0;
  for (int child = 1; child < order; child = subtree_end(levels, order, child))
  {
    children++;
  }

  if (children == 0) { return strbuf_append(buf, "f"); }

  if (children == 1)
  {
    return strbuf_append(buf, "f'") && diff_render(&levels[1], order - 1, buf);
  }

  if (children == 2 || children == 3)
  {
    if (!strbuf_append(buf, (children == 2) ? "f''(" : "f'''("))
    {
      return SUNFALSE;
    }
  }
  else
  {
    if (!strbuf_append(buf, "f^(") || !strbuf_append_long(buf, children) ||
        !strbuf_append(buf, ")("))
    {
      return SUNFALSE;
    }
  }

  int child = 1;
  for (int i = 0; i < children; i++)
  {
    int end = subtree_end(levels, order, child);
    if (i > 0 && !strbuf_append(buf, ",")) { return SUNFALSE; }
    if (!diff_render(&levels[child], end - child, buf)) { return SUNFALSE; }
    child = end;
  }
  return strbuf_append(buf, ")");
}

char* arkButcherTreeDiffString(const int* levels, int order)
{
  if (!tree_valid(levels, order)) { return NULL; }
  tree_strbuf buf = {NULL, 0, 0};
  if (!diff_render(levels, order, &buf))
  {
    free(buf.str);
    return NULL;
  }
  return buf.str;
}

static sunbooleantype weight_vec_render(const int* levels, int order,
                                        tree_strbuf* buf);

/* Renders the elementwise product of the stage vectors of the children of a
 * node.  Groups of m identical children render as a single factor raised to
 * the elementwise power m, compound (non-leaf) factors are parenthesized,
 * multi-factor products are joined with ".*" and parenthesized, and an empty
 * product renders as the ones vector "e". */
static sunbooleantype product_render(const int* levels, int order,
                                     tree_strbuf* buf)
{
  int children = 0;
  int groups   = 0;
  int prev     = -1;
  for (int child = 1; child < order; child = subtree_end(levels, order, child))
  {
    int end = subtree_end(levels, order, child);
    if (prev < 0 || end - child != child - prev ||
        memcmp(&levels[child], &levels[prev], (end - child) * sizeof(*levels)) !=
          0)
    {
      groups++;
    }
    prev = child;
    children++;
  }

  if (children == 0) { return strbuf_append(buf, "e"); }

  if (groups > 1 && !strbuf_append(buf, "(")) { return SUNFALSE; }

  int child = 1;
  int index = 0;
  while (child < order)
  {
    int end = subtree_end(levels, order, child);

    /* count the multiplicity of this child subtree */
    int mult = 1;
    int next = end;
    while (
      next < order && subtree_end(levels, order, next) - next == end - child &&
      memcmp(&levels[next], &levels[child], (end - child) * sizeof(*levels)) == 0)
    {
      mult++;
      next = subtree_end(levels, order, next);
    }

    if (index > 0 && !strbuf_append(buf, ".*")) { return SUNFALSE; }

    sunbooleantype leaf = (end - child == 1);
    if (!leaf && !strbuf_append(buf, "(")) { return SUNFALSE; }
    if (!weight_vec_render(&levels[child], end - child, buf))
    {
      return SUNFALSE;
    }
    if (!leaf && !strbuf_append(buf, ")")) { return SUNFALSE; }
    if (mult > 1 && (!strbuf_append(buf, ".^") || !strbuf_append_long(buf, mult)))
    {
      return SUNFALSE;
    }

    child = next;
    index++;
  }

  if (groups > 1 && !strbuf_append(buf, ")")) { return SUNFALSE; }
  return SUNTRUE;
}

/* Renders the stage vector of a non-root subtree: a leaf is the abscissa
 * vector "c" while an internal node applies A to the product of its
 * children */
static sunbooleantype weight_vec_render(const int* levels, int order,
                                        tree_strbuf* buf)
{
  if (order == 1) { return strbuf_append(buf, "c"); }
  return strbuf_append(buf, "A*") && product_render(levels, order, buf);
}

char* arkButcherTreeWeightString(const int* levels, int order)
{
  if (!tree_valid(levels, order)) { return NULL; }
  tree_strbuf buf = {NULL, 0, 0};
  if (!strbuf_append(&buf, "b'*") || !product_render(levels, order, &buf))
  {
    free(buf.str);
    return NULL;
  }
  return buf.str;
}

/*---------------------------------------------------------------
  Elementary weight evaluation
  ---------------------------------------------------------------*/

static sunrealtype tree_dot(const sunrealtype* x, const sunrealtype* y, int n)
{
  sunrealtype total = ZERO;
  sunrealtype err   = ZERO;
  for (int i = 0; i < n; i++)
  {
    sunCompensatedSum(total, x[i] * y[i], &total, &err);
  }
  return total;
}

/* Computes the stage vector of a non-root subtree into val (length s).
 * Each recursion depth uses two length-s slices of work: one accumulating
 * the elementwise product of the children and one holding a child value. */
static void tree_stage_vec(const int* levels, int order, ARKodeButcherTable B,
                           sunrealtype* work, sunrealtype* val)
{
  int s = B->stages;

  if (order == 1)
  {
    memcpy(val, B->c, s * sizeof(*val));
    return;
  }

  sunrealtype* prod      = work;
  sunrealtype* child_val = &work[s];
  for (int i = 0; i < s; i++) { prod[i] = ONE; }
  for (int child = 1; child < order;)
  {
    int end = subtree_end(levels, order, child);
    tree_stage_vec(&levels[child], end - child, B, &work[2 * s], child_val);
    for (int i = 0; i < s; i++) { prod[i] *= child_val[i]; }
    child = end;
  }

  for (int i = 0; i < s; i++) { val[i] = tree_dot(B->A[i], prod, s); }
}

int arkButcherTreePhi(const int* levels, int order, ARKodeButcherTable B,
                      sunrealtype* phi, sunrealtype* phi_hat)
{
  if (!tree_valid(levels, order) || B == NULL || B->stages < 1 ||
      B->A == NULL || B->b == NULL || B->c == NULL || phi == NULL)
  {
    return ARK_ILL_INPUT;
  }

  int s             = B->stages;
  sunrealtype* work = (sunrealtype*)malloc(2 * order * s * sizeof(*work));
  if (work == NULL) { return ARK_MEM_FAIL; }

  /* accumulate the product over the root's children (the empty product is
   * the ones vector, so a single node yields phi = sum(b)) */
  sunrealtype* prod      = work;
  sunrealtype* child_val = &work[s];
  for (int i = 0; i < s; i++) { prod[i] = ONE; }
  for (int child = 1; child < order;)
  {
    int end = subtree_end(levels, order, child);
    tree_stage_vec(&levels[child], end - child, B, &work[2 * s], child_val);
    for (int i = 0; i < s; i++) { prod[i] *= child_val[i]; }
    child = end;
  }

  *phi = tree_dot(B->b, prod, s);
  if (phi_hat != NULL)
  {
    *phi_hat = (B->d != NULL) ? tree_dot(B->d, prod, s) : ZERO;
  }

  free(work);
  return ARK_SUCCESS;
}
