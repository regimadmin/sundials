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
 * This is the implementation file for Butcher table structure
 * for the ARKODE infrastructure.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>

#include "arkode_butcher_trees_impl.h"
#include "arkode_impl.h"

/* tolerance for checking order conditions */
#define TOL (SUNRsqrt(SUN_UNIT_ROUNDOFF))

/* Maximum order for which analytical order conditions (one per
   rooted tree, cf. arkode_butcher_trees_impl.h) are checked; above
   this we revert to the Butcher simplifying assumptions.  The
   condition right-hand sides 1/gamma(t) decay like 1/order!, so the
   achievable order depends on the working precision. */
#if defined(SUNDIALS_SINGLE_PRECISION)
#define ARK_BUTCHER_MAX_CHECK_ORDER 6
#else
#define ARK_BUTCHER_MAX_CHECK_ORDER 9
#endif

/* Maximum order for which the analytical order conditions of an ARK
   pair (one per 2-colored rooted tree) are checked */
#define ARK_BUTCHER_MAX_ARK_CHECK_ORDER 6

/* Private utility functions for checking method order */
static int arkode_butcher_vp(sunrealtype* x, int l, int s, sunrealtype* z);
static int arkode_butcher_dot(sunrealtype* x, sunrealtype* y, int s,
                              sunrealtype* d);
static sunbooleantype arkode_butcher_rowsum(sunrealtype** A, sunrealtype* c,
                                            int s);
static int arkode_butcher_check_conditions(sunrealtype* b, sunrealtype** A,
                                           sunrealtype* c, int s, int order,
                                           sunrealtype* work, const char* name,
                                           FILE* outfile);
static int arkode_butcher_check_ark_conditions(sunrealtype* const* b,
                                               sunrealtype** const* A,
                                               sunrealtype* const* c, int s,
                                               int order, sunrealtype* work,
                                               const char* name, FILE* outfile);
static int __ButcherSimplifyingAssumptions(sunrealtype** A, sunrealtype* b,
                                           sunrealtype* c, int s);

/*---------------------------------------------------------------
  Routine to allocate an empty Butcher table structure
  ---------------------------------------------------------------*/
ARKodeButcherTable ARKodeButcherTable_Alloc(int stages, sunbooleantype embedded)
{
  int i;
  ARKodeButcherTable B;

  /* Check for legal 'stages' value */
  if (stages < 1) { return (NULL); }

  /* Allocate Butcher table structure */
  B = NULL;
  B = (ARKodeButcherTable)malloc(sizeof(struct ARKodeButcherTableMem));
  if (B == NULL) { return (NULL); }

  /* initialize pointers in B structure to NULL */
  B->A = NULL;
  B->b = NULL;
  B->c = NULL;
  B->d = NULL;

  /* set stages into B structure */
  B->stages = stages;

  /*
   * Allocate fields within Butcher table structure
   */

  /* allocate rows of A */
  B->A = (sunrealtype**)calloc(stages, sizeof(sunrealtype*));
  if (B->A == NULL)
  {
    ARKodeButcherTable_Free(B);
    return (NULL);
  }

  /* initialize each row of A to NULL */
  for (i = 0; i < stages; i++) { B->A[i] = NULL; }

  /* allocate columns of A */
  for (i = 0; i < stages; i++)
  {
    B->A[i] = (sunrealtype*)calloc(stages, sizeof(sunrealtype));
    if (B->A[i] == NULL)
    {
      ARKodeButcherTable_Free(B);
      return (NULL);
    }
  }

  B->b = (sunrealtype*)calloc(stages, sizeof(sunrealtype));
  if (B->b == NULL)
  {
    ARKodeButcherTable_Free(B);
    return (NULL);
  }

  B->c = (sunrealtype*)calloc(stages, sizeof(sunrealtype));
  if (B->c == NULL)
  {
    ARKodeButcherTable_Free(B);
    return (NULL);
  }

  if (embedded)
  {
    B->d = (sunrealtype*)calloc(stages, sizeof(sunrealtype));
    if (B->d == NULL)
    {
      ARKodeButcherTable_Free(B);
      return (NULL);
    }
  }

  /* initialize order parameters */
  B->q = 0;
  B->p = 0;

  return (B);
}

/*---------------------------------------------------------------
  Routine to allocate and fill a Butcher table structure
  ---------------------------------------------------------------*/
ARKodeButcherTable ARKodeButcherTable_Create(int s, int q, int p,
                                             sunrealtype* c, sunrealtype* A,
                                             sunrealtype* b, sunrealtype* d)
{
  int i, j;
  ARKodeButcherTable B;
  sunbooleantype embedded;

  /* Check for legal number of stages */
  if (s < 1) { return (NULL); }

  /* Does the table have an embedding? */
  embedded = (d != NULL) ? SUNTRUE : SUNFALSE;

  /* Allocate Butcher table structure */
  B = ARKodeButcherTable_Alloc(s, embedded);
  if (B == NULL) { return (NULL); }

  /* set the relevant parameters */
  B->stages = s;
  B->q      = q;
  B->p      = p;

  for (i = 0; i < s; i++)
  {
    B->c[i] = c[i];
    B->b[i] = b[i];
    for (j = 0; j < s; j++) { B->A[i][j] = A[i * s + j]; }
  }

  if (embedded)
  {
    for (i = 0; i < s; i++) { B->d[i] = d[i]; }
  }

  return (B);
}

/*---------------------------------------------------------------
  Routine to copy a Butcher table structure
  ---------------------------------------------------------------*/
ARKodeButcherTable ARKodeButcherTable_Copy(ARKodeButcherTable B)
{
  int i, j, s;
  ARKodeButcherTable Bcopy;
  sunbooleantype embedded;

  /* Check for legal input */
  if (B == NULL) { return (NULL); }

  /* Get the number of stages */
  s = B->stages;

  /* Does the table have an embedding? */
  embedded = (B->d != NULL) ? SUNTRUE : SUNFALSE;

  /* Allocate Butcher table structure */
  Bcopy = ARKodeButcherTable_Alloc(s, embedded);
  if (Bcopy == NULL) { return (NULL); }

  /* set the relevant parameters */
  Bcopy->stages = B->stages;
  Bcopy->q      = B->q;
  Bcopy->p      = B->p;

  /* Copy Butcher table */
  for (i = 0; i < s; i++)
  {
    Bcopy->c[i] = B->c[i];
    Bcopy->b[i] = B->b[i];
    for (j = 0; j < s; j++) { Bcopy->A[i][j] = B->A[i][j]; }
  }

  if (embedded)
  {
    for (i = 0; i < s; i++) { Bcopy->d[i] = B->d[i]; }
  }

  return (Bcopy);
}

/*---------------------------------------------------------------
  Routine to query the Butcher table structure workspace size
  ---------------------------------------------------------------*/
void ARKodeButcherTable_Space(ARKodeButcherTable B, sunindextype* liw,
                              sunindextype* lrw)
{
  /* initialize outputs and return if B is not allocated */
  *liw = 0;
  *lrw = 0;
  if (B == NULL) { return; }

  /* fill outputs based on B */
  *liw = 3;
  if (B->d != NULL) { *lrw = B->stages * (B->stages + 3); }
  else { *lrw = B->stages * (B->stages + 2); }
}

/*---------------------------------------------------------------
  Routine to free a Butcher table structure
  ---------------------------------------------------------------*/
void ARKodeButcherTable_Free(ARKodeButcherTable B)
{
  int i;

  /* Free each field within Butcher table structure, and then
     free structure itself */
  if (B != NULL)
  {
    if (B->d != NULL) { free(B->d); }
    if (B->c != NULL) { free(B->c); }
    if (B->b != NULL) { free(B->b); }
    if (B->A != NULL)
    {
      for (i = 0; i < B->stages; i++)
      {
        if (B->A[i] != NULL) { free(B->A[i]); }
      }
      free(B->A);
    }

    free(B);
  }
}

/*---------------------------------------------------------------
  Routine to print a Butcher table structure
  ---------------------------------------------------------------*/
void ARKodeButcherTable_Write(ARKodeButcherTable B, FILE* outfile)
{
  int i, j;

  /* check for valid table */
  if (B == NULL) { return; }
  if (B->A == NULL) { return; }
  for (i = 0; i < B->stages; i++)
  {
    if (B->A[i] == NULL) { return; }
  }
  if (B->c == NULL) { return; }
  if (B->b == NULL) { return; }

  fprintf(outfile, "  A = \n");
  for (i = 0; i < B->stages; i++)
  {
    fprintf(outfile, "      ");
    for (j = 0; j < B->stages; j++)
    {
      fprintf(outfile, SUN_FORMAT_E "  ", B->A[i][j]);
    }
    fprintf(outfile, "\n");
  }

  fprintf(outfile, "  c = ");
  for (i = 0; i < B->stages; i++)
  {
    fprintf(outfile, SUN_FORMAT_E "  ", B->c[i]);
  }
  fprintf(outfile, "\n");

  fprintf(outfile, "  b = ");
  for (i = 0; i < B->stages; i++)
  {
    fprintf(outfile, SUN_FORMAT_E "  ", B->b[i]);
  }
  fprintf(outfile, "\n");

  if (B->d != NULL)
  {
    fprintf(outfile, "  d = ");
    for (i = 0; i < B->stages; i++)
    {
      fprintf(outfile, SUN_FORMAT_E "  ", B->d[i]);
    }
    fprintf(outfile, "\n");
  }
}

sunbooleantype ARKodeButcherTable_IsStifflyAccurate(ARKodeButcherTable B)
{
  int i;
  for (i = 0; i < B->stages; i++)
  {
    if (SUNRabs(B->b[i] - B->A[B->stages - 1][i]) > 100 * SUN_UNIT_ROUNDOFF)
    {
      return SUNFALSE;
    }
  }
  return SUNTRUE;
}

/*---------------------------------------------------------------
  Routine to determine the analytical order of accuracy for a
  specified Butcher table.  We check the analytical [necessary]
  order conditions, generated from the rooted trees associated
  with the elementary differentials of the ODE right-hand side
  (see arkode_butcher_trees_impl.h), up through order
  ARK_BUTCHER_MAX_CHECK_ORDER.  After that, we revert to the
  [sufficient] Butcher simplifying assumptions.

  Inputs:
     B: Butcher table to check
     outfile: file pointer to print results; if NULL then no
        outputs are printed

  Outputs:
     q: measured order of accuracy for method
     p: measured order of accuracy for embedding [0 if not present]

  Return values:
     0 (success): internal {q,p} values match analytical order
     1 (warning): internal {q,p} values are lower than analytical
        order, or method achieves maximum order possible with this
        routine and internal {q,p} are higher.
    -1 (failure): internal p and q values are higher than analytical
         order
    -2 (failure): NULL-valued B (or critical contents)

  Note: for embedded methods, if the return flags for p and q would
  differ, failure takes precedence over warning, which takes
  precedence over success.
  ---------------------------------------------------------------*/
int ARKodeButcherTable_CheckOrder(ARKodeButcherTable B, int* q, int* p,
                                  FILE* outfile)
{
  /* local variables */
  int q_SA, p_SA, i, k, s, retval;
  sunrealtype **A, *b, *c, *d;
  sunrealtype* work;
  (*q) = (*p) = 0;

  /* verify non-NULL Butcher table structure and contents */
  if (B == NULL) { return (-2); }
  if (B->stages < 1) { return (-2); }
  if (B->A == NULL) { return (-2); }
  for (i = 0; i < B->stages; i++)
  {
    if (B->A[i] == NULL) { return (-2); }
  }
  if (B->c == NULL) { return (-2); }
  if (B->b == NULL) { return (-2); }

  /* set shortcuts for Butcher table components */
  A = B->A;
  b = B->b;
  c = B->c;
  d = B->d;
  s = B->stages;

  /* allocate workspace for elementary-weight evaluations */
  work = (sunrealtype*)calloc(2 * (ARK_BUTCHER_MAX_CHECK_ORDER + 1) * s,
                              sizeof(sunrealtype));
  if (work == NULL) { return (-2); }

  /* check method order */
  if (outfile) { fprintf(outfile, "ARKodeButcherTable_CheckOrder:\n"); }

  /*    row sum condition */
  if (arkode_butcher_rowsum(A, c, s)) { (*q) = 0; }
  else
  {
    (*q) = -1;
    if (outfile) { fprintf(outfile, "  method fails row sum condition\n"); }
  }
  /*    order conditions, one per rooted tree of each order */
  for (k = 1; k <= ARK_BUTCHER_MAX_CHECK_ORDER; k++)
  {
    if ((*q) != k - 1) { break; }
    retval = arkode_butcher_check_conditions(b, A, c, s, k, work, "method",
                                             outfile);
    if (retval < 0)
    {
      free(work);
      return (-2);
    }
    if (retval == 1) { (*q) = k; }
  }
  /*    higher order conditions (via simplifying assumptions) */
  if ((*q) == ARK_BUTCHER_MAX_CHECK_ORDER)
  {
    if (outfile)
    {
      fprintf(outfile,
              "  method order >= %i; reverting to simplifying assumptions\n",
              ARK_BUTCHER_MAX_CHECK_ORDER);
    }
    q_SA = __ButcherSimplifyingAssumptions(A, b, c, s);
    (*q) = SUNMAX((*q), q_SA);
    if (outfile) { fprintf(outfile, "  method order = %i\n", (*q)); }
  }

  /* check embedding order */
  if (d)
  {
    if (outfile) { fprintf(outfile, "\n"); }

    /*    row sum condition */
    if (arkode_butcher_rowsum(A, c, s)) { (*p) = 0; }
    else
    {
      (*p) = -1;
      if (outfile)
      {
        fprintf(outfile, "  embedding fails row sum condition\n");
      }
    }
    /*    order conditions, one per rooted tree of each order */
    for (k = 1; k <= ARK_BUTCHER_MAX_CHECK_ORDER; k++)
    {
      if ((*p) != k - 1) { break; }
      retval = arkode_butcher_check_conditions(d, A, c, s, k, work, "embedding",
                                               outfile);
      if (retval < 0)
      {
        free(work);
        return (-2);
      }
      if (retval == 1) { (*p) = k; }
    }
    /*    higher order conditions (via simplifying assumptions) */
    if ((*p) == ARK_BUTCHER_MAX_CHECK_ORDER)
    {
      if (outfile)
      {
        fprintf(outfile,
                "  embedding order >= %i; reverting to simplifying "
                "assumptions\n",
                ARK_BUTCHER_MAX_CHECK_ORDER);
      }
      p_SA = __ButcherSimplifyingAssumptions(A, d, c, s);
      (*p) = SUNMAX((*p), p_SA);
      if (outfile) { fprintf(outfile, "  embedding order = %i\n", (*p)); }
    }
  }

  /* clean up */
  free(work);

  /* compare results against stored values and return */

  /*    check failure modes first */
  if (((*q) < B->q) && ((*q) < ARK_BUTCHER_MAX_CHECK_ORDER)) { return (-1); }
  if (d)
  {
    if (((*p) < B->p) && ((*p) < ARK_BUTCHER_MAX_CHECK_ORDER)) { return (-1); }
  }

  /*    check warning modes */
  if ((*q) > B->q) { return (1); }
  if (d)
  {
    if ((*p) > B->p) { return (1); }
  }
  if (((*q) < B->q) && ((*q) >= ARK_BUTCHER_MAX_CHECK_ORDER)) { return (1); }
  if (d)
  {
    if (((*p) < B->p) && ((*p) >= ARK_BUTCHER_MAX_CHECK_ORDER)) { return (1); }
  }

  /*    return success */
  return (0);
}

/*---------------------------------------------------------------
  Routine to determine the analytical order of accuracy for a
  specified pair of Butcher tables in an ARK pair.  We check the
  analytical order conditions, generated from the 2-colored rooted
  trees associated with the elementary differentials of the
  additively-partitioned ODE right-hand side, up through order
  ARK_BUTCHER_MAX_ARK_CHECK_ORDER.

  Inputs:
     B1, B2: Butcher tables to check
     outfile: file pointer to print results; if NULL then no
        outputs are printed

  Outputs:
     q: measured order of accuracy for method
     p: measured order of accuracy for embedding [0 if not present]

  Return values:
     0 (success): completed checks
     1 (warning): internal {q,p} values are lower than analytical
        order, or method achieves maximum order possible with this
        routine and internal {q,p} are higher.
    -1 (failure): NULL-valued B1, B2 (or critical contents)

  Note: for embedded methods, if the return flags for p and q would
  differ, warning takes precedence over success.
  ---------------------------------------------------------------*/
int ARKodeButcherTable_CheckARKOrder(ARKodeButcherTable B1, ARKodeButcherTable B2,
                                     int* q, int* p, FILE* outfile)
{
  /* local variables */
  int i, k, s, retval;
  sunrealtype **A[2], *b[2], *c[2], *d[2];
  sunrealtype* work;
  (*q) = (*p) = 0;

  /* verify non-NULL Butcher table structure and contents */
  if (B1 == NULL) { return (-1); }
  if (B1->stages < 1) { return (-1); }
  if (B1->A == NULL) { return (-1); }
  for (i = 0; i < B1->stages; i++)
  {
    if (B1->A[i] == NULL) { return (-1); }
  }
  if (B1->c == NULL) { return (-1); }
  if (B1->b == NULL) { return (-1); }
  if (B2 == NULL) { return (-1); }
  if (B2->stages < 1) { return (-1); }
  if (B2->A == NULL) { return (-1); }
  for (i = 0; i < B2->stages; i++)
  {
    if (B2->A[i] == NULL) { return (-1); }
  }
  if (B2->c == NULL) { return (-1); }
  if (B2->b == NULL) { return (-1); }
  if (B1->stages != B2->stages) { return (-1); }

  /* set shortcuts for Butcher table components */
  A[0] = B1->A;
  b[0] = B1->b;
  c[0] = B1->c;
  d[0] = B1->d;
  A[1] = B2->A;
  b[1] = B2->b;
  c[1] = B2->c;
  d[1] = B2->d;
  s    = B1->stages;

  /* allocate workspace for elementary-weight evaluations */
  work = (sunrealtype*)calloc(2 * (ARK_BUTCHER_MAX_ARK_CHECK_ORDER + 1) * s,
                              sizeof(sunrealtype));
  if (work == NULL) { return (-1); }

  /* check method order */
  if (outfile) { fprintf(outfile, "ARKodeButcherTable_CheckARKOrder:\n"); }

  /*    row sum conditions */
  if (arkode_butcher_rowsum(A[0], c[0], s) && arkode_butcher_rowsum(A[1], c[1], s))
  {
    (*q) = 0;
  }
  else
  {
    (*q) = -1;
    if (outfile) { fprintf(outfile, "  method fails row sum conditions\n"); }
  }
  /*    order conditions, one per 2-colored rooted tree of each order */
  for (k = 1; k <= ARK_BUTCHER_MAX_ARK_CHECK_ORDER; k++)
  {
    if ((*q) != k - 1) { break; }
    retval = arkode_butcher_check_ark_conditions(b, A, c, s, k, work, "method",
                                                 outfile);
    if (retval < 0)
    {
      free(work);
      return (-1);
    }
    if (retval == 1) { (*q) = k; }
  }

  /* check embedding order */
  if (d[0] && d[1])
  {
    if (outfile) { fprintf(outfile, "\n"); }

    /*    row sum conditions */
    if (arkode_butcher_rowsum(A[0], c[0], s) &&
        arkode_butcher_rowsum(A[1], c[1], s))
    {
      (*p) = 0;
    }
    else
    {
      (*p) = -1;
      if (outfile)
      {
        fprintf(outfile, "  embedding fails row sum conditions\n");
      }
    }
    /*    order conditions, one per 2-colored rooted tree of each order */
    for (k = 1; k <= ARK_BUTCHER_MAX_ARK_CHECK_ORDER; k++)
    {
      if ((*p) != k - 1) { break; }
      retval = arkode_butcher_check_ark_conditions((sunrealtype* const*)d, A, c,
                                                   s, k, work, "embedding",
                                                   outfile);
      if (retval < 0)
      {
        free(work);
        return (-1);
      }
      if (retval == 1) { (*p) = k; }
    }
  }

  /* clean up */
  free(work);

  /* compare results against stored values and return */

  /*    check warning modes */
  if ((*q) > B1->q) { return (1); }
  if ((*q) > B2->q) { return (1); }
  if (d[0] && d[1])
  {
    if ((*p) > B1->p) { return (1); }
    if ((*p) > B2->p) { return (1); }
  }
  if (((*q) < B1->q) && ((*q) == ARK_BUTCHER_MAX_ARK_CHECK_ORDER))
  {
    return (1);
  }
  if (((*q) < B2->q) && ((*q) == ARK_BUTCHER_MAX_ARK_CHECK_ORDER))
  {
    return (1);
  }
  if (d[0] && d[1])
  {
    if (((*p) < B1->p) && ((*p) == ARK_BUTCHER_MAX_ARK_CHECK_ORDER))
    {
      return (1);
    }
    if (((*p) < B2->p) && ((*p) == ARK_BUTCHER_MAX_ARK_CHECK_ORDER))
    {
      return (1);
    }
  }

  /*    return success */
  return (0);
}

/*---------------------------------------------------------------
  Private utility routines for checking method order
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  Utility routine to compute small vector .^ int
       z = x.^l   [Matlab notation]
  Here all vectors are (s x 1).   Returns 0 on success,
  nonzero on failure.
  ---------------------------------------------------------------*/
static int arkode_butcher_vp(sunrealtype* x, int l, int s, sunrealtype* z)
{
  int i;
  if ((x == NULL) || (z == NULL) || (s < 1)) { return (1); }
  for (i = 0; i < s; i++) { z[i] = SUNRpowerI(x[i], l); }
  return (0);
}

/*---------------------------------------------------------------
  Utility routine to compute small vector dot product:
       d = dot(x,y)
  Here x and y are (s x 1), and d is scalar.   Returns 0 on success,
  nonzero on failure.
  ---------------------------------------------------------------*/
static int arkode_butcher_dot(sunrealtype* x, sunrealtype* y, int s,
                              sunrealtype* d)
{
  int i;
  if ((x == NULL) || (y == NULL) || (d == NULL) || (s < 1)) { return (1); }
  (*d) = SUN_RCONST(0.0);
  for (i = 0; i < s; i++) { (*d) += x[i] * y[i]; }
  return (0);
}

/*---------------------------------------------------------------
  Utility routine to check the row sum condition, c(i) = sum(A(i,:)).
  Returns SUNTRUE on success, SUNFALSE on failure.
  ---------------------------------------------------------------*/
static sunbooleantype arkode_butcher_rowsum(sunrealtype** A, sunrealtype* c, int s)
{
  int i, j;
  sunrealtype rsum;
  for (i = 0; i < s; i++)
  {
    rsum = SUN_RCONST(0.0);
    for (j = 0; j < s; j++) { rsum += A[i][j]; }
    if (SUNRabs(rsum - c[i]) > TOL) { return (SUNFALSE); }
  }
  return (SUNTRUE);
}

/*---------------------------------------------------------------
  Utility routine to check all analytical order conditions of a
  given order for a single Butcher table (b or d, A, c).  The
  conditions Phi(t) = 1/gamma(t) are generated from the rooted
  trees of the requested order.  Any failed condition is reported
  to outfile (if non-NULL) using the supplied name ("method" or
  "embedding").  The work array must hold at least
  2*(order+1)*s entries.

  Returns 1 if all conditions hold, 0 if any condition fails, and
  -1 on an internal error.
  ---------------------------------------------------------------*/
static int arkode_butcher_check_conditions(sunrealtype* b, sunrealtype** A,
                                           sunrealtype* c, int s, int order,
                                           sunrealtype* work, const char* name,
                                           FILE* outfile)
{
  ARKodeButcherTreeIter iter;
  long int gamma;
  sunrealtype phi;
  sunbooleantype alltrue, more;
  char expr[128];
  sunrealtype* b_[1];
  sunrealtype** A_[1];
  sunrealtype* c_[1];

  b_[0] = b;
  A_[0] = A;
  c_[0] = c;

  if (arkodeButcherTrees_IterInit(&iter, order)) { return (-1); }
  alltrue = SUNTRUE;
  more    = SUNTRUE;
  while (more)
  {
    gamma = arkodeButcherTrees_Density(iter.levels, order);
    if ((gamma < 1) || arkodeButcherTrees_Weight(iter.levels, NULL, order, b_,
                                                 A_, c_, s, work, &phi))
    {
      arkodeButcherTrees_IterFree(&iter);
      return (-1);
    }
    if (SUNRabs(phi - SUN_RCONST(1.0) / ((sunrealtype)gamma)) > TOL)
    {
      alltrue = SUNFALSE;
      if (outfile)
      {
        if (arkodeButcherTrees_WeightString(iter.levels, order, expr,
                                            sizeof(expr)))
        {
          expr[0] = '\0';
        }
        fprintf(outfile, "  %s fails order %i condition %s = 1/%ld\n", name,
                order, expr, gamma);
      }
    }
    more = arkodeButcherTrees_IterNext(&iter);
  }
  arkodeButcherTrees_IterFree(&iter);
  return (alltrue ? 1 : 0);
}

/*---------------------------------------------------------------
  Utility routine to check all analytical order conditions of a
  given order for an ARK pair of Butcher tables.  For each rooted
  tree of the requested order every assignment of the two tables
  ("colors") to the tree nodes is checked: the root color selects
  the b (or d) vector while every other node selects the A matrix
  (or, at a leaf, the c vector) of the corresponding table.  Any
  tree with a failed coloring is reported to outfile (if non-NULL)
  using the supplied name ("method" or "embedding").  The work
  array must hold at least 2*(order+1)*s entries.

  Returns 1 if all conditions hold, 0 if any condition fails, and
  -1 on an internal error.
  ---------------------------------------------------------------*/
static int arkode_butcher_check_ark_conditions(sunrealtype* const* b,
                                               sunrealtype** const* A,
                                               sunrealtype* const* c, int s,
                                               int order, sunrealtype* work,
                                               const char* name, FILE* outfile)
{
  ARKodeButcherTreeIter iter;
  long int gamma, mask, ncolorings;
  int i, colors[ARK_BUTCHER_MAX_ARK_CHECK_ORDER];
  sunrealtype phi;
  sunbooleantype alltrue, treetrue, more;
  char expr[128];

  if (order > ARK_BUTCHER_MAX_ARK_CHECK_ORDER) { return (-1); }
  if (arkodeButcherTrees_IterInit(&iter, order)) { return (-1); }
  alltrue = SUNTRUE;
  more    = SUNTRUE;
  while (more)
  {
    gamma = arkodeButcherTrees_Density(iter.levels, order);
    if (gamma < 1)
    {
      arkodeButcherTrees_IterFree(&iter);
      return (-1);
    }
    treetrue   = SUNTRUE;
    ncolorings = 1L << order;
    for (mask = 0; mask < ncolorings; mask++)
    {
      for (i = 0; i < order; i++) { colors[i] = (int)((mask >> i) & 1L); }
      if (arkodeButcherTrees_Weight(iter.levels, colors, order, b, A, c, s,
                                    work, &phi))
      {
        arkodeButcherTrees_IterFree(&iter);
        return (-1);
      }
      if (SUNRabs(phi - SUN_RCONST(1.0) / ((sunrealtype)gamma)) > TOL)
      {
        treetrue = SUNFALSE;
        break;
      }
    }
    if (!treetrue)
    {
      alltrue = SUNFALSE;
      if (outfile)
      {
        if (arkodeButcherTrees_WeightString(iter.levels, order, expr,
                                            sizeof(expr)))
        {
          expr[0] = '\0';
        }
        fprintf(outfile, "  %s fails order %i conditions %s = 1/%ld\n", name,
                order, expr, gamma);
      }
    }
    more = arkodeButcherTrees_IterNext(&iter);
  }
  arkodeButcherTrees_IterFree(&iter);
  return (alltrue ? 1 : 0);
}

/*---------------------------------------------------------------
  Utility routine to check Butcher's simplifying assumptions.
  Returns the maximum predicted order.
  ---------------------------------------------------------------*/
static int __ButcherSimplifyingAssumptions(sunrealtype** A, sunrealtype* b,
                                           sunrealtype* c, int s)
{
  int P, Q, R, i, j, k, q;
  sunrealtype RHS, LHS;
  sunbooleantype alltrue;
  sunrealtype* tmp = calloc(s, sizeof(sunrealtype));

  /* B(P) */
  P = 0;
  for (i = 1; i < 1000; i++)
  {
    if (arkode_butcher_vp(c, i - 1, s, tmp))
    {
      free(tmp);
      return (0);
    }
    if (arkode_butcher_dot(b, tmp, s, &LHS))
    {
      free(tmp);
      return (0);
    }
    RHS = SUN_RCONST(1.0) / i;
    if (SUNRabs(RHS - LHS) > TOL) { break; }
    P++;
  }

  /* C(Q) */
  Q = 0;
  for (k = 1; k < 1000; k++)
  {
    alltrue = SUNTRUE;
    for (i = 0; i < s; i++)
    {
      if (arkode_butcher_vp(c, k - 1, s, tmp))
      {
        free(tmp);
        return (0);
      }
      if (arkode_butcher_dot(A[i], tmp, s, &LHS))
      {
        free(tmp);
        return (0);
      }
      RHS = SUNRpowerI(c[i], k) / k;
      if (SUNRabs(RHS - LHS) > TOL)
      {
        alltrue = SUNFALSE;
        break;
      }
    }
    if (alltrue) { Q++; }
    else { break; }
  }

  /* D(R) */
  R = 0;
  for (k = 1; k < 1000; k++)
  {
    alltrue = SUNTRUE;
    for (j = 0; j < s; j++)
    {
      LHS = SUN_RCONST(0.0);
      for (i = 0; i < s; i++)
      {
        LHS += A[i][j] * b[i] * SUNRpowerI(c[i], k - 1);
      }
      RHS = b[j] / k * (SUN_RCONST(1.0) - SUNRpowerI(c[j], k));
      if (SUNRabs(RHS - LHS) > TOL)
      {
        alltrue = SUNFALSE;
        break;
      }
    }
    if (alltrue) { R++; }
    else { break; }
  }

  /* determine q, clean up and return */
  q = 0;
  for (i = 1; i <= P; i++)
  {
    if ((q > Q + R + 1) || (q > 2 * Q + 2)) { break; }
    q++;
  }
  free(tmp);
  return (q);
}

/*---------------------------------------------------------------
  EOF
  ---------------------------------------------------------------*/