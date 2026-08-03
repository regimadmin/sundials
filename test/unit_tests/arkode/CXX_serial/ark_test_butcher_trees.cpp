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
 *----------------------------------------------------------------
 * Routine to test the ARKODE rooted-tree utilities that generate
 * the Runge--Kutta order conditions from the elementary
 * differentials of the ODE right-hand side:
 * - tree counts per order match OEIS A000081,
 * - the per-order identity sum_t 1/(sigma(t)*gamma(t)) = 1/n holds,
 * - elementary weights Phi(t) certify/reject known Butcher tables,
 * - rendered expression strings match the literature forms.
 *-----------------------------------------------------------------*/

// Header files
#include <arkode/arkode_butcher.h>
#include <arkode/arkode_butcher_erk.h>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include "arkode/arkode_butcher_trees_impl.h"

// expected number of rooted trees per order (OEIS A000081)
static const long int A000081[] = {1, 1, 2, 4, 9, 20, 48, 115, 286, 719};
static const int NORDERS        = 10;

// evaluate all order conditions of a given order for a Butcher table,
// returning the number of failed conditions
static int check_conditions(ARKodeButcherTable B, int order, sunrealtype tol)
{
  ARKodeButcherTreeIter iter;
  sunrealtype* b[1];
  sunrealtype** A[1];
  sunrealtype* c[1];
  int nfailed = 0;

  b[0] = B->b;
  A[0] = B->A;
  c[0] = B->c;

  sunrealtype* work = new sunrealtype[2 * (order + 1) * B->stages];
  if (arkodeButcherTrees_IterInit(&iter, order))
  {
    delete[] work;
    return -1;
  }
  sunbooleantype more = SUNTRUE;
  while (more)
  {
    long int gamma  = arkodeButcherTrees_Density(iter.levels, order);
    sunrealtype phi = SUN_RCONST(0.0);
    if ((gamma < 1) || arkodeButcherTrees_Weight(iter.levels, NULL, order, b,
                                                 A, c, B->stages, work, &phi))
    {
      arkodeButcherTrees_IterFree(&iter);
      delete[] work;
      return -1;
    }
    if (SUNRabs(phi - SUN_RCONST(1.0) / ((sunrealtype)gamma)) > tol)
    {
      nfailed++;
    }
    more = arkodeButcherTrees_IterNext(&iter);
  }
  arkodeButcherTrees_IterFree(&iter);
  delete[] work;
  return nfailed;
}

// Main Program
int main()
{
  int numfails = 0;

  // Test 1: enumeration counts match OEIS A000081 and the
  // recurrence-based counter
  printf("Rooted-tree enumeration counts (OEIS A000081):\n");
  for (int n = 1; n <= NORDERS; n++)
  {
    ARKodeButcherTreeIter iter;
    long int count = 0;
    if (arkodeButcherTrees_IterInit(&iter, n)) { count = -1; }
    else
    {
      sunbooleantype more = SUNTRUE;
      while (more)
      {
        count++;
        more = arkodeButcherTrees_IterNext(&iter);
      }
      arkodeButcherTrees_IterFree(&iter);
    }
    const bool pass = (count == A000081[n - 1]) &&
                      (arkodeButcherTrees_Count(n) == A000081[n - 1]);
    printf("  order %2i: %4ld trees  %s\n", n, count, pass ? "PASS" : "FAIL");
    if (!pass) { numfails++; }
  }

  // Test 2: per-order identity sum_t n!/(sigma(t)*gamma(t)) = (n-1)!,
  // i.e. sum_t 1/(sigma(t)*gamma(t)) = 1/n, with every summand
  // alpha(t) = n!/(sigma(t)*gamma(t)) an exact integer
  printf("\nPer-order identity sum(1/(sigma*gamma)) = 1/n:\n");
  for (int n = 1; n <= NORDERS; n++)
  {
    long int nfact = 1;
    for (int i = 2; i <= n; i++) { nfact *= i; }
    ARKodeButcherTreeIter iter;
    long int sum = 0;
    bool pass    = true;
    if (arkodeButcherTrees_IterInit(&iter, n)) { pass = false; }
    else
    {
      sunbooleantype more = SUNTRUE;
      while (more)
      {
        const long int gamma = arkodeButcherTrees_Density(iter.levels, n);
        const long int sigma = arkodeButcherTrees_Symmetry(iter.levels, n);
        if ((gamma < 1) || (sigma < 1) || (nfact % (sigma * gamma) != 0))
        {
          pass = false;
        }
        else { sum += nfact / (sigma * gamma); }
        more = arkodeButcherTrees_IterNext(&iter);
      }
      arkodeButcherTrees_IterFree(&iter);
    }
    pass = pass && (sum == nfact / n);
    printf("  order %2i: %s\n", n, pass ? "PASS" : "FAIL");
    if (!pass) { numfails++; }
  }

  // Test 3: expression strings for all trees of orders 1-4
  printf("\nElementary differentials and weights (orders 1-4):\n");
  for (int n = 1; n <= 4; n++)
  {
    ARKodeButcherTreeIter iter;
    if (arkodeButcherTrees_IterInit(&iter, n))
    {
      numfails++;
      continue;
    }
    sunbooleantype more = SUNTRUE;
    while (more)
    {
      char diff[128], weight[128];
      if (arkodeButcherTrees_ElementaryDiff(iter.levels, n, diff,
                                            sizeof(diff)) ||
          arkodeButcherTrees_WeightString(iter.levels, n, weight,
                                          sizeof(weight)))
      {
        numfails++;
      }
      else
      {
        printf("  order %i: F(t) = %-16s Phi(t) = %-24s gamma(t) = %ld\n", n,
               diff, weight, arkodeButcherTrees_Density(iter.levels, n));
      }
      more = arkodeButcherTrees_IterNext(&iter);
    }
    arkodeButcherTrees_IterFree(&iter);
  }

  // Test 4: the previously missing order-6 condition
  // b'*((A*(c.*c)).*(A*c)) = 1/36 is enumerated
  {
    ARKodeButcherTreeIter iter;
    bool found = false;
    if (!arkodeButcherTrees_IterInit(&iter, 6))
    {
      sunbooleantype more = SUNTRUE;
      while (more)
      {
        char weight[128];
        if (!arkodeButcherTrees_WeightString(iter.levels, 6, weight,
                                             sizeof(weight)) &&
            !strcmp(weight, "b'*((A*(c.*c)).*(A*c))") &&
            (arkodeButcherTrees_Density(iter.levels, 6) == 36))
        {
          found = true;
        }
        more = arkodeButcherTrees_IterNext(&iter);
      }
      arkodeButcherTrees_IterFree(&iter);
    }
    printf("\nOrder-6 condition b'*((A*(c.*c)).*(A*c)) = 1/36 enumerated: %s\n",
           found ? "PASS" : "FAIL");
    if (!found) { numfails++; }
  }

  // Test 5: elementary weights certify classical RK4 exactly through
  // order 4 and reject it at order 5
  {
    ARKodeButcherTable B = ARKodeButcherTable_Alloc(4, SUNFALSE);
    B->q                 = 4;
    B->A[1][0]           = SUN_RCONST(0.5);
    B->A[2][1]           = SUN_RCONST(0.5);
    B->A[3][2]           = SUN_RCONST(1.0);
    B->b[0]              = SUN_RCONST(1.0) / SUN_RCONST(6.0);
    B->b[1]              = SUN_RCONST(1.0) / SUN_RCONST(3.0);
    B->b[2]              = SUN_RCONST(1.0) / SUN_RCONST(3.0);
    B->b[3]              = SUN_RCONST(1.0) / SUN_RCONST(6.0);
    B->c[1]              = SUN_RCONST(0.5);
    B->c[2]              = SUN_RCONST(0.5);
    B->c[3]              = SUN_RCONST(1.0);
    const sunrealtype tol = SUNRsqrt(SUN_UNIT_ROUNDOFF);
    bool pass             = true;
    for (int n = 1; n <= 4; n++)
    {
      if (check_conditions(B, n, tol) != 0) { pass = false; }
    }
    printf("\nClassical RK4 satisfies all conditions of orders 1-4: %s\n",
           pass ? "PASS" : "FAIL");
    if (!pass) { numfails++; }
    const int nfail5 = check_conditions(B, 5, tol);
    printf("Classical RK4 fails %s order-5 condition: %s\n",
           (nfail5 > 0) ? "at least one" : "no",
           (nfail5 > 0) ? "PASS" : "FAIL");
    if (nfail5 <= 0) { numfails++; }
    ARKodeButcherTable_Free(B);
  }

  // Test 6: the 3-stage Gauss--Legendre table satisfies all conditions
  // through order 6 and fails some order-7 condition
  {
    const sunrealtype r15 = SUNRsqrt(SUN_RCONST(15.0));
    ARKodeButcherTable B  = ARKodeButcherTable_Alloc(3, SUNFALSE);
    B->q                  = 6;
    B->A[0][0]            = SUN_RCONST(5.0) / SUN_RCONST(36.0);
    B->A[0][1] = SUN_RCONST(2.0) / SUN_RCONST(9.0) - r15 / SUN_RCONST(15.0);
    B->A[0][2] = SUN_RCONST(5.0) / SUN_RCONST(36.0) - r15 / SUN_RCONST(30.0);
    B->A[1][0] = SUN_RCONST(5.0) / SUN_RCONST(36.0) + r15 / SUN_RCONST(24.0);
    B->A[1][1] = SUN_RCONST(2.0) / SUN_RCONST(9.0);
    B->A[1][2] = SUN_RCONST(5.0) / SUN_RCONST(36.0) - r15 / SUN_RCONST(24.0);
    B->A[2][0] = SUN_RCONST(5.0) / SUN_RCONST(36.0) + r15 / SUN_RCONST(30.0);
    B->A[2][1] = SUN_RCONST(2.0) / SUN_RCONST(9.0) + r15 / SUN_RCONST(15.0);
    B->A[2][2] = SUN_RCONST(5.0) / SUN_RCONST(36.0);
    B->b[0]    = SUN_RCONST(5.0) / SUN_RCONST(18.0);
    B->b[1]    = SUN_RCONST(4.0) / SUN_RCONST(9.0);
    B->b[2]    = SUN_RCONST(5.0) / SUN_RCONST(18.0);
    B->c[0]    = SUN_RCONST(0.5) - r15 / SUN_RCONST(10.0);
    B->c[1]    = SUN_RCONST(0.5);
    B->c[2]    = SUN_RCONST(0.5) + r15 / SUN_RCONST(10.0);
    const sunrealtype tol = SUNRsqrt(SUN_UNIT_ROUNDOFF);
    bool pass             = true;
    for (int n = 1; n <= 6; n++)
    {
      if (check_conditions(B, n, tol) != 0) { pass = false; }
    }
    printf("\nGauss-Legendre (3 stage) satisfies all conditions of orders "
           "1-6: %s\n",
           pass ? "PASS" : "FAIL");
    if (!pass) { numfails++; }
    const int nfail7 = check_conditions(B, 7, tol);
    printf("Gauss-Legendre (3 stage) fails %s order-7 condition: %s\n",
           (nfail7 > 0) ? "at least one" : "no",
           (nfail7 > 0) ? "PASS" : "FAIL");
    if (nfail7 <= 0) { numfails++; }
    ARKodeButcherTable_Free(B);
  }

  if (numfails) { printf("\nFailed %i tests\n", numfails); }
  else { printf("\nAll rooted-tree tests passed!\n"); }
  return numfails;
}

/*---- end of file ----*/
