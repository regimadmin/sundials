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
 * Routine to test the ARKODE rooted-tree module that expresses
 * the elementary differentials of an autonomous ODE y' = f(y):
 *   * enumeration counts per order match OEIS A000081 and the
 *     divisor-sum recurrence,
 *   * the per-tree functionals satisfy the labeling identity
 *     sum over trees of order n of n!/(sigma(t)*gamma(t)) = (n-1)!,
 *   * rendered elementary differential and weight expressions
 *     match their literature forms through order 4, and
 *   * elementary weights Phi(t) reproduce the analytic orders of
 *     classical RK4, the 3-stage Gauss method, and the Heun-Euler
 *     embedded pair.
 *-----------------------------------------------------------------*/

// Header files
#include <arkode/arkode.h>
#include <arkode/arkode_butcher.h>
#include <cstring>
#include <iostream>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include "arkode/arkode_butcher_trees_impl.h"

#define MAX_TREE_ORDER 10

// Reference counts of rooted trees per order, OEIS A000081
static const long int A000081[MAX_TREE_ORDER + 1] = {0,  1,  1,   2,   4,  9,
                                                     20, 48, 115, 286, 719};

static long int factorial(int n)
{
  long int nfact = 1;
  for (int i = 2; i <= n; i++) { nfact *= i; }
  return nfact;
}

// Verify the enumeration visits A000081(n) trees for each order n, in
// agreement with the recurrence-based counter
static int test_counts()
{
  int numfails = 0;
  std::cout << "\nRooted tree enumeration and A000081 counts:\n";
  for (int n = 1; n <= MAX_TREE_ORDER; n++)
  {
    ARKodeButcherTreeIter iter;
    if (arkButcherTreeIterInit(n, &iter) != ARK_SUCCESS)
    {
      std::cout << "  order " << n << ": iterator initialization failed\n";
      numfails++;
      continue;
    }
    long int enumerated = 0;
    while (arkButcherTreeIterNext(&iter)) { enumerated++; }
    arkButcherTreeIterFree(&iter);

    long int counted = arkButcherTreeCount(n);
    std::cout << "  order " << n << ": " << enumerated << " tree(s)";
    if (enumerated != A000081[n] || counted != A000081[n])
    {
      std::cout << " but A000081 gives " << A000081[n] << " and recurrence "
                << counted << "  FAIL";
      numfails++;
    }
    std::cout << "\n";
  }
  if (numfails == 0)
  {
    std::cout << "  enumeration matches A000081 and the recurrence counter\n";
  }
  return numfails;
}

// Verify the joint identity on the enumeration, symmetry, and density:
// alpha(t) = n!/(sigma(t)*gamma(t)) is a positive integer for every tree,
// and the alpha values of the trees of order n sum to (n-1)!
static int test_alpha_identity()
{
  int numfails = 0;
  std::cout << "\nLabeling identity, sum of n!/(sigma*gamma) = (n-1)!:\n";
  for (int n = 1; n <= MAX_TREE_ORDER; n++)
  {
    ARKodeButcherTreeIter iter;
    if (arkButcherTreeIterInit(n, &iter) != ARK_SUCCESS)
    {
      std::cout << "  order " << n << ": iterator initialization failed\n";
      numfails++;
      continue;
    }
    long int nfact = factorial(n);
    long int sum   = 0;
    while (arkButcherTreeIterNext(&iter))
    {
      long int sigma = arkButcherTreeSigma(iter.levels, n);
      long int gamma = arkButcherTreeGamma(iter.levels, n);
      if (sigma < 1 || gamma < 1 || nfact % (sigma * gamma) != 0)
      {
        std::cout << "  order " << n << ": invalid sigma/gamma pair: " << sigma
                  << "/" << gamma << "  FAIL\n";
        numfails++;
      }
      else { sum += nfact / (sigma * gamma); }
    }
    arkButcherTreeIterFree(&iter);

    std::cout << "  order " << n << ": sum = " << sum;
    if (sum != factorial(n - 1))
    {
      std::cout << " but (n-1)! = " << factorial(n - 1) << "  FAIL";
      numfails++;
    }
    std::cout << "\n";
  }
  if (numfails == 0) { std::cout << "  identity satisfied for all orders\n"; }
  return numfails;
}

// Verify the rendered elementary differentials and weights against their
// literature forms through order 4 (trees in Beyer-Hedetniemi order)
static int test_expressions()
{
  struct expected_tree
  {
    int order;
    const char* diff;
    const char* weight;
    long int sigma;
    long int gamma;
  };

  static const expected_tree expected[] = {{1, "f", "b'*e", 1, 1},
                                           {2, "f'f", "b'*c", 1, 2},
                                           {3, "f'f'f", "b'*(A*c)", 1, 6},
                                           {3, "f''(f,f)", "b'*c.^2", 2, 3},
                                           {4, "f'f'f'f", "b'*(A*(A*c))", 1, 24},
                                           {4, "f'f''(f,f)", "b'*(A*c.^2)", 2, 12},
                                           {4, "f''(f'f,f)", "b'*((A*c).*c)", 1,
                                            8},
                                           {4, "f'''(f,f,f)", "b'*c.^3", 6, 4}};

  int numfails = 0;
  int index    = 0;
  std::cout << "\nElementary differentials and weights through order 4:\n";
  for (int n = 1; n <= 4; n++)
  {
    ARKodeButcherTreeIter iter;
    if (arkButcherTreeIterInit(n, &iter) != ARK_SUCCESS)
    {
      std::cout << "  order " << n << ": iterator initialization failed\n";
      numfails++;
      continue;
    }
    while (arkButcherTreeIterNext(&iter))
    {
      char* diff     = arkButcherTreeDiffString(iter.levels, n);
      char* weight   = arkButcherTreeWeightString(iter.levels, n);
      long int sigma = arkButcherTreeSigma(iter.levels, n);
      long int gamma = arkButcherTreeGamma(iter.levels, n);
      if (diff == NULL || weight == NULL)
      {
        std::cout << "  order " << n << ": rendering failed  FAIL\n";
        numfails++;
      }
      else
      {
        std::cout << "  order " << n << ": F(t) = " << diff
                  << ", Phi(t) = " << weight << " = 1/" << gamma
                  << ", sigma = " << sigma << "\n";
        const expected_tree& ref = expected[index];
        if (n != ref.order || strcmp(diff, ref.diff) != 0 ||
            strcmp(weight, ref.weight) != 0 || sigma != ref.sigma ||
            gamma != ref.gamma)
        {
          std::cout << "    expected F(t) = " << ref.diff
                    << ", Phi(t) = " << ref.weight << " = 1/" << ref.gamma
                    << ", sigma = " << ref.sigma << "  FAIL\n";
          numfails++;
        }
      }
      free(diff);
      free(weight);
      index++;
    }
    arkButcherTreeIterFree(&iter);
  }
  if (numfails == 0)
  {
    std::cout << "  rendered expressions match reference forms\n";
  }
  return numfails;
}

// Count the order-n conditions Phi(t) = 1/gamma(t) that the method (or, with
// embedding = SUNTRUE, its embedding) fails to satisfy; -1 signals an error
static long int failed_conditions(ARKodeButcherTable B, int n,
                                  sunbooleantype embedding)
{
  sunrealtype tol = SUN_RCONST(100.0) * B->stages * SUN_UNIT_ROUNDOFF;

  ARKodeButcherTreeIter iter;
  if (arkButcherTreeIterInit(n, &iter) != ARK_SUCCESS) { return -1; }
  long int fails = 0;
  while (arkButcherTreeIterNext(&iter))
  {
    sunrealtype phi, phi_hat;
    long int gamma = arkButcherTreeGamma(iter.levels, n);
    if (gamma < 1 ||
        arkButcherTreePhi(iter.levels, n, B, &phi, &phi_hat) != ARK_SUCCESS)
    {
      arkButcherTreeIterFree(&iter);
      return -1;
    }
    sunrealtype value = embedding ? phi_hat : phi;
    if (SUNRabs(value - SUN_RCONST(1.0) / gamma) > tol) { fails++; }
  }
  arkButcherTreeIterFree(&iter);
  return fails;
}

// Verify a method (embedding = SUNFALSE) or its embedding (embedding =
// SUNTRUE) satisfies all conditions through the expected order and fails at
// least one condition of the next order
static int check_expected_order(ARKodeButcherTable B, const char* name,
                                sunbooleantype embedding, int order)
{
  int numfails    = 0;
  const char* who = embedding ? " embedding" : "";

  for (int n = 1; n <= order; n++)
  {
    long int fails = failed_conditions(B, n, embedding);
    if (fails != 0)
    {
      std::cout << "  " << name << who << " fails " << fails << " order " << n
                << " condition(s)  FAIL\n";
      numfails++;
    }
  }
  if (numfails == 0)
  {
    std::cout << "  " << name << who << " satisfies all conditions through order "
              << order << ": PASS\n";
  }

  long int fails = failed_conditions(B, order + 1, embedding);
  if (fails > 0)
  {
    std::cout << "  " << name << who << " fails at least one order "
              << order + 1 << " condition: PASS\n";
  }
  else
  {
    std::cout << "  " << name << who << " unexpectedly satisfies all order "
              << order + 1 << " conditions  FAIL\n";
    numfails++;
  }
  return numfails;
}

// Verify elementary weights reproduce the analytic orders of classical
// tables: RK4 (order 4), the 3-stage Gauss method (order 6), and the
// Heun-Euler pair (order 2 with an order 1 embedding)
static int test_elementary_weights()
{
  int numfails = 0;
  std::cout << "\nElementary weight checks:\n";

  // Classical 4th-order Runge-Kutta method
  sunrealtype rk4_c[]    = {SUN_RCONST(0.0), SUN_RCONST(0.5), SUN_RCONST(0.5),
                            SUN_RCONST(1.0)};
  sunrealtype rk4_A[]    = {SUN_RCONST(0.0), SUN_RCONST(0.0), SUN_RCONST(0.0),
                            SUN_RCONST(0.0), SUN_RCONST(0.5), SUN_RCONST(0.0),
                            SUN_RCONST(0.0), SUN_RCONST(0.0), SUN_RCONST(0.0),
                            SUN_RCONST(0.5), SUN_RCONST(0.0), SUN_RCONST(0.0),
                            SUN_RCONST(0.0), SUN_RCONST(0.0), SUN_RCONST(1.0),
                            SUN_RCONST(0.0)};
  sunrealtype rk4_b[]    = {SUN_RCONST(1.0) / SUN_RCONST(6.0),
                            SUN_RCONST(1.0) / SUN_RCONST(3.0),
                            SUN_RCONST(1.0) / SUN_RCONST(3.0),
                            SUN_RCONST(1.0) / SUN_RCONST(6.0)};
  ARKodeButcherTable rk4 = ARKodeButcherTable_Create(4, 4, 0, rk4_c, rk4_A,
                                                     rk4_b, NULL);
  if (rk4 == NULL)
  {
    std::cout << "  error creating RK4 table\n";
    return 1;
  }
  numfails += check_expected_order(rk4, "RK4", SUNFALSE, 4);
  ARKodeButcherTable_Free(rk4);

  // 3-stage Gauss-Legendre method of order 6
  sunrealtype w         = SUNRsqrt(SUN_RCONST(15.0));
  sunrealtype gauss_c[] = {SUN_RCONST(0.5) - w / SUN_RCONST(10.0),
                           SUN_RCONST(0.5),
                           SUN_RCONST(0.5) + w / SUN_RCONST(10.0)};
  sunrealtype gauss_A[] =
    {SUN_RCONST(5.0) / SUN_RCONST(36.0),
     SUN_RCONST(2.0) / SUN_RCONST(9.0) - w / SUN_RCONST(15.0),
     SUN_RCONST(5.0) / SUN_RCONST(36.0) - w / SUN_RCONST(30.0),
     SUN_RCONST(5.0) / SUN_RCONST(36.0) + w / SUN_RCONST(24.0),
     SUN_RCONST(2.0) / SUN_RCONST(9.0),
     SUN_RCONST(5.0) / SUN_RCONST(36.0) - w / SUN_RCONST(24.0),
     SUN_RCONST(5.0) / SUN_RCONST(36.0) + w / SUN_RCONST(30.0),
     SUN_RCONST(2.0) / SUN_RCONST(9.0) + w / SUN_RCONST(15.0),
     SUN_RCONST(5.0) / SUN_RCONST(36.0)};
  sunrealtype gauss_b[]    = {SUN_RCONST(5.0) / SUN_RCONST(18.0),
                              SUN_RCONST(4.0) / SUN_RCONST(9.0),
                              SUN_RCONST(5.0) / SUN_RCONST(18.0)};
  ARKodeButcherTable gauss = ARKodeButcherTable_Create(3, 6, 0, gauss_c,
                                                       gauss_A, gauss_b, NULL);
  if (gauss == NULL)
  {
    std::cout << "  error creating Gauss table\n";
    return numfails + 1;
  }
  numfails += check_expected_order(gauss, "Gauss-3", SUNFALSE, 6);
  ARKodeButcherTable_Free(gauss);

  // Heun-Euler embedded pair: order 2 method with order 1 embedding
  sunrealtype he_c[]    = {SUN_RCONST(0.0), SUN_RCONST(1.0)};
  sunrealtype he_A[]    = {SUN_RCONST(0.0), SUN_RCONST(0.0), SUN_RCONST(1.0),
                           SUN_RCONST(0.0)};
  sunrealtype he_b[]    = {SUN_RCONST(0.5), SUN_RCONST(0.5)};
  sunrealtype he_d[]    = {SUN_RCONST(1.0), SUN_RCONST(0.0)};
  ARKodeButcherTable he = ARKodeButcherTable_Create(2, 2, 1, he_c, he_A, he_b,
                                                    he_d);
  if (he == NULL)
  {
    std::cout << "  error creating Heun-Euler table\n";
    return numfails + 1;
  }
  numfails += check_expected_order(he, "Heun-Euler", SUNFALSE, 2);
  numfails += check_expected_order(he, "Heun-Euler", SUNTRUE, 1);
  ARKodeButcherTable_Free(he);

  return numfails;
}

// Verify invalid inputs are rejected
static int test_invalid_inputs()
{
  int numfails = 0;
  std::cout << "\nInvalid input checks:\n";

  ARKodeButcherTreeIter iter;
  if (arkButcherTreeIterInit(0, &iter) != ARK_ILL_INPUT) { numfails++; }
  arkButcherTreeIterFree(&iter);
  if (arkButcherTreeCount(0) != 0) { numfails++; }

  const int not_a_tree[] = {1, 3};
  if (arkButcherTreeOrder(not_a_tree, 2) != 0) { numfails++; }
  if (arkButcherTreeSigma(not_a_tree, 2) != 0) { numfails++; }
  if (arkButcherTreeGamma(not_a_tree, 2) != 0) { numfails++; }
  if (arkButcherTreeDiffString(not_a_tree, 2) != NULL) { numfails++; }
  if (arkButcherTreeWeightString(not_a_tree, 2) != NULL) { numfails++; }

  const int leaf[] = {1};
  sunrealtype phi;
  if (arkButcherTreePhi(leaf, 1, NULL, &phi, NULL) != ARK_ILL_INPUT)
  {
    numfails++;
  }

  if (numfails == 0) { std::cout << "  invalid inputs rejected\n"; }
  else { std::cout << "  " << numfails << " invalid input check(s)  FAIL\n"; }
  return numfails;
}

// Main Program
int main()
{
  int numfails = 0;

  numfails += test_counts();
  numfails += test_alpha_identity();
  numfails += test_expressions();
  numfails += test_elementary_weights();
  numfails += test_invalid_inputs();

  // determine overall success/failure and return
  if (numfails == 0)
  {
    std::cout << "\nAll rooted-tree tests passed\n";
    return 0;
  }
  else
  {
    std::cout << "\n" << numfails << " rooted-tree test(s) failed\n";
    return 1;
  }
}

/*---- end of file ----*/
