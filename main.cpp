#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <ostream>
#include <sstream>
#include <vector>

#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#include <objscip/objscipdefplugins.h>

#include "cons_labsnogood.hpp"

struct Quad
{
  int indices[4];

  Quad(const Quad& other)
  {
    indices[0] = other.indices[0];
    indices[1] = other.indices[1];
    indices[2] = other.indices[2];
    indices[3] = other.indices[3];
  }

  Quad(int* inds)
  {
    indices[0] = inds[0];
    indices[1] = inds[1];
    indices[2] = inds[2];
    indices[3] = inds[3];
    std::sort(&indices[0], &indices[4]);
  }

  Quad(int a, int b, int c, int d)
  {
    indices[0] = a;
    indices[1] = b;
    indices[2] = c;
    indices[3] = d;
    std::sort(&indices[0], &indices[4]);
  }

  bool operator==(const Quad& other) const
  {
    if (indices[0] != other.indices[0])
      return false;
    if (indices[1] != other.indices[1])
      return false;
    if (indices[2] != other.indices[2])
      return false;
    if (indices[3] != other.indices[3])
      return false;

    return true;
  }
};

namespace std
{
  template<> struct hash<Quad>
  {
    std::size_t operator()(const Quad& quad) const
    {
      size_t hash0 = std::hash<int>()(quad.indices[0]);
      size_t hash1 = std::hash<int>()(quad.indices[1]);
      size_t hash2 = std::hash<int>()(quad.indices[2]);
      size_t hash3 = std::hash<int>()(quad.indices[3]);
      return hash0 ^ (hash1-1) ^ (hash2-2) ^ (hash3-3);
    }
  };
};

std::ostream& operator<<(std::ostream& stream, const Quad& quad)
{
  return stream << "(" << quad.indices[0] << "," << quad.indices[1] << "," << quad.indices[2] << "," << quad.indices[3] << ")";
}

SCIP_RETCODE solveMonomials(int N, int R, std::vector<int>& optimalSequence, double timeLimit,
  int& numVariables, int& numFixedConstraints, int& numLazyConstraints, long& numNodes, double& finalDualBound,
  double& firstLPDualboundRoot, double &time)
{
  std::unordered_map<Quad, int> polynomial;

  /* Create the polynomial. */
  for (int i = 0; i <= N - R; ++i)
  {
    for (int d = 1; d <= R-1; ++d)
    {
      /* j1 and j2 go through the inner sum of C(d,i,R) to build C^2(d,i,R). */
      for (int j1 = i; j1 <= i + R - 1 - d; ++j1)
      {
        for (int j2 = i; j2 <= i + R - 1 - d; ++j2)
        {
          /* We add (2x_{j1}-x_0) * (2x_{j1+d}-x_0) * (2x_{j2}-x_0) * (2x_{j2+d}-x_0) with x_0 = 1 */
          for (int choice = 0; choice < 16; ++choice)
          {
            int coefficient = 1;
            int indices[4];
            if (choice & 1)
            {
              indices[0] = j1;
              coefficient *= 2;
            }
            else
            {
              indices[0] = -1;
              coefficient *= -1;
            }
            if (choice & 2)
            {
              indices[1] = j1+d;
              coefficient *= 2;
            }
            else
            {
              indices[1] = -1;
              coefficient *= -1;
            }
            if (choice & 4)
            {
              indices[2] = j2;
              coefficient *= 2;
            }
            else
            {
              indices[2] = -1;
              coefficient *= -1;
            }
            if (choice & 8)
            {
              indices[3] = j2+d;
              coefficient *= 2;
            }
            else
            {
              indices[3] = -1;
              coefficient *= -1;
            }
            Quad quad(indices);
            auto iter = polynomial.find(quad);
            if (iter == polynomial.end())
              polynomial.insert(std::make_pair(quad, coefficient));
            else
              polynomial[quad] += coefficient;
          }
        }
      }
      std::vector<std::pair<int,int>> terms;
      for (int j = i; j <= i + R - 1 - d; ++j)
        terms.push_back(std::make_pair(j, j+d));
    }
  }

  std::vector<Quad> zeroQuads;
  for (auto iter : polynomial)
  {
    if (iter.second == 0)
      zeroQuads.push_back(iter.first);
  }
  for (auto quad : zeroQuads)
    polynomial.erase(quad);

  std::cout << "After deleting " << zeroQuads.size() << ", the polynomial has " << polynomial.size() << " terms."
    << std::endl;

  /* Create the model. */
  SCIP* scip = NULL;
  SCIP_CALL( SCIPcreate(&scip) );
  SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
  SCIP_CALL( SCIPcreateProbBasic(scip, "LABS mono") );
  SCIP_CALL( SCIPsetRealParam(scip, "limits/time", timeLimit) );

  std::vector<SCIP_VAR*> xVars(N, NULL);
  for (int i = 0; i < N; ++i)
  {
    char name[256];
    snprintf(name, 256, "x_%d", i);
    SCIP_CALL( SCIPcreateVarBasic(scip, &xVars[i], name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
    SCIP_CALL( SCIPaddVar(scip, xVars[i]) );
  }

  std::unordered_map<Quad, SCIP_VAR*> monomialVars;
  for (auto iter : polynomial)
  {
    const Quad& quad = iter.first;
    SCIP_VAR* var = NULL;
    SCIP_CONS* cons = NULL;
    std::stringstream stream;
    int degree = 0;
    for (int i = 0; i < 4; ++i)
    {
      if (quad.indices[i] >= 0)
      {
        stream << "_" << quad.indices[i];
        ++degree;
      }
    }
    std::string name = "z" + stream.str();
    SCIP_CALL( SCIPcreateVarBasic(scip, &var, name.c_str(), 0.0, 1.0, iter.second, SCIP_VARTYPE_BINARY) );
    SCIP_CALL( SCIPaddVar(scip, var) );
    monomialVars[iter.first] = var;

    /* Add long constraint: z_product + (1-x_1) + (1-x_2) + ... >= 1 */
    name = "long" + stream.str();
    SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name.c_str(), 0, 0, 0, 1.0 - degree, 1.0) );
    SCIP_CALL( SCIPaddCoefLinear(scip, cons, var, 1.0) );
    for (int i = 0; i < 4; ++i)
    {
      if (quad.indices[i] >= 0)
        SCIP_CALL( SCIPaddCoefLinear(scip, cons, xVars[quad.indices[i]], -1.0) );
    }
    SCIP_CALL( SCIPaddCons(scip, cons) );
    SCIP_CALL( SCIPreleaseCons(scip, &cons) );

    /* Add short constraints: z_product <= x_i:  0 <= -z_product + x_i <= 1 */
    name = stream.str();
    for (int i = 0; i < 4; ++i)
    {
      if (quad.indices[i] >= 0)
      {
        stream.clear();
        stream << "short_" << name << "_" << quad.indices[i];
        SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, stream.str().c_str(), 0, 0, 0, 0.0, 1.0) );
        SCIP_CALL( SCIPaddCoefLinear(scip, cons, var, -1.0) );
        SCIP_CALL( SCIPaddCoefLinear(scip, cons, xVars[quad.indices[i]], +1.0) );
        SCIP_CALL( SCIPaddCons(scip, cons) );
        SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
    }

    SCIP_CALL( SCIPreleaseVar(scip, &var) );
  }

  optimalSequence.clear();
  numVariables = SCIPgetNOrigVars(scip);
  numFixedConstraints = SCIPgetNOrigConss(scip);
  numLazyConstraints = 0;

  SCIP_CALL( SCIPsolve(scip) );

  time = SCIPgetTotalTime(scip);
  numNodes = SCIPgetNNodes(scip);
  finalDualBound = SCIPgetDualbound(scip);
  firstLPDualboundRoot = SCIPgetFirstLPDualboundRoot(scip);
  SCIP_SOL* bestSol = SCIPgetBestSol(scip);

  optimalSequence.resize(N);
  for (int i = 0; i < N; ++i)
    optimalSequence[i] = (SCIPgetSolVal(scip, bestSol, xVars[i]) > 0.5) ? -1 : +1;

  for (int i = 0; i < N; ++i)
    SCIP_CALL( SCIPreleaseVar(scip, &xVars[i]) );

  SCIP_CALL( SCIPfreeProb(scip) );
  SCIP_CALL( SCIPfree(&scip) );

  return SCIP_OKAY;
}

SCIP_RETCODE solveValue(int N, int R, std::vector<int>& optimalSequence, bool modelProduct, double timeLimit,
  int& numVariables, int& numFixedConstraints, int& numLazyConstraints, long& numNodes, double& finalDualBound,
  double& firstLPDualboundRoot, double &time)
{
  /* Create the model. */
  SCIP* scip = NULL;
  SCIP_CALL( SCIPcreate(&scip) );
  SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
  SCIP_CALL( SCIPcreateProbBasic(scip, modelProduct ? "LABS value-prod" : "LABS value-same") );
  SCIP_CALL( SCIPsetRealParam(scip, "limits/time", timeLimit) );

  /* x-variables: x_i = 0  <=> s_i = -1 */

  std::vector<SCIP_VAR*> xVars(N, NULL);
  for (int i = 0; i < N; ++i)
  {
    char name[256];
    snprintf(name, 256, "x_%d", i);
    SCIP_CALL( SCIPcreateVarBasic(scip, &xVars[i], name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
    SCIP_CALL( SCIPaddVar(scip, xVars[i]) );
  }

  /* y-variables:
   * modelProduct: y_ij = x_i * x_j
   * modelProduct: y_ij = 1 <=>  x_i == x_j
   */

  std::vector<std::vector<SCIP_VAR*>> yVars(N);
  for (int i = 0; i < N; ++i)
  {
    yVars[i].resize(N);
    for (int j = i + 1; j < N; ++j)
    {
      char name[256];
      snprintf(name, 256, "y_%d_%d", i, j);
      SCIP_CALL( SCIPcreateVarBasic(scip, &yVars[i][j], name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
      SCIP_CALL( SCIPaddVar(scip, yVars[i][j]) );
    }
  }

  /* y_ii = x_i and y_ji = y_ij */

  for (int i = 0; i < N; ++i)
  {
    yVars[i][i] = xVars[i];
    SCIP_CALL( SCIPcaptureVar(scip, yVars[i][i]) );
    for (int j = 0; j < i; ++j)
    {
      yVars[i][j] = yVars[j][i];
      SCIP_CALL( SCIPcaptureVar(scip, yVars[i][j]) );
    }
  }

  /* z-variables: z_i_d_l = 1 if and only if C(d,i,R) = l */

  std::vector<std::vector<std::unordered_map<int, SCIP_VAR*>>> zVars(N-R+1);
  for (int i = 0; i <= N-R; ++i)
  {
    zVars[i].resize(R);
    for (int d = 1; d <= R-1; ++d)
    {
      int maxValue = R - d;
      for (int l = -maxValue; l <= maxValue; l += 2)
      {
        char name[256];
        snprintf(name, 256, "z_%d_%d_%d", i, d, l);
        SCIP_VAR* var = NULL;
        SCIP_CALL( SCIPcreateVarBasic(scip, &var, name, 0.0, 1.0, l * l, SCIP_VARTYPE_BINARY) );
        SCIP_CALL( SCIPaddVar(scip, var) );
        zVars[i][d][l] = var;
      }
    }
  }

  /* Quadratic constraints. */
  for (int i = 0; i < N; ++i)
  {
    for (int j = i + 1; j < N; ++j)
    {
      SCIP_CONS* cons = NULL;
      char name[SCIP_MAXSTRLEN];

      if (modelProduct)
      {
        /* y_ij - x_i <= 0 */
        SCIPsnprintf(name, SCIP_MAXSTRLEN, "short_%d_%d_%d", i, j, i);
        SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, 0, 0, -SCIPinfinity(scip), 0.0) );
        SCIP_CALL( SCIPaddCoefLinear(scip, cons, yVars[i][j], +1.0) );
        SCIP_CALL( SCIPaddCoefLinear(scip, cons, xVars[i], -1.0) );
        SCIP_CALL( SCIPaddCons(scip, cons) );
        SCIP_CALL( SCIPreleaseCons(scip, &cons) );

        /* y_ij - x_j <= 0 */
        SCIPsnprintf(name, SCIP_MAXSTRLEN, "short_%d_%d_%d", i, j, j);
        SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, 0, 0, -SCIPinfinity(scip), 0.0) );
        SCIP_CALL( SCIPaddCoefLinear(scip, cons, yVars[i][j], +1.0) );
        SCIP_CALL( SCIPaddCoefLinear(scip, cons, xVars[j], -1.0) );
        SCIP_CALL( SCIPaddCons(scip, cons) );
        SCIP_CALL( SCIPreleaseCons(scip, &cons) );

        /* 0 <= -y_ij + x_i + x_i <= 1 */
        SCIPsnprintf(name, SCIP_MAXSTRLEN, "short_%d_%d_%d", i, j, j);
        SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, 0, 0, -SCIPinfinity(scip), 1.0) );
        SCIP_CALL( SCIPaddCoefLinear(scip, cons, yVars[i][j], -1.0) );
        SCIP_CALL( SCIPaddCoefLinear(scip, cons, xVars[i], +1.0) );
        SCIP_CALL( SCIPaddCoefLinear(scip, cons, xVars[j], +1.0) );
        SCIP_CALL( SCIPaddCons(scip, cons) );
        SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
      else
      {
        /* x_i + x_j + y_ij >= 1 */
        SCIPsnprintf(name, SCIP_MAXSTRLEN, "same1_%d_%d", i, j);
        SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, 0, 0, 1.0, SCIPinfinity(scip)) );
        SCIP_CALL( SCIPaddCoefLinear(scip, cons, xVars[i], +1.0) );
        SCIP_CALL( SCIPaddCoefLinear(scip, cons, xVars[j], +1.0) );
        SCIP_CALL( SCIPaddCoefLinear(scip, cons, yVars[i][j], +1.0) );
        SCIP_CALL( SCIPaddCons(scip, cons) );
        SCIP_CALL( SCIPreleaseCons(scip, &cons) );

        /* x_i - x_j - y_ij >= -1 */
        SCIPsnprintf(name, SCIP_MAXSTRLEN, "same2_%d_%d", i, j);
        SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, 0, 0, -1.0, SCIPinfinity(scip)) );
        SCIP_CALL( SCIPaddCoefLinear(scip, cons, xVars[i], +1.0) );
        SCIP_CALL( SCIPaddCoefLinear(scip, cons, xVars[j], -1.0) );
        SCIP_CALL( SCIPaddCoefLinear(scip, cons, yVars[i][j], -1.0) );
        SCIP_CALL( SCIPaddCons(scip, cons) );
        SCIP_CALL( SCIPreleaseCons(scip, &cons) );

        /* -x_i + x_j - y_ij >= -1 */
        SCIPsnprintf(name, SCIP_MAXSTRLEN, "same3_%d_%d", i, j);
        SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, 0, 0, -1.0, SCIPinfinity(scip)) );
        SCIP_CALL( SCIPaddCoefLinear(scip, cons, xVars[i], -1.0) );
        SCIP_CALL( SCIPaddCoefLinear(scip, cons, xVars[j], +1.0) );
        SCIP_CALL( SCIPaddCoefLinear(scip, cons, yVars[i][j], -1.0) );
        SCIP_CALL( SCIPaddCons(scip, cons) );
        SCIP_CALL( SCIPreleaseCons(scip, &cons) );

        /* x_i + x_j - y_ij <= 1 */
        SCIPsnprintf(name, SCIP_MAXSTRLEN, "same4_%d_%d", i, j);
        SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, 0, 0, -SCIPinfinity(scip), 1.0) );
        SCIP_CALL( SCIPaddCoefLinear(scip, cons, xVars[i], +1.0) );
        SCIP_CALL( SCIPaddCoefLinear(scip, cons, xVars[j], +1.0) );
        SCIP_CALL( SCIPaddCoefLinear(scip, cons, yVars[i][j], -1.0) );
        SCIP_CALL( SCIPaddCons(scip, cons) );
        SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
    }
  }

  /* Cardinality constraints. */
  for (int i = 0; i <= N-R; ++i)
  {
    for (int d = 1; d <= R-1; ++d)
    {
      SCIP_CONS* cons = NULL;
      char name[SCIP_MAXSTRLEN];
      SCIPsnprintf(name, SCIP_MAXSTRLEN, "card_%d_%d", i, d);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, 0, 0, 1.0, 1.0) );
      for (auto iter : zVars[i][d])
        SCIP_CALL( SCIPaddCoefLinear(scip, cons, iter.second, 1.0) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
    }
  }

  /* Selection constraints:
   *
   * modelProduct: \sum_j (4y_{j,j+d} - 2x_j - 2x_{j+d} + 1) - \sum_l l*z_{i,d,l} == 0
   * modelProduct: \sum_j (4y_{j,j+d} - 2x_j - 2x_{j+d} + 1) - \sum_l l*z_{i,d,l} == 0
   */
  for (int i = 0; i <= N-R; ++i)
  {
    for (int d = 1; d <= R-1; ++d)
    {
      SCIP_CONS* cons = NULL;
      char name[SCIP_MAXSTRLEN];
      SCIPsnprintf(name, SCIP_MAXSTRLEN, "select_%d_%d", i, d);

      if (modelProduct)
      {
        SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, 0, 0, -(R-d), -(R-d)) );
        for (int j = i; j <= i + R - 1 - d; ++j)
        {
          SCIP_CALL( SCIPaddCoefLinear(scip, cons, yVars[j][j+d], 4.0) );
          SCIP_CALL( SCIPaddCoefLinear(scip, cons, xVars[j], -2.0) );
          SCIP_CALL( SCIPaddCoefLinear(scip, cons, xVars[j+d], -2.0) );
          /* constant +1 is taken into account by rhs R-d. */
        }
      }
      else
      {
        SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, 0, 0, R-d, R-d) );
        for (int j = i; j <= i + R - 1 - d; ++j)
        {
          SCIP_CALL( SCIPaddCoefLinear(scip, cons, yVars[j][j+d], 2.0) );
          /* constant -1 is taken into account by rhs R-d. */
        }
      }
      for (auto iter : zVars[i][d])
        SCIP_CALL( SCIPaddCoefLinear(scip, cons, iter.second, -iter.first) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
    }
  }

  optimalSequence.clear();
  numVariables = SCIPgetNOrigVars(scip);
  numFixedConstraints = SCIPgetNOrigConss(scip);
  numLazyConstraints = 0;

  SCIP_CALL( SCIPsolve(scip) );

  time = SCIPgetTotalTime(scip);
  numNodes = SCIPgetNNodes(scip);
  finalDualBound = SCIPgetDualbound(scip);
  firstLPDualboundRoot = SCIPgetFirstLPDualboundRoot(scip);
  SCIP_SOL* bestSol = SCIPgetBestSol(scip);

  optimalSequence.resize(N);
  for (int i = 0; i < N; ++i)
    optimalSequence[i] = (SCIPgetSolVal(scip, bestSol, xVars[i]) > 0.5) ? -1 : +1;

  /* Release variables. */
  for (int i = 0; i < N; ++i)
  {
    SCIP_CALL( SCIPreleaseVar(scip, &xVars[i]) );
    for (int j = 0; j < N; ++j)
      SCIP_CALL( SCIPreleaseVar(scip, &yVars[i][j]) );
  }
  for (int i = 0; i <= N-R; ++i)
  {
    for (int d = 1; d <= R-1; ++d)
    {
      for (auto iter : zVars[i][d])
        SCIP_CALL( SCIPreleaseVar(scip, &iter.second) );
    }
  }

  SCIP_CALL( SCIPfreeProb(scip) );
  SCIP_CALL( SCIPfree(&scip) );

  return SCIP_OKAY;
}

SCIP_RETCODE solveNogood(int N, int R, std::vector<int>& optimalSequence, bool cardinalityEquation, double timeLimit,
  int& numVariables, int& numFixedConstraints, int& numLazyConstraints, long& numNodes, double& finalDualBound,
  double& firstLPDualboundRoot, double &time)
{
  /* Create the model. */
  SCIP* scip = NULL;
  SCIP_CALL( SCIPcreate(&scip) );
  SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
  SCIP_CALL( SCIPincludeObjConshdlr(scip, new ConshdlrLABSnogood(scip), TRUE) );
  SCIP_CALL( SCIPcreateProbBasic(scip, cardinalityEquation ? "LABS nogood+card" : "LABS nogood") );
  SCIP_CALL( SCIPsetRealParam(scip, "limits/time", timeLimit) );

  /* x-variables: x_i = 0  <=> s_i = -1 */

  std::vector<SCIP_VAR*> xVars(N, NULL);
  for (int i = 0; i < N; ++i)
  {
    char name[256];
    snprintf(name, 256, "x_%d", i);
    SCIP_CALL( SCIPcreateVarBasic(scip, &xVars[i], name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
    SCIP_CALL( SCIPaddVar(scip, xVars[i]) );
  }

  /* z-variables: z_i_d_l = 1 if and only if C(d,i,R) = l */

  std::vector<std::vector<std::unordered_map<int, SCIP_VAR*>>> zVars(N-R+1);
  for (int i = 0; i <= N-R; ++i)
  {
    zVars[i].resize(R);
    for (int d = 1; d <= R-1; ++d)
    {
      int maxValue = R - d;
      for (int l = -maxValue; l <= maxValue; l += 2)
      {
        char name[256];
        snprintf(name, 256, "z_%d_%d_%d", i, d, l);
        SCIP_VAR* var = NULL;
        SCIP_CALL( SCIPcreateVarBasic(scip, &var, name, 0.0, 1.0, l * l, SCIP_VARTYPE_BINARY) );
        SCIP_CALL( SCIPaddVar(scip, var) );
        zVars[i][d][l] = var;
      }
    }
  }

  /* Cardinality constraints. */
  numFixedConstraints = 0;
  if (cardinalityEquation)
  {
    for (int i = 0; i <= N-R; ++i)
    {
      for (int d = 1; d <= R-1; ++d)
      {
        SCIP_CONS* cons = NULL;
        char name[SCIP_MAXSTRLEN];
        SCIPsnprintf(name, SCIP_MAXSTRLEN, "card_%d_%d", i, d);
        SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, 0, 0, 1.0, 1.0) );
        for (auto iter : zVars[i][d])
          SCIP_CALL( SCIPaddCoefLinear(scip, cons, iter.second, 1.0) );
        SCIP_CALL( SCIPaddCons(scip, cons) );
        SCIP_CALL( SCIPreleaseCons(scip, &cons) );
        numFixedConstraints++;
      }
    }
  }

  /* Add lazy constraints. */
  for (int i = 0; i <= N-R; ++i)
  {
    for (int d = 1; d <= R-1; ++d)
    {
      SCIP_CONS* cons = NULL;
      char name[SCIP_MAXSTRLEN];
      SCIPsnprintf(name, SCIP_MAXSTRLEN, "labs-nogood-%d-%d", i, d);
      SCIP_CALL( SCIPcreateConsLABSnogood(scip, &cons, "LABSnogood", R, d, &xVars[i], zVars[i][d], TRUE, FALSE, TRUE,
        TRUE, FALSE, FALSE, FALSE, TRUE, TRUE) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
    }
  }

  optimalSequence.clear();
  numVariables = SCIPgetNOrigVars(scip);
  numLazyConstraints = 0;

  SCIP_CALL( SCIPsolve(scip) );

  time = SCIPgetTotalTime(scip);
  numNodes = SCIPgetNNodes(scip);
  finalDualBound = SCIPgetDualbound(scip);
  firstLPDualboundRoot = SCIPgetFirstLPDualboundRoot(scip);
  numLazyConstraints = SCIPconshdlrGetNCutsApplied(SCIPfindConshdlr(scip, "LABSnogood"));
  SCIP_SOL* bestSol = SCIPgetBestSol(scip);

  optimalSequence.resize(N);
  for (int i = 0; i < N; ++i)
    optimalSequence[i] = (SCIPgetSolVal(scip, bestSol, xVars[i]) > 0.5) ? -1 : +1;

  /* Release variables. */
  for (int i = 0; i < N; ++i)
    SCIP_CALL( SCIPreleaseVar(scip, &xVars[i]) );
  for (int i = 0; i <= N-R; ++i)
  {
    for (int d = 1; d <= R-1; ++d)
    {
      for (auto iter : zVars[i][d])
        SCIP_CALL( SCIPreleaseVar(scip, &iter.second) );
    }
  }

  SCIP_CALL( SCIPfreeProb(scip) );
  SCIP_CALL( SCIPfree(&scip) );

  return SCIP_OKAY;
}

void printUsage(const char* binaryName)
{
  std::cerr << "Usage: " << binaryName << " N R METHOD TIME-LIMIT\n";
  std::cerr << "Parameters:\n";
  std::cerr << "  N          Length of sequence.\n";
  std::cerr << "  R          Interaction range.\n";
  std::cerr << "  METHOD     MIP model to use among {mono,nogood,nogood+card,value-prod,value-same}.\n";
  std::cerr << "  TIME-LIMIT Time limit\n";
  std::cerr << std::flush;
}

int main(int argc, char** argv)
{
  if (argc < 4)
  {
    printUsage(argv[0]);
    return EXIT_FAILURE;
  }

  int N, R;
  if (sscanf(argv[1], "%d", &N) != 1 || N <= 0)
  {
    std::cerr << "Error: invalid sequence length <" << argv[1] << ">.\n";
    printUsage(argv[0]);
    return EXIT_FAILURE;
  }
  if (sscanf(argv[2], "%d", &R) != 1 || R <= 0 || R > N)
  {
    std::cerr << "Error: invalid interaction range <" << argv[2] << ">.\n";
    printUsage(argv[0]);
    return EXIT_FAILURE;
  }
  double timeLimit = 1e20;
  if (argc > 4 && (sscanf(argv[4], "%lf", &timeLimit) != 1 || timeLimit <= 0))
  {
    std::cerr << "Error: invalid time limit <" << argv[4] << ">.\n";
    printUsage(argv[0]);
    return EXIT_FAILURE;
  }

  const std::string model = argv[3];
  SCIP_RETCODE retcode;
  std::vector<int> optimalSequence;
  int numVariables;
  int numFixedConstraints;
  int numLazyConstraints;
  double firstLPDualboundRoot;
  double finalDualBound;
  long numNodes;
  double time;
  if (model == "mono")
  {
    retcode = solveMonomials(N, R, optimalSequence, timeLimit, numVariables, numFixedConstraints, numLazyConstraints,
      numNodes, finalDualBound, firstLPDualboundRoot, time);
  }
  else if (model == "value-prod")
    retcode = solveValue(N, R, optimalSequence, true, timeLimit, numVariables, numFixedConstraints, numLazyConstraints,
      numNodes, finalDualBound, firstLPDualboundRoot, time);
  else if (model == "value-same")
    retcode = solveValue(N, R, optimalSequence, false, timeLimit, numVariables, numFixedConstraints, numLazyConstraints,
      numNodes, finalDualBound, firstLPDualboundRoot, time);
  else if (model == "nogood")
    retcode = solveNogood(N, R, optimalSequence, false, timeLimit, numVariables, numFixedConstraints, numLazyConstraints,
      numNodes, finalDualBound, firstLPDualboundRoot, time);
  else if (model == "nogood+card")
    retcode = solveNogood(N, R, optimalSequence, true, timeLimit, numVariables, numFixedConstraints, numLazyConstraints,
      numNodes, finalDualBound, firstLPDualboundRoot, time);
  else
  {
    std::cerr << "Error: invalid method <" << argv[3] << ">.\n";
    printUsage(argv[0]);
    return EXIT_FAILURE;
  }

  if (retcode != SCIP_OKAY)
  {
    SCIPprintError(retcode);
    return EXIT_FAILURE;
  }

  std::cout << "\nBest found sequence:";
  for (int i = 0; i < N; ++i)
  {
    std::cout << " " << (optimalSequence[i] > 0 ? "+" : "-");
  }
  std::cout << std::endl;
  int objective = 0;
  for (int i = 0; i <= N - R; ++i)
  {
    for (int d = 1; d <= R-1; ++d)
    {
      int C = 0;
      for (int j = i; j <= i + R - 1 - d; ++j)
        C += optimalSequence[j] * optimalSequence[j+d];
      std::cout << "Contribution from i=" << i << ", d=" << d << ": " << C*C << std::endl;
      objective += C*C;
    }
  }
  std::cout << "Sequence length: " << N << std::endl;
  std::cout << "Interaction range: " << N << std::endl;
  std::cout << "Model: " << model << std::endl;
  std::cout << "Primal bound: " << objective << std::endl;
  std::cout << "Dual bound: " << finalDualBound << std::endl;
  std::cout << "LP dual bound: " << firstLPDualboundRoot << std::endl;
  std::cout << "Time: " << time << std::endl;
  std::cout << "#vars: " << numVariables << std::endl;
  std::cout << "#fixed conss: " << numFixedConstraints << std::endl;
  std::cout << "#lazy conss: " << numLazyConstraints << std::endl;

  return EXIT_SUCCESS;
}
