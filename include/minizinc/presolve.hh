/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */

/*
 *  Main authors:
 *     Jip J. Dekker <jip@dekker.li>
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef __MINIZINC_PRESOLVE_HH__
#define __MINIZINC_PRESOLVE_HH__

#include <minizinc/model.hh>
#include <minizinc/copy.hh>
#include <minizinc/solver_instance_base.hh>

namespace MiniZinc {
  class Presolver {
  public:
    struct Options;
    class Solns2Vector;
    class TableBuilder;
  protected:
    class Subproblem;
    class GlobalSubproblem;
    class ModelSubproblem;
    class CallsSubproblem;

    Env& env;
    Model* model;
    std::vector<Subproblem*> subproblems;

    Options& options;

    void findPresolveAnnotations();

    void findPresolvedCalls();

  public:
    Presolver(Env& env, Model* m, Options& options)
            : env(env), model(m), options(options) {}

    virtual ~Presolver();

    void presolve();
  };

//  TODO: Group Flattener options and use those, to avoid copying.
  struct Presolver::Options {
    std::vector<std::string> includePaths;
    std::string stdLibDir;
    std::string globalsDir;
    std::string modelOutput;

    bool verbose = false;
    bool newfzn = false;
    bool optimize = true;
    bool onlyRangeDomains = false;
  };


  class Presolver::Subproblem{
  protected:
    Model* origin;
    EnvI& origin_env;
    // predicate around which the submodel revolves
    FunctionI* predicate;
    // calls to the predicate in the model
    std::vector<Call*> calls;
    // Flattener for compiling flag access.
    Options& options;
    // FZN solver to be used.
    std::string solver;

    // Constructed Model & Env to solve subproblem.
    Model* m = nullptr;
    Env* e = nullptr;
    // Solver used to solve constructed model.
    SolverInstanceBase* si = nullptr;
    Solns2Vector* solns = nullptr;

    enum Constraint { BoolTable, IntTable, Element };

  public:

    Subproblem(Model* origin, EnvI& origin_env, FunctionI* predicate, Options& options, std::string solver="");
    virtual ~Subproblem();

    void addCall(Call* c) { calls.push_back(c); }

    FunctionI* getPredicate() const { return predicate; }

    virtual void solve();

  protected:
    virtual void constructModel() = 0;

    virtual void solveModel();

    virtual void replaceUsage() = 0;
  };

  class Presolver::GlobalSubproblem : public Presolver::Subproblem {
  public:
    GlobalSubproblem(Model* origin, EnvI& origin_env, FunctionI* predicate, Options& options, std::string solver="") : Subproblem(
            origin, origin_env, predicate, options, solver) { }

  protected:
    virtual void constructModel();

    virtual void replaceUsage();
  };

  class Presolver::ModelSubproblem : public Presolver::GlobalSubproblem {
  public:
    ModelSubproblem(Model* origin, EnvI& origin_env, FunctionI* predicate, Options& options, std::string solver="") : GlobalSubproblem(
            origin, origin_env, predicate, options, solver) { }

  protected:
    virtual void constructModel();
  };

  class Presolver::CallsSubproblem : public Presolver::Subproblem {
  protected:
    Call* currentCall = nullptr;

  public:
    CallsSubproblem(Model* origin, EnvI& origin_env, FunctionI* predicate, Options& options, std::string solver="") : Subproblem(
            origin, origin_env, predicate, options, solver) { }

    virtual void solve();

  protected:
    virtual void constructModel();

    virtual void replaceUsage();
  };

}

#endif
