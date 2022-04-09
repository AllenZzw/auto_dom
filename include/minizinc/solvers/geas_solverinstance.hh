/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */

/*
 *  Main authors:
 *     Jip J. Dekker <jip.dekker@monash.edu>
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include <minizinc/flattener.hh>
#include <minizinc/solver.hh>

#include <geas/solver/solver.h>

namespace MiniZinc {
enum Monotonicity { VAR_NONE, VAR_INC, VAR_DEC, VAR_EQL };

typedef void (*analyzer)(SolverInstanceBase&, Call* call);
class Preprocessor {
protected:
  ManagedASTStringMap<analyzer> _preprocessor;
  SolverInstanceBase& _base;

public:
  Preprocessor(SolverInstanceBase& base) : _base(base) {}
  void add(ASTString name, analyzer p);
  void add(const std::string& name, analyzer p);
  void post(Call* c);
  void cleanup() { _preprocessor.clear(); }
};

typedef void (*dominator)(SolverInstanceBase&, const Call* call, const std::unordered_set<Id*>& fixedVars);
class DomPoster {
protected:
  ManagedASTStringMap<dominator> _domposter;
  SolverInstanceBase& _base;

public:
  DomPoster(SolverInstanceBase& base) : _base(base) {}
  void add(ASTString name, dominator p);
  void add(const std::string& name, dominator p);
  void post(Call* c, const std::unordered_set<Id*>& fixedVars);
  void cleanup() { _domposter.clear(); }
};

class GeasOptions : public SolverInstanceBase::Options {
public:
  bool allSolutions = false;
  int conflicts = 0;
  bool freeSearch = false;
  int nrSolutions = 1;
  int objProbeLimit = 0;
  bool statistics = false;
  std::chrono::milliseconds time = std::chrono::milliseconds(0);
  int minNogoodLength = 2; 
  int maxNogoodLength = 2; 
};

class GeasVariable {
public:
  enum Type { BOOL_TYPE, FLOAT_TYPE, INT_TYPE };

protected:
  Type _t;  // Type of the variable
  union {
    geas::patom_t bv;
    geas::fp::fpvar fv;
    geas::intvar iv;
  };

public:
  explicit GeasVariable(const geas::patom_t& bv0) : _t(BOOL_TYPE), bv(bv0){};
  explicit GeasVariable(const geas::fp::fpvar& fv0) : _t(FLOAT_TYPE), fv(fv0){};
  explicit GeasVariable(const geas::intvar& iv0) : _t(INT_TYPE), iv(iv0){};

  GeasVariable(const GeasVariable& gv) : _t(gv._t) {
    switch (_t) {
      case BOOL_TYPE:
        bv = gv.bv;
        break;
      case FLOAT_TYPE:
        fv = gv.fv;
        break;
      case INT_TYPE:
        iv = gv.iv;
        break;
    }
  }

  bool isBool() const { return _t == BOOL_TYPE; }
  bool isFloat() const { return _t == FLOAT_TYPE; }
  bool isInt() const { return _t == INT_TYPE; }

  geas::patom_t boolVar() const { return bv; }
  geas::fp::fpvar floatVar() const { return fv; }
  geas::intvar intVar() const { return iv; }
};

class GeasTypes {
public:
  typedef GeasVariable Variable;
  typedef MiniZinc::Statistics Statistics;
};

class GeasSolverInstance : public SolverInstanceImpl<GeasTypes> {
public:
  GeasSolverInstance(Env& env, std::ostream& log, SolverInstanceBase::Options* opt);
  ~GeasSolverInstance() override = default;
  void processFlatZinc() override;
  geas::solver_data* solverData() const { return _domsolver->data; }
  geas::solver& solver() { return *_domsolver; }

  Status solve() override;

  Status next() override { return SolverInstance::ERROR; }  // TODO: Implement
  void resetSolver() override; // assert(false)
  Expression* getSolutionValue(Id* id) override; // assert(false)
  void printStatistics() override; // assert(false)

  // MiniZinc to value conversions for posting constraints 
  bool asBool(Expression* e) { return eval_bool(env().envi(), e); }
  vec<bool> asBool(ArrayLit* al);
  int asInt(Expression* e) { return static_cast<int>(eval_int(env().envi(), e).toInt()); }
  vec<int> asInt(ArrayLit* al);

  /// information from the analysis of FlatZinc model 
  std::vector<Id*> _indepVars; // vector of indepent variables 
  std::map<Id*, Expression*> _alias; // map variable id to its alias (IntLit or BoolLit or another variable)
  std::map<Id*, std::set<Call*> > _argsToConstraint; // map variable id to the list of constraint that has the variable as arguments
  std::map<Call*, int> _constraintArgNr; // map experssion pointer to the number of unique variables that it depends on
  std::map<Id*, std::string> _variableName; // map variable id to variable name in the MiniZinc Model 

  std::map<Id*, Monotonicity> _variableMono; // map variable id to its monotonicity 
  std::map<Id*, bool> _variableSensitive; // map variable id to its sensitivity  

  // MiniZinc to Id pointers conversion 
  geas::intvar asIntVarDom(Expression* e, bool dominating);
  vec<geas::intvar> asIntVarDom(ArrayLit* al, bool dominating);
  geas::patom_t asBoolVarDom(Expression* e, bool dominating);
  vec<geas::patom_t> asBoolVarDom(ArrayLit* al, bool dominating);

  void createVar(geas::solver& s, VarDecl* vd); 

  // MiniZinc to value conversions 
  Id* resolveVarId(Expression* e);
  Id* asVarId(Expression* e);
  std::vector<Id*> asVarId(ArrayLit* al);

  // utilities 
  void replaceCallAlias(Call* c); 
  void replaceArrayAlias(ArrayLit* c); 
  std::set<Id*> getArguments(Call* c);

  IdMap<VarId> _variableMap0, _variableMap1; // map variable id to variable instance 
  std::vector<geas::intvar> _sensvar; // sensitive variable which makes the objective strictly better
  Id* _objVar; 
  // TODO: create only when necessary or use Geas internal
  geas::intvar zero;
protected:
  geas::solver _solver;
  Model* _flat;

  Preprocessor _constraintPreprocessor;
  DomPoster _domPoster; 

  geas::solver * _domsolver;
  SolveI::SolveType _objType = SolveI::ST_SAT;

  GeasTypes::Variable& resolveVarDom(Expression* e, bool dominating);

  void registerConstraint(const std::string& name, analyzer a, poster p, dominator d);
  void registerConstraints();
};

class GeasSolverFactory : public SolverFactory {
public:
  GeasSolverFactory();
  SolverInstanceBase::Options* createOptions() override;
  SolverInstanceBase* doCreateSI(Env& env, std::ostream& log,
                                 SolverInstanceBase::Options* opt) override;

  std::string getDescription(SolverInstanceBase::Options* opt) override {
    return "Elsie Geas - Another Lazy Clause Generation Solver";
  };
  std::string getVersion(SolverInstanceBase::Options* opt) override { return "0.0.1"; }
  std::string getId() override { return "org.minizinc.geas"; }

  bool processOption(SolverInstanceBase::Options* opt, int& i, std::vector<std::string>& argv,
                    const std::string& workingDir = std::string()) override;
  void printHelp(std::ostream& os) override;
};

}  // namespace MiniZinc
