/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */

/*
 *  Main authors:
 *     Jip J. Dekker <jip.dekker@monash.edu>
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-pro-type-static-cast-downcast"

#include <minizinc/solvers/geas/geas_constraints.hh>
#include <minizinc/solvers/geas_solverinstance.hh>

#include <geas/constraints/builtins.h>
#include <geas/constraints/flow/flow.h>

namespace MiniZinc {
namespace GeasConstraints {

#define SI (static_cast<GeasSolverInstance&>(s))
#define SD SI.solverData()
#define SOL SI.solver()
#define EXPR(X) call->arg(X)
#define BOOL(X) SI.asBool(EXPR(X))
#define BOOLARRAY(X) SI.asBool(ARRAY(X))
#define BOOLVAR0(X) SI.asBoolVarDom(EXPR(X), true)
#define BOOLVAR1(X) SI.asBoolVarDom(EXPR(X), false)
#define INT(X) SI.asInt(EXPR(X))
#define INTARRAY(X) SI.asInt(ARRAY(X))
#define PAR(X) call->arg(X)->type().isPar()
#define ARRAY(X) eval_array_lit(s.env().envi(), call->arg(X))
#define VARIDARRAY(X) SI.asVarId(ARRAY(X))
#define VARID(X) SI.asVarId(EXPR(X))
#define INTVAR0(X) SI.asIntVarDom(EXPR(X), true)
#define INTVAR1(X) SI.asIntVarDom(EXPR(X), false)
#define FUNCTIONNOTIMPLEMENT if ( call->ann().containsCall(std::string("defines_var")) ) throw InternalError(std::string("Constraint defining a variable not implemented: ") + std::string(call->id().c_str()) )
#define CONSTRAINTNOTIMPLEMENT if ( !call->ann().containsCall(std::string("defines_var")) ) throw InternalError(std::string("Constraint without defining a variable not implemented: ") + std::string(call->id().c_str()) )
#define CONSTANT(X) SOL.new_intvar( static_cast<geas::intvar::val_t>(X), static_cast<geas::intvar::val_t>(X) )
#define CONTEXT SI._variableMono[SI.asVarId(call->ann().getCall(std::string("defines_var"))->arg(0))]
#define SENSITIVE SI._variableSensitive[SI.asVarId(call->ann().getCall(std::string("defines_var"))->arg(0))]

void p_int_eq(SolverInstanceBase& s, const Call* call) { 
  FUNCTIONNOTIMPLEMENT; 
  geas::int_eq(SD, INTVAR0(0), INTVAR0(1)); 
  geas::int_eq(SD, INTVAR1(0), INTVAR1(1)); 
}

void p_int_ne(SolverInstanceBase& s, const Call* call) { 
  FUNCTIONNOTIMPLEMENT; 
  geas::int_ne(SD, INTVAR0(0), INTVAR0(1)); 
  geas::int_ne(SD, INTVAR1(0), INTVAR1(1)); 
}

void p_int_le(SolverInstanceBase& s, const Call* call) {
  FUNCTIONNOTIMPLEMENT; 
  geas::int_le(SD, INTVAR0(0), INTVAR0(1), 0);
  geas::int_le(SD, INTVAR1(0), INTVAR1(1), 0);
}

void p_int_lt(SolverInstanceBase& s, const Call* call) {
  FUNCTIONNOTIMPLEMENT; 
  geas::int_le(SD, INTVAR0(0), INTVAR0(1), -1);
  geas::int_le(SD, INTVAR1(0), INTVAR1(1), -1);
}

void p_int_eq_reif(SolverInstanceBase& s, const Call* call) {
  if (PAR(2)) {
    if (BOOL(2)) {
      p_int_eq(s, call);
    } else {
      p_int_ne(s, call);
    }
  } else {
    geas::int_eq(SD, INTVAR0(0), INTVAR0(1), BOOLVAR0(2));
    geas::int_ne(SD, INTVAR0(0), INTVAR0(1), ~BOOLVAR0(2));
    geas::int_eq(SD, INTVAR1(0), INTVAR1(1), BOOLVAR1(2));
    geas::int_ne(SD, INTVAR1(0), INTVAR1(1), ~BOOLVAR1(2));
  }
}

void p_int_ne_reif(SolverInstanceBase& s, const Call* call) {
  if (PAR(2)) {
    if (BOOL(2)) {
      p_int_ne(s, call);
    } else {
      p_int_eq(s, call);
    }
  } else {
    geas::int_ne(SD, INTVAR0(0), INTVAR0(1), BOOLVAR0(2));
    geas::int_eq(SD, INTVAR0(0), INTVAR0(1), ~BOOLVAR0(2));
    geas::int_ne(SD, INTVAR1(0), INTVAR1(1), BOOLVAR1(2));
    geas::int_eq(SD, INTVAR1(0), INTVAR1(1), ~BOOLVAR1(2));
  }
}

void p_int_le_reif(SolverInstanceBase& s, const Call* call) {
  if (PAR(2)) {
    if (BOOL(2)) {
      p_int_le(s, call);
    } else {
      auto* nc = new Call(Location().introduce(), call->id(), {call->arg(1), call->arg(0)});
      p_int_lt(s, nc);
    }
  } else {
    geas::int_le(SD, INTVAR0(0), INTVAR0(1), 0, BOOLVAR0(2));
    geas::int_le(SD, INTVAR0(1), INTVAR0(0), -1, ~BOOLVAR0(2));
    geas::int_le(SD, INTVAR1(0), INTVAR1(1), 0, BOOLVAR1(2));
    geas::int_le(SD, INTVAR1(1), INTVAR1(0), -1, ~BOOLVAR1(2));
  }
}

void p_int_lt_reif(SolverInstanceBase& s, const Call* call) {
  if (PAR(2)) {
    if (BOOL(2)) {
      p_int_lt(s, call);
    } else {
      auto* nc = new Call(Location().introduce(), call->id(), {call->arg(1), call->arg(0)});
      p_int_le(s, nc);
    }
  } else {
    geas::int_le(SD, INTVAR0(0), INTVAR0(1), -1, BOOLVAR0(2));
    geas::int_le(SD, INTVAR0(1), INTVAR0(0), 0, ~BOOLVAR0(2));
    geas::int_le(SD, INTVAR1(0), INTVAR1(1), -1, BOOLVAR1(2));
    geas::int_le(SD, INTVAR1(1), INTVAR1(0), 0, ~BOOLVAR1(2));
  }
}

void p_int_abs(SolverInstanceBase& s, const Call* call) { 
  // CONSTRAINTNOTIMPLEMENT; 
  geas::int_abs(SD, INTVAR0(1), INTVAR0(0)); 
  geas::int_abs(SD, INTVAR1(1), INTVAR1(0)); 
}

void p_int_times(SolverInstanceBase& s, const Call* call) {
  CONSTRAINTNOTIMPLEMENT; 
  geas::int_mul(SD, INTVAR0(2), INTVAR0(0), INTVAR0(1));
  geas::int_mul(SD, INTVAR1(2), INTVAR1(0), INTVAR1(1));
}

void p_int_div(SolverInstanceBase& s, const Call* call) {
  CONSTRAINTNOTIMPLEMENT; 
  geas::int_div(SD, INTVAR0(2), INTVAR0(0), INTVAR0(1));
  geas::int_div(SD, INTVAR1(2), INTVAR1(0), INTVAR1(1));
}

void p_int_max(SolverInstanceBase& s, const Call* call) {
  vec<geas::intvar> vars0 = {INTVAR0(0), INTVAR0(1)}; 
  vec<geas::intvar> vars1 = {INTVAR1(0), INTVAR1(1)}; 
  geas::int_max(SD, INTVAR0(2), vars0); 
  geas::int_max(SD, INTVAR1(2), vars1); 
}

void p_int_min(SolverInstanceBase& s, const Call* call) {
  vec<geas::intvar> vars0 = {-INTVAR0(0), -INTVAR0(1)}; 
  vec<geas::intvar> vars1 = {-INTVAR1(0), -INTVAR1(1)}; 
  geas::int_max(SD, -INTVAR0(2), vars0); 
  geas::int_max(SD, -INTVAR1(2), vars1); 
}

void p_int_lin_eq(SolverInstanceBase& s, const Call* call) {
  vec<int> pos = INTARRAY(0);
  vec<int> neg(pos.size());
  for (int i = 0; i < neg.size(); ++i)
    neg[i] = -pos[i];
  vec<geas::intvar> vars0 = SI.asIntVarDom(ARRAY(1), true);
  vec<geas::intvar> vars1 = SI.asIntVarDom(ARRAY(1), false);
  geas::linear_le(SD, pos, vars0, INT(2)); 
  geas::linear_le(SD, pos, vars1, INT(2)); 
  geas::linear_le(SD, neg, vars0, -INT(2));
  geas::linear_le(SD, neg, vars1, -INT(2));
}

void p_int_lin_ne(SolverInstanceBase& s, const Call* call) {
  vec<int> pos = INTARRAY(0);
  vec<geas::intvar> vars0 = SI.asIntVarDom(ARRAY(1), true);
  vec<geas::intvar> vars1 = SI.asIntVarDom(ARRAY(1), false);
  geas::linear_ne(SD, pos, vars0, INT(2)); 
  geas::linear_ne(SD, pos, vars1, INT(2)); 
}

void p_int_lin_le(SolverInstanceBase& s, const Call* call) {
  vec<int> cons = INTARRAY(0);
  vec<geas::intvar> vars0 = SI.asIntVarDom(ARRAY(1), true);
  vec<geas::intvar> vars1 = SI.asIntVarDom(ARRAY(1), false);
  geas::linear_le(SD, cons, vars0, INT(2));
  geas::linear_le(SD, cons, vars1, INT(2));
}

void p_int_lin_eq_reif(SolverInstanceBase& s, const Call* call) {
  vec<int> pos = INTARRAY(0);
  vec<int> neg(pos.size());
  for (int i = 0; i < neg.size(); ++i)
    neg[i] = -pos[i];
  vec<geas::intvar> vars0 = SI.asIntVarDom(ARRAY(1), true);
  geas::linear_le(SD, pos, vars0, INT(2), BOOLVAR0(3));
  geas::linear_le(SD, neg, vars0, -INT(2), BOOLVAR0(3));
  geas::linear_ne(SD, pos, vars0, INT(2), ~BOOLVAR0(3));

  vec<geas::intvar> vars1 = SI.asIntVarDom(ARRAY(1), false);
  geas::linear_le(SD, pos, vars1, INT(2), BOOLVAR1(3));
  geas::linear_le(SD, neg, vars1, -INT(2), BOOLVAR1(3));
  geas::linear_ne(SD, pos, vars1, INT(2), ~BOOLVAR1(3));
}

void p_int_lin_ne_reif(SolverInstanceBase& s, const Call* call) {
  vec<int> pos = INTARRAY(0);
  vec<int> neg(pos.size());
  for (int i = 0; i < neg.size(); ++i)
    neg[i] = -pos[i];
  vec<geas::intvar> vars0 = SI.asIntVarDom(ARRAY(1), true);
  geas::linear_ne(SD, pos, vars0, INT(2), BOOLVAR0(3));
  geas::linear_le(SD, pos, vars0, INT(2), ~BOOLVAR0(3));
  geas::linear_le(SD, neg, vars0, -INT(2), ~BOOLVAR0(3));

  vec<geas::intvar> vars1 = SI.asIntVarDom(ARRAY(1), false);
  geas::linear_ne(SD, pos, vars1, INT(2), BOOLVAR1(3));
  geas::linear_le(SD, pos, vars1, INT(2), ~BOOLVAR1(3));
  geas::linear_le(SD, neg, vars1, -INT(2), ~BOOLVAR1(3));
}

void p_int_lin_le_reif(SolverInstanceBase& s, const Call* call) {
  vec<int> pos = INTARRAY(0);
  vec<int> neg(pos.size());
  for (int i = 0; i < neg.size(); ++i)
    neg[i] = -pos[i];
  vec<geas::intvar> vars0 = SI.asIntVarDom(ARRAY(1), true);
  vec<geas::intvar> vars1 = SI.asIntVarDom(ARRAY(1), false);
  geas::linear_le(SD, pos, vars0, INT(2), BOOLVAR0(3));
  geas::linear_le(SD, neg, vars0, -INT(2) - 1, ~BOOLVAR0(3));
  geas::linear_le(SD, pos, vars1, INT(2), BOOLVAR1(3));
  geas::linear_le(SD, neg, vars1, -INT(2) - 1, ~BOOLVAR1(3));
}

void p_bool_eq(SolverInstanceBase& s, const Call* call) {
  if (PAR(0)) {
    SOL.post(BOOL(0) ? BOOLVAR0(1) : ~BOOLVAR0(1));
    SOL.post(BOOL(0) ? BOOLVAR1(1) : ~BOOLVAR1(1));
  } else if (PAR(1)) {
    SOL.post(BOOL(1) ? BOOLVAR0(0) : ~BOOLVAR0(0));
    SOL.post(BOOL(1) ? BOOLVAR1(0) : ~BOOLVAR1(0));
  } else {
    geas::add_clause(SD, BOOLVAR0(0), ~BOOLVAR0(1));
    geas::add_clause(SD, ~BOOLVAR0(0), BOOLVAR0(1));
    geas::add_clause(SD, BOOLVAR1(0), ~BOOLVAR1(1));
    geas::add_clause(SD, ~BOOLVAR1(0), BOOLVAR1(1));
  }
}

void p_bool_ne(SolverInstanceBase& s, const Call* call) {
  if (PAR(0)) {
    SOL.post(BOOL(0) ? ~BOOLVAR0(1) : BOOLVAR0(1));
    SOL.post(BOOL(0) ? ~BOOLVAR1(1) : BOOLVAR1(1));
  } else if (PAR(1)) {
    SOL.post(BOOL(1) ? ~BOOLVAR0(0) : BOOLVAR0(0));
    SOL.post(BOOL(1) ? ~BOOLVAR1(0) : BOOLVAR1(0));
  } else {
    geas::add_clause(SD, BOOLVAR0(0), BOOLVAR0(1));
    geas::add_clause(SD, ~BOOLVAR0(0), ~BOOLVAR0(1));
    geas::add_clause(SD, BOOLVAR1(0), BOOLVAR1(1));
    geas::add_clause(SD, ~BOOLVAR1(0), ~BOOLVAR1(1));
  }
}

void p_bool_le(SolverInstanceBase& s, const Call* call) {
  if (PAR(0)) {
    if (BOOL(0)) {
      SOL.post(BOOLVAR0(1));
      SOL.post(BOOLVAR1(1));
    }
  } else if (PAR(1)) {
    if (!BOOL(1)) {
      SOL.post(~BOOLVAR0(0));
      SOL.post(~BOOLVAR1(0));
    }
  } else {
    geas::add_clause(SD, ~BOOLVAR0(0), BOOLVAR0(1));
    geas::add_clause(SD, ~BOOLVAR1(0), BOOLVAR1(1));
  }
}

void p_bool_lt(SolverInstanceBase& s, const Call* call) {
  SOL.post(~BOOLVAR0(0));
  SOL.post(BOOLVAR0(1));
  SOL.post(~BOOLVAR1(0));
  SOL.post(BOOLVAR1(1));
}

void p_bool_eq_reif(SolverInstanceBase& s, const Call* call) {
  if (PAR(2)) {
    if (BOOL(2))
      p_bool_eq(s, call);
    else
      p_bool_ne(s, call);
  } else {
    geas::add_clause(SD, BOOLVAR0(2), BOOLVAR0(0), BOOLVAR0(1));
    geas::add_clause(SD, BOOLVAR0(2), ~BOOLVAR0(0), ~BOOLVAR0(1));
    geas::add_clause(SD, ~BOOLVAR0(2), ~BOOLVAR0(0), BOOLVAR0(1));
    geas::add_clause(SD, ~BOOLVAR0(2), BOOLVAR0(0), ~BOOLVAR0(1));

    geas::add_clause(SD, BOOLVAR1(2), BOOLVAR1(0), BOOLVAR1(1));
    geas::add_clause(SD, BOOLVAR1(2), ~BOOLVAR1(0), ~BOOLVAR1(1));
    geas::add_clause(SD, ~BOOLVAR1(2), ~BOOLVAR1(0), BOOLVAR1(1));
    geas::add_clause(SD, ~BOOLVAR1(2), BOOLVAR1(0), ~BOOLVAR1(1));
  }
}

void p_bool_ne_reif(SolverInstanceBase& s, const Call* call) {
  if (PAR(2)) {
    if (BOOL(2)) {
      p_bool_ne(s, call);
    } else {
      p_bool_eq(s, call);
    }
  } else {
    geas::add_clause(SD, BOOLVAR0(2), ~BOOLVAR0(0), BOOLVAR0(1));
    geas::add_clause(SD, BOOLVAR0(2), BOOLVAR0(0), ~BOOLVAR0(1));
    geas::add_clause(SD, ~BOOLVAR0(2), BOOLVAR0(0), BOOLVAR0(1));
    geas::add_clause(SD, ~BOOLVAR0(2), ~BOOLVAR0(0), ~BOOLVAR0(1));

    geas::add_clause(SD, BOOLVAR1(2), ~BOOLVAR1(0), BOOLVAR1(1));
    geas::add_clause(SD, BOOLVAR1(2), BOOLVAR1(0), ~BOOLVAR1(1));
    geas::add_clause(SD, ~BOOLVAR1(2), BOOLVAR1(0), BOOLVAR1(1));
    geas::add_clause(SD, ~BOOLVAR1(2), ~BOOLVAR1(0), ~BOOLVAR1(1));
  }
}

void p_bool_le_reif(SolverInstanceBase& s, const Call* call) {
  if (PAR(2)) {
    if (BOOL(2)) {
      p_bool_le(s, call);
    } else {
      auto* nc = new Call(Location().introduce(), call->id(), {call->arg(1), call->arg(0)});
      p_bool_lt(s, nc);
    }
  } else {
    geas::add_clause(SD, BOOLVAR0(2), ~BOOLVAR0(1));
    geas::add_clause(SD, BOOLVAR0(2), BOOLVAR0(0));
    geas::add_clause(SD, ~BOOLVAR0(2), ~BOOLVAR0(0), BOOLVAR0(1));

    geas::add_clause(SD, BOOLVAR1(2), ~BOOLVAR1(1));
    geas::add_clause(SD, BOOLVAR1(2), BOOLVAR1(0));
    geas::add_clause(SD, ~BOOLVAR1(2), ~BOOLVAR1(0), BOOLVAR1(1));
  }
}

void p_bool_lt_reif(SolverInstanceBase& s, const Call* call) {
  if (PAR(2)) {
    if (BOOL(2)) {
      p_int_lt(s, call);
    } else {
      auto* nc = new Call(Location().introduce(), call->id(), {call->arg(1), call->arg(0)});
      p_int_le(s, nc);
    }
  } else {
    geas::add_clause(SD, ~BOOLVAR0(2), ~BOOLVAR0(0));
    geas::add_clause(SD, ~BOOLVAR0(2), BOOLVAR0(1));
    geas::add_clause(SD, BOOLVAR0(2), BOOLVAR0(0), ~BOOLVAR0(1));

    geas::add_clause(SD, ~BOOLVAR1(2), ~BOOLVAR1(0));
    geas::add_clause(SD, ~BOOLVAR1(2), BOOLVAR1(1));
    geas::add_clause(SD, BOOLVAR1(2), BOOLVAR1(0), ~BOOLVAR1(1));
  }
}

void p_bool_or(SolverInstanceBase& s, const Call* call) {
  geas::add_clause(SD, BOOLVAR0(2), ~BOOLVAR0(0));
  geas::add_clause(SD, BOOLVAR0(2), ~BOOLVAR0(1));
  geas::add_clause(SD, ~BOOLVAR0(2), BOOLVAR0(0), BOOLVAR0(1));

  geas::add_clause(SD, BOOLVAR1(2), ~BOOLVAR1(0));
  geas::add_clause(SD, BOOLVAR1(2), ~BOOLVAR1(1));
  geas::add_clause(SD, ~BOOLVAR1(2), BOOLVAR1(0), BOOLVAR1(1));
}

void p_bool_and(SolverInstanceBase& s, const Call* call) {
  geas::add_clause(SD, ~BOOLVAR0(2), BOOLVAR0(0));
  geas::add_clause(SD, ~BOOLVAR0(2), BOOLVAR0(1));
  geas::add_clause(SD, BOOLVAR0(2), ~BOOLVAR0(0), ~BOOLVAR0(1));

  geas::add_clause(SD, ~BOOLVAR1(2), BOOLVAR1(0));
  geas::add_clause(SD, ~BOOLVAR1(2), BOOLVAR1(1));
  geas::add_clause(SD, BOOLVAR1(2), ~BOOLVAR1(0), ~BOOLVAR1(1));
}

void p_bool_xor(SolverInstanceBase& s, const Call* call) {
  if (call->argCount() == 2) {
    p_bool_ne(s, call);
  } else {
    p_bool_ne_reif(s, call);
  }
}

void p_bool_not(SolverInstanceBase& s, const Call* call) { p_bool_ne(s, call); }

void p_bool_clause(SolverInstanceBase& s, const Call* call) {
  auto& gi = static_cast<GeasSolverInstance&>(s);
  auto* pos = ARRAY(0);
  auto* neg = ARRAY(1);
  vec<geas::clause_elt> cl0, cl1; 
  for (unsigned int i = 0; i != pos->size(); i++) {
    cl0.push(SI.asBoolVarDom((*pos)[i], true));
    cl1.push(SI.asBoolVarDom((*pos)[i], false));
  }
  for (unsigned int j = 0; j != neg->size(); j++) {
    cl0.push(~SI.asBoolVarDom((*neg)[j], true));
    cl1.push(~SI.asBoolVarDom((*neg)[j], false));
  }
  geas::add_clause(*SD, cl0);
  geas::add_clause(*SD, cl1);
}

void p_array_bool_or(SolverInstanceBase& s, const Call* call) {
  auto* arr = ARRAY(0);
  vec<geas::clause_elt> clause0, clause1;
  clause0.push(~BOOLVAR0(1)); 
  clause1.push(~BOOLVAR1(1)); 
  for (unsigned int i = 0; i < arr->size(); i++) {
    geas::patom_t elem0 = SI.asBoolVarDom((*arr)[i], true);
    geas::patom_t elem1 = SI.asBoolVarDom((*arr)[i], false);
    geas::add_clause(SD, BOOLVAR0(1), ~elem0);
    geas::add_clause(SD, BOOLVAR1(1), ~elem1);
    clause0.push(elem0);
    clause1.push(elem1);
  }
  geas::add_clause(*SD, clause0);
  geas::add_clause(*SD, clause1);
}

void p_array_bool_and(SolverInstanceBase& s, const Call* call) {
  auto* arr = ARRAY(0);
  vec<geas::clause_elt> clause0, clause1;
  clause0.push(BOOLVAR0(1));
  clause1.push(BOOLVAR1(1));
  for (int i = 0; i < arr->size(); ++i) {
    geas::patom_t elem0 = SI.asBoolVarDom((*arr)[i], true);
    geas::patom_t elem1 = SI.asBoolVarDom((*arr)[i], false);
    geas::add_clause(SD, ~BOOLVAR0(1), elem0);
    geas::add_clause(SD, ~BOOLVAR1(1), elem1);
    clause0.push(~elem0);
    clause1.push(~elem1);
  }
  geas::add_clause(*SD, clause0);
  geas::add_clause(*SD, clause1);
}

void p_bool_clause_reif(SolverInstanceBase& s, const Call* call) {
  auto* pos = ARRAY(0); 
  auto* neg = ARRAY(1);
  vec<geas::clause_elt> cl0, cl1;
  cl0.push(~BOOLVAR0(2));
  cl1.push(~BOOLVAR1(2));
  for (int i = 0; i < pos->size(); ++i) {
    geas::patom_t elem0 = SI.asBoolVarDom((*pos)[i], true); 
    geas::patom_t elem1 = SI.asBoolVarDom((*pos)[i], false); 
    geas::add_clause(SD, BOOLVAR0(2), ~elem0);
    geas::add_clause(SD, BOOLVAR1(2), ~elem1);
    cl0.push(elem0); 
    cl1.push(elem1); 
  }
  for (int j = 0; j < neg->size(); ++j) {
    geas::patom_t elem0 = SI.asBoolVarDom((*neg)[j], true);
    geas::patom_t elem1 = SI.asBoolVarDom((*neg)[j], false);
    geas::add_clause(SD, BOOLVAR0(2), elem0);
    geas::add_clause(SD, BOOLVAR1(2), elem1);
    cl0.push(~elem0); 
    cl1.push(~elem1); 
  }
  geas::add_clause(*SD, cl0);
  geas::add_clause(*SD, cl1);
}

void p_bool_lin_eq(SolverInstanceBase& s, const Call* call) {
  vec<int> cons = INTARRAY(0);
  vec<geas::patom_t> vars0 = SI.asBoolVarDom(ARRAY(1), true);
  geas::bool_linear_le(SD, geas::at_True, SI.zero, cons, vars0, -INT(2));
  geas::bool_linear_ge(SD, geas::at_True, SI.zero, cons, vars0, -INT(2));

  vec<geas::patom_t> vars1 = SI.asBoolVarDom(ARRAY(1), true);
  geas::bool_linear_le(SD, geas::at_True, SI.zero, cons, vars1, -INT(2));
  geas::bool_linear_ge(SD, geas::at_True, SI.zero, cons, vars1, -INT(2));
}

void p_bool_lin_ne(SolverInstanceBase& s, const Call* call) {
  vec<int> cons = INTARRAY(0);
  vec<geas::patom_t> vars0 = SI.asBoolVarDom(ARRAY(1), true);
  geas::bool_linear_ne(SD, cons, vars0, INT(2));

  vec<geas::patom_t> vars1 = SI.asBoolVarDom(ARRAY(1), true);
  geas::bool_linear_ne(SD, cons, vars1, INT(2));
}

void p_bool_lin_le(SolverInstanceBase& s, const Call* call) {
  vec<int> cons = INTARRAY(0);
  vec<geas::patom_t> vars0 = SI.asBoolVarDom(ARRAY(1), true);
  geas::bool_linear_ge(SD, geas::at_True, SI.zero, cons, vars0, -INT(2));

  vec<geas::patom_t> vars1 = SI.asBoolVarDom(ARRAY(1), true);
  geas::bool_linear_ge(SD, geas::at_True, SI.zero, cons, vars1, -INT(2));
}

void p_bool_lin_eq_reif(SolverInstanceBase& s, const Call* call) {
  vec<int> cons = INTARRAY(0);
  vec<geas::patom_t> vars0 = SI.asBoolVarDom(ARRAY(1), true);
  geas::bool_linear_le(SD, BOOLVAR0(3), SI.zero, cons, vars0, -INT(2));
  geas::bool_linear_ge(SD, BOOLVAR0(3), SI.zero, cons, vars0, -INT(2));
  geas::bool_linear_ne(SD, cons, vars0, INT(2), ~BOOLVAR0(3));

  vec<geas::patom_t> vars1 = SI.asBoolVarDom(ARRAY(1), true);
  geas::bool_linear_le(SD, BOOLVAR1(3), SI.zero, cons, vars1, -INT(2));
  geas::bool_linear_ge(SD, BOOLVAR1(3), SI.zero, cons, vars1, -INT(2));
  geas::bool_linear_ne(SD, cons, vars1, INT(2), ~BOOLVAR1(3));
}

void p_bool_lin_ne_reif(SolverInstanceBase& s, const Call* call) {
  vec<int> cons = INTARRAY(0);
  vec<geas::patom_t> vars0 = SI.asBoolVarDom(ARRAY(1), true);
  geas::bool_linear_ne(SD, cons, vars0, INT(2), BOOLVAR0(3));
  geas::bool_linear_le(SD, ~BOOLVAR0(3), SI.zero, cons, vars0, -INT(2));
  geas::bool_linear_ge(SD, ~BOOLVAR0(3), SI.zero, cons, vars0, -INT(2));

  vec<geas::patom_t> vars1 = SI.asBoolVarDom(ARRAY(1), false);
  geas::bool_linear_ne(SD, cons, vars0, INT(2), BOOLVAR1(3));
  geas::bool_linear_le(SD, ~BOOLVAR1(3), SI.zero, cons, vars1, -INT(2));
  geas::bool_linear_ge(SD, ~BOOLVAR1(3), SI.zero, cons, vars1, -INT(2));
}

void p_bool_lin_le_reif(SolverInstanceBase& s, const Call* call) {
  vec<int> cons = INTARRAY(0);
  vec<geas::patom_t> vars0 = SI.asBoolVarDom(ARRAY(1), true);
  geas::bool_linear_ge(SD, BOOLVAR0(3), SI.zero, cons, vars0, -INT(2));
  geas::bool_linear_le(SD, ~BOOLVAR0(3), SI.zero, cons, vars0, -INT(2) - 1);

  vec<geas::patom_t> vars1 = SI.asBoolVarDom(ARRAY(1), true);
  geas::bool_linear_ge(SD, BOOLVAR1(3), SI.zero, cons, vars1, -INT(2));
  geas::bool_linear_le(SD, ~BOOLVAR1(3), SI.zero, cons, vars1, -INT(2) - 1);
}

void p_bool2int(SolverInstanceBase& s, const Call* call) {
  CONSTRAINTNOTIMPLEMENT; 
  geas::add_clause(SD, BOOLVAR0(0), INTVAR0(1) <= 0); 
  geas::add_clause(SD, ~BOOLVAR0(0), INTVAR0(1) >= 1); 
  geas::add_clause(SD, BOOLVAR1(0), INTVAR1(1) <= 0); 
  geas::add_clause(SD, ~BOOLVAR1(0), INTVAR1(1) >= 1); 
}

void p_array_int_element(SolverInstanceBase& s, const Call* call) {
  if ( call->ann().containsCall(std::string("defines_var")) || PAR(2) ) {
    assert(ARRAY(1)->min(0) == 1 && ARRAY(1)->max(0) == ARRAY(1)->size() + 1);
    vec<int> vals = INTARRAY(1);
    geas::int_element(SD, INTVAR0(2), INTVAR0(0), vals);
    geas::int_element(SD, INTVAR1(2), INTVAR1(0), vals);
  } else 
    throw InternalError(std::string(call->id().c_str()) + std::string(" constraint with non-variablve second argument"));
}

void p_array_bool_element(SolverInstanceBase& s, const Call* call) {
  if ( call->ann().containsCall(std::string("defines_var")) || PAR(2) ) {
    assert(ARRAY(1)->min(0) == 1 && ARRAY(1)->max(0) == ARRAY(1)->size() + 1);
    vec<bool> vals = BOOLARRAY(1);
    for (int j = 0; j < vals.size(); ++j) {
      geas::add_clause(SD, INTVAR0(0) != j + 1, vals[j] ? BOOLVAR0(2) : ~BOOLVAR0(2)); 
      geas::add_clause(SD, INTVAR1(0) != j + 1, vals[j] ? BOOLVAR1(2) : ~BOOLVAR1(2)); 
    }
  } else 
    throw InternalError(std::string(call->id().c_str()) + std::string(" constraint with non-variablve second argument"));
}

void p_all_different(SolverInstanceBase& s, const Call* call) {
  FUNCTIONNOTIMPLEMENT; 
  vec<geas::intvar> vars0 = SI.asIntVarDom(ARRAY(0), true);
  vec<geas::intvar> vars1 = SI.asIntVarDom(ARRAY(0), false);
  geas::all_different_int(SD, vars0);
  geas::all_different_int(SD, vars1);
}

void p_all_different_except_0(SolverInstanceBase& s, const Call* call) {
  FUNCTIONNOTIMPLEMENT; 
  vec<geas::intvar> vars0 = SI.asIntVarDom(ARRAY(0), true);
  vec<geas::intvar> vars1 = SI.asIntVarDom(ARRAY(0), false);
  geas::all_different_except_0(SD, vars0);
  geas::all_different_except_0(SD, vars1);
}

void p_array_int_minimum(SolverInstanceBase& s, const Call* call) {
  auto arr = ARRAY(1);
  vec<geas::intvar> aux0 = SI.asIntVarDom(ARRAY(1), true);
  vec<geas::intvar> aux1 = SI.asIntVarDom(ARRAY(1), false);
  vec<geas::intvar> vars0(aux0.size()), vars1(aux1.size()); 
  for (unsigned int i = 0; i != aux0.size(); i++) {
    vars0[i] = -aux0[i];
    vars1[i] = -aux1[i];
  }
  geas::int_max(SD, -INTVAR0(0), vars0);
  geas::int_max(SD, -INTVAR1(0), vars1);
}

void p_array_int_maximum(SolverInstanceBase& s, const Call* call) {
  vec<geas::intvar> vars0 = SI.asIntVarDom(ARRAY(1), true);
  vec<geas::intvar> vars1 = SI.asIntVarDom(ARRAY(1), false);
  geas::int_max(SD, INTVAR0(0), vars0); 
  geas::int_max(SD, INTVAR1(0), vars1);
}

// void p_at_most(SolverInstanceBase& s, const Call* call) {
//   vec<geas::intvar> ivars = INTVARARRAY(1);
//   vec<geas::patom_t> bvars;
//   for (auto& ivar : ivars) {
//     bvars.push(ivar == INT(2));
//   }

//   if (INT(0) == 1) {
//     geas::atmost_1(SD, bvars);
//   } else {
//     geas::atmost_k(SD, bvars, INT(0));
//   }
// }

// void p_at_most1(SolverInstanceBase& s, const Call* call) {
//   vec<geas::intvar> ivars = INTVARARRAY(0);
//   vec<geas::patom_t> bvars;
//   for (auto& ivar : ivars) {
//     bvars.push(ivar == INT(1));
//   }
//   geas::atmost_1(SD, bvars);
// }

// void p_cumulative(SolverInstanceBase& s, const Call* call) {
//   vec<geas::intvar> st = INTVARARRAY(0);
//   if (PAR(1) && PAR(2) && PAR(3)) {
//     vec<int> d = INTARRAY(1);
//     vec<int> r = INTARRAY(2);
//     geas::cumulative(SD, st, d, r, INT(3));
//   } else {
//     vec<geas::intvar> d = INTVARARRAY(1);
//     vec<geas::intvar> r = INTVARARRAY(2);
//     geas::cumulative_var(SD, st, d, r, INTVAR(3));
//   }
// }

// void p_disjunctive(SolverInstanceBase& s, const Call* call) {
//   vec<geas::intvar> st = INTVARARRAY(0);
//   if (PAR(1)) {
//     vec<int> d = INTARRAY(1);
//     geas::disjunctive_int(SD, st, d);
//   } else {
//     vec<geas::intvar> d = INTVARARRAY(1);
//     geas::disjunctive_var(SD, st, d);
//   }
// }

// void p_global_cardinality(SolverInstanceBase& s, const Call* call) {
//   vec<geas::intvar> x = INTVARARRAY(0);
//   vec<int> cover = INTARRAY(1);
//   vec<int> count = INTARRAY(2);

//   vec<int> srcs(x.size(), 1);
//   vec<geas::bflow> flows;
//   for (int i = 0; i < x.size(); ++i) {
//     for (int j = 0; j < cover.size(); ++j) {
//       if (x[i].lb(SD) <= cover[j] && cover[j] <= x[i].ub(SD)) {
//         flows.push({i, j, x[i] == cover[j]});
//       }
//     }
//   }
//   geas::bipartite_flow(SD, srcs, count, flows);
// }

void p_table_int(SolverInstanceBase& s, const Call* call) {
  vec<geas::intvar> vars0 = SI.asIntVarDom(ARRAY(0), true);
  vec<geas::intvar> vars1 = SI.asIntVarDom(ARRAY(0), false);
  
  vec<int> tmp = INTARRAY(1);
  assert(tmp.size() % vars0.size() == 0);
  vec<vec<int>> table(tmp.size() == 0 ? 0 : tmp.size() / vars0.size());
  for (int i = 0; i < table.size(); ++i) {
    table[i].growTo(vars0.size());
    for (int j = 0; j < vars0.size(); ++j) {
      table[i][j] = tmp[i * vars0.size() + j];
    }
  }
  geas::table_id id = geas::table::build(SD, table);
  geas::table::post(SD, id, vars0); 
  geas::table::post(SD, id, vars1); 
}

void p_set_in(SolverInstanceBase& s, const Call* call) {
  if (PAR(1) && !PAR(0)) {
    VarDecl* vd = VARID(0)->decl();
    geas::intvar v0 = INTVAR0(0);
    geas::intvar v1 = INTVAR1(0);
    IntSetVal* isv = eval_intset(SI.env().envi(), call->arg(1));
    if (isv->size() == 1) {
      geas::int_le(SD, INTVAR0(0), CONSTANT(isv->max().toInt()), 0); 
      geas::int_le(SD, CONSTANT(isv->min().toInt()), INTVAR0(0), 0); 
      geas::int_le(SD, INTVAR1(0), CONSTANT(isv->max().toInt()), 0); 
      geas::int_le(SD, CONSTANT(isv->min().toInt()), INTVAR1(0), 0); 
    } else {
      vec<geas::clause_elt> cl0, cl1; 
      if (vd->ti()->domain() != nullptr) {
        IntSetVal* dom = eval_intset(SI.env().envi(), vd->ti()->domain());
        for (IntSetRanges isr(isv); isr(); ++isr) {
          for (IntVal v = isr.min(); v < isr.max(); ++v)
            if (dom->contains(v)) {
              geas::patom_t b0 = SOL.new_boolvar(); 
              geas::patom_t b1 = SOL.new_boolvar(); 

              geas::int_eq(SD, INTVAR0(0), CONSTANT(v.toInt()), b0);
              geas::int_ne(SD, INTVAR0(0), CONSTANT(v.toInt()), ~b0);
              geas::int_eq(SD, INTVAR1(0), CONSTANT(v.toInt()), b1);
              geas::int_ne(SD, INTVAR1(0), CONSTANT(v.toInt()), ~b1);
              cl0.push(b0); 
              cl1.push(b1); 
            }
        }
      }
      geas::add_clause(*SD, cl0);
      geas::add_clause(*SD, cl1);
    }
  }
  
}

// analyzer for variables' monotonicity (whether increment or decrement is prefer)
void mono_inc(Monotonicity &child_mono, const Monotonicity &parent_mono)
{
  if (child_mono == VAR_NONE) 
    child_mono = parent_mono; 
  else if (parent_mono == VAR_EQL || child_mono != parent_mono) 
    child_mono = VAR_EQL; 
}

void mono_dec(Monotonicity &child_mono, const Monotonicity &parent_mono)
{
  if (child_mono == VAR_NONE) {
    if (parent_mono == VAR_INC)
      child_mono = VAR_DEC; 
    else if (parent_mono == VAR_DEC)
      child_mono = VAR_INC; 
    else 
      child_mono = VAR_EQL; 
  }
  else if (parent_mono == VAR_EQL || child_mono == parent_mono)
    child_mono = VAR_EQL; 
}

void a_set_in(SolverInstanceBase& s, Call* call) { 
  FUNCTIONNOTIMPLEMENT; 
  if (!PAR(0) && PAR(1) ) 
    SI._variableMono[VARID(0)] = VAR_EQL; 
  else 
    throw InternalError(std::string("Unknown structure for constraint ") + std::string(call->id().c_str()) );
}

/* constraint analyzer implementation */ 
void a_eql(SolverInstanceBase& s, Call* call) { 
  FUNCTIONNOTIMPLEMENT;
  Id* var0 = VARID(0); 
  Id* var1 = VARID(1); 
  if ( var0 != nullptr )
    SI._variableMono[var0] = VAR_EQL; 
  if ( var1 != nullptr )
    SI._variableMono[var1] = VAR_EQL; 
}

void a_get(SolverInstanceBase& s, Call* call) {
  FUNCTIONNOTIMPLEMENT;
  Id* var0 = VARID(0); 
  Id* var1 = VARID(1); 
  if (var0 != nullptr) 
    SI._variableMono[var0] = VAR_INC; 
  if (var1 != nullptr) 
    SI._variableMono[var1] = VAR_DEC; 
}

void a_let(SolverInstanceBase& s, Call* call) {
  FUNCTIONNOTIMPLEMENT;
  Id* var0 = VARID(0); 
  Id* var1 = VARID(1); 
  if (var0 != nullptr) 
    SI._variableMono[var0] = VAR_DEC; 
  if (var1 != nullptr) 
    SI._variableMono[var1] = VAR_INC; 
}

void a_eql_reif(SolverInstanceBase& s, Call* call) { 
  if (PAR(2)) 
    a_eql(s, call); 
  else {
    if (!call->ann().containsCall(std::string("defines_var")) && VARID(2) != nullptr) 
      SI._variableMono[VARID(2)] = VAR_EQL; 
    if ( VARID(0) != nullptr ) SI._variableMono[VARID(0)] = VAR_EQL; 
    if ( VARID(1) != nullptr ) SI._variableMono[VARID(1)] = VAR_EQL; 
  }
}

void a_let_reif(SolverInstanceBase& s, Call* call) {
  if (PAR(2)) {
    if (BOOL(2)) 
      a_let(s, call); 
    else 
      a_get(s, call); 
  }
  else if (!PAR(2) && call->ann().containsCall(std::string("defines_var"))) {
    Id* var0 = VARID(0); 
    Id* var1 = VARID(1); 
    if ( var0 != nullptr ) mono_dec(SI._variableMono[var0], CONTEXT); 
    if ( var1 != nullptr ) mono_inc(SI._variableMono[var1], CONTEXT); 
  } else if (!PAR(2) && !call->ann().containsCall(std::string("defines_var"))) {
    if (VARID(0)) SI._variableMono[VARID(0)] = VAR_EQL; 
    if (VARID(1)) SI._variableMono[VARID(1)] = VAR_EQL;
    if (VARID(2)) SI._variableMono[VARID(2)] = VAR_EQL; 
  } else 
    throw InternalError(std::string("Unknown structure for constraint ") + std::string(call->id().c_str()) );
}

void a_int_abs(SolverInstanceBase& s, Call* call) { 
  if (VARID(0) != nullptr) SI._variableMono[VARID(0)] = VAR_EQL; 
  if (VARID(1) != nullptr) SI._variableMono[VARID(1)] = VAR_EQL; 
}

// todo: inspect the domain of variables to determine the monotonicity for arguments of int_tims/int_div
void a_eql_binary_op(SolverInstanceBase& s, Call* call) {
  CONSTRAINTNOTIMPLEMENT; 
  if (VARID(0) != nullptr)
    SI._variableMono[VARID(0)] = VAR_EQL; 
  if (VARID(1) != nullptr)
    SI._variableMono[VARID(1)] = VAR_EQL;  
}

void a_inc_binary_op(SolverInstanceBase& s, Call* call) {
  if ( call->ann().containsCall(std::string("defines_var")) ) {
    if ( VARID(0) != nullptr )
      mono_inc(SI._variableMono[VARID(0)], CONTEXT); 
    if ( VARID(1) != nullptr )
      mono_inc(SI._variableMono[VARID(1)], CONTEXT); 
  } else {
    if ( VARID(0) != nullptr )
      SI._variableMono[VARID(0)] = VAR_EQL; 
    if ( VARID(1) != nullptr )
      SI._variableMono[VARID(1)] = VAR_EQL; 
    if ( VARID(2) != nullptr )
      SI._variableMono[VARID(2)] = VAR_EQL; 
  }
}

void a_lin_eq(SolverInstanceBase& s, Call* call) {
  if ( call->ann().containsCall(std::string("defines_var")) ) {
    // defines a variable 
    Id* defVar = SI.asVarId(call->ann().getCall(std::string("defines_var"))->arg(0)); 
    vec<int> cons = INTARRAY(0);
    std::vector<Id*> vars = VARIDARRAY(1);
    assert(cons.size() == vars.size());
    int a = 1; 
    for (unsigned int i = 0; i != vars.size(); i++) 
      if (vars[i] == defVar) {
        a = -cons[i]; 
        break;
      }
    for (unsigned int i = 0; i != vars.size(); i++) {
      if (vars[i] != nullptr && vars[i] != defVar) {
        if (cons[i]/a > 0) 
          mono_inc(SI._variableMono[vars[i]], CONTEXT); 
        else if (cons[i]/a < 0)
          mono_dec(SI._variableMono[vars[i]], CONTEXT); 
        SI._variableSensitive[vars[i]] = SENSITIVE; 
      }
    }
  }
  else {
    // a pure constraint 
    std::vector<Id*> vars = VARIDARRAY(1);
    for (unsigned int i = 0; i != vars.size(); i++) {
      if (vars[i] != nullptr)
        SI._variableMono[vars[i]] = VAR_EQL; 
    }
  }
}

void a_lin_ne(SolverInstanceBase& s, Call* call) {
  FUNCTIONNOTIMPLEMENT; 
  std::vector<Id*> vars = VARIDARRAY(1); 
  for (unsigned int i = 0; i != vars.size(); i++) {
    if (vars[i] != nullptr) 
      SI._variableMono[vars[i]] = VAR_EQL; 
  }
}

void a_lin_le(SolverInstanceBase& s, Call* call) {
  FUNCTIONNOTIMPLEMENT; 
  vec<int> cons = INTARRAY(0);
  std::vector<Id*> vars = VARIDARRAY(1);
  assert(cons.size() == vars.size());
  for (unsigned int i = 0; i != vars.size(); i++) {
    if (vars[i] != nullptr) {
      if (cons[i] > 0)
        SI._variableMono[vars[i]] = VAR_DEC; 
      else if (cons[i] < 0)
        SI._variableMono[vars[i]] = VAR_INC; 
    }
  }
}

void a_lin_le_reif(SolverInstanceBase& s, Call* call) {
  if (PAR(3) && BOOL(3))
    a_lin_le(s, call); 
  else if (PAR(3) && !BOOL(3))
    a_lin_ne(s, call); 
  else if (!PAR(3) && call->ann().containsCall(std::string("defines_var")) ) {
    vec<int> cons = INTARRAY(0);
    std::vector<Id*> vars = VARIDARRAY(1);
    assert(cons.size() == vars.size());
    for (unsigned int i = 0; i != vars.size(); i++) {
      if (vars[i] != nullptr) {
        if (cons[i] > 0)
          mono_dec(SI._variableMono[vars[i]], CONTEXT); 
        else if (cons[i] < 0)
          mono_inc(SI._variableMono[vars[i]], CONTEXT); 
      }
    }
  } else if (!PAR(3) && !call->ann().containsCall(std::string("defines_var")) ) {
    std::vector<Id*> vars = VARIDARRAY(1);
    for (unsigned int i = 0; i != vars.size(); i++) {
      if (vars[i] != nullptr)
        SI._variableMono[vars[i]] = VAR_EQL; 
    }
  } else 
    throw InternalError(std::string("Unknown structure for constraint ") + std::string(call->id().c_str()) );
}

void a_lin_eql_reif(SolverInstanceBase& s, Call* call) {
  CONSTRAINTNOTIMPLEMENT; 
  std::vector<Id*> vars = VARIDARRAY(1);
  assert(cons.size() == vars.size());
  for (unsigned int i = 0; i != vars.size(); i++) {
    if (vars[i] != nullptr)
      SI._variableMono[vars[i]] = VAR_EQL; 
  }
}

void a_bool_not(SolverInstanceBase& s, Call* call) {
  CONSTRAINTNOTIMPLEMENT; 
  if (VARID(1) != nullptr) 
    mono_dec(SI._variableMono[VARID(0)], CONTEXT); 
  else 
    throw InternalError(std::string("bool2int constraint with non-variablve integer argument"));
}

void a_bool_clause(SolverInstanceBase& s, Call* call) {
  FUNCTIONNOTIMPLEMENT; 
  std::vector<Id*> pos = VARIDARRAY(0);
  std::vector<Id*> neg = VARIDARRAY(1);
  for (unsigned int i = 0; i != pos.size(); i++) {
    if (pos[i] != nullptr)
      mono_inc(SI._variableMono[pos[i]], VAR_INC);
  }
  for (unsigned int i = 0; i != neg.size(); i++) {
    if (neg[i] != nullptr)
      mono_dec(SI._variableMono[neg[i]], VAR_INC);
  }
}

void a_array_bool_and_or(SolverInstanceBase& s, Call* call) {
  Monotonicity context0; 
  if ( call->ann().containsCall(std::string("defines_var")) )
    context0 = CONTEXT; 
  else if ( PAR(1) ) {
    if (BOOL(1)) 
      context0 = VAR_INC; 
    else 
      context0 = VAR_DEC; 
  }
  else if ( !PAR(1) ) {
    context0 = VAR_EQL; 
    SI._variableMono[VARID(1)] = VAR_EQL; 
  }
  else 
    throw InternalError(std::string("Constraint with exceptional case not implemented: ") + std::string(call->id().c_str()) );
  
  std::vector<Id*> arr = VARIDARRAY(0); 
  for (unsigned int i = 0; i != arr.size(); i++) {
    if (arr[i] != nullptr) 
      mono_inc(SI._variableMono[arr[i]], context0); 
  }
}

void a_array_bool_xor(SolverInstanceBase& s, Call* call) {
  FUNCTIONNOTIMPLEMENT; 
  std::vector<Id*> vars = VARIDARRAY(0);
  for (unsigned int i = 0; i != vars.size(); i++) {
    if (vars[i] != nullptr) {
      SI._variableMono[vars[i]] = VAR_EQL; 
    }
  }
}

void a_bool_clause_reif(SolverInstanceBase& s, Call* call) {
  CONSTRAINTNOTIMPLEMENT; 
  std::vector<Id*> pos = VARIDARRAY(0);
  std::vector<Id*> neg = VARIDARRAY(1);
  for (unsigned int i = 0; i != pos.size(); i++) {
    if (pos[i] != nullptr)
      mono_inc(SI._variableMono[pos[i]], CONTEXT);
  }
  for (unsigned int i = 0; i != neg.size(); i++) {
    if (neg[i] != nullptr)
      mono_dec(SI._variableMono[neg[i]], CONTEXT);
  }
}

void a_bool2int(SolverInstanceBase& s, Call* call) {
  CONSTRAINTNOTIMPLEMENT; 
  Id* bvid = VARID(0); 
  Id* ivid = VARID(1);
  if (ivid != nullptr) 
    mono_inc(SI._variableMono[bvid], CONTEXT); 
  else 
    throw InternalError(std::string("bool2int constraint with non-variablve integer argument"));
}

void a_array_lit_element(SolverInstanceBase& s, Call* call) {
  if ( call->ann().containsCall(std::string("defines_var")) || PAR(2) ) {
    assert(ARRAY(1)->min(0) == 1 && ARRAY(1)->max(0) == ARRAY(1)->size() + 1);
    Id* id = VARID(0);
    if (id != nullptr) 
      SI._variableMono[id] = VAR_EQL;
  } else 
    throw InternalError(std::string(call->id().c_str()) + std::string(" constraint with non-variablve second argument"));
}

void a_all_different(SolverInstanceBase& s, Call* call) {
  FUNCTIONNOTIMPLEMENT; 
  std::vector<Id*> vars = VARIDARRAY(0);
  for (unsigned int i = 0; i != vars.size(); i++) {
    if (vars[i] != nullptr) {
      SI._variableMono[vars[i]] = VAR_EQL; 
    }
  }
}

void a_array_int_min_max(SolverInstanceBase& s, Call* call) {
  if ( call->ann().containsCall(std::string("defines_var")) ) {
    std::vector<Id*> vars = VARIDARRAY(1); 
    for (unsigned int i = 0; i != vars.size(); i++) {
      if (vars[i] != nullptr) 
        mono_inc(SI._variableMono[vars[i]], CONTEXT); 
    }
  } else {
    std::vector<Id*> vars = VARIDARRAY(1); 
    for (unsigned int i = 0; i != vars.size(); i++) {
      if (vars[i] != nullptr) 
        SI._variableMono[vars[i]] = VAR_EQL; 
    }
    if (VARID(0) != nullptr)
      SI._variableMono[VARID(0)] = VAR_EQL; 
  }
}

void a_table_int(SolverInstanceBase& s, Call* call) {
  FUNCTIONNOTIMPLEMENT; 
  std::vector<Id*> vars = VARIDARRAY(0); 
  for (unsigned int i = 0; i != vars.size(); i++)
    if (vars[i] != nullptr) 
      SI._variableMono[vars[i]] = VAR_EQL; 
}

/* dominance condition poster */ 
void int_lin_partial_sum(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars, geas::intvar& sum0, geas::intvar& sum1) {
  vec<int> cons = INTARRAY(0);
  std::vector<Id*> varIds = VARIDARRAY(1);
  assert(varIds.size() == cons.size()); 
  vec<geas::intvar> vars0, vars1; 
  vec<int> pos, neg; 
  for (unsigned int i = 0; i != varIds.size(); i++) {
    if (varIds[i] != nullptr && fixedVars.find(varIds[i]) != fixedVars.end()) {
      vars0.push(SI.asIntVarDom(varIds[i], true)); 
      vars1.push(SI.asIntVarDom(varIds[i], false)); 
      pos.push(cons[i]);
      neg.push(-cons[i]);
    }
  }

  vars0.push(sum0); 
  vars1.push(sum1); 
  pos.push(-1); 
  neg.push(1);

  geas::linear_le(SD, pos, vars0, 0); 
  geas::linear_le(SD, pos, vars1, 0); 
  geas::linear_le(SD, neg, vars0, 0);
  geas::linear_le(SD, neg, vars1, 0);
}

void d_int_eql(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) { 
  FUNCTIONNOTIMPLEMENT; 
  Id* x = VARID(0); 
  if (x == nullptr || fixedVars.find(x) == fixedVars.end())
    x = VARID(1);
  geas::intvar x0 = SI.asIntVarDom(x, true);
  geas::intvar x1 = SI.asIntVarDom(x, false);
  geas::int_eq(SD, x0, x1);
}

void d_int_get(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  FUNCTIONNOTIMPLEMENT; 
  Id* x = VARID(0); 
  Id* y = VARID(1); 
  if (x != nullptr && fixedVars.find(x) != fixedVars.end()) {
    geas::intvar x0 = SI.asIntVarDom(x, true);
    geas::intvar x1 = SI.asIntVarDom(x, false);
    geas::int_le(SD, x1, x0, 0); 
  } else if (y != nullptr && fixedVars.find(y) != fixedVars.end()) {
    geas::intvar y0 = SI.asIntVarDom(y, true);
    geas::intvar y1 = SI.asIntVarDom(y, false);
    geas::int_le(SD, y0, y1, 0); 
  }
}

void d_int_let(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  FUNCTIONNOTIMPLEMENT; 
  Id* x = VARID(0); 
  Id* y = VARID(1); 
  if (x != nullptr && fixedVars.find(x) != fixedVars.end()) {
    geas::intvar x0 = SI.asIntVarDom(x, true);
    geas::intvar x1 = SI.asIntVarDom(x, false);
    geas::int_le(SD, x0, x1, 0); 
  } else if (y != nullptr && fixedVars.find(y) != fixedVars.end()) {
    geas::intvar y0 = SI.asIntVarDom(y, true);
    geas::intvar y1 = SI.asIntVarDom(y, false);
    geas::int_le(SD, y1, y0, 0); 
  }
}

void d_int_eql_reif(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) { 
  // CONSTRAINTNOTIMPLEMENT;
  if (PAR(2)) 
    d_int_eql(s, call, fixedVars); 
  else {
    if ( VARID(0) != nullptr && fixedVars.find(VARID(0)) != fixedVars.end()) 
      geas::int_eq(SD, INTVAR0(0), INTVAR1(0));
    if ( VARID(1) != nullptr && fixedVars.find(VARID(1)) != fixedVars.end()) 
      geas::int_eq(SD, INTVAR0(1), INTVAR1(1));
    if ( VARID(2) != nullptr && fixedVars.find(VARID(2)) != fixedVars.end()) {
      geas::add_clause(SD, SI.asBoolVarDom(VARID(2), true), ~SI.asBoolVarDom(VARID(2), false));
      geas::add_clause(SD, ~SI.asBoolVarDom(VARID(2), true), SI.asBoolVarDom(VARID(2), false)); 
    }
  }
}

// x ≤ y <-> b, there are three cases: 1) x is an int; 2) y is an int; 3) x, y are int vars
void d_int_let_reif(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  if (PAR(2)) {
    if (BOOL(2)) 
      d_int_let(s, call, fixedVars); 
    else 
      d_int_get(s, call, fixedVars); 
  } else if (!PAR(2) && call->ann().containsCall(std::string("defines_var"))) {
    Monotonicity mono = SI._variableMono[SI.asVarId(call->ann().getCall(std::string("defines_var"))->arg(0))]; 
    if (VARID(0) == nullptr || VARID(1) == nullptr) // if there is only one variable, throw exception 
      throw InternalError(std::string("no dominance condition for unary functional constraint: ") + std::string(call->id().c_str()) ); 
    else {
      // if there are two variables and only one is fixed, enforce ≤ or ≥
      if (fixedVars.find(VARID(0)) != fixedVars.end() && fixedVars.find(VARID(1)) == fixedVars.end()) {
        if (mono == VAR_DEC) 
          geas::int_le(SD, SI.asIntVarDom(VARID(0), false), SI.asIntVarDom(VARID(0), true), 0); 
        else 
          geas::int_le(SD, SI.asIntVarDom(VARID(0), true), SI.asIntVarDom(VARID(0), false), 0); 
      } else if (fixedVars.find(VARID(1)) != fixedVars.end() && fixedVars.find(VARID(0)) == fixedVars.end()) {
        if (mono == VAR_DEC)
          geas::int_le(SD, SI.asIntVarDom(VARID(1), true), SI.asIntVarDom(VARID(1), false), 0); 
        else 
          geas::int_le(SD, SI.asIntVarDom(VARID(1), false), SI.asIntVarDom(VARID(1), true), 0); 
      } else 
        throw InternalError(std::string("Unknown structure for constraint: ") + std::string(call->id().c_str()) ); 
    }
  } else if (!PAR(2) && !call->ann().containsCall(std::string("defines_var"))) {
    // enforce all equal constraint for all input parameters 
    if (VARID(0) != nullptr && fixedVars.find(VARID(0)) != fixedVars.end()) 
      geas::int_eq(SD, SI.asIntVarDom(VARID(0), true), SI.asIntVarDom(VARID(0), false));
    if (VARID(1) != nullptr && fixedVars.find(VARID(1)) != fixedVars.end()) 
      geas::int_eq(SD, SI.asIntVarDom(VARID(1), true), SI.asIntVarDom(VARID(1), false));
    if (VARID(2) != nullptr && fixedVars.find(VARID(2)) != fixedVars.end()) {
      geas::add_clause(SD, SI.asBoolVarDom(VARID(2), true), ~SI.asBoolVarDom(VARID(2), false));
      geas::add_clause(SD, ~SI.asBoolVarDom(VARID(2), true), SI.asBoolVarDom(VARID(2), false)); 
    }
  } else 
    throw InternalError(std::string("Unknown structure for constraint ") + std::string(call->id().c_str()) );
}

// x op y = z, there are three cases: 1) x is an int; 2) y is an int; 3) x, y are int vars
void d_int_binary_op(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  if ( call->ann().containsCall(std::string("defines_var")) ) {
    Id* x = VARID(0); 
    if (x == nullptr || fixedVars.find(x) == fixedVars.end()) 
      x = VARID(1);
    if (SI._variableMono[x] == VAR_INC) 
      geas::int_le(SD, SI.asIntVarDom(x, false), SI.asIntVarDom(x, true), 0); 
    else if (SI._variableMono[x] == VAR_DEC) 
      geas::int_le(SD, SI.asIntVarDom(x, true), SI.asIntVarDom(x, false), 0);
    else 
      geas::int_eq(SD, SI.asIntVarDom(x, true), SI.asIntVarDom(x, false));  
  } else {
    if (VARID(0) == nullptr || fixedVars.find(VARID(0)) == fixedVars.end()) 
      geas::int_eq(SD, SI.asIntVarDom(VARID(0), true), SI.asIntVarDom(VARID(0), false));  
    if (VARID(1) == nullptr || fixedVars.find(VARID(1)) == fixedVars.end()) 
      geas::int_eq(SD, SI.asIntVarDom(VARID(1), true), SI.asIntVarDom(VARID(1), false));  
    if (VARID(2) == nullptr || fixedVars.find(VARID(2)) == fixedVars.end()) 
      geas::int_eq(SD, SI.asIntVarDom(VARID(2), true), SI.asIntVarDom(VARID(2), false));  
  }
}

// todo: need to compute the upper bound and lower bound of the sum 
void d_int_lin_eq_func(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  int a = 1; 
  Id* defVar = SI.asVarId(call->ann().getCall(std::string("defines_var"))->arg(0)); ; 
  Monotonicity mono = SI._variableMono[defVar]; 

  vec<int> cons = INTARRAY(0);
  std::vector<Id*> varIds = VARIDARRAY(1);
  assert(varIds.size() == cons.size()); 
  for (unsigned int i = 0; i != varIds.size(); i++) {
    if (varIds[i] != nullptr && varIds[i] == defVar) {
      a = -cons[i]; 
      break;
    }
  }
  
  vec<geas::intvar> vars0, vars1; 
  vec<int> pos, neg; 
  for (unsigned int i = 0; i != varIds.size(); i++) {
    if (varIds[i] != nullptr && fixedVars.find(varIds[i]) != fixedVars.end()) {
      vars0.push(SI.asIntVarDom(varIds[i], true)); 
      vars1.push(SI.asIntVarDom(varIds[i], false)); 
      pos.push(cons[i]/a);
      neg.push(-cons[i]/a);
    }
  }
  geas::intvar sum0 = SOL.new_intvar(SHRT_MIN, SHRT_MAX); 
  geas::intvar sum1 = SOL.new_intvar(SHRT_MIN, SHRT_MAX); 
  vars0.push(sum0); 
  vars1.push(sum1); 
  pos.push(-1); 
  neg.push(1); 

  if (mono == VAR_INC) {
    geas::linear_le(SD, pos, vars1, 0); 
    geas::linear_le(SD, neg, vars0, 0);
    geas::int_le(SD, sum1, sum0, 0);
    if (defVar != nullptr && SI._variableSensitive[defVar]) {
      SI._sensvar.push_back(sum1);
      SI._sensvar.push_back(sum0);
    }
  }
  else if (mono == VAR_DEC) {
    geas::linear_le(SD, pos, vars0, 0); 
    geas::linear_le(SD, neg, vars1, 0);
    geas::int_le(SD, sum0, sum1, 0);
    if (defVar != nullptr && SI._variableSensitive[defVar]) {
      SI._sensvar.push_back(sum0);
      SI._sensvar.push_back(sum1);
    }
  }
  else if (mono == VAR_EQL) {
    geas::linear_le(SD, pos, vars0, 0); 
    geas::linear_le(SD, pos, vars1, 0); 
    geas::linear_le(SD, neg, vars0, 0);
    geas::linear_le(SD, neg, vars1, 0);
    geas::int_eq(SD, sum0, sum1);
  }
}

void d_int_lin_eq(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  if (call->ann().containsCall(std::string("defines_var"))) 
    d_int_lin_eq_func(s, call, fixedVars); 
  else 
    d_int_lin_eql(s, call, fixedVars); 
}

void d_int_lin_eql(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  FUNCTIONNOTIMPLEMENT; 
  geas::intvar sum0 = SOL.new_intvar(SHRT_MIN, SHRT_MAX); 
  geas::intvar sum1 = SOL.new_intvar(SHRT_MIN, SHRT_MAX);

  vec<int> cons = INTARRAY(0);
  std::vector<Id*> varIds = VARIDARRAY(1);
  assert(varIds.size() == cons.size()); 
  vec<geas::intvar> vars0, vars1; 
  bool allpos = true, allneg = true; 
  vec<int> pos, neg; 
  for (unsigned int i = 0; i != varIds.size(); i++) {
    if (varIds[i] != nullptr && fixedVars.find(varIds[i]) != fixedVars.end()) {
      vars0.push(SI.asIntVarDom(varIds[i], true)); 
      vars1.push(SI.asIntVarDom(varIds[i], false)); 
      pos.push(cons[i]);
      neg.push(-cons[i]);
      allpos = allpos && (cons[i] >= 0);
      allneg = allneg && (cons[i] <= 0);
    }
  }
  vars0.push(sum0); 
  vars1.push(sum1); 
  pos.push(-1); 
  neg.push(1);

  geas::linear_le(SD, pos, vars0, 0); 
  geas::linear_le(SD, pos, vars1, 0); 
  geas::linear_le(SD, neg, vars0, 0);
  geas::linear_le(SD, neg, vars1, 0);

  if (allpos) {
    SOL.post(sum0 <= INT(2)); 
    SOL.post(sum1 <= INT(2)); 
  } else if (allneg) {
    SOL.post(sum0 >= INT(2)); 
    SOL.post(sum1 >= INT(2)); 
  }
  geas::int_eq(SD, sum0, sum1);
}

void d_int_lin_le(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  FUNCTIONNOTIMPLEMENT; 
  vec<int> cons = INTARRAY(0);
  std::vector<Id*> varIds = VARIDARRAY(1);
  assert(varIds.size() == cons.size()); 
  vec<geas::intvar> vars0, vars1;
  vec<int> pos, neg;  
  vec<int> coeff; 
  for (unsigned int i = 0; i != varIds.size(); i++) {
    if (varIds[i] != nullptr && fixedVars.find(varIds[i]) != fixedVars.end()) {
      vars0.push(SI.asIntVarDom(varIds[i], true)); 
      vars1.push(SI.asIntVarDom(varIds[i], false)); 
      pos.push(cons[i]);
      neg.push(-cons[i]);
    }
  }
  geas::intvar sum0 = SOL.new_intvar(SHRT_MIN, SHRT_MAX); 
  geas::intvar sum1 = SOL.new_intvar(SHRT_MIN, SHRT_MAX); 
  vars0.push(sum0); 
  vars1.push(sum1); 
  pos.push(-1); 
  neg.push(1); 

  geas::linear_le(SD, pos, vars0, 0); 
  geas::linear_le(SD, neg, vars1, 0);

  geas::int_le(SD, sum0, sum1, 0);

  SI._sensvar.push_back(sum0);
  SI._sensvar.push_back(sum1);
}

void d_int_lin_eql_reif(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  CONSTRAINTNOTIMPLEMENT; 
  geas::intvar v0 = SOL.new_intvar(SHRT_MIN, SHRT_MAX); 
  geas::intvar v1 = SOL.new_intvar(SHRT_MIN, SHRT_MAX);
  int_lin_partial_sum(s, call, fixedVars, v0, v1); 
  geas::int_eq(SD, v0, v1);
}

void d_int_lin_le_reif(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  if (PAR(3) && BOOL(3))
    d_int_lin_le(s, call, fixedVars); 
  else if (PAR(3) && !BOOL(3))
    d_int_lin_eql(s, call, fixedVars); 
  else if (!PAR(3)) {
    Monotonicity mono = VAR_EQL; 
    if (call->ann().containsCall(std::string("defines_var")))
      mono = SI._variableMono[SI.asVarId(call->ann().getCall(std::string("defines_var"))->arg(0))]; 
    else if ( VARID(3) != nullptr && fixedVars.find(VARID(3)) != fixedVars.end() ) {
      // enforce equality of boolean variable 
      geas::add_clause(SD, SI.asBoolVarDom(VARID(3), true), ~SI.asBoolVarDom(VARID(3), false));
      geas::add_clause(SD, ~SI.asBoolVarDom(VARID(3), true), SI.asBoolVarDom(VARID(3), false)); 
    }
    geas::intvar v0 = SOL.new_intvar(SHRT_MIN, SHRT_MAX); 
    geas::intvar v1 = SOL.new_intvar(SHRT_MIN, SHRT_MAX);
    int_lin_partial_sum(s, call, fixedVars, v0, v1); 

    if (mono == VAR_INC)
      geas::int_le(SD, v0, v1, 0);
    else if (mono == VAR_DEC)
      geas::int_le(SD, v1, v0, 0);
    else if (mono == VAR_EQL)
      geas::int_eq(SD, v0, v1); 
  } else 
    throw InternalError(std::string("Unknown structure for constraint ") + std::string(call->id().c_str()) );
}

void d_bool_eql(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  FUNCTIONNOTIMPLEMENT; 
  Id* x = VARID(0); 
  if (x == nullptr || fixedVars.find(x) == fixedVars.end())
    x = VARID(1);
  geas::add_clause(SD, SI.asBoolVarDom(x, true), ~SI.asBoolVarDom(x, false));
  geas::add_clause(SD, ~SI.asBoolVarDom(x, true), SI.asBoolVarDom(x, false)); 
}

void d_bool_let(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  FUNCTIONNOTIMPLEMENT; 
  Id* x = VARID(0); 
  if (x != nullptr && fixedVars.find(x) == fixedVars.end()) 
    geas::add_clause(SD, SI.asBoolVarDom(x, true), ~SI.asBoolVarDom(x, false)); 
  else {
    x = VARID(1);
    geas::add_clause(SD, SI.asBoolVarDom(x, false), ~SI.asBoolVarDom(x, true)); 
  }
}

void d_bool_get(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  FUNCTIONNOTIMPLEMENT; 
  Id* x = VARID(0); 
  if (x != nullptr && fixedVars.find(x) == fixedVars.end()) 
    geas::add_clause(SD, SI.asBoolVarDom(x, false), ~SI.asBoolVarDom(x, true)); 
  else {
    x = VARID(1);
    geas::add_clause(SD, SI.asBoolVarDom(x, true), ~SI.asBoolVarDom(x, false)); 
  }
}

void d_bool_eql_reif(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  if (PAR(2)) {
      d_bool_eql(s, call, fixedVars);
  } else {
    if (VARID(0) != nullptr && fixedVars.find(VARID(0)) != fixedVars.end() ) { 
      geas::add_clause(SD, BOOLVAR0(0), ~BOOLVAR1(0));
      geas::add_clause(SD, ~BOOLVAR0(0), BOOLVAR0(0));
    } 
    if (VARID(1) != nullptr && fixedVars.find(VARID(1)) != fixedVars.end() ) {
      geas::add_clause(SD, BOOLVAR0(1), ~BOOLVAR1(1));
      geas::add_clause(SD, ~BOOLVAR0(1), BOOLVAR0(1));
    } 
    if (VARID(2) != nullptr && fixedVars.find(VARID(2)) != fixedVars.end() ) {
      geas::add_clause(SD, BOOLVAR0(2), ~BOOLVAR1(2));
      geas::add_clause(SD, ~BOOLVAR0(2), BOOLVAR0(2));
    }
  }
}

void d_bool_let_reif(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  if (PAR(2)) {
    if (BOOL(2)) 
      d_int_let(s, call, fixedVars); 
    else 
      d_int_get(s, call, fixedVars); 
  } else if (!PAR(2) && call->ann().containsCall(std::string("defines_var"))) {
    Monotonicity mono = SI._variableMono[SI.asVarId(call->ann().getCall(std::string("defines_var"))->arg(0))]; 
    if (VARID(0) == nullptr || VARID(1) == nullptr) // if there is only one variable, throw exception 
      throw InternalError(std::string("no dominance condition for unary functional constraint: ") + std::string(call->id().c_str()) ); 
    else {
      if (fixedVars.find(VARID(0)) != fixedVars.end() && fixedVars.find(VARID(1)) == fixedVars.end()) {
        if (mono == VAR_DEC) 
          geas::add_clause(SD, SI.asBoolVarDom(VARID(0), true), ~SI.asBoolVarDom(VARID(0), false)); 
        else 
          geas::add_clause(SD, SI.asBoolVarDom(VARID(0), false), ~SI.asBoolVarDom(VARID(0), true)); 
      } else if (fixedVars.find(VARID(1)) != fixedVars.end() && fixedVars.find(VARID(0)) == fixedVars.end()) {
        if (mono == VAR_DEC) 
          geas::add_clause(SD, ~SI.asBoolVarDom(VARID(1), true), SI.asBoolVarDom(VARID(1), false)); 
        else 
          geas::add_clause(SD, ~SI.asBoolVarDom(VARID(1), false), SI.asBoolVarDom(VARID(1), true)); 
      } else 
        throw InternalError(std::string("Unknown structure for constraint: ") + std::string(call->id().c_str()) ); 
    }
  } else if (!PAR(2) && !call->ann().containsCall(std::string("defines_var"))) {
    // enforce all equal constraint for all input parameters 
    if (VARID(0) != nullptr && fixedVars.find(VARID(0)) != fixedVars.end()) {
      geas::add_clause(SD, SI.asBoolVarDom(VARID(0), true), ~SI.asBoolVarDom(VARID(0), false));
      geas::add_clause(SD, ~SI.asBoolVarDom(VARID(0), true), SI.asBoolVarDom(VARID(0), false)); 
    }
    if (VARID(1) != nullptr && fixedVars.find(VARID(1)) != fixedVars.end()) {
      geas::add_clause(SD, SI.asBoolVarDom(VARID(1), true), ~SI.asBoolVarDom(VARID(1), false));
      geas::add_clause(SD, ~SI.asBoolVarDom(VARID(1), true), SI.asBoolVarDom(VARID(1), false)); 
    }
    if (VARID(2) != nullptr && fixedVars.find(VARID(2)) != fixedVars.end()) {
      geas::add_clause(SD, SI.asBoolVarDom(VARID(2), true), ~SI.asBoolVarDom(VARID(2), false));
      geas::add_clause(SD, ~SI.asBoolVarDom(VARID(2), true), SI.asBoolVarDom(VARID(2), false)); 
    }
  } else 
    throw InternalError(std::string("Unknown structure for constraint ") + std::string(call->id().c_str()) );
}

// x op y = b, there are three cases: 1) x is a bool; 2) y is a bool; 3) x, y are bool vars
void d_bool_binary_op(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  if ( call->ann().containsCall(std::string("defines_var")) ) {
    Id* x = VARID(0); 
    if (x == nullptr || fixedVars.find(x) == fixedVars.end()) 
      x = VARID(1);
    if (SI._variableMono[x] == VAR_INC) 
      geas::add_clause(SD, SI.asBoolVarDom(x, false), ~SI.asBoolVarDom(x, true));
    else if (SI._variableMono[x] == VAR_DEC) 
      geas::add_clause(SD, SI.asBoolVarDom(x, true), ~SI.asBoolVarDom(x, false));
    else {
      geas::add_clause(SD, SI.asBoolVarDom(x, false), ~SI.asBoolVarDom(x, true));
      geas::add_clause(SD, SI.asBoolVarDom(x, true), ~SI.asBoolVarDom(x, false));
    }
  } else {
    if (VARID(0) == nullptr || fixedVars.find(VARID(0)) == fixedVars.end()) {
      geas::add_clause(SD, SI.asBoolVarDom(VARID(0), false), ~SI.asBoolVarDom(VARID(0), true));
      geas::add_clause(SD, SI.asBoolVarDom(VARID(0), true), ~SI.asBoolVarDom(VARID(0), false));
    }
    if (VARID(1) == nullptr || fixedVars.find(VARID(1)) == fixedVars.end()) {
      geas::add_clause(SD, SI.asBoolVarDom(VARID(1), false), ~SI.asBoolVarDom(VARID(1), true));
      geas::add_clause(SD, SI.asBoolVarDom(VARID(1), true), ~SI.asBoolVarDom(VARID(1), false));
    }
    if (VARID(2) == nullptr || fixedVars.find(VARID(2)) == fixedVars.end()) {
      geas::add_clause(SD, SI.asBoolVarDom(VARID(2), false), ~SI.asBoolVarDom(VARID(2), true));
      geas::add_clause(SD, SI.asBoolVarDom(VARID(2), true), ~SI.asBoolVarDom(VARID(2), false));
    }
  }
}

void d_bool_clause(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  FUNCTIONNOTIMPLEMENT; 
  std::vector<Id*> pos = VARIDARRAY(0); 
  std::vector<Id*> neg = VARIDARRAY(1); 
  vec<geas::clause_elt> cl0, cl1; 
  geas::patom_t bv0 = SOL.new_boolvar(); 
  geas::patom_t bv1 = SOL.new_boolvar(); 
  cl0.push(~bv0);
  cl1.push(~bv1);
  for (unsigned int i = 0; i != pos.size(); i++) {
    if (pos[i] != nullptr && fixedVars.find(pos[i]) != fixedVars.end()) {
      geas::patom_t elem0 = SI.asBoolVarDom(pos[i], true);
      geas::add_clause(SD, bv0, ~elem0); 
      cl0.push(elem0);

      geas::patom_t elem1 = SI.asBoolVarDom(pos[i], false);
      geas::add_clause(SD, bv1, ~elem1); 
      cl1.push(elem1);
    }
  }
  for (unsigned int j = 0; j != neg.size(); j++) {
    if (neg[j] != nullptr && fixedVars.find(neg[j]) != fixedVars.end()) {
      geas::patom_t elem0 = SI.asBoolVarDom(neg[j], true);
      geas::add_clause(SD, bv0, elem0);
      cl0.push(~elem0);

      geas::patom_t elem1 = SI.asBoolVarDom(neg[j], false);
      geas::add_clause(SD, bv1, elem1);
      cl1.push(~elem1);
    }
  }
  geas::add_clause(*SD, cl0);
  geas::add_clause(*SD, cl1);
  geas::add_clause(SD, ~bv1, bv0); 
}

void d_array_bool_or(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  Monotonicity mono = VAR_DEC; 
  if ( call->ann().containsCall(std::string("defines_var")) )
    mono = SI._variableMono[SI.asVarId(call->ann().getCall(std::string("defines_var"))->arg(0))]; 
  else if ( PAR(1) ) { 
    if (BOOL(1)) 
      mono = VAR_INC; 
    else 
      mono = VAR_DEC; 
  }
  else if ( !PAR(1) ) {
    mono = VAR_EQL; 
    if (VARID(1) != nullptr && fixedVars.find(VARID(1)) != fixedVars.end()) {
      geas::add_clause(SD, BOOLVAR0(1), ~BOOLVAR1(1));
      geas::add_clause(SD, ~BOOLVAR0(1), BOOLVAR1(1));
    }
  }
  else 
    throw InternalError(std::string("Constraint with exceptional case not implemented: ") + std::string(call->id().c_str()) ); 
  
  std::vector<Id*> vars = VARIDARRAY(0); 
  vec<geas::patom_t> vars0, vars1; 
  for (unsigned int i = 0; i != vars.size(); i++) {
    if (vars[i] != nullptr && fixedVars.find(vars[i]) != fixedVars.end()) {
      vars0.push(SI.asBoolVarDom(vars[i], true)); 
      vars1.push(SI.asBoolVarDom(vars[i], false)); 
    }
  }
  // post array or constraint 
  if (vars0.size() != 0) {
    geas::patom_t bv0 = SOL.new_boolvar(); 
    geas::patom_t bv1 = SOL.new_boolvar(); 
    vec<geas::clause_elt> clause0, clause1;
    clause0.push(~bv0); 
    clause1.push(~bv1); 
    for (unsigned int i = 0; i != vars0.size(); i++) {
      geas::add_clause(SD, bv0, ~vars0[i]);
      geas::add_clause(SD, bv1, ~vars1[i]);
      clause0.push(vars0[i]);
      clause1.push(vars1[i]);
    }
    geas::add_clause(*SD, clause0);
    geas::add_clause(*SD, clause1);

    if (mono == VAR_INC)
      geas::add_clause(SD, ~bv1, bv0);
    else if (mono == VAR_DEC)
      geas::add_clause(SD, ~bv0, bv1);
    else if (mono == VAR_EQL) {
      geas::add_clause(SD, bv0, ~bv1);
      geas::add_clause(SD, ~bv0, bv1);
    }
  } 
}

void d_array_bool_and(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  Monotonicity mono = VAR_DEC; 
  if ( call->ann().containsCall(std::string("defines_var")) )
    mono = SI._variableMono[SI.asVarId(call->ann().getCall(std::string("defines_var"))->arg(0))]; 
  else if ( PAR(1) ) {
    if (BOOL(1)) 
      mono = VAR_INC; 
    else 
      mono = VAR_DEC; 
  } else if (!PAR(1)) {
    mono = VAR_EQL; 
    if (VARID(1) != nullptr && fixedVars.find(VARID(1)) != fixedVars.end()) {
      geas::add_clause(SD, BOOLVAR0(1), ~BOOLVAR1(1));
      geas::add_clause(SD, ~BOOLVAR0(1), BOOLVAR1(1));
    }
  }
  else 
    throw InternalError(std::string("Constraint with exceptional case not implemented: ") + std::string(call->id().c_str()) ); 
  
  std::vector<Id*> vars = VARIDARRAY(0); 
  vec<geas::patom_t> vars0, vars1; 
  for (unsigned int i = 0; i != vars.size(); i++) {
    if (vars[i] != nullptr && fixedVars.find(vars[i]) != fixedVars.end()) {
      vars0.push(SI.asBoolVarDom(vars[i], true)); 
      vars1.push(SI.asBoolVarDom(vars[i], false)); 
    }
  }
  // post array and constraint 
  if (vars0.size() != 0) {
    geas::patom_t bv0 = SOL.new_boolvar(); 
    geas::patom_t bv1 = SOL.new_boolvar(); 
    vec<geas::clause_elt> clause0, clause1;
    clause0.push(bv0); 
    clause1.push(bv1); 
    for (unsigned int i = 0; i != vars0.size(); i++) {
      geas::add_clause(SD, ~bv0, vars0[i]);
      geas::add_clause(SD, ~bv1, vars1[i]);
      clause0.push(~vars0[i]);
      clause1.push(~vars1[i]);
    }
    geas::add_clause(*SD, clause0);
    geas::add_clause(*SD, clause1);

    if (mono == VAR_INC)
      geas::add_clause(SD, ~bv1, bv0);
    else if (mono == VAR_DEC)
      geas::add_clause(SD, ~bv0, bv1);
    else if (mono == VAR_EQL) {
      geas::add_clause(SD, ~bv0, bv1);
      geas::add_clause(SD, bv0, ~bv1);
    }
  } 
}

void d_bool_clause_reif(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  CONSTRAINTNOTIMPLEMENT; 
  Monotonicity mono = SI._variableMono[SI.asVarId(call->ann().getCall(std::string("defines_var"))->arg(0))]; 
  std::vector<Id*> pos = VARIDARRAY(0); 
  std::vector<Id*> neg = VARIDARRAY(1); 
  vec<geas::clause_elt> cl0, cl1; 
  geas::patom_t bv0 = SOL.new_boolvar(); 
  geas::patom_t bv1 = SOL.new_boolvar(); 
  cl0.push(~bv0);
  cl1.push(~bv1);
  for (unsigned int i = 0; i != pos.size(); i++) {
    if (pos[i] != nullptr && fixedVars.find(pos[i]) != fixedVars.end()) {
      geas::patom_t elem0 = SI.asBoolVarDom(pos[i], true);
      geas::add_clause(SD, bv0, ~elem0); 
      cl0.push(elem0);

      geas::patom_t elem1 = SI.asBoolVarDom(pos[i], false);
      geas::add_clause(SD, bv1, ~elem1); 
      cl1.push(elem1);
    }
  }
  for (unsigned int j = 0; j != neg.size(); j++) {
    if (neg[j] != nullptr && fixedVars.find(neg[j]) != fixedVars.end()) {
      geas::patom_t elem0 = SI.asBoolVarDom(neg[j], true);
      geas::add_clause(SD, bv0, elem0);
      cl0.push(~elem0);

      geas::patom_t elem1 = SI.asBoolVarDom(neg[j], false);
      geas::add_clause(SD, bv1, elem1);
      cl1.push(~elem1);
    }
  }
  geas::add_clause(*SD, cl0);
  geas::add_clause(*SD, cl1);
  if (mono == VAR_INC || mono == VAR_EQL) 
    geas::add_clause(SD, ~bv1, bv0);
  if (mono == VAR_DEC || mono == VAR_EQL) 
    geas::add_clause(SD, ~bv0, bv1);
}

// todo: implement the dominance condition for boolean linear constraints 
void bool_lin_partial_sum(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars, geas::intvar& sum0, geas::intvar& sum1) {
  vec<int> cons = INTARRAY(0);
  std::vector<Id*> varIds = VARIDARRAY(1);
  assert(varIds.size() == cons.size()); 
  vec<geas::patom_t> vars0, vars1; 
  vec<int> parcons; 
  for (unsigned int i = 0; i != varIds.size(); i++) {
    if (varIds[i] != nullptr && fixedVars.find(varIds[i]) != fixedVars.end()) {
      vars0.push(SI.asBoolVarDom(varIds[i], true)); 
      vars1.push(SI.asBoolVarDom(varIds[i], false)); 
      parcons.push(cons[i]);
    }
  }

  geas::bool_linear_le(SD, geas::at_True, sum0, parcons, vars0, 0);
  geas::bool_linear_ge(SD, geas::at_True, sum0, parcons, vars0, 0);
  geas::bool_linear_le(SD, geas::at_True, sum1, parcons, vars1, 0);
  geas::bool_linear_ge(SD, geas::at_True, sum1, parcons, vars1, 0); 
}

void d_bool_lin_eql(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  FUNCTIONNOTIMPLEMENT; 
  geas::intvar sum0 = SOL.new_intvar(SHRT_MIN, SHRT_MAX); 
  geas::intvar sum1 = SOL.new_intvar(SHRT_MIN, SHRT_MAX);
  
  vec<int> cons = INTARRAY(0);
  std::vector<Id*> varIds = VARIDARRAY(1);
  assert(varIds.size() == cons.size()); 
  vec<geas::patom_t> vars0, vars1; 
  bool allpos = true, allneg = true; 
  vec<int> parcons; 
  for (unsigned int i = 0; i != varIds.size(); i++) {
    if (varIds[i] != nullptr && fixedVars.find(varIds[i]) != fixedVars.end()) {
      vars0.push(SI.asBoolVarDom(varIds[i], true)); 
      vars1.push(SI.asBoolVarDom(varIds[i], false)); 
      parcons.push(cons[i]);
      allpos = allpos && (cons[i] >= 0);
      allneg = allneg && (cons[i] <= 0);
    }
  }

  geas::bool_linear_le(SD, geas::at_True, sum0, parcons, vars0, 0);
  geas::bool_linear_ge(SD, geas::at_True, sum0, parcons, vars0, 0);
  geas::bool_linear_le(SD, geas::at_True, sum1, parcons, vars1, 0);
  geas::bool_linear_ge(SD, geas::at_True, sum1, parcons, vars1, 0); 
  
  if (allpos) {
    SOL.post(sum0 <= INT(2)); 
    SOL.post(sum1 <= INT(2)); 
  } else if (allneg) {
    SOL.post(sum0 >= INT(2)); 
    SOL.post(sum1 >= INT(2)); 
  }
  geas::int_eq(SD, sum0, sum1);
}

void d_bool_lin_le(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  FUNCTIONNOTIMPLEMENT; 
  vec<int> cons = INTARRAY(0);
  std::vector<Id*> varIds = VARIDARRAY(1);
  assert(varIds.size() == cons.size()); 
  vec<geas::patom_t> vars0, vars1; 
  bool allpos = true, allneg = true; 
  vec<int> parcons; 
  for (unsigned int i = 0; i != varIds.size(); i++) {
    if (varIds[i] != nullptr && fixedVars.find(varIds[i]) != fixedVars.end()) {
      vars0.push(SI.asBoolVarDom(varIds[i], true)); 
      vars1.push(SI.asBoolVarDom(varIds[i], false)); 
      parcons.push(cons[i]);
      allpos = allpos && (cons[i] >= 0);
      allneg = allneg && (cons[i] <= 0);
    }
  }

  geas::intvar sum0 = SOL.new_intvar(SHRT_MIN, SHRT_MAX); 
  geas::intvar sum1 = SOL.new_intvar(SHRT_MIN, SHRT_MAX);
  
  geas::bool_linear_ge(SD, geas::at_True, sum0, parcons, vars0, 0);
  geas::bool_linear_le(SD, geas::at_True, sum1, parcons, vars1, 0);
  geas::int_le(SD, sum0, sum1, 0); 

  SI._sensvar.push_back(sum0);
  SI._sensvar.push_back(sum1);
}

void d_bool_lin_eql_reif(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  CONSTRAINTNOTIMPLEMENT; 
  geas::intvar v0 = SOL.new_intvar(SHRT_MIN, SHRT_MAX); 
  geas::intvar v1 = SOL.new_intvar(SHRT_MIN, SHRT_MAX);
  bool_lin_partial_sum(s, call, fixedVars, v0, v1); 
  geas::int_eq(SD, v0, v1);
}

void d_bool_lin_le_reif(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  if (PAR(3) && BOOL(3))
    d_bool_lin_le(s, call, fixedVars); 
  else if (PAR(3) && !BOOL(3))
    d_bool_lin_eql(s, call, fixedVars); 
  else if (!PAR(3)) {
    Monotonicity mono = VAR_EQL; 
    if (call->ann().containsCall(std::string("defines_var")))
      mono = SI._variableMono[SI.asVarId(call->ann().getCall(std::string("defines_var"))->arg(0))]; 
    else if ( VARID(3) != nullptr && fixedVars.find(VARID(3)) != fixedVars.end() ) {
      // enforce equality of boolean variable 
      geas::add_clause(SD, SI.asBoolVarDom(VARID(3), true), ~SI.asBoolVarDom(VARID(3), false));
      geas::add_clause(SD, ~SI.asBoolVarDom(VARID(3), true), SI.asBoolVarDom(VARID(3), false)); 
    }
    geas::intvar v0 = SOL.new_intvar(SHRT_MIN, SHRT_MAX); 
    geas::intvar v1 = SOL.new_intvar(SHRT_MIN, SHRT_MAX);
    bool_lin_partial_sum(s, call, fixedVars, v0, v1); 

    if (mono == VAR_INC)
      geas::int_le(SD, v0, v1, 0);
    else if (mono == VAR_DEC)
      geas::int_le(SD, v1, v0, 0);
    else if (mono == VAR_EQL)
      geas::int_eq(SD, v0, v1); 
  } else 
    throw InternalError(std::string("Unknown structure for constraint ") + std::string(call->id().c_str()) );
}

void d_no_dominance(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  throw InternalError(std::string("no dominance condition for unary functional constraint: ") + std::string(call->id().c_str()) ); 
}

void d_all_different(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  FUNCTIONNOTIMPLEMENT; 
  std::vector<Id*> args; 
  std::vector<Id*> vars = VARIDARRAY(0);
  std::unordered_set<int> vals; 
  for (unsigned int i = 0; i != vars.size(); i++) {
    if (vars[i] != nullptr) {
      if (fixedVars.find(vars[i]) == fixedVars.end()) {
        Expression* domain = vars[i]->decl()->ti()->domain(); 
        IntSetVal* isv = eval_intset(SI.env().envi(), domain);
        if (isv->size() >= 1) {
          for (int j = 0; j < isv->size(); ++j)
            for (auto k = isv->min(j).toInt(); k <= isv->max(j).toInt(); ++k)
              vals.insert(static_cast<int>(k)); 
        }
      } else
        args.push_back(vars[i]); 
    }
  }

  if (args.size() != 0) {
    for (auto& val: vals) {
      geas::intvar c0 = SOL.new_intvar(0,1);
      geas::intvar c1 = SOL.new_intvar(0,1);
      vec<int> ks(args.size(), 1); 
      vec<geas::patom_t> bv0, bv1; 
      for (auto& id: args) {
        bv0.push(SI.asIntVarDom(id, true) == val);
        bv1.push(SI.asIntVarDom(id, false) == val);
      }
      geas::bool_linear_le(SD, geas::at_True, c0, ks, bv0, 0); 
      geas::bool_linear_le(SD, geas::at_True, c1, ks, bv1, 0); 
      geas::bool_linear_ge(SD, geas::at_True, c0, ks, bv0, 0); 
      geas::bool_linear_ge(SD, geas::at_True, c1, ks, bv1, 0); 
      geas::int_le(SD, c0, c1, 0); 
    }
  }
}

void d_all_different_except_0(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  FUNCTIONNOTIMPLEMENT; 
  std::vector<Id*> args; 
  std::vector<Id*> vars = VARIDARRAY(0);
  std::unordered_set<int> vals; 
  for (unsigned int i = 0; i != vars.size(); i++) {
    if (vars[i] != nullptr) {
      if (fixedVars.find(vars[i]) == fixedVars.end()) {
        Expression* domain = vars[i]->decl()->ti()->domain(); 
        IntSetVal* isv = eval_intset(SI.env().envi(), domain);
        if (isv->size() >= 1) {
          for (int j = 0; j < isv->size(); ++j)
            for (auto k = isv->min(j).toInt(); k <= isv->max(j).toInt(); ++k)
              if (k != 0)
                vals.insert(static_cast<int>(k)); 
        }
      } else
        args.push_back(vars[i]); 
    }
  }

  if (args.size() != 0) {
    for (auto& val: vals) {
      geas::intvar c0 = SOL.new_intvar(0,1);
      geas::intvar c1 = SOL.new_intvar(0,1);
      vec<int> ks(args.size(), 1); 
      vec<geas::patom_t> bv0, bv1; 
      for (auto& id: args) {
        bv0.push(SI.asIntVarDom(id, true) == val);
        bv1.push(SI.asIntVarDom(id, false) == val);
      }
      geas::bool_linear_le(SD, geas::at_True, c0, ks, bv0, 0); 
      geas::bool_linear_le(SD, geas::at_True, c1, ks, bv1, 0); 
      geas::bool_linear_ge(SD, geas::at_True, c0, ks, bv0, 0); 
      geas::bool_linear_ge(SD, geas::at_True, c1, ks, bv1, 0); 
      geas::int_le(SD, c0, c1, 0); 
    }
  }
}

void d_array_int_minimum(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  Monotonicity mono = VAR_EQL; 
  if ( call->ann().containsCall(std::string("defines_var")) )
    mono = SI._variableMono[SI.asVarId(call->ann().getCall(std::string("defines_var"))->arg(0))]; 
  else if (VARID(0) != nullptr && fixedVars.find(VARID(0)) != fixedVars.end() ) // enforce equality for the first argument 
    geas::int_eq(SD, SI.asIntVarDom(VARID(0), true), SI.asIntVarDom(VARID(0), false)); 

  std::vector<Id*> vars = VARIDARRAY(1); 
  vec<geas::intvar> vars0, vars1; 
  for (unsigned int i = 0; i != vars.size(); i++) {
    if (vars[i] != nullptr && fixedVars.find(vars[i]) != fixedVars.end()) {
      vars0.push(-SI.asIntVarDom(vars[i], true)); 
      vars1.push(-SI.asIntVarDom(vars[i], false));
    }
  }

  if (vars0.size() != 0) {
    geas::intvar r0 = SOL.new_intvar(SHRT_MIN, SHRT_MAX);
    geas::intvar r1 = SOL.new_intvar(SHRT_MIN, SHRT_MAX); 
    geas::int_max(SD, -r0, vars0); 
    geas::int_max(SD, -r1, vars1);
    
    if (mono == VAR_INC)
      geas::int_le(SD, r1, r0, 0);
    else if (mono == VAR_DEC)
      geas::int_le(SD, r0, r1, 0);
    else if (mono == VAR_EQL) 
      geas::int_eq(SD, r0, r1);  
  } 
}

void d_array_int_maximum(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  Monotonicity mono = VAR_EQL; 
  if ( call->ann().containsCall(std::string("defines_var")) )
    mono = SI._variableMono[SI.asVarId(call->ann().getCall(std::string("defines_var"))->arg(0))]; 
  else if (VARID(0) != nullptr && fixedVars.find(VARID(0)) != fixedVars.end()) // enforce equality for the first argument 
    geas::int_eq(SD, SI.asIntVarDom(VARID(0), true), SI.asIntVarDom(VARID(0), false)); 

  std::vector<Id*> vars = VARIDARRAY(1); 
  vec<geas::intvar> vars0, vars1; 
  for (unsigned int i = 0; i != vars.size(); i++) {
    if (vars[i] != nullptr && fixedVars.find(vars[i]) != fixedVars.end()) {
      vars0.push(SI.asIntVarDom(vars[i], true)); 
      vars1.push(SI.asIntVarDom(vars[i], false));
    }
  }
  
  if (vars0.size() != 0) {
    // todo: need to compute the bound 
    geas::intvar r0 = SOL.new_intvar(SHRT_MIN, SHRT_MAX);
    geas::intvar r1 = SOL.new_intvar(SHRT_MIN, SHRT_MAX); 
    geas::int_max(SD, r0, vars0); 
    geas::int_max(SD, r1, vars1);
    
    if (mono == VAR_INC)
      geas::int_le(SD, r1, r0, 0);
    else if (mono == VAR_DEC)
      geas::int_le(SD, r0, r1, 0);
    else if (mono == VAR_EQL) 
      geas::int_eq(SD, r0, r1); 
  }   
}

void d_table_int(SolverInstanceBase& s, const Call* call, const std::unordered_set<Id*>& fixedVars) {
  FUNCTIONNOTIMPLEMENT; 
  std::vector<Id*> vars = VARIDARRAY(0); 
  for (unsigned int i = 0; i != vars.size(); i++) {
    if (vars[i] != nullptr && fixedVars.find(vars[i]) != fixedVars.end()) 
      geas::int_eq(SD, SI.asIntVarDom(vars[i], true), SI.asIntVarDom(vars[i], false)); 
  }
}

}  // namespace GeasConstraints
}  // namespace MiniZinc

#pragma clang diagnostic pop
