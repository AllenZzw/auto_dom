/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */

/*
 *  Main authors:
 *     Jip J. Dekker <jip.dekker@monash.edu>
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include <minizinc/ast.hh>
#include <minizinc/solver_instance_base.hh>

namespace MiniZinc {
namespace GeasConstraints {

#define PosterImpl(X) void X(SolverInstanceBase& s, const Call* ce)
#define AnalyzerImpl(X) void X(SolverInstanceBase& s, Call* ce)
#define DominanceImpl(X) void X(SolverInstanceBase& s, const Call* ce, const std::unordered_set<Id*>& fixedVars)

AnalyzerImpl(a_set_in); 

PosterImpl(p_set_in); 

/* Integer Comparison Constraints */
PosterImpl(p_int_eq);
PosterImpl(p_int_ne);
PosterImpl(p_int_le);
PosterImpl(p_int_lt);
PosterImpl(p_int_eq_imp);
PosterImpl(p_int_ne_imp);
PosterImpl(p_int_le_imp);
PosterImpl(p_int_lt_imp);
PosterImpl(p_int_eq_reif);
PosterImpl(p_int_ne_reif);
PosterImpl(p_int_le_reif);
PosterImpl(p_int_lt_reif);

AnalyzerImpl(a_eql); 
AnalyzerImpl(a_let);  
AnalyzerImpl(a_eql_reif); 
AnalyzerImpl(a_let_reif); 

DominanceImpl(d_int_eql); 
DominanceImpl(d_int_let); 
DominanceImpl(d_int_eql_reif); 
DominanceImpl(d_int_let_reif); 

/* Integer Arithmetic Constraints */
PosterImpl(p_int_abs);
PosterImpl(p_int_times);
PosterImpl(p_int_div);
PosterImpl(p_int_max);
PosterImpl(p_int_min);

AnalyzerImpl(a_int_abs); 
AnalyzerImpl(a_eql_binary_op); 
AnalyzerImpl(a_inc_binary_op); 

DominanceImpl(d_int_binary_op); 

/* Integer Linear Constraints */
PosterImpl(p_int_lin_eq);
PosterImpl(p_int_lin_ne);
PosterImpl(p_int_lin_le);
PosterImpl(p_int_lin_eq_imp);
PosterImpl(p_int_lin_ne_imp);
PosterImpl(p_int_lin_le_imp);
PosterImpl(p_int_lin_eq_reif);
PosterImpl(p_int_lin_ne_reif);
PosterImpl(p_int_lin_le_reif);
AnalyzerImpl(a_lin_eq);
AnalyzerImpl(a_lin_ne);
AnalyzerImpl(a_lin_le);
AnalyzerImpl(a_lin_eql_reif);
// AnalyzerImpl(a_lin_le_imp);
AnalyzerImpl(a_lin_le_reif);

DominanceImpl(d_int_lin_eq); 
DominanceImpl(d_int_lin_eql); 
DominanceImpl(d_int_lin_le);
DominanceImpl(d_int_lin_eql_reif);
// DominanceImpl(a_lin_le_imp);
DominanceImpl(d_int_lin_le_reif);


/* Boolean Comparison Constraints */
PosterImpl(p_bool_eq);
PosterImpl(p_bool_ne);
PosterImpl(p_bool_le);
PosterImpl(p_bool_lt);
PosterImpl(p_bool_eq_imp);
PosterImpl(p_bool_ne_imp);
PosterImpl(p_bool_le_imp);
PosterImpl(p_bool_lt_imp);
PosterImpl(p_bool_eq_reif);
PosterImpl(p_bool_ne_reif);
PosterImpl(p_bool_le_reif);
PosterImpl(p_bool_lt_reif);

DominanceImpl(d_bool_eql);
DominanceImpl(d_bool_let);
DominanceImpl(d_bool_eql_reif);
DominanceImpl(d_bool_let_reif);

/* Boolean Arithmetic Constraints */
PosterImpl(p_bool_or);
PosterImpl(p_bool_and);
PosterImpl(p_bool_xor);
PosterImpl(p_bool_not);
PosterImpl(p_bool_or_imp);
PosterImpl(p_bool_and_imp);
PosterImpl(p_bool_xor_imp);

DominanceImpl(d_bool_binary_op); 

PosterImpl(p_bool_clause);
PosterImpl(p_array_bool_or);
PosterImpl(p_array_bool_and);
PosterImpl(p_bool_clause_imp);
PosterImpl(p_array_bool_or_imp);
PosterImpl(p_array_bool_and_imp);
PosterImpl(p_bool_clause_reif);

AnalyzerImpl(a_bool_not);
AnalyzerImpl(a_bool_clause);
AnalyzerImpl(a_array_bool_and_or);
AnalyzerImpl(a_bool_clause_reif);

DominanceImpl(d_bool_clause);
DominanceImpl(d_array_bool_and);
DominanceImpl(d_array_bool_or);
DominanceImpl(d_bool_clause_reif);

/* Boolean Linear Constraints */
PosterImpl(p_bool_lin_eq);
PosterImpl(p_bool_lin_ne);
PosterImpl(p_bool_lin_le);
PosterImpl(p_bool_lin_eq_imp);
PosterImpl(p_bool_lin_ne_imp);
PosterImpl(p_bool_lin_le_imp);
PosterImpl(p_bool_lin_lt_imp);
PosterImpl(p_bool_lin_eq_reif);
PosterImpl(p_bool_lin_ne_reif);
PosterImpl(p_bool_lin_le_reif);

DominanceImpl(d_bool_lin_eq); 
DominanceImpl(d_bool_lin_eql); 
DominanceImpl(d_bool_lin_le);
DominanceImpl(d_bool_lin_eql_reif);
DominanceImpl(d_bool_lin_le_reif);

/* Coercion Constraints */
PosterImpl(p_bool2int);
AnalyzerImpl(a_bool2int); 

/* Element Constraints */
PosterImpl(p_array_int_element);
PosterImpl(p_array_bool_element);
PosterImpl(p_array_var_int_element);
PosterImpl(p_array_var_bool_element);

AnalyzerImpl(a_array_lit_element);

/* Global Constraints */
PosterImpl(p_all_different);
PosterImpl(p_all_different_except_0);
PosterImpl(p_array_int_maximum);
PosterImpl(p_array_int_minimum);
PosterImpl(p_at_most);
PosterImpl(p_at_most1);
PosterImpl(p_cumulative);
PosterImpl(p_disjunctive);
PosterImpl(p_global_cardinality);
PosterImpl(p_table_int);

AnalyzerImpl(a_all_different);
AnalyzerImpl(a_array_int_min_max); 
AnalyzerImpl(a_table_int);

DominanceImpl(d_all_different);
DominanceImpl(d_all_different_except_0);
DominanceImpl(d_array_int_minimum);
DominanceImpl(d_array_int_maximum);
DominanceImpl(d_table_int);

DominanceImpl(d_no_dominance);

}  // namespace GeasConstraints
}  // namespace MiniZinc