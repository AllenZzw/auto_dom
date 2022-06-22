/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */

/*
 *  Main authors:
 *     Jip J. Dekker <jip.dekker@monash.edu>
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <minizinc/solvers/geas/geas_constraints.hh>
#include <minizinc/solvers/geas_solverinstance.hh>

#include <geas/constraints/builtins.h>
#include <list>

namespace MiniZinc {
void Preprocessor::add(const ASTString name, analyzer p) { _preprocessor.insert(std::make_pair(name, p)); }
void Preprocessor::add(const std::string& name, analyzer p) {
  GCLock lock;
  ASTString str(name);
  return add(str, p);
}
void Preprocessor::post(Call* c) {
  auto it = _preprocessor.find(c->id());
  if (it == _preprocessor.end()) {
    std::ostringstream ss;
    ss << "Error: preprocessor cannot handle constraint: " << c->id();
    throw InternalError(ss.str());
  }
  return it->second(_base, c);
}

void DomPoster::add(const ASTString name, dominator p) { _domposter.insert(std::make_pair(name, p)); }
void DomPoster::add(const std::string& name, dominator p) {
  GCLock lock;
  ASTString str(name);
  return add(str, p);
}
void DomPoster::post(Call* c, const std::unordered_set<Id*>& fixedVars) {
  auto it = _domposter.find(c->id());
  if (it == _domposter.end()) {
    std::ostringstream ss;
    ss << "Error: preprocessor cannot handle constraint: " << c->id();
    throw InternalError(ss.str());
  }
  it->second(_base, c, fixedVars);
}

GeasSolverInstance::GeasSolverInstance(Env& env, std::ostream& log,
                                       SolverInstanceBase::Options* opt)
    : SolverInstanceImpl<GeasTypes>(env, log, opt), _flat(env.flat()), _constraintPreprocessor(*this), _domPoster(*this) {
  registerConstraints();
}

void GeasSolverInstance::registerConstraint(const std::string& name, analyzer a, poster p, dominator d) {
  _constraintPreprocessor.add("dom_" + name, a);
  _constraintPreprocessor.add(name, a);
  
  _constraintRegistry.add("dom_" + name, p);
  _constraintRegistry.add(name, p);
  
  _domPoster.add("dom_" + name, d);
  _domPoster.add(name, d);
}

void GeasSolverInstance::registerConstraints() {
  GCLock lock;
  registerConstraint("set_in", GeasConstraints::a_set_in, GeasConstraints::p_set_in, GeasConstraints::d_no_dominance);

  /* Integer Comparison Constraints */
  registerConstraint("int_eq", GeasConstraints::a_eql, GeasConstraints::p_int_eq, GeasConstraints::d_int_eql);
  registerConstraint("int_ne", GeasConstraints::a_eql, GeasConstraints::p_int_ne, GeasConstraints::d_int_eql);
  registerConstraint("int_le", GeasConstraints::a_let, GeasConstraints::p_int_le, GeasConstraints::d_int_let);
  registerConstraint("int_lt", GeasConstraints::a_let, GeasConstraints::p_int_lt, GeasConstraints::d_int_let);
  registerConstraint("int_eq_reif", GeasConstraints::a_eql_reif, GeasConstraints::p_int_eq_reif, GeasConstraints::d_int_eql_reif);
  registerConstraint("int_ne_reif", GeasConstraints::a_eql_reif, GeasConstraints::p_int_ne_reif, GeasConstraints::d_int_eql_reif);
  registerConstraint("int_le_reif", GeasConstraints::a_let_reif, GeasConstraints::p_int_le_reif, GeasConstraints::d_int_let_reif);
  registerConstraint("int_lt_reif", GeasConstraints::a_let_reif, GeasConstraints::p_int_lt_reif, GeasConstraints::d_int_let_reif);

  /* Integer Arithmetic Constraints */
  registerConstraint("int_abs", GeasConstraints::a_int_abs, GeasConstraints::p_int_abs, GeasConstraints::d_no_dominance);
  registerConstraint("int_times", GeasConstraints::a_eql_binary_op, GeasConstraints::p_int_times, GeasConstraints::d_int_binary_op);
  // registerConstraint("int_pow", GeasConstraints::a_eql_binary_op, GeasConstraints::p_int_pow, GeasConstraints::d_int_binary_op);
  // registerConstraint("int_plus", GeasConstraints::a_eql_binary_op, GeasConstraints::p_int_plus, GeasConstraints::d_int_binary_op);
  // registerConstraint("int_mod", GeasConstraints::a_eql_binary_op, GeasConstraints::p_int_mod, GeasConstraints::d_int_binary_op);
  registerConstraint("int_div", GeasConstraints::a_eql_binary_op, GeasConstraints::p_int_div, GeasConstraints::d_int_binary_op);
  registerConstraint("int_min", GeasConstraints::a_inc_binary_op, GeasConstraints::p_int_min, GeasConstraints::d_int_binary_op);
  registerConstraint("int_max", GeasConstraints::a_inc_binary_op, GeasConstraints::p_int_max, GeasConstraints::d_int_binary_op);

  /* Integer Linear Constraints */
  registerConstraint("int_lin_eq", GeasConstraints::a_lin_eq, GeasConstraints::p_int_lin_eq, GeasConstraints::d_int_lin_eq);
  registerConstraint("int_lin_ne", GeasConstraints::a_lin_ne, GeasConstraints::p_int_lin_ne, GeasConstraints::d_int_lin_eql);
  registerConstraint("int_lin_le", GeasConstraints::a_lin_le, GeasConstraints::p_int_lin_le, GeasConstraints::d_int_lin_le);
  registerConstraint("int_lin_eq_reif", GeasConstraints::a_lin_eql_reif, GeasConstraints::p_int_lin_eq_reif, GeasConstraints::d_int_lin_eql_reif);
  registerConstraint("int_lin_ne_reif", GeasConstraints::a_lin_eql_reif, GeasConstraints::p_int_lin_ne_reif, GeasConstraints::d_int_lin_eql_reif);
  registerConstraint("int_lin_le_reif", GeasConstraints::a_lin_le_reif, GeasConstraints::p_int_lin_le_reif, GeasConstraints::d_int_lin_le_reif);

  /* Boolean Comparison Constraints todo: handle Boolean reification constraint */ 
  registerConstraint("bool_eq", GeasConstraints::a_eql, GeasConstraints::p_bool_eq, GeasConstraints::d_bool_eql);
  registerConstraint("bool_ne", GeasConstraints::a_eql, GeasConstraints::p_bool_ne, GeasConstraints::d_bool_eql);
  registerConstraint("bool_le", GeasConstraints::a_let, GeasConstraints::p_bool_le, GeasConstraints::d_bool_let);
  registerConstraint("bool_lt", GeasConstraints::a_let, GeasConstraints::p_bool_lt, GeasConstraints::d_bool_let);
  registerConstraint("bool_eq_reif", GeasConstraints::a_eql_reif, GeasConstraints::p_bool_eq_reif, GeasConstraints::d_bool_eql_reif);
  registerConstraint("bool_ne_reif", GeasConstraints::a_eql_reif, GeasConstraints::p_bool_ne_reif, GeasConstraints::d_bool_eql_reif);
  registerConstraint("bool_le_reif", GeasConstraints::a_let_reif, GeasConstraints::p_bool_le_reif, GeasConstraints::d_bool_let_reif);
  registerConstraint("bool_lt_reif", GeasConstraints::a_let_reif, GeasConstraints::p_bool_lt_reif, GeasConstraints::d_bool_let_reif);

  /* Boolean Arithmetic Constraints todo: handle Boolean arithmetic */ 
  registerConstraint("bool_or", GeasConstraints::a_inc_binary_op, GeasConstraints::p_bool_or, GeasConstraints::d_bool_binary_op);
  registerConstraint("bool_and", GeasConstraints::a_inc_binary_op, GeasConstraints::p_bool_and, GeasConstraints::d_bool_binary_op);
  registerConstraint("bool_xor", GeasConstraints::a_eql_binary_op, GeasConstraints::p_bool_xor, GeasConstraints::d_bool_binary_op);
  registerConstraint("bool_not", GeasConstraints::a_bool_not, GeasConstraints::p_bool_not, GeasConstraints::d_no_dominance);
  
  registerConstraint("bool_clause", GeasConstraints::a_bool_clause, GeasConstraints::p_bool_clause, GeasConstraints::d_bool_clause);
  registerConstraint("array_bool_or", GeasConstraints::a_array_bool_and_or, GeasConstraints::p_array_bool_or, GeasConstraints::d_array_bool_or);
  registerConstraint("array_bool_and", GeasConstraints::a_array_bool_and_or, GeasConstraints::p_array_bool_and, GeasConstraints::d_array_bool_and);
  registerConstraint("bool_clause_reif", GeasConstraints::a_bool_clause_reif, GeasConstraints::p_bool_clause_reif, GeasConstraints::d_bool_clause_reif);

  // /* Boolean Linear Constraints todo: handle Boolean linear constraints */
  registerConstraint("bool_lin_eq", GeasConstraints::a_lin_eq, GeasConstraints::p_bool_lin_eq, GeasConstraints::d_bool_lin_eql);
  registerConstraint("bool_lin_ne", GeasConstraints::a_lin_ne, GeasConstraints::p_bool_lin_ne, GeasConstraints::d_bool_lin_eql);
  registerConstraint("bool_lin_le", GeasConstraints::a_lin_le, GeasConstraints::p_bool_lin_le, GeasConstraints::d_bool_lin_le);
  registerConstraint("bool_lin_eq_reif", GeasConstraints::a_lin_eql_reif, GeasConstraints::p_bool_lin_eq_reif, GeasConstraints::d_bool_lin_eql_reif);
  registerConstraint("bool_lin_ne_reif", GeasConstraints::a_lin_eql_reif, GeasConstraints::p_bool_lin_ne_reif, GeasConstraints::d_bool_lin_eql_reif);
  registerConstraint("bool_lin_le_reif", GeasConstraints::a_lin_le_reif, GeasConstraints::p_bool_lin_le_reif, GeasConstraints::d_bool_lin_le_reif);

  /* Coercion Constraints */
  registerConstraint("bool2int", GeasConstraints::a_bool2int, GeasConstraints::p_bool2int, GeasConstraints::d_no_dominance);

  /* Element Constraints */
  registerConstraint("array_int_element", GeasConstraints::a_array_lit_element, GeasConstraints::p_array_int_element, GeasConstraints::d_no_dominance);
  registerConstraint("array_bool_element", GeasConstraints::a_array_lit_element, GeasConstraints::p_array_bool_element, GeasConstraints::d_no_dominance);
  // registerConstraint("array_var_int_element", GeasConstraints::a_array_lit_element, GeasConstraints::p_array_var_int_element, GeasConstraints::d_var_int_element);
  // registerConstraint("array_var_bool_element", GeasConstraints::a_array_lit_element, GeasConstraints::p_array_var_bool_element, GeasConstraints::d_var_bool_element);

  // /* Global Constraints */
  registerConstraint("all_different_int", GeasConstraints::a_all_different, GeasConstraints::p_all_different, GeasConstraints::d_all_different);
  registerConstraint("all_different_except_0", GeasConstraints::a_all_different, GeasConstraints::p_all_different_except_0, GeasConstraints::d_all_different_except_0);
  registerConstraint("array_int_minimum", GeasConstraints::a_array_int_min_max, GeasConstraints::p_array_int_minimum, GeasConstraints::d_array_int_minimum);
  registerConstraint("array_int_maximum", GeasConstraints::a_array_int_min_max, GeasConstraints::p_array_int_maximum, GeasConstraints::d_array_int_maximum);
  registerConstraint("table_int", GeasConstraints::a_table_int, GeasConstraints::p_table_int, GeasConstraints::d_table_int);
}


// todo: factor out the analyzer to solver_instance_base 
// compute the monotonicity of each internal variables and the dependency of variables 
void GeasSolverInstance::processFlatZinc() {
  // constraints or the constraint defining the objective variable 
  std::vector<Expression*> _constriant; 
  std::unordered_set<Id*> _varSet;
  std::map<Id*, Call*> _variableDefiningExp; // map variable id to experssion that defines the variable 
   
  // initialize the Monotonicity of variables 
  for (auto it = _flat->vardecls().begin(); it != _flat->vardecls().end(); ++it)
    if (!it->removed() && it->e()->type().isvar() && it->e()->type().dim() == 0 ) {
      VarDecl* vd = it->e(); 
      Id* id = resolveVarId(it->e()); 
      // test whether the variable is a constant 
      if (vd->type().isbool()) {
        if (vd->e() == nullptr)
          _varSet.insert(id); 
        else {
          Expression* init = vd->e();
          if (init->isa<Id>() || init->isa<ArrayAccess>()) {
            // alias boolean variable 
            Id* aliasId = resolveVarId(init); 
            _alias[id] = aliasId; 
          } else 
            _alias[id] = init->cast<BoolLit>(); 
        }
      } else if (vd->type().isint()) {
        if (vd->e() == nullptr) {
          Expression* domain = vd->ti()->domain();
          if (domain != nullptr) {
            IntSetVal* isv = eval_intset(_env.envi(), domain);
            if (isv->min().toInt() == isv->max().toInt())
              _alias[id] = IntLit::a(isv->min().toInt()); 
            else
              _varSet.insert(id); 
          } else {
            _varSet.insert(id); 
            // std::ostringstream ss;
            // ss << "GeasSolverInstance::processFlatZinc: Error: Unbounded variable: " << vd->id()->str();
            // throw Error(ss.str());
          }
        } else {
          // alias integer variable 
          Expression* init = vd->e();
          if (init->isa<Id>() || init->isa<ArrayAccess>())
            _alias[id] = resolveVarId(init); 
          else
            _alias[id] = init->cast<IntLit>();
        }
      }
    }

  // record the objective variable
  SolveI* si = _flat->solveItem();
  _objVar = nullptr; 
  if (si->e() != nullptr) {
    _objVar = resolveVarId(si->e()); 
    if (_alias.count(_objVar) != 0) _objVar = resolveVarId(_alias[_objVar]); 
    _variableMono[_objVar] = si->st()==SolveI::ST_MIN?VAR_DEC:VAR_INC; 
    _variableSensitive[_objVar] = true; 
  } 

  // put constriants into different category 
  for (ConstraintIterator it = _flat->constraints().begin(); it != _flat->constraints().end(); ++it) {
    if (!it->removed()) {
      if (auto* c = it->e()->dynamicCast<Call>()) {
        replaceCallAlias(c); 
        if (c->ann().containsCall(std::string("defines_var"))) {
          replaceCallAlias(c->ann().getCall(std::string("defines_var"))); 
          Id* defVar = asVarId(c->ann().getCall(std::string("defines_var"))->arg(0)); 
          if (defVar == nullptr) {
            c->ann().removeCall(std::string("defines_var")); 
            _constriant.push_back(c); 
          } else {
            _variableDefiningExp[defVar] = c;
            _varSet.erase(defVar); 
          }
        } else 
          _constriant.push_back(c); 
        
        std::set<Id*> args = getArguments(c); 
        _constraintArgNr[c] = args.size();
        for (auto vit = args.begin(); vit != args.end(); vit++)
          _argsToConstraint[*vit].insert(c); 
      }
    }
  }

  // iterate over set of ids that have an output annotation 
  // obtain their right hand side from the flat model
  if (_varsWithOutput.empty()) {
    for (VarDeclIterator it = getEnv()->flat()->vardecls().begin(); it != getEnv()->flat()->vardecls().end(); ++it) {
      if (!it->removed()) {
        VarDecl* vd = it->e();
        if ( Call* output_array_ann = Expression::dynamicCast<Call>(get_annotation(vd->ann(), std::string("output_array"))) ) {
          auto* al = vd->e()->dynamicCast<ArrayLit>(); 
          std::vector<Expression*> array_elems;
          ArrayLit& array = *al; 
          std::vector<std::pair<int, int> > dims_v;
          std::vector<int> cur_idx;
          ArrayLit* dims;
          Expression* e = output_array_ann->arg(0);
          if (auto* al = e->dynamicCast<ArrayLit>())
            dims = al;
          else if (Id* id = e->dynamicCast<Id>())
            dims = id->decl()->e()->cast<ArrayLit>();
          else
            throw -1;
          for (int i = 0; i < dims->length(); i++) {
            IntSetVal* isv = eval_intset(getEnv()->envi(), (*dims)[i]);
            if (isv->size() == 0)
              throw InternalError("zero dimension for array variable");
            else {
              cur_idx.emplace_back(static_cast<int>(isv->min().toInt())); 
              dims_v.emplace_back(static_cast<int>(isv->min().toInt()),static_cast<int>(isv->max().toInt()));
            }
          }
          // iterate over the array and record the name 
          replaceArrayAlias(al); 
          for (unsigned int j = 0; j < array.size(); j++) {
            if (array[j]->dynamicCast<Id>() != nullptr) {
              Id* id = resolveVarId(array[j]->dynamicCast<Id>());
              if (_varSet.count(id) != 0 && _variableName.find(id) == _variableName.end()) {
                std::ostringstream ss;
                ss << vd->id()->str() << "["; 
                for (unsigned int i = 0; i != cur_idx.size(); i++) {
                  ss << cur_idx[i]; 
                  if (i+1 != cur_idx.size()) 
                    ss << ","; 
                }
                ss << "]";
                _variableName[id] = ss.str();
                _variableMono[id] = VAR_NONE; 
                _variableSensitive[id] = false;
                _indepVars.push_back(id);
              }
            }
            // increment cur_idx 
            if (j+1 != array.size()) {
              for (unsigned int i = cur_idx.size()-1; i >= 0; i--) {
                cur_idx[i]++; 
                if (cur_idx[i] > dims_v[i].second)
                  cur_idx[i] = dims_v[i].first; 
                else 
                  break; 
              }
            }
          }
        } else if ( get_annotation(vd->ann(), std::string("output_var")) != nullptr ) {
          Id* id = resolveVarId(vd->id()); 
          if (_varSet.count(id) != 0 && _variableName.find(id) == _variableName.end()) {
            _variableName[id] = std::string(id->str().c_str());
            _variableMono[id] = VAR_NONE; 
            _variableSensitive[id] = false;
            _indepVars.push_back(id);
          }
            
        }
      }
    }
  } else 
    throw InternalError(std::string("Not Implemented: Cannot handle MiniZinc Model with output items"));
  
  // analyze the monotonicity of auxiliary variables 
  std::list<Id*> depVars; 
  if (_objVar != nullptr)
    depVars.push_back(_objVar); 
  
  // set monotonicity for variable with domain as constraint 
  for (auto it = _variableDefiningExp.begin(); it != _variableDefiningExp.end(); it++)
    if (_argsToConstraint[it->first].size() == 0 && it->first != _objVar) {
      _variableMono[it->first] = VAR_EQL; 
      depVars.push_back(it->first); 
    }
  
  for (auto it = _constriant.begin(); it != _constriant.end(); it++ ) {
    Call* c = (*it)->dynamicCast<Call>(); 
    std::set<Id*> args = getArguments(c); 
    _constraintPreprocessor.post(c); 
    depVars.insert(depVars.end(), args.begin(), args.end()); 
  }

  // set the monotonicity of variable recursively for dependent variables 
  while (!depVars.empty()) {
    Id* curVar = depVars.front(); 
    if (_variableDefiningExp.count(curVar) != 0) {
      std::set<Id*> args = getArguments(_variableDefiningExp[curVar]);
      _constraintPreprocessor.post(_variableDefiningExp[curVar]); 
      _variableDefiningExp.erase(curVar); 
      depVars.insert(depVars.end(), args.begin(), args.end());
    }
    depVars.pop_front();
  }

  // for (auto mit = _variableName.begin(); mit != _variableName.end(); mit++)
  //   std::cout << mit->first->str() << ": " << mit->second << std::endl; 
  // for (auto it = _indepVars.begin(); it != _indepVars.end(); it++ ) 
  //   std::cout << _variableName[*it] << " "; 
  // std::cout << std::endl; 
  // for (auto mit = _variableMono.begin(); mit != _variableMono.end(); mit++)
  //   std::cout << mit->first->str() << ": " << mit->second << std::endl; 
  // for (auto it = _indepVars.begin(); it != _indepVars.end(); it++ ) 
  //   std::cout << (*it)->str() << ": " << _variableName[*it] << std::endl;
  // for (auto mit = _alias.begin(); mit != _alias.end(); mit++) {
  //   std::cout << mit->first->str() << " : "; 
  //   if (mit->second != nullptr && mit->second->dynamicCast<Id>() != nullptr) 
  //     std::cout << mit->second->dynamicCast<Id>()->str(); 
  //   std::cout << "\\" << std::endl;  
  // }
  // for (auto mit = _argsToConstraint.begin(); mit != _argsToConstraint.end(); mit++) {
  //   std::cout << mit->first->str() << " : "; 
  //   for (auto lit = mit->second.begin(); lit != mit->second.end(); lit++)
  //     std::cout << (*lit)->id() << ", "; 
  //   std::cout << std::endl; 
  // }
}

void MiniZinc::GeasSolverInstance::createVar(geas::solver& s, VarDecl* vd) {
  if (vd->type().isbool()) {
    if (vd->e() == nullptr) {
      Expression* domain = vd->ti()->domain();
      long long int lb = 0, ub = 1;
      if (domain != nullptr) {
        IntBounds ib = compute_int_bounds(_env.envi(), domain);
        lb = ib.l.toInt();
        ub = ib.u.toInt();
      } 
      if (lb == ub) {
        geas::patom_t val = (lb == 0) ? geas::at_False : geas::at_True;
        _variableMap0.insert(vd->id(), GeasVariable(val));
        _variableMap1.insert(vd->id(), GeasVariable(val));
      } else {
        auto var0 = s.new_boolvar();
        auto var1 = s.new_boolvar();
        _variableMap0.insert(vd->id(), GeasVariable(var0));
        _variableMap1.insert(vd->id(), GeasVariable(var1));
      }
    } else {
      std::stringstream ssm;
      ssm << "Independent variable " << vd->id()->str() << " has initlization expression.";
      throw InternalError(ssm.str());
    }
  } else if (vd->type().isint()) {
    if (vd->e() == nullptr) {
      Expression* domain = vd->ti()->domain();
      if (domain != nullptr) {
        IntSetVal* isv = eval_intset(_env.envi(), domain);
        auto var0 = s.new_intvar(static_cast<geas::intvar::val_t>(isv->min().toInt()), static_cast<geas::intvar::val_t>(isv->max().toInt()));
        auto var1 = s.new_intvar(static_cast<geas::intvar::val_t>(isv->min().toInt()), static_cast<geas::intvar::val_t>(isv->max().toInt()));
        if (isv->size() > 1) {
          vec<int> vals(static_cast<int>(isv->card().toInt()));
          int i = 0;
          for (int j = 0; j < isv->size(); ++j)
            for (auto k = isv->min(j).toInt(); k <= isv->max(j).toInt(); ++k)
              vals[i++] = static_cast<int>(k);
          assert(i == isv->card().toInt());
          auto res = geas::make_sparse(var0, vals) && geas::make_sparse(var1, vals);
          assert(res);
        }
        _variableMap0.insert(vd->id(), GeasVariable(var0)); 
        _variableMap1.insert(vd->id(), GeasVariable(var1));
      } else {
        auto var0 = s.new_intvar(SHRT_MIN, SHRT_MAX);
        auto var1 = s.new_intvar(SHRT_MIN, SHRT_MAX);
        _variableMap0.insert(vd->id(), GeasVariable(var0)); 
        _variableMap1.insert(vd->id(), GeasVariable(var1));
        // std::ostringstream ss;
        // ss << "GeasSolverInstance::processFlatZinc: Error: Unbounded variable: " << vd->id()->str();
        // throw Error(ss.str());
      }
    } else {
      std::stringstream ssm;
      ssm << "Independent variable " << vd->id()->str() << " has initlization expression.";
      throw InternalError(ssm.str());
    }
  } else {
    std::stringstream ssm;
    ssm << "Independent variable of type " << *vd->ti() << " is currently not supported by Geas.";
    throw InternalError(ssm.str());
  }
}

// Solving dominance breaking nogood generation problems 
SolverInstanceBase::Status MiniZinc::GeasSolverInstance::solve() { 
  int n = _indepVars.size(); 
  auto _opt = static_cast<GeasOptions&>(*_options);
  assert(_opt.minNogoodLength <= _opt.maxNogoodLength); 
  auto remaining_time = [_opt] {
    if (_opt.time == std::chrono::milliseconds(0)) {
      return 0.0;
    }
    using geas_time = std::chrono::duration<double>;
    static auto timeout = std::chrono::high_resolution_clock::now() + _opt.time;
    return geas_time(timeout - std::chrono::high_resolution_clock::now()).count();
  };
  std::map< std::vector<Id*>, std::list<std::vector<int>> > nogoodLists; 
  for (unsigned int l = _opt.minNogoodLength; l <= _opt.maxNogoodLength; l++) {
    if (l > n)
      break; 
    std::vector<bool> v(n);
    std::fill(v.begin(), v.begin() + l, true);
    do {
      if (remaining_time() < 0.0) 
        return SolverInstance::SAT;
      std::vector<Id*> srcVars;
      for (int i = 0; i < n; ++i) { if (v[i]) { srcVars.push_back(_indepVars[i]); } }
      // for (int i = 0; i != srcVars.size(); i++) { std::cout << srcVars[i]->str() << ": " << _variableName[srcVars[i]] << " "; }
      // std::cout << std::endl;
      _domsolver = new geas::solver(); 
      _variableMap0.clear(); 
      _variableMap1.clear(); 
      _sensvar.clear();
      // initialize independent variables 
      for (auto it = srcVars.begin(); it != srcVars.end(); it++)
        createVar(*_domsolver, (*it)->decl()); 

      std::unordered_set<Id*> fixedVars(srcVars.begin(), srcVars.end()); 
      std::list<Id*> activeVars(srcVars.begin(), srcVars.end());
      std::map<Call*, int> constraintArgNr(_constraintArgNr); 
      std::set<Call*> activeCons; 
      
      while (!activeVars.empty()) {
        Id* curVar = activeVars.front(); 
        // decrease the number of dependent variable count 
        for (auto mit = _argsToConstraint[curVar].begin(); mit != _argsToConstraint[curVar].end(); mit++) {
          Call * c = *mit; 
          constraintArgNr[c]--; 
          if (constraintArgNr[c] == 0) {
            // std::cout << "post constraint " << c->id() << std::endl; 
            // if the constriant defines a variable, create variable and post constraint 
            if (c->ann().containsCall(std::string("defines_var"))) {
              Id* defVar = resolveVarId(c->ann().getCall(std::string("defines_var"))->arg(0));
              createVar(*_domsolver, defVar->decl());
              activeVars.push_back(defVar);
              fixedVars.insert(defVar);
            }
            // post the constraint directly 
            _constraintRegistry.post(c); 
            activeCons.erase(c); 
          } else 
            activeCons.insert(c);
        }
        activeVars.pop_front(); 
      }

      // post dominance condition for the constraint defining the objective variable 
      for (auto sit = activeCons.begin(); sit != activeCons.end(); sit++) {
        if ((*sit)->ann().containsCall(std::string("defines_var"))) {
            Id* defVar = resolveVarId((*sit)->ann().getCall(std::string("defines_var"))->arg(0));
            if (defVar == _objVar)
              _domPoster.post(*sit, fixedVars); 
        }
      }
      
      // post dominance condition for constriants that are not covered 
      for (auto sit = activeCons.begin(); sit != activeCons.end(); sit++) {
        if ((*sit)->ann().containsCall(std::string("defines_var"))) {
          Id* defVar = resolveVarId((*sit)->ann().getCall(std::string("defines_var"))->arg(0));
          if (defVar != _objVar) {
            // std::cout << "post dominance condition for " << (*sit)->id() << std::endl; 
            _domPoster.post(*sit, fixedVars); 
          }
        }
        else {
          // std::cout << "post dominance condition for " << (*sit)->id() << std::endl; 
          _domPoster.post(*sit, fixedVars);   
        }
      }
      
      // add clauses to prevent redundant solving 
      for (unsigned int l2 = _opt.minNogoodLength; l2 <= srcVars.size(); l2++) {
        std::vector<bool> v2(srcVars.size());
        std::fill(v2.begin(), v2.begin() + l2, true);
        do {
          std::vector<Id*> subsetVars; 
          for (int i = 0; i < srcVars.size(); ++i) { if (v2[i]) { subsetVars.push_back(srcVars[i]); } }
          std::list<std::vector<int>> nglst = nogoodLists[subsetVars]; 
          for (auto it = nglst.begin(); it != nglst.end(); it++) {
            vec<geas::clause_elt> clause;
            for (unsigned int j = 0; j != l2; j++) {
              auto geas_var = _variableMap1.get(subsetVars[j]); 
              if (geas_var.isBool()) {
                geas::patom_t bv = geas_var.boolVar();
                clause.push( ((*it)[clause.size()] == 1)? ~bv : bv);
              } else {
                geas::intvar iv = geas_var.intVar();
                clause.push(iv != (*it)[clause.size()]);
              }
            }
            if (!geas::add_clause(*(_domsolver->data), clause))
              break;
          }
        } while (std::prev_permutation(v2.begin(), v2.end()));
      }
      
      // lexicographical ordering  
      {
        vec<geas::intvar> x, y;
        for(unsigned int i=0; i<_sensvar.size(); i=i+2) {
          x.push(_sensvar[i]);
          y.push(_sensvar[i+1]);
        }
        for (unsigned int i=0; i!=srcVars.size(); i++) {
          GeasVariable& var0 = resolveVarDom(follow_id_to_decl(srcVars[i]), true); 
          GeasVariable& var1 = resolveVarDom(follow_id_to_decl(srcVars[i]), false); 
          if (var0.isInt()) {
            x.push(var0.intVar()); 
            y.push(var1.intVar()); 
          } else if (var0.isBool()) {
            geas::intvar iv0 = _domsolver->new_intvar(0, 1);
            geas::intvar iv1 = _domsolver->new_intvar(0, 1);
            geas::add_clause(_domsolver->data, var0.boolVar(), iv0 <= 0);
            geas::add_clause(_domsolver->data, ~var0.boolVar(), iv0 >= 1);
            geas::add_clause(_domsolver->data, var1.boolVar(), iv1 <= 0);
            geas::add_clause(_domsolver->data, ~var1.boolVar(), iv1 >= 1);
            x.push(iv0);
            y.push(iv1);
          }
        }
        vec<geas::patom_t> b(x.size()+1);
        vec<geas::patom_t> r(x.size());
        b[0] = geas::at_True;  
        for (int l = 1; l < x.size(); l++) 
          b[l] = _domsolver->new_boolvar();
        b[x.size()] = geas::at_False; 
        for (int l = 0; l < x.size(); l++)
          r[l] = _domsolver->new_boolvar();
        for (int l = 0; l < x.size(); l++) {
          geas::int_le(_domsolver->data, x[l], y[l], -1, r[l]);
          geas::add_clause(_domsolver->data, b[l], ~r[l]);
          geas::add_clause(_domsolver->data, b[l], ~b[l+1]);
          geas::add_clause(_domsolver->data, ~b[l], r[l], b[l+1]);
        }
        for (int l = 0; l < x.size(); l++)
          geas::int_le(_domsolver->data, x[l], y[l], 0, b[l]); 
      }

      // common assignment elimination 
      // {
      //   for (unsigned int i=0; i!=srcVars.size(); i++) {
      //     GeasVariable& var0 = resolveVarDom(follow_id_to_decl(srcVars[i]), true); 
      //     GeasVariable& var1 = resolveVarDom(follow_id_to_decl(srcVars[i]), false); 
      //     geas::int_ne(_domsolver->data, var0.intVar(), var1.intVar()); 
      //   }
      // }

      // solve and store the dominance breaking nogoods 
      geas::solver::result res = geas::solver::SAT;
      int solutionNr = 0; 
      while (res == geas::solver::SAT)
      {
        // res = _domsolver->solve({0.5, 0});
        res = _domsolver->solve({remaining_time(), 0});
        if(res == geas::solver::SAT) {
          solutionNr++; 
          _domsolver->restart(); 
          // get dominance breaking nogoods 
          geas::model solution = _domsolver->get_model(); 
          vec<geas::clause_elt> clause;
          std::vector<int> nogoods(srcVars.size()); 
          for (unsigned int i=0; i!=srcVars.size(); i++) {
            auto geas_var = _variableMap1.get(srcVars[i]);
            if (geas_var.isBool()) {
              geas::patom_t bv = geas_var.boolVar();
              clause.push(solution.value(bv)? ~bv:bv);
              nogoods[i] = solution.value(bv)? 1:0; 
            } else {
              geas::intvar iv = geas_var.intVar();
              clause.push(iv != solution[iv]);
              nogoods[i] = solution[iv]; 
            }
          }
          nogoodLists[srcVars].push_back(nogoods); 

          // output dominance breaking nogood 
          std::cout << "constraint "; 
          for (unsigned int i = 0; i != srcVars.size(); i++) {
            auto geas_var = _variableMap1.get(srcVars[i]); 
            std::cout <<  _variableName[srcVars[i]] << " != "; 
            if (geas_var.isBool())
              std::cout << (solution.value(geas_var.boolVar())?"true":"false"); 
            else 
              std::cout << solution[geas_var.intVar()]; 
            if (i+1 != srcVars.size())
              std::cout << " \\/ "; 
            else 
              std::cout << ";" << std::endl;  
          }

          // std::cout << "constraint "; 
          // for (unsigned int i = 0; i != srcVars.size(); i++) {
          //   auto geas_var = _variableMap0.get(srcVars[i]); 
          //   std::cout << _variableName[srcVars[i]] << " = "; 
          //   if (geas_var.isBool())
          //     std::cout << (solution.value(geas_var.boolVar())?"true":"false"); 
          //   else 
          //     std::cout << solution[geas_var.intVar()]; 
          //   if (i+1 != srcVars.size())
          //     std::cout << " \\/ "; 
          //   else 
          //     std::cout << ";" << std::endl;  
          // }

          // // output all variables 
          // std::cout << "variableMap 0: "; 
          // for (auto id: fixedVars) {
          //   auto geas_var = _variableMap0.get(id); 
          //   std::cout << " " << ((_variableName.find(id) == _variableName.end())?id->str():_variableName[id]) << ": ";
          //   if (geas_var.isBool())
          //     std::cout << (solution.value(geas_var.boolVar())?"true":"false");
          //   else 
          //     std::cout << solution[geas_var.intVar()]; 
          //   std::cout << "(" << _variableMono[id] << ")" << ","; 
          // } 
          // std::cout << std::endl; 
          
          // std::cout << "variableMap 1: "; 
          // for (auto id: fixedVars) {
          //   auto geas_var = _variableMap1.get(id); 
          //   std::cout << " " << ((_variableName.find(id) == _variableName.end())?id->str():_variableName[id]) << ": ";
          //   if (geas_var.isBool())
          //     std::cout << (solution.value(geas_var.boolVar())?"true":"false");
          //   else 
          //     std::cout << solution[geas_var.intVar()]; 
          //   std::cout << "(" << _variableMono[id] << ")" << ","; 
          // } 
          // std::cout << std::endl; 
          
          if (!geas::add_clause(*(_domsolver->data), clause))
            res = geas::solver::UNSAT;
        }
      }
      delete _domsolver;
    } while (std::prev_permutation(v.begin(), v.end()));
  }
  // output nogoodLists 
  return SolverInstance::SAT;
}

Expression* GeasSolverInstance::getSolutionValue(Id* id) { assert(false); return nullptr; } 

void GeasSolverInstance::resetSolver() { assert(false); }


GeasTypes::Variable& GeasSolverInstance::resolveVarDom(Expression* e, bool dominating) {
  if (auto* id = e->dynamicCast<Id>()) {
    return dominating?_variableMap0.get(id->decl()->id()):_variableMap1.get(id->decl()->id());
  }
  if (auto* vd = e->dynamicCast<VarDecl>()) {
    return dominating?_variableMap0.get(vd->id()->decl()->id()):_variableMap1.get(vd->id()->decl()->id());
  }
  if (auto* aa = e->dynamicCast<ArrayAccess>()) {
    auto* ad = aa->v()->cast<Id>()->decl();
    auto idx = aa->idx()[0]->cast<IntLit>()->v().toInt();
    auto* al = eval_array_lit(_env.envi(), ad->e());
    return dominating?_variableMap0.get((*al)[idx]->cast<Id>()):_variableMap1.get((*al)[idx]->cast<Id>());
  }
  std::stringstream ssm;
  ssm << "Expected Id, VarDecl or ArrayAccess instead of \"" << *e << "\"";
  throw InternalError(ssm.str());
}

vec<bool> GeasSolverInstance::asBool(ArrayLit* al) {
  vec<bool> vec(static_cast<int>(al->size()));
  for (int i = 0; i < al->size(); ++i) {
    vec[i] = asBool((*al)[i]);
  }
  return vec;
}

geas::patom_t GeasSolverInstance::asBoolVarDom(Expression * e, bool dominating) {
  if (e->type().isvar()) {
    GeasVariable& var = resolveVarDom(follow_id_to_decl(e), dominating);
    assert(var.isBool());
    return var.boolVar();
  }
  if (auto* bl = e->dynamicCast<BoolLit>()) {
    return bl->v() ? geas::at_True : geas::at_False;
  }
  std::stringstream ssm;
  ssm << "Expected bool literal instead of: " << *e;
  throw InternalError(ssm.str());
}

vec<geas::patom_t> GeasSolverInstance::asBoolVarDom(ArrayLit* al, bool dominating) {
  vec<geas::patom_t> vec(static_cast<int>(al->size()));
  for (int i = 0; i < al->size(); ++i)
    vec[i] = this->asBoolVarDom((*al)[i], dominating);
  return vec;
}

vec<int> GeasSolverInstance::asInt(ArrayLit* al) {
  vec<int> vec(static_cast<int>(al->size()));
  for (int i = 0; i < al->size(); ++i) {
    vec[i] = this->asInt((*al)[i]);
  }
  return vec;
}

geas::intvar GeasSolverInstance::asIntVarDom(Expression * e, bool dominating) {
  if (e->type().isvar()) {
    GeasVariable& var = resolveVarDom(follow_id_to_decl(e), dominating);
    assert(var.isInt());
    return var.intVar();
  }
  IntVal i;
  if (auto* il = e->dynamicCast<IntLit>()) {
    i = il->v().toInt();
  } else if (auto* bl = e->dynamicCast<BoolLit>()) {
    i = static_cast<long long>(bl->v());
  } else {
    std::stringstream ssm;
    ssm << "Expected bool or int literal instead of: " << *e;
    throw InternalError(ssm.str());
  }
  
  return _domsolver->new_intvar(static_cast<geas::intvar::val_t>(i.toInt()),
                            static_cast<geas::intvar::val_t>(i.toInt()));
}



vec<geas::intvar> GeasSolverInstance::asIntVarDom(ArrayLit* al, bool dominating) {
  vec<geas::intvar> vec(static_cast<int>(al->size()));
  for (int i = 0; i < al->size(); ++i)
    vec[i] = this->asIntVarDom((*al)[i], dominating);
  return vec;
}

Id* GeasSolverInstance::resolveVarId(Expression* e) {
  if (auto* id = e->dynamicCast<Id>()) {
    return id->decl()->id(); 
  }
  if (auto* vd = e->dynamicCast<VarDecl>()) {
    return vd->id()->decl()->id(); 
  }
  if (auto* aa = e->dynamicCast<ArrayAccess>()) {
    auto* ad = aa->v()->cast<Id>()->decl();
    auto idx = aa->idx()[0]->cast<IntLit>()->v().toInt();
    auto* al = eval_array_lit(_env.envi(), ad->e());
    return (*al)[idx]->cast<Id>(); 
  }
  std::stringstream ssm;
  ssm << "Expected Id, VarDecl or ArrayAccess instead of \"" << *e << "\"";
  throw InternalError(ssm.str());
}

Id* GeasSolverInstance::asVarId(Expression* e) {
  if (e->type().isvar())
    return resolveVarId(follow_id_to_decl(e)); 
  else 
    return nullptr; 
}

std::vector<Id*> GeasSolverInstance::asVarId(ArrayLit* al) {
  std::vector<Id*> vec(static_cast<int>(al->size()));
  for (int i = 0; i < al->size(); ++i) {
    vec[i] = this->asVarId((*al)[i]);
  }
  return vec;
}

void GeasSolverInstance::replaceCallAlias(Call* c) { 
  for (unsigned int i = 0; i != c->argCount(); i++) {
    if (Id* id = c->arg(i)->dynamicCast<Id>()) {
      if (id != nullptr) {
        id = resolveVarId(id);
        if(_alias[id] != nullptr) 
          c->arg(i, _alias[id]);
      }
    } else if (ArrayLit* a = c->arg(i)->dynamicCast<ArrayLit>())
      replaceArrayAlias(a);
  }
}

void GeasSolverInstance::replaceArrayAlias(ArrayLit* a) { 
  for (unsigned int i = 0; i != a->size(); i++) {
    Id* id = (*a)[i]->dynamicCast<Id>(); 
    if (id != nullptr) {
      id = resolveVarId(id);
      if(_alias[id] != nullptr) 
        a->set(i, _alias[id]);
    }
  }
}

std::set<Id*> GeasSolverInstance::getArguments(Call* c) {
  int sz = c->argCount();
  std::set<Id*> args; 
  std::list<Expression*> exps; 
  for (unsigned int i = 0; i != sz; i++) 
    exps.push_back(c->arg(i)); 
  while (!exps.empty()) {
    Expression* e = exps.front(); 
    if (Id* id = e->dynamicCast<Id>()) {
      if ( id->decl()->ti()->isarray() ) {
        ArrayLit* a = id->decl()->e()->dynamicCast<ArrayLit>(); 
        for (unsigned int j = 0; j != a->size(); j++) 
          exps.push_back(a->getSlice(j)); 
      }
      else
        args.insert(resolveVarId(id));
    } 
    else if (ArrayLit* a = e->dynamicCast<ArrayLit>()) {
      for (unsigned int j = 0; j != a->size(); j++) 
        exps.push_back(a->getSlice(j)); 
    }
    exps.pop_front();
  }
  if (c->ann().containsCall(std::string("defines_var")))
    args.erase(asVarId(c->ann().getCall(std::string("defines_var"))->arg(0))); 
  
  return args; 
}

void GeasSolverInstance::printStatistics() {
  assert(false);
  // todo: print number of nogoods, total time 
}

GeasSolverFactory::GeasSolverFactory() {
  SolverConfig sc("org.minizinc.geas", getVersion(nullptr));
  sc.name("Geas");
  sc.mznlib("-Gdominance");
  sc.mznlibVersion(1);
  sc.supportsMzn(false);
  sc.description(getDescription(nullptr));
  sc.tags({
      "api",
      "cp",
      "float",
      "int",
      "lcg",
  });
  sc.stdFlags({"-a", "-f", "-n", "-s", "-t"});
  sc.extraFlags({
      SolverConfig::ExtraFlag("--conflicts",
                              "Limit the maximum number of conflicts to be used during solving.",
                              SolverConfig::ExtraFlag::FlagType::T_INT, {}, "0"),
      SolverConfig::ExtraFlag(
          "--obj-probe",
          "Number of conflicts to use to probe for better solutions after a new solution is found.",
          SolverConfig::ExtraFlag::FlagType::T_INT, {}, "0"),
  });
  SolverConfigs::registerBuiltinSolver(sc);
};

SolverInstanceBase::Options* GeasSolverFactory::createOptions() { return new GeasOptions; }

SolverInstanceBase* GeasSolverFactory::doCreateSI(Env& env, std::ostream& log,
                                                  SolverInstanceBase::Options* opt) {
  return new GeasSolverInstance(env, log, opt);
}

bool GeasSolverFactory::processOption(SolverInstanceBase::Options* opt, int& i,
                                      std::vector<std::string>& argv,
                                      const std::string& workingDir) {
  auto* _opt = static_cast<GeasOptions*>(opt);
  if (argv[i] == "-a" || argv[i] == "--all-solutions") {
    _opt->allSolutions = true;
  } else if (argv[i] == "--conflicts") {
    if (++i == argv.size()) {
      return false;
    }
    int nodes = atoi(argv[i].c_str());
    if (nodes >= 0) {
      _opt->conflicts = nodes;
    }
  } else if (argv[i] == "-f") {
    _opt->freeSearch = true;
  } else if (argv[i] == "-n") {
    if (++i == argv.size()) {
      return false;
    }
    int n = atoi(argv[i].c_str());
    if (n >= 0) {
      _opt->nrSolutions = n;
    }
  } else if (argv[i] == "--minnogoodlen") {
    if (++i == argv.size()) {
      return false;
    }
    int n = atoi(argv[i].c_str());
    if (n >= 1) {
      _opt->minNogoodLength = n;
    }
  } else if (argv[i] == "--maxnogoodlen") {
    if (++i == argv.size()) {
      return false;
    }
    int n = atoi(argv[i].c_str());
    if (n >= 1) {
      _opt->maxNogoodLength = n;
    }
  }
  else if (argv[i] == "--obj-probe") {
    if (++i == argv.size()) {
      return false;
    }
    int limit = atoi(argv[i].c_str());
    if (limit >= 0) {
      _opt->objProbeLimit = limit;
    }
  } else if (argv[i] == "--solver-statistics" || argv[i] == "-s") {
    _opt->statistics = true;
  } else if (argv[i] == "--solver-time-limit" || argv[i] == "-t") {
    if (++i == argv.size()) {
      return false;
    }
    int time = atoi(argv[i].c_str());
    if (time >= 0) {
      _opt->time = std::chrono::milliseconds(time);
    }
  } else {
    return false;
  }
  return true;
}

void GeasSolverFactory::printHelp(std::ostream& os) {
  os << "Geas solver plugin options:" << std::endl
     << "  --minnogoodlen <int>" << std::endl 
     << "    Minimum length of generated dominance breaking nogoods." << std::endl
     << "  --maxnogoodlen <int>" << std::endl 
     << "    Maximum length of generated dominance breaking nogoods." << std::endl
     << std::endl;
}
}  // namespace MiniZinc
