// * -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */

/*
 *  Main authors:
 *     Gleb Belov <gleb.belov@monash.edu>
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifdef _MSC_VER 
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <ctime>
#include <cmath>
#include <stdexcept>

using namespace std;

#include <minizinc/solvers/MIP/MIP_gurobi_wrap.hh>
#include <minizinc/utils.hh>

/// Linking this module provides these functions:
MIP_wrapper* MIP_WrapperFactory::GetDefaultMIPWrapper() {
  return new MIP_gurobi_wrapper;
}

string MIP_WrapperFactory::getVersion( ) {
  ostringstream oss;
  oss << "  MIP wrapper for Gurobi library ";
  int major, minor, technical;
  GRBversion(&major, &minor, &technical);
  oss << major << '.' << minor << '.' << technical;
  oss << ".  Compiled  " __DATE__ "  " __TIME__;
  return oss.str();
}

void MIP_WrapperFactory::printHelp(ostream& os) {
  os
  << "GUROBI MIP wrapper options:" << std::endl
  // -s                  print statistics
  //            << "  --readParam <file>  read GUROBI parameters from file
  //               << "--writeParam <file> write GUROBI parameters to file
  //               << "--tuneParam         instruct GUROBI to tune parameters instead of solving
  << "--writeModel <file> write model to <file> (.lp, .mps, .sav, ...)" << std::endl
  << "-a                  print intermediate solutions (use for optimization problems only TODO)" << std::endl
  << "-p <N>              use N threads, default: 1." << std::endl
//   << "--nomippresolve     disable MIP presolving   NOT IMPL" << std::endl
  << "--timeout <N>       stop search after N seconds" << std::endl
//   << "--workmem <N>       maximal amount of RAM used, MB" << std::endl
  << "--readParam <file>  read GUROBI parameters from file" << std::endl
  << "--writeParam <file> write GUROBI parameters to file" << std::endl
//   << "--tuneParam         instruct GUROBI to tune parameters instead of solving   NOT IMPL"

  << "--absGap <n>        absolute gap |primal-dual| to stop. Default 0.99" << std::endl
  << "--relGap <n>        relative gap |primal-dual|/<solver-dep> to stop. Default 1e-8" << std::endl
  << "--intTol <n>        integrality tolerance for a variable. Default 1e-6" << std::endl
//   << "--objDiff <n>       objective function discretization. Default 1.0" << std::endl

  << std::endl;
}

  static inline bool beginswith(string s, string t) {
    return s.compare(0, t.length(), t)==0;
  }

            /// SOLVER PARAMS ????
 static   int nThreads=1;
 static   string sExportModel;
 static   double nTimeout=-1;
 static   double nWorkMemLimit=-1;
 static   string sReadParams;
 static   string sWriteParams;
 static   bool flag_all_solutions = false;

 static   double absGap=0.99;
 static   double relGap=1e-8;
 static   double intTol=1e-6;
 static   double objDiff=1.0;

bool MIP_WrapperFactory::processOption(int& i, int argc, const char** argv) {
  MiniZinc::CLOParser cop( i, argc, argv );
  if ( string(argv[i])=="-a"
      || string(argv[i])=="--all"
      || string(argv[i])=="--all-solutions" ) {
    flag_all_solutions = true;
  } else if (string(argv[i])=="-f") {
//     std::cerr << "  Flag -f: ignoring fixed strategy anyway." << std::endl;
  } else if ( cop.get( "--writeModel", &sExportModel ) ) {
  } else if ( cop.get( "-p", &nThreads ) ) {
  } else if ( cop.get( "--timeout", &nTimeout ) ) {
  } else if ( cop.get( "--workmem", &nWorkMemLimit ) ) {
  } else if ( cop.get( "--readParam", &sReadParams ) ) {
  } else if ( cop.get( "--writeParam", &sWriteParams ) ) {
  } else if ( cop.get( "--absGap", &absGap ) ) {
  } else if ( cop.get( "--relGap", &relGap ) ) {
  } else if ( cop.get( "--intTol", &intTol ) ) {
//   } else if ( cop.get( "--objDiff", &objDiff ) ) {
  } else
    return false;
  return true;
error:
  return false;
}

void MIP_gurobi_wrapper::wrap_assert(bool cond, string msg, bool fTerm)
{
   if ( !cond ) {
      gurobi_buffer = "[NO ERROR STRING GIVEN]";
      if (error) {
         gurobi_buffer = GRBgeterrormsg(env);
      }
      string msgAll = ("  MIP_gurobi_wrapper runtime error:  " + msg + "  " + gurobi_buffer);
      cerr << msgAll << endl;
      if (fTerm) {
        cerr << "TERMINATING." << endl;
        throw runtime_error(msgAll);
      }
   }
}

void MIP_gurobi_wrapper::openGUROBI()
{
   /* Initialize the GUROBI environment */
   error = GRBloadenv (&env, "mzn-gurobi.log");
   wrap_assert ( !error, "Could not open GUROBI environment." );
   error = GRBsetintparam(env, "OutputFlag", 0);  // Switch off output
//    error = GRBsetintparam(env, "LogToConsole", 
//                             fVerbose ? 1 : 0);  // also when flag_all_solutions?  TODO
  /* Create the problem. */
   error = GRBnewmodel(env, &model, "mzn_gurobi", 0, NULL, NULL, NULL, NULL, NULL);
   wrap_assert ( model, "Failed to create LP." );
}

void MIP_gurobi_wrapper::closeGUROBI()
{
  /// Freeing the problem can be slow both in C and C++, see IBM forums. Skipping.
     /* Free up the problem as allocated by GRB_createprob, if necessary */
  /* Free model */

  GRBfreemodel(model);      
  model = 0;

  /* Free environment */

  if (env)
    GRBfreeenv(env);
  /// and at last:
//   MIP_wrapper::cleanup();
}

void MIP_gurobi_wrapper::doAddVars
  (size_t n, double* obj, double* lb, double* ub, MIP_wrapper::VarType* vt, string *names)
{
  /// Convert var types:
  vector<char> ctype(n);
  vector<char*> pcNames(n);
  for (size_t i=0; i<n; ++i) {
    pcNames[i] = (char*)names[i].c_str();
    switch (vt[i]) {
      case REAL:
        ctype[i] = GRB_CONTINUOUS;
        break;
      case INT:
        ctype[i] = GRB_INTEGER;
        break;
      case BINARY:
        ctype[i] = GRB_BINARY;
        break;
      default:
        throw runtime_error("  MIP_wrapper: unknown variable type");
    }
  }
  error = GRBaddvars(model, n, 0, NULL, NULL, NULL, obj, lb, ub, &ctype[0], &pcNames[0]);
  wrap_assert( !error,  "Failed to declare variables." );
  error = GRBupdatemodel(model);
  wrap_assert( !error,  "Failed to update model." );
}

static char getGRBSense( MIP_wrapper::LinConType s ) {
    switch (s) {
      case MIP_wrapper::LQ:
        return GRB_LESS_EQUAL;
      case MIP_wrapper::EQ:
        return GRB_EQUAL;
      case MIP_wrapper::GQ:
        return GRB_GREATER_EQUAL;
      default:
        throw runtime_error("  MIP_gurobi_wrapper: unknown constraint sense");
    }
}

void MIP_gurobi_wrapper::addRow
  (int nnz, int* rmatind, double* rmatval, MIP_wrapper::LinConType sense,
   double rhs, int mask, string rowName)
{
  /// Convert var types:
  char ssense=getGRBSense(sense);
  const int ccnt=0;
  const int rcnt=1;
  const int rmatbeg[] = { 0 };
  const char * pRName = rowName.c_str();
//   if (MaskConsType_Normal & mask) {
  /// User cuts & lazy only by callback in Gurobi:
    error = GRBaddconstr(model, nnz, rmatind, rmatval, ssense, rhs, pRName);
    wrap_assert( !error,  "Failed to add constraint." );
//   }
}

/// SolutionCallback ------------------------------------------------------------------------
/// Gurobi ensures thread-safety
static int __stdcall
solcallback(GRBmodel *model,
           void     *cbdata,
           int       where,
           void     *usrdata)
{
  MIP_wrapper::CBUserInfo *info = (MIP_wrapper::CBUserInfo*) usrdata;
  double nodecnt=0.0, actnodes=0.0, objVal=0.0;
  int    solcnt=0;
  int    newincumbent=0;

  if ( GRB_CB_MIP==where ) {
      /* General MIP callback */
      GRBcbget(cbdata, where, GRB_CB_MIP_OBJBND, &info->pOutput->bestBound);
        GRBcbget(cbdata, where, GRB_CB_MIP_NODLFT, &actnodes);
      info->pOutput->nOpenNodes = actnodes;
  } else if ( GRB_CB_MESSAGE==where ) {
    /* Message callback */
    if ( info->fVerb ) {
      char *msg;
      GRBcbget(cbdata, where, GRB_CB_MSG_STRING, &msg);
      cerr << msg << flush;
    }
  } else if ( GRB_CB_MIPSOL==where ) {
      /* MIP solution callback */
      GRBcbget(cbdata, where, GRB_CB_MIPSOL_NODCNT, &nodecnt);
      info->pOutput->nNodes = nodecnt;
      GRBcbget(cbdata, where, GRB_CB_MIPSOL_OBJ, &objVal);
      GRBcbget(cbdata, where, GRB_CB_MIPSOL_SOLCNT, &solcnt);

      if ( solcnt ) {
        
        if ( fabs(info->pOutput->objVal - objVal) > 1e-12*(1.0 + fabs(objVal)) ) {
          newincumbent = 1;
          info->pOutput->objVal = objVal;
          info->pOutput->status = MIP_wrapper::SAT;
          info->pOutput->statusName = "feasible from a callback";
        }
      }
    if ( newincumbent ) {
        assert(info->pOutput->x);
        GRBcbget(cbdata, where, GRB_CB_MIPSOL_SOL, (void*)info->pOutput->x);
        
        info->pOutput->dCPUTime = -1;

        /// Call the user function:
        if (info->solcbfn)
            (*info->solcbfn)(*info->pOutput, info->ppp);
    }
    /// Callback for lazy cuts
    if ( info->cutcbfn && info->cutMask&MIP_wrapper::MaskConsType_Lazy ) {
      MIP_wrapper::CutInput cutInput;
      info->cutcbfn( *info->pOutput, cutInput, info->ppp, true );
      for ( auto& cd : cutInput ) {
//         assert( cd.mask & MIP_wrapper::MaskConsType_Lazy );
        if ( cd.mask & MIP_wrapper::MaskConsType_Lazy ) {
          int error = GRBcblazy(cbdata, cd.rmatind.size(),
                  cd.rmatind.data(), cd.rmatval.data(), 
                  getGRBSense(cd.sense), cd.rhs);
          if (error)
            cerr << "  GRB_wrapper: failed to add lazy cut. " << endl;
        }
      }
    }
  } else if ( GRB_CB_MIPNODE==where  ) {
    int status;
    GRBcbget(cbdata, where, GRB_CB_MIPNODE_STATUS, &status);
    if ( status == GRB_OPTIMAL && info->cutcbfn ) {    // if cut handler given
      MIP_wrapper::Output outpRlx;
      outpRlx.x = info->pOutput->x;  // using the sol output storage  TODO?
      outpRlx.nCols = info->pOutput->nCols;
      assert( outpRlx.x && outpRlx.nCols );
//       GRBcbget(cbdata, where, GRB_CB_MIPNODE_RELOBJ, outpRlx.objVal);
      GRBcbget(cbdata, where, GRB_CB_MIPNODE_REL, (void*)outpRlx.x);
      MIP_wrapper::CutInput cutInput;
      info->cutcbfn( outpRlx, cutInput, info->ppp, false );
//       static int nCuts=0;
//       nCuts += cutInput.size();
//       if ( cutInput.size() )
//         cerr << "\n   N CUTS:  " << nCuts << endl;
      for ( auto& cd : cutInput ) {
        assert( cd.mask &
          (MIP_wrapper::MaskConsType_Usercut|MIP_wrapper::MaskConsType_Lazy) );
        if ( cd.mask & MIP_wrapper::MaskConsType_Usercut ) {
          int error = GRBcbcut(cbdata, cd.rmatind.size(),
                  cd.rmatind.data(), cd.rmatval.data(), 
                  getGRBSense(cd.sense), cd.rhs);
          if (error)
            cerr << "  GRB_wrapper: failed to add user cut. " << endl;
        }
        if ( cd.mask & MIP_wrapper::MaskConsType_Lazy ) {
          int error = GRBcblazy(cbdata, cd.rmatind.size(),
                  cd.rmatind.data(), cd.rmatval.data(), 
                  getGRBSense(cd.sense), cd.rhs);
          if (error)
            cerr << "  GRB_wrapper: failed to add lazy cut. " << endl;
        }
      }
    }
  }  
  return 0;
} /* END logcallback */
// end SolutionCallback ---------------------------------------------------------------------


MIP_gurobi_wrapper::Status MIP_gurobi_wrapper::convertStatus(int gurobiStatus)
{
  Status s = Status::UNKNOWN;
  ostringstream oss;
   /* Converting the status. */
  if (gurobiStatus == GRB_OPTIMAL) {
    s = Status::OPT;
    oss << "Optimal";
  } else if (gurobiStatus == GRB_INF_OR_UNBD) {
    s = Status::UNSATorUNBND;
    oss << "Infeasible or unbounded";
  } else if (gurobiStatus == GRB_INFEASIBLE) {
    s = Status::UNSAT;
    oss << "Infeasible";
  } else if (gurobiStatus == GRB_UNBOUNDED) {
    oss << "Unbounded";
    s = Status::UNBND;
  } else {
    int solcount=0;
    error = GRBgetintattr(model, "SolCount", &solcount);
    wrap_assert(!error, "  Failure to access solution count.", false);
    if (solcount)
      s = Status::SAT;
    oss << "Gurobi stopped with status " << gurobiStatus;
  }
  output.statusName = gurobi_status_buffer = oss.str();
  return s;
}


void MIP_gurobi_wrapper::solve() {  // Move into ancestor?
   error = GRBupdatemodel(model);                  // for model export
   wrap_assert( !error,  "Failed to update model." );

  /////////////// Last-minute solver options //////////////////
  /* Turn on output to file */
   error = GRBsetstrparam(GRBgetenv(model), "LogFile", "");  // FAILS to switch off in Ubuntu 15.04
  /* Turn on output to the screen */
   error = GRBsetintparam(GRBgetenv(model), "OutputFlag", 
                             /*fVerbose ? 1 :*/ 0);  // switch off, redirect in callback
//    error = GRBsetintparam(GRBgetenv(model), "LogToConsole", 
//                             fVerbose ? 1 : 0);  // also when flag_all_solutions?  TODO
   wrap_assert(!error, "  GUROBI Warning: Failure to switch screen indicator.", false);
//    error =  GRB_setintparam (env, GRB_PARAM_ClockType, 1);            // CPU time
//    error =  GRB_setintparam (env, GRB_PARAM_MIP_Strategy_CallbackReducedLP, GRB__OFF);    // Access original model
   if (sExportModel.size()) {
     error = GRBwrite(model, sExportModel.c_str());
     wrap_assert(!error, "Failed to write LP to disk.", false);
   }

   /// TODO
//     if(all_solutions && obj.getImpl()) {
//       IloNum lastObjVal = (obj.getSense() == IloObjective::Minimize ) ?
//       _ilogurobi->use(SolutionCallback(_iloenv, lastObjVal, *this));
      // Turn off GUROBI logging

   if (nThreads>0) {
     error = GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_THREADS, nThreads);
//      int nn;    // THE SETTING FAILS TO WORK IN 6.0.5.
//      error = GRBgetintparam(env, GRB_INT_PAR_THREADS, &nn);
//      cerr << "Set " << nThreads << " threads, reported " << nn << endl;
     wrap_assert(!error, "Failed to set GRB_INT_PAR_THREADS.", false);
   }

    if (nTimeout>0) {
     error = GRBsetdblparam(GRBgetenv(model), GRB_DBL_PAR_TIMELIMIT, nTimeout);
     wrap_assert(!error, "Failed to set GRB_PARAM_TimeLimit.", false);
    }

//     if (nWorkMemLimit>0) {
//      error =  GRB_setdblparam (env, GRB_PARAM_MIP_Limits_TreeMemory, nWorkMemLimit);
//      wrap_assert(!error, "Failed to set GRB_PARAM_MIP_Limits_TreeMemory.", false);
//     }

   if ( true ) {
     error = GRBsetdblparam( GRBgetenv(model),  "MIPGapAbs", absGap );
     wrap_assert(!error, "Failed to set  MIPGapAbs.", false);
   }
   if ( true ) {
     error = GRBsetdblparam( GRBgetenv(model),  "MIPGap", relGap );
     wrap_assert(!error, "Failed to set  MIPGap.", false);
   }
   if ( true ) {
     error = GRBsetdblparam( GRBgetenv(model),  "IntFeasTol", intTol );
     wrap_assert(!error, "Failed to set   IntFeasTol.", false);
   }

    
       /// Solution callback
   output.nCols = colObj.size();
   x.resize(output.nCols);
   output.x = &x[0];
   if (true) {                 // Need for logging
      cbui.fVerb = fVerbose;
      if ( !flag_all_solutions )
        cbui.solcbfn = 0;
      if ( cbui.cutcbfn ) {
        assert( cbui.cutMask & (MaskConsType_Usercut|MaskConsType_Lazy) );
        if ( cbui.cutMask & MaskConsType_Usercut ) {
          // For user cuts, needs to keep some info after presolve
          if ( fVerbose )
            cerr << "  MIP_gurobi_wrapper: user cut callback enabled, setting PreCrush=1" << endl;
          error = GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_PRECRUSH, 1);
          wrap_assert(!error, "Failed to set GRB_INT_PAR_PRECRUSH.", false);
        }
        if ( cbui.cutMask & MaskConsType_Lazy ) {
          // For lazy cuts, Gurobi disables some presolves
          if ( fVerbose )
            cerr << "  MIP_gurobi_wrapper: lazy cut callback enabled, setting LazyConstraints=1" << endl;
          error = GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_LAZYCONSTRAINTS, 1);
          wrap_assert(!error, "Failed to set GRB_INT_PAR_LAZYCONSTRAINTS.", false);
        }
      }
      error = GRBsetcallbackfunc(model, solcallback, (void *) &cbui);
      wrap_assert(!error, "Failed to set callback", false);
   }

   /// after all modifs
    if (sReadParams.size()) {
     error = GRBreadparams (GRBgetenv(model), sReadParams.c_str());
     wrap_assert(!error, "Failed to read GUROBI parameters.", false);
    }
    
    if (sWriteParams.size()) {
     error = GRBwriteparams (GRBgetenv(model), sWriteParams.c_str());
     wrap_assert(!error, "Failed to write GUROBI parameters.", false);
    }

   output.dCPUTime = std::clock();

   /* Optimize the problem and obtain solution. */
   error = GRBoptimize(model);
   wrap_assert( !error,  "Failed to optimize MIP." );

   output.dCPUTime = (std::clock() - output.dCPUTime) / CLOCKS_PER_SEC;

   int solstat;
   error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &solstat);
   wrap_assert(!error, "Failed to get MIP status.", false);
   output.status = convertStatus(solstat);

   /// Continuing to fill the output object:
   if (Status::OPT == output.status || Status::SAT ==output.status) {
      error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &output.objVal);
      wrap_assert( !error, "No MIP objective value available." );

    //    int cur_numrows = GRB_getnumrows (env, lp);
      int cur_numcols = getNCols();
      assert(cur_numcols == colObj.size());
      
      x.resize(cur_numcols);
      output.x = &x[0];
      error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, cur_numcols, (double*)output.x);
      wrap_assert(!error, "Failed to get variable values.");
   }
   output.bestBound = 1e308;
   error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJBOUNDC, &output.bestBound);
   wrap_assert(!error, "Failed to get the best bound.", false);
   double nNodes=-1;
   error = GRBgetdblattr(model, GRB_DBL_ATTR_NODECOUNT, &nNodes);
   output.nNodes = nNodes;
   output.nOpenNodes = 0;
}

void MIP_gurobi_wrapper::setObjSense(int s)
{
  error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE,
                        s>0 ? GRB_MAXIMIZE : GRB_MINIMIZE);
  wrap_assert(!error, "Failed to set obj sense.");
}

