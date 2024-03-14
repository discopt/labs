// #define SCIP_DEBUG

#include <cassert>
#include <string>
#include <iostream>
#include <vector>

#include "cons_labsnogood.hpp"

#include "objscip/objscip.h"

#include "scip/cons_linear.h"

ConshdlrLABSnogood::ConshdlrLABSnogood(SCIP* scip)
  : ObjConshdlr(scip, "LABSnogood", "LABS nogood constraints", 1000000, -2000000, -2000000, 1, -1, 1, 0, FALSE, FALSE,
  TRUE, SCIP_PROPTIMING_BEFORELP, SCIP_PRESOLTIMING_FAST)
{

}

ConshdlrLABSnogood::~ConshdlrLABSnogood()
{

}

struct SCIP_ConsData
{
  int shift;
  std::vector<SCIP_VAR*> varSequence;
  std::unordered_map<int, SCIP_VAR*> varMap;
};

/** frees specific constraint data */
SCIP_DECL_CONSDELETE(ConshdlrLABSnogood::scip_delete)
{
  assert(consdata != NULL);

  for (SCIP_VAR* var : (*consdata)->varSequence)
    SCIP_CALL( SCIPreleaseVar(scip, &var) );
  for (auto iter : (*consdata)->varMap)
    SCIP_CALL( SCIPreleaseVar(scip, &iter.second) );

  delete (*consdata);
  *consdata = NULL;

  return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
SCIP_DECL_CONSTRANS(ConshdlrLABSnogood::scip_trans)
{
  SCIP_CONSDATA* sourcedata = SCIPconsGetData(sourcecons);
  assert(sourcedata);

  SCIP_CONSDATA* targetdata = new SCIP_CONSDATA;
  targetdata->shift = sourcedata->shift;

  for (SCIP_VAR* var : sourcedata->varSequence)
  {
    SCIP_VAR* transvar = NULL;
    SCIP_CALL( SCIPgetTransformedVar(scip, var, &transvar) );
    SCIP_CALL( SCIPcaptureVar(scip, transvar) );
    targetdata->varSequence.push_back(transvar);
  }

  for (auto iter : sourcedata->varMap)
  {
    SCIP_VAR* transvar = NULL;
    SCIP_CALL( SCIPgetTransformedVar(scip, iter.second, &transvar) );
    SCIP_CALL( SCIPcaptureVar(scip, transvar) );
    targetdata->varMap[iter.first] = transvar;
  }

  SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
    SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
    SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),  SCIPconsIsLocal(sourcecons),
    SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons),
    SCIPconsIsStickingAtNode(sourcecons)) );

  return SCIP_OKAY;
}

SCIP_DECL_CONSSEPALP(ConshdlrLABSnogood::scip_sepalp)
{
  assert(result);

  for (int c = 0; c < nconss; ++c)
    SCIP_CALL( separate(scip, conshdlr, conss[c], NULL, *result) );

  if( *result == SCIP_DIDNOTFIND )
    *result = SCIP_FEASIBLE;

  return SCIP_OKAY;
}

SCIP_DECL_CONSSEPASOL(ConshdlrLABSnogood::scip_sepasol)
{
  assert(result);

  for (int c = 0; c < nconss; ++c)
    SCIP_CALL( separate(scip, conshdlr, conss[c], sol, *result) );

  if( *result == SCIP_DIDNOTFIND )
    *result = SCIP_FEASIBLE;

  return SCIP_OKAY;
}

SCIP_DECL_CONSENFOLP(ConshdlrLABSnogood::scip_enfolp)
{
  assert(result);

  for (int c = 0; c < nconss; ++c)
    SCIP_CALL( separate(scip, conshdlr, conss[c], NULL, *result) );

  if( *result == SCIP_DIDNOTFIND )
    *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

SCIP_DECL_CONSENFOPS(ConshdlrLABSnogood::scip_enfops)
{
  assert(result != NULL);

  *result = SCIP_INFEASIBLE;

  return SCIP_OKAY;
}

SCIP_DECL_CONSCHECK(ConshdlrLABSnogood::scip_check)
{
  assert(result != NULL);
  *result = SCIP_FEASIBLE;

  for (int c = 0; c < nconss; ++c)
  {
    SCIP_CALL( check(scip, conss[c], sol, *result) );
    if (*result == SCIP_INFEASIBLE)
      break;
  }

  return SCIP_OKAY;
}

SCIP_DECL_CONSPROP(ConshdlrLABSnogood::scip_prop)
{
  assert(result != NULL);

  *result = SCIP_DIDNOTRUN;
  return SCIP_OKAY;
}

SCIP_DECL_CONSLOCK(ConshdlrLABSnogood::scip_lock)
{
  SCIP_CONSDATA* consdata = SCIPconsGetData(cons);
  assert(consdata != NULL);

  for (auto var : consdata->varSequence)
    SCIP_CALL( SCIPaddVarLocksType(scip, var, SCIP_LOCKTYPE_MODEL, nlocksneg + nlockspos, nlocksneg + nlockspos) );
  for (auto iter : consdata->varMap)
    SCIP_CALL( SCIPaddVarLocksType(scip, iter.second, SCIP_LOCKTYPE_MODEL, nlocksneg + nlockspos, nlocksneg + nlockspos) );

  return SCIP_OKAY;
}

SCIP_DECL_CONSDELVARS(ConshdlrLABSnogood::scip_delvars)
{
  return SCIP_OKAY;
}

SCIP_DECL_CONSPRINT(ConshdlrLABSnogood::scip_print)
{
  SCIP_CONSDATA* consdata = SCIPconsGetData(cons);
  assert(consdata != NULL);

  SCIPinfoMessage(scip, file, "LABSnogood for length %zu and shift %d.\n", consdata->varSequence.size(),
    consdata->shift);

  return SCIP_OKAY;
}

SCIP_DECL_CONSHDLRCLONE(scip::ObjProbCloneable* ConshdlrLABSnogood::clone)
{
   assert(valid != NULL);

   *valid = TRUE;
   return new ConshdlrLABSnogood(scip);
}

SCIP_DECL_CONSCOPY(ConshdlrLABSnogood::scip_copy)
{
  assert(valid != NULL);

#if 0

  /* Find the constraint handler */
  SCIP_CONSHDLR* conshdlr = SCIPfindConshdlr(scip, "LABSnogood");
  if(!conshdlr)
  {
    SCIPerrorMessage("LABSnogood constraint handler not found.\n");
    return SCIP_PLUGINNOTFOUND;
  }

  SCIP_CONSDATA* sourceconsdata = SCIPconsGetData(sourcecons);

  SCIP_CONSDATA* consdata = new SCIP_CONSDATA;
  consdata->shift = sourceconsdata->shift;
  consdata->varSequence.resize(sourceconsdata->varSequence.size());
  for (int i = 0; i < sourceconsdata->varSequence.size(); ++i)
  {
    SCIPdebugMsg(scip, "Copy of <%s>...\n", SCIPvarGetName(sourceconsdata->varSequence[i]));
    fflush(stdout);
    SCIP_VAR* copiedvar = static_cast<SCIP_VAR*>(SCIPhashmapGetImage(varmap, sourceconsdata->varSequence[i]));
    SCIPdebugMsg(scip, "... is <%s>.\n", SCIPvarGetName(copiedvar));
    fflush(stdout);

    consdata->varSequence[i] = copiedvar;
    SCIP_CALL( SCIPcaptureVar(scip, copiedvar) );
  }
  for (auto iter : sourceconsdata->varMap)
  {
    SCIPdebugMsg(scip, "Copy of <%s>...\n", SCIPvarGetName(iter.second));
    fflush(stdout);
    SCIP_VAR* copiedvar = static_cast<SCIP_VAR*>(SCIPhashmapGetImage(varmap, iter.second));
    SCIPdebugMsg(scip, "... is <%s>.\n", SCIPvarGetName(copiedvar));
    fflush(stdout);
    consdata->varMap[iter.first] = copiedvar;
    SCIP_CALL( SCIPcaptureVar(scip, copiedvar) );
  }

  SCIP_CALL( SCIPcreateCons(scip, cons, (name == NULL) ? SCIPconsGetName(sourcecons) : name, conshdlr, consdata,
    initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE) );

  *valid = TRUE;

#endif

  return SCIP_OKAY;
}

SCIP_RETCODE ConshdlrLABSnogood::check(SCIP* scip, SCIP_CONS* cons, SCIP_SOL* sol, SCIP_RESULT& result)
{
  SCIP_CONSDATA* consdata = SCIPconsGetData(cons);
  assert(consdata);

  std::vector<int> valsSequence(consdata->varSequence.size());
//   SCIPdebugMsg(scip, "Checking for sequence of length %zu with shift %d.\n", valsSequence.size(), consdata->shift);
  for (int i = 0; i < valsSequence.size(); ++i)
  {
    double x = SCIPgetSolVal(scip, sol, consdata->varSequence[i]);
    valsSequence[i] = x > 0.5;
//     SCIPdebugMsg(scip, "  x[%d] = %d\n", i, valsSequence[i]);
  }

  int correct_value = 0;
  for (int i = 0; i + consdata->shift < consdata->varSequence.size(); ++i)
  {
    int s1 = 2 * valsSequence[i] - 1;
    int s2 = 2 * valsSequence[i + consdata->shift] - 1;
    correct_value += s1 * s2;
  }

  for (auto iter : consdata->varMap)
  {
    if (((iter.first == correct_value) && (SCIPgetSolVal(scip, sol, iter.second) < 0.5))
      || ((iter.first != correct_value) && (SCIPgetSolVal(scip, sol, iter.second) > 0.5)))
    {
      result = SCIP_INFEASIBLE;
    }
  }

  return SCIP_OKAY;
}

SCIP_RETCODE ConshdlrLABSnogood::separate(SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS* cons, SCIP_SOL* sol,
  SCIP_RESULT& result)
{
  SCIP_CONSDATA* consdata = SCIPconsGetData(cons);
  assert(consdata);

  std::vector<int> valsSequence(consdata->varSequence.size());
  double distanceToBinary = 0.0;
  SCIPdebugMsg(scip, "Separating for sequence of length %zu with shift %d.\n", valsSequence.size(), consdata->shift);
  int num1s = 0;
  for (int i = 0; i < valsSequence.size(); ++i)
  {
    double x = SCIPgetSolVal(scip, sol, consdata->varSequence[i]);
    valsSequence[i] = x > 0.5;
    distanceToBinary += fabs(x - valsSequence[i]);
    SCIPdebugMsg(scip, "  x[%d] = %d\n", i, valsSequence[i]);
    num1s += valsSequence[i];
  }

  int correct_value = 0;
  for (int i = 0; i + consdata->shift < consdata->varSequence.size(); ++i)
  {
    int s1 = 2 * valsSequence[i] - 1;
    int s2 = 2 * valsSequence[i + consdata->shift] - 1;
    correct_value += s1 * s2;
  }

  for (auto iter : consdata->varMap)
  {
    double slack = distanceToBinary;
    bool match;
    if (correct_value == iter.first)
    {
      slack += SCIPgetSolVal(scip, sol, iter.second);
      match = true;
    }
    else
    {
      slack += 1.0 - SCIPgetSolVal(scip, sol, iter.second);
      match = false;
    }

    if (SCIPisLT(scip, slack, 1.0))
    {
      SCIP_ROW* row = NULL;
      char name[SCIP_MAXSTRLEN];
      SCIPsnprintf(name, SCIP_MAXSTRLEN, "nogood");

      SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &row, conshdlr, name, (match ? 1.0 : 0.0) - num1s,
        SCIPinfinity(scip), FALSE, FALSE, TRUE) );
      SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
      for (int i = 0; i < valsSequence.size(); ++i)
        SCIP_CALL( SCIPaddVarToRow(scip, row, consdata->varSequence[i], valsSequence[i] ? -1.0 : +1.0) );
      SCIP_CALL( SCIPaddVarToRow(scip, row, iter.second, match ? +1.0 : -1.0) );
      SCIP_CALL( SCIPflushRowExtensions(scip, row) );

      SCIP_Bool cutoff;
      SCIP_CALL( SCIPaddRow(scip, row, TRUE, &cutoff) );
      SCIP_CALL( SCIPreleaseRow(scip, &row));

      if (cutoff)
        result = SCIP_CUTOFF;
      else
        result = SCIP_SEPARATED;
    }
  }

  return SCIP_OKAY;
}

SCIP_RETCODE SCIPcreateConsLABSnogood(
  SCIP* scip,                                 /**< SCIP data structure */
  SCIP_CONS** cons,                           /**< pointer to hold the created constraint */
  const char* name,                           /**< name of constraint */
  int length,                                 /**< Length of sequences. */
  int shift,                                  /**< Shift distance. */
  SCIP_VAR** varSequence,                     /**< Sequence variables. */
  std::unordered_map<int, SCIP_VAR*>& varMap, /**< Mapping of possible values to variables. */
  SCIP_Bool initial,                          /**< should the LP relaxation of constraint be in the initial LP? */
  SCIP_Bool separate,                         /**< should the constraint be separated during LP processing? */
  SCIP_Bool enforce,                          /**< should the constraint be enforced during node processing? */
  SCIP_Bool check,                            /**< should the constraint be checked for feasibility? */
  SCIP_Bool propagate,                        /**< should the constraint be propagated during node processing? */
  SCIP_Bool local,                            /**< is constraint only valid locally? */
  SCIP_Bool modifiable,                       /**< is constraint modifiable (subject to column generation)? */
  SCIP_Bool dynamic,                          /**< is constraint dynamic? */
  SCIP_Bool removable                         /**< should the constraint be removed from the LP due to aging or cleanup? */
  )
{
  SCIP_CONSHDLR* conshdlr = SCIPfindConshdlr(scip, "LABSnogood");
  if (!conshdlr)
  {
    SCIPerrorMessage("LABSnogood constraint handler not found.\n");
    return SCIP_PLUGINNOTFOUND;
  }

  SCIP_CONSDATA* consdata = new SCIP_CONSDATA;
  consdata->shift = shift;
  for (int i = 0; i < length; ++i)
  {
    SCIP_VAR* var = varSequence[i];
    consdata->varSequence.push_back(var);
    SCIP_CALL( SCIPcaptureVar(scip, var) );
  }
  for (auto iter : varMap)
  {
    consdata->varMap[iter.first] = iter.second;
    SCIP_CALL( SCIPcaptureVar(scip, iter.second) );
  }

  SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
    local, modifiable, dynamic, removable, FALSE) );

  return SCIP_OKAY;
}

