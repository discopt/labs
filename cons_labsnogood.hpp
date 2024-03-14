#ifndef CONS_LABSNOGOOD_H_
#define CONS_LABSNOGOOD_H_

#include <unordered_map>

#include "objscip/objscip.h"

class ConshdlrLABSnogood : public scip::ObjConshdlr
{
public:
   /** default constructor */
   ConshdlrLABSnogood(SCIP* scip);

   /** destructor */
   virtual ~ConshdlrLABSnogood();

   /** frees specific constraint data */
   virtual SCIP_DECL_CONSDELETE(scip_delete);

   /** transforms constraint data into data belonging to the transformed problem */
   virtual SCIP_DECL_CONSTRANS(scip_trans);

   /** separation method of constraint handler for LP solution */
   virtual SCIP_DECL_CONSSEPALP(scip_sepalp);

   /** separation method of constraint handler for arbitrary primal solution */
   virtual SCIP_DECL_CONSSEPASOL(scip_sepasol);

   /** constraint enforcing method of constraint handler for LP solutions */
   virtual SCIP_DECL_CONSENFOLP(scip_enfolp);

   /** constraint enforcing method of constraint handler for pseudo solutions */
   virtual SCIP_DECL_CONSENFOPS(scip_enfops);

   /** feasibility check method of constraint handler for primal solution */
   virtual SCIP_DECL_CONSCHECK(scip_check);

   /** domain propagation method of constraint handler */
   virtual SCIP_DECL_CONSPROP(scip_prop);

   /** variable rounding lock method of constraint handler */
   virtual SCIP_DECL_CONSLOCK(scip_lock);

   /** variable deletion method of constraint handler */
   virtual SCIP_DECL_CONSDELVARS(scip_delvars);

   /** constraint display method of constraint handler */
   virtual SCIP_DECL_CONSPRINT(scip_print);

   /** returns whether the objective plugin is copyable */
   virtual SCIP_DECL_CONSHDLRISCLONEABLE(iscloneable)
   {
      return TRUE;
   }

   /** clone method which will be used to copy a objective plugin */
   virtual SCIP_DECL_CONSHDLRCLONE(scip::ObjProbCloneable* clone); /*lint !e665*/

   /** constraint copying method of constraint handler */
   virtual SCIP_DECL_CONSCOPY(scip_copy);

private:
   /** Actual check routine. */
   SCIP_RETCODE check(
      SCIP* scip,
      SCIP_CONS* cons,
      SCIP_SOL* sol,
      SCIP_RESULT& result
   );

   /** Actual separation routine. */
   SCIP_RETCODE separate(
      SCIP* scip,
      SCIP_CONSHDLR* conshdlr,
      SCIP_CONS* cons,
      SCIP_SOL* sol,
      SCIP_RESULT& result
   );
};

/** creates and captures a LABS nogood constraint. */
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
  );

#endif

