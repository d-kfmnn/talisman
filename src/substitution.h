/*------------------------------------------------------------------------*/
/*! \file substitution.h
    \brief Used to identify special subcircuits

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_SUBSTITUTION_H_
#define TALISMAN_SRC_SUBSTITUTION_H_
/*------------------------------------------------------------------------*/
#include <vector>

#include "gate.h"
/*------------------------------------------------------------------------*/


/**
    Routine for identifying the final stage adder

    @return True if all checks succeeded
*/
bool identify_final_stage_adder();

void mark_bottom_of_circuit(Gate * g);

void unmark_fsa();

/*------------------------------------------------------------------------*/





#endif  // TALISMAN_SRC_SUBSTITUTION_H_
