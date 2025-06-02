/*------------------------------------------------------------------------*/
/*! \file extensions.h
    \brief contains functions to add extension variables

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_EXT_H_
#define TALISMAN_SRC_EXT_H_
/*------------------------------------------------------------------------*/
#include "gate.h"
#include "pac.h"

Term *extend_var_gates(Term *t);

void adjust_level_of_extended_gates();

#endif // TALISMAN_SRC_EXT_H_