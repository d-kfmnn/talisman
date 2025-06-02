/*------------------------------------------------------------------------*/
/*! \file vanishing_constraints.h
    \brief contains functions used in the polynomial solver

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_VANISHING_H_
#define TALISMAN_SRC_VANISHING_H_
/*------------------------------------------------------------------------*/
#include <vector>

#include "gate.h"
#include "pac.h"
#include "reductionmethods.h"
/*------------------------------------------------------------------------*/


void find_vanishing_constraints();
void find_vanishing_constraints_light();
/*------------------------------------------------------------------------*/

#endif // TALISMAN_SRC_VANISHING_H_
