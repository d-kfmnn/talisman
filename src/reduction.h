/*------------------------------------------------------------------------*/
/*! \file reduction.h
    \brief contains functions used in the polynomial solver

  This file contains the main reduction loop

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_REDUCTION_H_ 
#define TALISMAN_SRC_REDUCTION_H_
/*------------------------------------------------------------------------*/
#include <string.h>
#include <vector>

#include "fglm.h"
#include "substitution.h"
/*------------------------------------------------------------------------*/
// / 1 for pac
// / 2 for hybrid
// / 3 for nss
extern int proof;



/*------------------------------------------------------------------------*/

Polynomial *reduce(Polynomial *spec);

#endif // TALISMAN_SRC_REDUCTION_H_
