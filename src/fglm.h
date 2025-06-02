/*------------------------------------------------------------------------*/
/*! \file fglm.h
    \brief contains functions used to linearize via fglm or gap


  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
  
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_GB_H_
#define TALISMAN_SRC_GB_H_
/*------------------------------------------------------------------------*/
#include <set>
#include <sstream>
#include <string.h>
#include <vector>

#include "gate.h"
#include "polynomial.h"
#include "propagate.h"
#include "subcircuit.h"
#include "specpoly.h"

/*------------------------------------------------------------------------*/

bool linearize_via_fglm_or_gap(Gate *g);

#endif // TALISMAN_SRC_GB_H_
