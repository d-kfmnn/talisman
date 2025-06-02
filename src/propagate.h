/*------------------------------------------------------------------------*/
/*! \file propagate.h
    \brief contains functions to check whether a var is constant 

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_PROPAGATE_H_ 
#define TALISMAN_SRC_PROPAGATE_H_
/*------------------------------------------------------------------------*/
#include "gate.h"
#include "polynomial.h"
#include "reductionmethods.h"
/*------------------------------------------------------------------------*/

bool check_if_propagate(Polynomial *p); 
/*------------------------------------------------------------------------*/
#endif // TALISMAN_SRC_PROPAGATE_H_
