/*------------------------------------------------------------------------*/
/*! \file specpoly.h
    \brief Used to generate specification polynomials

 Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_SPECPOLY_H_
#define TALISMAN_SRC_SPECPOLY_H_
/*------------------------------------------------------------------------*/
#include <stdarg.h>
#include <vector>

#include "gate.h"
#include "polynomial.h"
/*------------------------------------------------------------------------*/
Polynomial *mult_spec_poly();
Polynomial *miter_spec_poly();
Polynomial *assertion_spec_poly();

/**
    Reads the input spec given in the file called i

    @return Polynomial *
*/
Polynomial *parse_specification_polynomial(const char *file_name);

#endif // TALISMAN_SRC_SPECPOLY_H_
