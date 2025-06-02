/*------------------------------------------------------------------------*/
/*! \file witness.h
    \brief contains the polynomial solving routine

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_WITNESS_H_
#define TALISMAN_SRC_WITNESS_H_
/*------------------------------------------------------------------------*/
#include "gate.h"
/*------------------------------------------------------------------------*/
// / If final remainder is not equal to zero a counter example is generated and
// / printed to file <input_name>.wit, default is true, can be turned of
//  using command line input
extern bool gen_witness;
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
// Functions used to find counter examples
/**
    Checks whether the polynomial contains only input variables

    @param p Polynomial* to be checked

    @return true if p contains only input variables
*/
bool check_inputs_only(const Polynomial *p);

/**
    Writes a vector to the file such that all input variables in the term t
    are set to 1, all other to 0.

    @param t Term
    @param file output file for the counter example
*/
void write_witness_vector(const Term *t);

/**
    Search for the smallest term and calls write_witness_vector on that

    @param p Polynomial* for which counter examples are generated
    @param file output file for the counter example
*/
void write_witnesses(const Polynomial *p);

/**
    Generates a witness for the remainder polynomial p, by identifying
    the smallest term in the polynomial and setting all its variables
    to 1.

    @param p Polynomial* for which counter examples are generated
    @param name prefix name for the output file, suffix is '.cex'
*/
void generate_witness(const Polynomial *p, const char *name);

#endif // TALISMAN_SRC_WITNESS_H_
