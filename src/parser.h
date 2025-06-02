/*------------------------------------------------------------------------*/
/*! \file parser.h
    \brief contains functions necessary to parse the AIG

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_PARSER_H_
#define TALISMAN_SRC_PARSER_H_
/*------------------------------------------------------------------------*/
#include "aig.h"
/*------------------------------------------------------------------------*/

/**
    Reads the input aiger given in the file called input_name to the aiger 'model'
    using the parserer function of aiger.h

    @param input_name char * ame of input file
*/
void parse_aig(const char *input_name);

#endif // TALISMAN_SRC_PARSER_H_
