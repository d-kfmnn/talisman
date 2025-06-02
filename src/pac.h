/*------------------------------------------------------------------------*/
/*! \file pac.h
    \brief contains functions necessary to generate PAC proofs

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_PAC_H_
#define TALISMAN_SRC_PAC_H_
/*------------------------------------------------------------------------*/

#include "gate.h"
/*------------------------------------------------------------------------*/

extern Var * var_x;
extern Term * term_x;

/**
    Prints all initial gate constraints to the file(with indices)

    @param file output file
*/
void print_circuit_poly(FILE * file);

void init_proof_logging(size_t level);

void print_refutation_spec(FILE * file);

int print_pac_pattern_out_rules(FILE* file, std::vector<Polynomial*> lin_poly, int i);

void pac_add_circuit_poly(FILE * file, Polynomial * p);
/**
    Prints a deletion rule

    @param p1 Polynomial*
    @param file output file
*/
void print_pac_del_rule(FILE * file, const Polynomial *p1);

/**
    Prints the modulo rule

    @param file output file
    @param p1   Polynomial*, factor
    @param p    Polynomial*, conclusion
*/
void print_pac_mod_rule(FILE * file, const Polynomial *p1, Polynomial *p);

/**
    Prints an addition rule for pac

    @param file output file
    @param p1   Polynomial*, First summand
    @param p2   Polynomial*, Second summand
    @param p    Polynomial*, Conclusion
*/
void print_pac_add_rule(
  FILE * file, const Polynomial *p1, const Polynomial *p2, Polynomial *p);

  /**
      Prints a combination rule for pac sum(indices) = p

      @param file output file
      @param indices std::vector<int>
      @param p    Polynomial*, Conclusion
  */
void print_pac_vector_add_rule(
  FILE * file, std::vector<int> indices, Polynomial * p);


/**
    Prints a combination rule for pac p1*p2+p3=p

    @param file output file
    @param p    Polynomial*, Conclusion
    @param p1   Polynomial*, first poly
    @param p2   Polynomial*, first factor
    @param p3   Polynomial*, second poly
    @param p3   Polynomial*, second factor
*/
void print_pac_combi_rule(
  FILE * file, const Polynomial *p1, const Polynomial *p2,
  const Polynomial *p3, const Polynomial *p4, Polynomial *p);

  void print_pac_combi_monomial_rule(
    FILE * file, const Polynomial *p1, const Monomial *m2,
    const Polynomial *p3, const Polynomial *p4, Polynomial *p);

  /**
      Prints a combination rule for pac <indices>*<co_factors> = p

      @param file output file
      @param indices std::vector<int>
      @param co_factors std::vector<const Polynomial*>
      @param p    Polynomial*, Conclusion
  */
void print_pac_vector_combi_rule(
  FILE * file, std::vector<int> indices,
  std::vector<const Polynomial*> co_factors, Polynomial * p);

/**
    Prints the multiplication rule of pac

    @param file output file
    @param p1   Polynomial*, First factor
    @param p2   Polynomial*, Second factor
    @param p    Polynomial*, Conclusion
*/
void print_pac_mul_rule(
  FILE * file, const Polynomial *p1, const Polynomial *p2, Polynomial *p);

  void print_pac_mul_const_rule(
    FILE * file, const Polynomial *p1, int n, Polynomial *p);

void print_pac_extension_rule_for_mon(FILE * file,Gate * g, const Term * t, Polynomial *p);

void print_dual_constraints(FILE * file);

#endif  // TALISMAN_SRC_PAC_H_
