/*------------------------------------------------------------------------*/
/*! \file monomial.h
    \brief contains the class Monomial and further functions to
    manipulate monomials

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_MONOMIAL_H_
#define TALISMAN_SRC_MONOMIAL_H_
/*------------------------------------------------------------------------*/
#include <gmp.h>

#include "term.h"
/*------------------------------------------------------------------------*/

/** \class Monomial
    This class is used to represent monomials in a polynomial.
    A monomial consist of a coefficient and a term.
*/

class Monomial {
  Term* term;  // term
  unsigned ref;// reference counter

  unsigned get_ref() const { return ref; }
  unsigned dec_ref() { return --ref; }

  public:
  mpz_t coeff;// coefficient

  /** Constructor

      @param c mpz_t coefficient
      @param t Term*
  */
  Monomial(mpz_t _c, Term* _t);

  /** Getter for member term

      @return Term*
  */
  Term* get_term() const { return term; }

  /** Getter for member term, calls copy routine of Term

      @return a copy of Term* term
  */
  Term* get_term_copy() const { return term->copy(); }// TODO remove

  /** Returns the size fo the term

      @return unsigned, the size of the term
  */
  unsigned get_term_size() const { return term->degree(); }// TODO remove

  /**
  Copy routine

  @return Same monomial with increased reference counter
  */
  Monomial* copy();

  /**
  Printing routine

  @param file Output file
  @param lm bool indicating whether we print a leading monomial
  */
  void print(FILE* file, bool lm = 0) const;

  /** Destructor */
  ~Monomial();

  int evalute() {
    int res = this->term->evalute();
    return static_cast<int>(mpz_get_si(coeff)) * res;
  }

  friend void deallocate_monomial(Monomial* m);
};
/*------------------------------------------------------------------------*/

/**
  Multiplies two monomials

  @param m1 Monomial*
  @param m2 Monomial*

  @return Product of m1*m2
*/
Monomial*
multiply_monomial(const Monomial* m1, const Monomial* m2);

/**
  Wrapper for deconstructor, reduces the references of m until 0.
  Then deconstructor is called.

  @param m Monomial*
*/
void
deallocate_monomial(Monomial* m);

#endif// TALISMAN_SRC_MONOMIAL_H_
