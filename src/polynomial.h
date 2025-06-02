/*------------------------------------------------------------------------*/
/*! \file polynomial.h
    \brief contains arithmetic operations for polynomials

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_POLYNOMIAL_H_
#define TALISMAN_SRC_POLYNOMIAL_H_
/*------------------------------------------------------------------------*/
#include <cstring>
#include <deque>
#include <list>

#include "monomial.h"
/*------------------------------------------------------------------------*/
extern size_t running_idx;

/** \class Polynomial
    This class is used to polynomials.
*/
class Polynomial {

  Monomial** mon;

  size_t num_mon = 0;
  size_t idx;
  int ref = 1;

  int deg = 0;

  public:
  /** Constructor */
  Polynomial();

  Polynomial(Monomial** m, size_t len, int d);

  size_t len() const { return num_mon; }

  size_t degree() const { return deg; }

  Monomial* get_lm() const { return mon[0]; }
  Monomial* get_mon(size_t i) const;

  Term* get_lt() const { return mon[0]->get_term(); }
  Term* get_tail_term() const { return mon[1]->get_term(); }

  Term* get_largest_term() const;
  Monomial* get_largest_mon() const;

  Polynomial* get_tail_poly() const;

  Var* contains_dual_var() const;

  Polynomial* copy();


  size_t get_idx() const { return idx; }

  void set_idx(size_t i) { idx = i; }

  /**
      Printing routine

      @param file Output file
      @param end if true we print trailing ";"
  */
  void print(FILE* file, bool end = 1) const;

  /**
      Returns whether the polynomial is the constant zero polynomial

      @return bool
  */
  bool is_constant_zero_poly() const;

  /**
      Returns whether the polynomial is the constant one polynomial

      @return bool
  */
  bool is_constant_one_poly() const;

  /**
      Returns the size of the smallest term

      @return integer
  */
  int min_term_size() const;

  /** Destructor */
  ~Polynomial();

  int evaluate() {
    int res = 0;
    for(size_t i = 1; i < this->len(); i++) {
      res += this->get_mon(i)->evalute();
    }
    return res;
  }
};
/*------------------------------------------------------------------------*/
// Polynomials are generated using a sorted array "mstack"

/**
    Deallocates mstack
*/
void
deallocate_mstack();

/**
    Pushes a monomial to the end of the stack

    @param t monomial to be added to the mstack
*/
void
push_mstack_end(Monomial* t);

/**
    Pushes a monomial to the stack such that mstack remains sorted

    @param t monomial to be added to the mstack
*/
void
push_mstack(Monomial* t);

/*------------------------------------------------------------------------*/
/**
    Generates a polynomial from mstack and clears mstack

    @return Polynomial*
*/
Polynomial*
build_poly();// TODO maybe from vector of monomials instead of global stack

// Generates the constraint -v-v_+1
Polynomial*
gen_dual_constraint(Var* v);

/*------------------------------------------------------------------------*/
bool
equal_poly(Polynomial* p1, Polynomial* p2);

int
cmp_poly(Polynomial* p1, Polynomial* p2);

/*---------------------------------------------------------------------------*/

/**
    Add two polynomials p1+p2

    @param p1 Polynomial*
    @param p2 Polynomial*

    @return Polynomial*, sum of p1+p2
*/
Polynomial*
add_poly(Polynomial* p1, Polynomial* p2);

/**
    Add two polynomials p1+p2

    @param p1 Polynomial*
    @param p2 Polynomial*

    @return Polynomial*, sum of p1-p2
*/
Polynomial*
sub_poly(Polynomial* p1, Polynomial* p2);

/**
    Multiplies two polynomials p1*p2

    @param p1 Polynomial*
    @param p2 Polynomial*

    @return Polynomial*, product of p1*p2
*/
Polynomial*
multiply_poly(Polynomial* p1, Polynomial* p2);

/**
    Multiplies a polynomial p1 with a constant c

    @param p1: Polynomial*
    @param c:  mpz_t object

    @return Polynomial*, product of c*p1
*/
Polynomial*
multiply_poly_with_constant(Polynomial* p1, mpz_t c);

/**
    Multiplies a polynomial p1 with a constant c

    @param p1: Polynomial*
    @param c:  mpz_t object

    @return Polynomial*, product of c*p1
*/
Polynomial*
multiply_poly_with_term(Polynomial* p1, Term* t);

/**
    Multiplies a polynomial p1 with a constant c

    @param p1: Polynomial*
    @param c:  mpz_t object

    @return Polynomial*, product of c*p1
*/
Polynomial*
multiply_poly_with_monomial(Polynomial* p1, Monomial* m);

/**
    Returns the quotient of dividing a polynomial p1 by a term t

    @param p1 Polynomial*
    @param t  Term *

    @return Polynomial defining the quotient of p1/t
*/
Polynomial*
divide_poly_by_term(Polynomial* p1, const Term* t);
Polynomial*
divide_by_var(Polynomial* p1, const Term* t);

/*---------------------------------------------------------------------------*/

// / gmp for 1
extern mpz_t one;

// / gmp for -1
extern mpz_t minus_one;

// / gmp for -2
extern mpz_t minus_two;

// / gmp for 2
extern mpz_t base;

// / gmp for 2^NN
extern mpz_t mod_coeff;
/*---------------------------------------------------------------------------*/
/**
    Initializes all global gmp objects

    @param exp unsigned exponent for mod coeff
*/
void
init_mpz(unsigned exp);

/**
    Clears all globally allocated gmp objects
*/
void
clear_mpz();

#endif// TALISMAN_SRC_POLYNOMIAL_H_
