/*------------------------------------------------------------------------*/
/*! \file term.h
    \brief contains the class Term and further functions to
    manipulate terms

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_TERM_H_
#define TALISMAN_SRC_TERM_H_
/*------------------------------------------------------------------------*/
#include <algorithm>
#include <list>
#include <vector>

#include "signal_statistics.h"
#include "variable.h"
/*------------------------------------------------------------------------*/

/** \class Term
    This class is used to represent terms in a polynomial.
    Terms are represented as ordered linked lists of variables.
*/

class Term {
  // / head variable
  Var* variable;

  // / tail in linked list
  Term* rest;

  // / reference counter
  uint64_t ref;

  // / hash value
  const uint64_t hash;

  // / hash collision chain link
  Term* next;

  // / Total degree of term
  const size_t deg;

  public:
  /** Constructor

      @param _v Var*
      @param _r Term*
      @param _hash uint64_t
      @param _n Term*
  */
  Term(Var* _v, Term* _r, uint64_t _hash, Term* _n);

  /**
    Copy routine

    @return A copy of the current term
*/
  Term* copy();

  /**
      Printing routine

      @param file Output file
  */
  void print(FILE* file) const;
  void print_orig(FILE* file) const;

  // GETTER & SETTER
  /** Getter for member variable

      @return Var*
  */
  Var* get_var() const { return variable; }

  /** Getter for level of variable

      @return integer
  */
  int get_var_level() const { return variable->get_level(); }

  /** Getter for num of variable

      @return integer
  */
  int get_var_num() const { return variable->get_num(); }

  /** Getter for member rest

      @return Term*
  */
  Term* get_rest() const { return rest; }

  /** Getter for member hash

      @return uint64_t
  */
  uint64_t get_hash() const { return hash; }

  /** Getter for member next

      @return Term*
  */
  Term* get_next() const { return next; }

  /** Setter for member term

      @param t Term*
  */
  void set_next(Term* t) { next = t; }

  /** Getter for member ref

      @return uint64_t
  */
  uint64_t get_ref() const { return ref; }

  /** Increases ref

      @return uint64_t
  */
  uint64_t inc_ref() { return ++ref; }

  /** Decreases ref

      @return uint64_t
  */
  uint64_t dec_ref() { return --ref; }

  /** Getter for member deg

      @return size_t
  */
  size_t degree() const { return deg; }

  /**
      Checks whether v is contained in Term

      @param v Var*

      @return true if v is contained in term
  */
  bool contains(Var* v) const;

  /**
      Checks whether t is contained in Term

      @param t Term*

      @return true if t is contained in term
  */
  bool contains_subterm(const Term* t) const;

  /**
      Returns the first dual variable in term

      @return Var* or 0
  */
  Var* extract_first_dual_var() const;

  size_t count_dual() const;
  /*
  Evalute the term t.
  Returns -1 if the value of t cannot be determined
  (if some variables are not initalized)
*/
  int evalute() {
    int res = 1;
    Term* tmp = this;
    while(tmp and res != 0) {
      if(tmp->get_var()->get_value() == -1)
        die(5, "Trying to evaluate variable that was not set");
      res *= tmp->get_var()->get_value();
      tmp = tmp->get_rest();
    }
    return res;
  }
};

/*------------------------------------------------------------------------*/
// We organize terms in a hash table that is dynamically enlarged.
// Every time a new term is defined, we compute a hash value and insert
// the term. Terms are counted using a reference counter, which is incremented
// and decremented depending how often the term occurs in polynomials.

/**
    Builds a term, where variable is added at the front of rest

    @param variable Variable*
    @param rest     Term*

    @return Term*
*/
Term*
new_term(Var* variable, Term* rest = 0);

/**
    Builds a quadratic term v1*v2

    @param v1 Variable*
    @param v2 Variable*

    @return Term*
*/
Term*
new_quadratic_term(Var* v1, Var* v2);

/**
    adds Var* to vstack
*/
void
add_to_vstack(Var* v);
/**
    removes elements from vstack
*/
void
clear_vstack();

/**
    Generates a term from the variable stack

    @return Term* generated from the variable stack
*/
Term*
build_term_from_stack(bool sort = 0);

/**
    provides a vstack as parameter //TODO this should replace global vstack
*/
Term*
sort_and_build_term_from_vector(std::vector<Var*> v);

/**
    Decrements the reference count of a term, and actually deletes a
    term if its reference count goes to zero.  In this case it also
    removes it from the hash table and applies the same procedure to the
    suffix 'rest'.

    @param t Term*
*/
void
deallocate_term(Term* t);

/**
    Deallocates the hash table "term_table"
*/
void
deallocate_terms();

/*------------------------------------------------------------------------*/
// ORDERS and COMPARISON
/*------------------------------------------------------------------------*/
/**
    Compares this term to term t using the levels of the variables, based on lex
   ordering

    @param t Term*

    @return +1 if t1 > t2
            -1 if t1 < t2
            0  if t1 = t2
*/
int
cmp_term(const Term* t1, const Term* t2);

bool
equal_up_to_duality(const Term* t1, const Term* t2);




/*------------------------------------------------------------------------*/
// ARITHMETIC OPERATIONS
/*------------------------------------------------------------------------*/

Term*
multiply_term(Term* t1, const Term* t2);

Term*
multiply_term_by_var(Term* t1, Var* v);

Term*
divide_by_var(const Term* t, Var* v);

Term*
divide_by_term(Term* t, const Term* t1);

#endif// TALISMAN_SRC_TERM_H_
