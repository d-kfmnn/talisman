/*------------------------------------------------------------------------*/
/*! \file propagate.cpp
    \brief contains functions used in the polynomial solver

  This file contains all functions used for propagating equivalent polynomials
  and constants

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#include "propagate.h"
/*------------------------------------------------------------------------*/
static void rewrite_parents(Gate *g, Polynomial *p)
{
  for (auto &g_parent : g->get_parents())
  {

    Polynomial *flip = unflip_poly(g_parent->get_gate_constraint());
    Polynomial *tmp = reduce_by_one_poly(flip, p);
    delete (flip);
    g_parent->update_gate_poly(tmp, 0);
    g_parent->children_remove(g);
    check_if_propagate(tmp);
  }
}

/*------------------------------------------------------------------------*/
static bool try_propagate_constant_zero(Polynomial *p)
{
  // p:= ax for a in int
  assert(p->degree() == 1 && p->len() == 1);
  Gate *g = gate(p->get_lt()->get_var_num());
  if (verbose > 1)
    msg("found constant 0: %s", g->get_var_name());

  rewrite_parents(g, p);
  return 1;
}

/*------------------------------------------------------------------------*/
static bool try_propagate_constant_one(Polynomial *p)
{
  // p:= ax-a for a in int
  assert(p->degree() == 1 && p->len() == 2);
  assert(!p->get_mon(1)->get_term());

  mpz_t neg;
  mpz_init(neg);
  mpz_neg(neg, p->get_mon(1)->coeff);

  if (mpz_cmp(p->get_lm()->coeff, neg) != 0)
  {
    mpz_clear(neg);
    return 0;
  }
  mpz_clear(neg);

  Gate *g = gate(p->get_lt()->get_var_num());
  if (verbose > 1)
    msg("found constant 1: %s", g->get_var_name());

  rewrite_parents(g, p);
  return 1;
}

/*------------------------------------------------------------------------*/
static bool try_propagate_equality(Polynomial *p)
{
  // p:= ax-ay for a in int
  assert(p->degree() == 1 && p->len() == 2);

  mpz_t neg;
  mpz_init(neg);
  mpz_neg(neg, p->get_mon(1)->coeff);

  if (mpz_cmp(p->get_lm()->coeff, neg) != 0)
  {
    mpz_clear(neg);
    return 0;
  }
  mpz_clear(neg);

  Gate *g = gate(p->get_lt()->get_var_num());
  if (verbose > 1)
    msg("found equality: %s", g->get_var_name());

  rewrite_parents(g, p);
  return 1;
}
/*------------------------------------------------------------------------*/
static bool try_propagate_negated_equality(Polynomial *p)
{
  // p:= ax+ay-a for a in int

  assert(p->degree() == 1 && p->len() == 3);
  if (p->get_mon(2)->get_term())
    return 0;

  if (mpz_cmp(p->get_lm()->coeff, p->get_mon(1)->coeff) != 0)
    return 0;

  mpz_t neg;
  mpz_init(neg);
  mpz_neg(neg, p->get_mon(2)->coeff);

  if (mpz_cmp(p->get_mon(1)->coeff, neg) != 0)
  {
    mpz_clear(neg);
    return 0;
  }
  mpz_clear(neg);

  Gate *g = gate(p->get_lt()->get_var_num());
  if (verbose > 1)
    msg("found negated equality: %s", g->get_var_name());

  rewrite_parents(g, p);
  return 1;
}

/*------------------------------------------------------------------------*/
bool check_if_propagate(Polynomial *p)
{
  assert(p->len() > 0);

  if (p->degree() > 1)
    return 0;
  if (p->len() > 3)
    return 0;

  if (p->len() == 1)
    return try_propagate_constant_zero(p);
  else if (p->len() == 2)
  {
    if (!p->get_tail_term())
      return try_propagate_constant_one(p);
    else
      return try_propagate_equality(p);
  }
  else
  {
    return try_propagate_negated_equality(p);
  }
  return 0;
}