/*------------------------------------------------------------------------*/
/*! \file preprocessing.cpp
    \brief contains functions for preprocessing
    
  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#include "preprocessing.h"

#include <algorithm>

#include "reductionmethods.h"
#include "vanishing_constraints.h"
/*------------------------------------------------------------------------*/

std::vector<Gate *> resolved;

std::vector<Polynomial *> resolved_gc;

std::vector<Polynomial *> init_gc;

std::vector<Polynomial *> proof_poly;
Polynomial *target_proof;

std::vector<Term *> mult1v;
std::vector<Term *> remv;
std::vector<Term *> outerv;
std::vector<Polynomial *> target_tmp;
/*------------------------------------------------------------------------*/

static bool is_unit(Gate *g) {
  if (g->get_elim())
    return 0;
  Polynomial *gc = g->get_gate_constraint();

  if (gc->len() > 2)
    return 0;
  if (gc->len() == 1)
    return 1;

  Term *t = gc->get_tail_term();
  if (t->degree() == 1)
    return 1;
  return 0;
}
/*----------------------------------------------------------------------------*/

static void remove_var_from_stored_terms(std::vector<Term *> vec, Var *v) {
  for (auto it = vec.begin(); it != vec.end(); ++it) {
    Term *t = *it;
    if (!t->contains(v) && !t->contains(v->get_dual()))
      continue;
    while (t) {
      if (t->get_var() != v && t->get_var() != v->get_dual()) {
        add_to_vstack(t->get_var());
      } else if (t->get_var() == v) {
        clear_vstack();
        break;
      }
      t = t->get_rest();
    }
    Term *new_t = build_term_from_stack();
    *it = new_t;
  }
}
/*----------------------------------------------------------------------------*/
static void replace_var_from_stored_terms(std::vector<Term *> vec, Var *v, Var *w) {
  for (auto it = vec.begin(); it != vec.end(); ++it) {
    Term *t = *it;
    if (!t->contains(v) && !t->contains(v->get_dual()))
      continue;
    while (t) {
      if (t->get_var() != v && t->get_var() != v->get_dual()) {
        add_to_vstack(t->get_var());
      } else if (t->get_var() == v) {
        add_to_vstack(w);
      } else {
        add_to_vstack(w->get_dual());
      }
      t = t->get_rest();
    }
    Term *new_t = build_term_from_stack();
    *it = new_t;
  }
}
/*----------------------------------------------------------------------------*/
static void update_stored_terms(Gate *g) {
  Polynomial *gc = g->get_gate_constraint();
  if (gc->len() == 1) {
    remove_var_from_stored_terms(mult1v, g->get_var());
    remove_var_from_stored_terms(remv, g->get_var());
  } else if (is_unit(g)) {
    replace_var_from_stored_terms(mult1v, g->get_var(), gc->get_tail_term()->get_var());
    replace_var_from_stored_terms(remv, g->get_var(), gc->get_tail_term()->get_var());
  }
}
/*----------------------------------------------------------------------------*/
static void eliminate_by_one_gate(Gate *n1, Gate *n2) {
  Polynomial *tmp_p1 = n1->get_gate_constraint();
  Polynomial *flip = n2->get_dual_constraint();
  Polynomial *p1 = reduce_by_one_poly(tmp_p1, flip);
  Polynomial *p2 = n2->get_gate_constraint();
  if (!p1 || !p2) {
    delete (p1);
    return;
  }
  Polynomial *negfactor = divide_poly_by_term(p1, p2->get_lt());

  if (negfactor->is_constant_zero_poly()) {
    delete (p1);
    delete (negfactor);

    return;
  }

  Polynomial *mult = multiply_poly(negfactor, p2);
  Polynomial *rem = add_poly(p1, mult);

  n1->update_gate_poly(rem);
  // rem->print(stdout);

  delete (mult);
  delete (negfactor);
  delete (p1);
}
/*----------------------------------------------------------------------------*/
static void eliminate_unit_gate(Gate *n) {
  update_stored_terms(n);
  for (auto &n_child : n->get_children()) {
    n_child->parents_remove(n);
  }

  for (auto &n_parent : n->get_parents()) {
    msg("before her1");
    eliminate_by_one_gate(n_parent, n);
    n_parent->children_remove(n);

    for (auto &n_child : n->get_children()) {
      if (!n_parent->is_child(n_child))
        n_parent->children_push_back(n_child);
      if (!n_child->is_in_parents(n_parent))
        n_child->parents_push_back(n_parent);
    }

    if (is_unit(n_parent)) {
      eliminate_unit_gate(n_parent);
    } else if (n_parent->children_size() == 1 && n_parent->get_gate_constraint()->len() == 3) {
      Gate *tmp = n_parent->children_front();
      Polynomial *flip = tmp->get_dual_constraint();
      Polynomial *rem1 = reduce_by_one_poly(n_parent->get_gate_constraint(), flip);

      if (rem1->len() != 2) {
        delete (rem1);
        flip = gen_dual_constraint(tmp->get_var()->get_dual());
        rem1 = reduce_by_one_poly(n_parent->get_gate_constraint(), flip);
        delete (flip);
      }

      n_parent->update_gate_poly(rem1);
      eliminate_unit_gate(n_parent);
    }
  }

  if (verbose > 2)
    msg("removed unit %s", n->get_var_name());
}

/*----------------------------------------------------------------------------*/

//
static void remove_only_positives(size_t parent_limit = 0) {
  msg("remove only positives");

  int counter = 0;
  for (unsigned i = M - 1; i >= NN; i--) {  // do not even consider changing direction
    Gate *n = gates[i];
    if (parent_limit > 0 && n->parents_size() > parent_limit)
      continue;
    if (!parent_limit && n->parents_size() == 1)
      continue;
    if (n->get_pp())
      continue;
    if (n->get_input())
      continue;

    if (n->get_elim())
      continue;
    if (n->get_output() || n->get_aig_output())
      continue;

    if (n->get_gate_constraint()->len() > 2)
      continue;
    bool flag = 0;

    for (auto &n_parent : n->get_parents()) {
      if (n_parent->get_gate_constraint()->len() > 2) {
        flag = 1;
        break;
      }

      Monomial *m = n_parent->get_gate_constraint()->get_mon(1);
      if (!m->get_term()->contains(n->get_var())) {
        flag = 1;
        break;
      }
    }
    if (flag)
      continue;

    for (auto &n_child : n->get_children()) {
      n_child->parents_remove(n);
    }

    for (auto &n_parent : n->get_parents()) {
      Polynomial *rem = reduce_by_one_poly(n_parent->get_gate_constraint(), n->get_gate_constraint());
      delete (n_parent->get_gate_constraint());
      n_parent->set_gate_constraint(rem);
      // rem->print(stdout);

      for (auto &n_child : n->get_children()) {
        n_child->parents_push_back(n_parent);
        n_parent->children_push_back(n_child);
      }
      n_parent->children_remove(n);
    }



    counter++;
  }
  if (verbose >= 1)
    msg("removed %i positive gates", counter);
}

/*============================================================================*/
static void print_proof_van_constraint(Gate *g, Gate *andg) {
  Polynomial *g_tmp = unflip_poly(g->get_gate_constraint());
  Polynomial *and_tmp = unflip_poly(andg->get_gate_constraint());

  std::vector<int> indices;
  std::vector<const Polynomial *> co_factors;

  // Step 1:
  Polynomial *left0 = multiply_poly_with_term(and_tmp, g_tmp->get_lt());
  push_mstack(new Monomial(one, g_tmp->get_lt()));
  Polynomial *resp0 = build_poly();
  indices.push_back(and_tmp->get_idx());
  co_factors.push_back(resp0);

  Polynomial *right0 = multiply_poly(g_tmp, and_tmp->get_tail_poly());
  indices.push_back(g_tmp->get_idx());
  co_factors.push_back(and_tmp->get_tail_poly());

  indices.push_back(and_tmp->get_idx());
  push_mstack(new Monomial(minus_one, 0));
  Polynomial *fac = build_poly();
  co_factors.push_back(fac);

  Polynomial *result0 = add_poly(left0, right0);
  Polynomial *result1 = sub_poly(result0, and_tmp);

  left0->print(stdout);
  right0->print(stdout);
  result1->print(stdout);

  print_pac_vector_combi_rule(proof_file, indices, co_factors, result1);
  dual_van_poly.insert({result1->get_lt(), result1});
}
/*============================================================================*/
static void check_for_new_vanishing_combinations(Gate *repl, Gate *g) {
  if (g->children_size() != 2) return;

  Gate *ch1 = g->children_front();
  Gate *ch2 = g->children_back();
  if (ch2 == repl) std::swap(ch1, ch2);
  assert(ch1 == repl);

  for (auto &ch1_aig_p : ch1->get_aig_parents()) {
    if (!(ch1_aig_p & 1)) continue;

    Gate *candidate = gate(ch1_aig_p);
    for (auto &cand_p : candidate->get_parents()) {
      if (cand_p->get_gate_constraint()->len() != 2) continue;
      if (cand_p->is_child(ch2)) {
        Term *t = divide_by_var(cand_p->get_gate_constraint()->get_tail_term(), candidate->get_var()->get_dual());
        Term *t1 = divide_by_var(g->get_gate_constraint()->get_tail_term(), repl->get_var());
        if (t->degree() > 1) continue;
        if (t1->degree() > 1) continue;
        if (t == t1) {
          msg("dual twins push back (prerp) %s %s", g->get_var_name(), cand_p->get_var_name());
          g->print_gate_constraint(stdout);
          cand_p->print_gate_constraint(stdout);
          if (proof_logging) print_proof_van_constraint(g, cand_p);
          g->dual_twins_push_back(cand_p);
        }
      }
    }
  }
}
/*--------------------------------------------------------------------*/
static bool do_backward_substitution(Gate *outer) {

  Polynomial *outer_gc = outer->get_gate_constraint();
  if (outer_gc->len() != 2) return 0;
 
  Term *outer_t = outer_gc->get_tail_term();
  Term *res = outer_t;
  Gate *repl = 0;

  Term *outer_t_it = outer_t;
  while (outer_t_it) {
    Var *v = outer_t_it->get_var();
    for (auto &par : gate(v->get_num())->get_parents()) {
      if (par == outer) continue;
      if (par->get_output()) continue;
      Polynomial *p_par = par->get_gate_constraint();
      if (p_par->len() != 2) continue;
      if (!p_par->get_tail_term()->contains(v)) continue;
 

      Term *t = divide_by_term(outer_t, p_par->get_tail_term());
      if (t == outer_t) continue;

      if (t->degree() < res->degree()) {
        res = t;
        repl = par;
        if (res->degree() == 1) break;
      }
    }
    if (res->degree() == 1) break;
    outer_t_it = outer_t_it->get_rest();
  }
 
  if (!repl) return 0;

  Term *t0 = res;
  Term *t1 = new_term(repl->get_var());
  Term *t2;
  if (t0)
    t2 = multiply_term(t0, t1);
  else
    t2 = t1;

  push_mstack_end(outer_gc->get_mon(0)->copy());
  Monomial *tmp = new Monomial(one, t2);
  push_mstack_end(tmp);
  Polynomial *rewr = build_poly();

  if (proof_logging) {
    Monomial *tmp2 = new Monomial(minus_one, t0);
    print_pac_combi_monomial_rule(proof_file, repl->get_gate_constraint(), tmp2, outer_gc, 0, rewr);
    delete (tmp2);
  }

  outer->update_gate_poly(rewr);

  if (verbose > 3)
    msg("substituted %s in %s", repl->get_var_name(), outer->get_var_name());

  if (do_vanishing_constraints) check_for_new_vanishing_combinations(repl, outer);  // TODO maybe needed

  return 1;
}

/*------------------------------------------------------------------------*/
std::vector<Gate *> sub;

/*------------------------------------------------------------------------*/

static void backward_substitution() {  
  msg("backward substitution");

  int counter = 0;

  Gate *outer;
  Polynomial *outer_gc;
  Monomial *outer_m;
  Term *outer_t;

  for (unsigned i = M - 2; i >= NN; i--) {
    outer = gates[i];

    if (outer->get_elim())
      continue;
    if (outer->get_pp())
      continue;

    outer_gc = outer->get_gate_constraint();
    if (outer_gc->len() != 2)
      continue;

    outer_m = outer_gc->get_mon(1);
    outer_t = outer_m->get_term();
    if (!outer_t)
      continue;

    if (outer_t->degree() < 3)
      continue;

    do_backward_substitution(outer);

    //    if (is_unit(outer)) eliminate_unit_gate(outer);

    // outer->get_gate_constraint()->print(stdout);
    sub.push_back(outer);
    counter++;
  }

  if (verbose >= 1)
    msg("backwards substitution done", counter);
  // linearize_backward_sub();
}
/*============================================================================*/

static Var *search_for_tail(Term *t) {
  assert(t);
  Gate *g = gate(t->get_var_num());
  for (auto &parent : g->get_parents()) {
    Polynomial *gc = parent->get_gate_constraint();
    if (gc->len() != 2)
      continue;
    if (t == gc->get_tail_term())
      return parent->get_var();
  }
  return 0;
}

static Term *backward_rewrite_term(Term *t) {
  if (t->degree() == 1)
    return t->copy();
  Var *tail;
  std::vector<Var *> remainder;
  std::vector<Var *> remainder2;

  while (t) {
    if (t->get_ref() > 1) {
      tail = search_for_tail(t);
      if (tail) {
        remainder.push_back(tail);
        break;
      }
    }
    remainder.push_back(t->get_var());
    remainder2.push_back(t->get_var());
    t = t->get_rest();
  }
  Term *quo = sort_and_build_term_from_vector(remainder);
  return quo;
}

Term *backward_rewrite_term_until_completion(Term *t) {
  if (t->degree() == 1)
    return t;
  Term *tmp = backward_rewrite_term(t);
  while (tmp != t) {
    deallocate_term(t);
    t = tmp;
    tmp = backward_rewrite_term(t);
  }
  deallocate_term(t);
  return tmp;
}

/*-----------------------------------------------------------------------------*/

void preprocessing() {
  msg("starting preprocessing");

  remove_only_positives(1);
  remove_only_positives(0);
  
  
  if (!force_guessing) { // if we do not strictly force guessing, we switch on vanishing constraints
    Gate *g = gate(slit(NN - 1));
    if (g->get_xor_gate() == 1) {
      g = xor_left_child(g)->get_xor_gate() == 1 ? xor_right_child(g) : xor_left_child(g);
    }

    if (g->get_gate_constraint()->degree() > NN / 4) {
      msg("potential CLA %i, better solved with FGLM", g->get_gate_constraint()->degree());
      unmark_fsa();
      do_vanishing_constraints = true;
      if(!force_vanishing_off) find_vanishing_constraints();
    }
  }

  backward_substitution();
  


  msg("finished preprocessing");
}
