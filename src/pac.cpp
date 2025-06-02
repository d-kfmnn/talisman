/*------------------------------------------------------------------------*/
/*! \file pac.cpp
    \brief contains functions necessary to generate PAC proofs

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#include "pac.h"
/*------------------------------------------------------------------------*/
static int poly_idx;
/*------------------------------------------------------------------------*/
// ERROR CODES:
static int err_rule = 81;  // error in proof rule
/*------------------------------------------------------------------------*/

Var *var_x;
Term *term_x;

void init_proof_logging(size_t level) {
  var_x = new Var("x", level);
  term_x = new_term(var_x);
}

void print_refutation_spec(FILE *file) {
  fputs("1;", file);
}

void print_circuit_poly(FILE *file) {
  fputs("1 ", file);
  mpz_out_str(file, 10, mod_coeff);
  fputs(";\n", file);
  poly_idx = 2;
  for (unsigned i = NN; i < num_gates; i++) {
    Polynomial *p = gen_gate_constraint(i);
    assert(p);

    fprintf(file, "%i ", poly_idx);
    p->print(file);
    p->set_idx(poly_idx);
    Gate * g = gates[i];
    g->get_aig_poly()->set_idx(poly_idx);
    poly_idx++;
  }
}

void pac_add_circuit_poly(FILE *file, Polynomial *p) {
  fprintf(file, "%i ", poly_idx);
  p->print(file);
  p->set_idx(poly_idx++);
}

void print_dual_constraints(FILE *file) {
  for (unsigned i = 0; i < M - 1; i++) {
    Gate *g = gates[i];

    Polynomial *p = g->get_dual_constraint();
    assert(p);

    // fprintf(file, "%lu = %s, -%s+1;\n", poly_idx, g->get_var()->get_dual()->get_name(),g->get_var_name());

    fprintf(file, "%i ", poly_idx);
    p->print(file);
    p->set_idx(poly_idx++);
  }
}

void print_pac_extension_rule_for_mon(FILE *file, Gate *g, const Term *t, Polynomial *p) {
  fprintf(file, "%i = %s, ", poly_idx, g->get_var_name());
  t->print(file);
  fprintf(file, ";\n");
  p->set_idx(poly_idx++);
}

/*------------------------------------------------------------------------*/

void print_pac_del_rule(FILE *file, const Polynomial *p1) {
  assert(p1);

  fprintf(file, "%lu d;\n", p1->get_idx());
}

int print_pac_pattern_out_rules(FILE *file, std::vector<Polynomial *> lin_poly, int i) {
  for (auto &p : lin_poly) {
    p->set_idx(poly_idx++);
    fprintf(file, "out%i %lu ", i++, p->get_idx());
    p->print(file);
  }
  return i;
}

/*------------------------------------------------------------------------*/

void print_pac_mod_rule(FILE *file, const Polynomial *p1, Polynomial *p) {
  assert(p1 && !p1->is_constant_zero_poly());
  assert(p && !p->is_constant_zero_poly());

  fprintf(file, "%i ", poly_idx);
  fputs("% 1 *(", file);
  p1->print(file, 0);
  fputs("), ", file);
  p->print(file);

  p->set_idx(poly_idx++);
}

/*------------------------------------------------------------------------*/

void print_pac_add_rule(
    FILE *file, const Polynomial *p1, const Polynomial *p2, Polynomial *p) {
  assert(p1 && !p1->is_constant_zero_poly());
  assert(p2 && !p2->is_constant_zero_poly());
  assert(p && !p->is_constant_zero_poly());

  fprintf(file, "%u %% %lu + %lu, ", poly_idx, p1->get_idx(), p2->get_idx());
  p->print(file);

  p->set_idx(poly_idx++);
}

/*------------------------------------------------------------------------*/

void print_pac_vector_add_rule(
    FILE *file, std::vector<int> indices, Polynomial *p) {
  fprintf(file, "%u %% ", poly_idx);

  int ind;
  while (!indices.empty()) {
    ind = indices.back();
    indices.pop_back();

    fprintf(file, "%i", ind);
    if (indices.size() > 0) fputs(" + ", file);
  }
  fprintf(file, ", ");
  p->print(file);
  p->set_idx(poly_idx++);
}
/*------------------------------------------------------------------------*/

void print_pac_combi_rule(
    FILE *file, const Polynomial *p1, const Polynomial *p2,
    const Polynomial *p3, const Polynomial *p4, Polynomial *p) {
  assert(p1 && !p1->is_constant_zero_poly());
  assert(p2 && !p2->is_constant_zero_poly());
  assert(p3 && !p3->is_constant_zero_poly());

  fprintf(file, "%u %% %lu", poly_idx, p1->get_idx());
  if (p2) {
    fputs(" *(", file);
    p2->print(file, 0);
    fputs(") ", file);
  }

  fprintf(file, "+ %lu", p3->get_idx());
  if (p4) {
    fputs(" *(", file);
    p4->print(file, 0);
    fputs(") ", file);
  }

  fprintf(file, ", ");
  if (p) {
    p->print(file);
    p->set_idx(poly_idx++);
  } else
    fprintf(file, "0;\n");
}

void print_pac_combi_monomial_rule(
    FILE *file, const Polynomial *p1, const Monomial *m2,
    const Polynomial *p3, const Polynomial *p4, Polynomial *p) {
  assert(p1 && !p1->is_constant_zero_poly());
  assert(m2);
  assert(p3 && !p3->is_constant_zero_poly());
  assert(p && !p->is_constant_zero_poly());

  fprintf(file, "%u %% %lu", poly_idx, p1->get_idx());
  if (m2) {
    fputs(" *(", file);
    m2->print(file, 0);
    fputs(") ", file);
  }

  fprintf(file, "+ %lu", p3->get_idx());
  if (p4) {
    fputs(" *(", file);
    p4->print(file, 0);
    fputs(") ", file);
  }

  fprintf(file, ", ");
  p->print(file);
  p->set_idx(poly_idx++);
}
/*------------------------------------------------------------------------*/

void print_pac_vector_combi_rule(
    FILE *file, std::vector<int> indices,
    std::vector<const Polynomial *> co_factors, Polynomial *p) {
  if (co_factors.size() != indices.size()) die(err_rule, "combination rule receives invalid arguments;");

  fprintf(file, "%u %% ", poly_idx);

  const Polynomial *tmp;
  int ind;
  while (!co_factors.empty()) {
    tmp = co_factors.back();
    co_factors.pop_back();
    ind = indices.back();
    indices.pop_back();

    fprintf(file, "%i", ind);
    if (tmp && !tmp->is_constant_one_poly() && !tmp->is_constant_zero_poly()) {
      fputs(" *(", file);
      tmp->print(file, 0);
      fputs(")", file);
    }

    if (co_factors.size() > 0) fputs(" + ", file);
    // delete(tmp);
  }
  fprintf(file, ", ");
  p->print(file);
  p->set_idx(poly_idx++);
}
/*------------------------------------------------------------------------*/

void print_pac_mul_rule(
    FILE *file, const Polynomial *p1, const Polynomial *p2, Polynomial *p) {
  assert(p1 && !p1->is_constant_zero_poly());
  assert(p2 && !p2->is_constant_zero_poly());
  assert(p && !p->is_constant_zero_poly());

  fprintf(file, "%u %% %lu *(", poly_idx, p1->get_idx());
  p2->print(file, 0);
  fputs("), ", file);
  p->print(file);

  p->set_idx(poly_idx++);
}

void print_pac_mul_const_rule(
    FILE *file, const Polynomial *p1, int n, Polynomial *p) {
  fprintf(file, "%u %% %lu *(%i), ", poly_idx, p1->get_idx(), n);

  p->print(file);

  p->set_idx(poly_idx++);
}
/*------------------------------------------------------------------------*/
