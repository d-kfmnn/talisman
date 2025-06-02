/*------------------------------------------------------------------------*/
/*! \file polynomial_solver.cpp
    \brief contains the polynomial solving routine

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#include "polynomial_solver.h"
/*------------------------------------------------------------------------*/
// Global variable
bool gen_witness = 1;
/*------------------------------------------------------------------------*/
// ERROR CODES:
static int err_writing = 41;   // cannot write to
static int err_rem_poly = 42;  // remainder poly no witness
/*------------------------------------------------------------------------*/
static Polynomial *linearize_spec(Polynomial *spec) {
  msg("");
  print_hline();
  msg("started reducing non linear terms in spec");
  if (spec->degree() == 1)
    return spec;

  std::vector<int> factor_idx;
  std::vector<const Polynomial *> factor_polys;

  Polynomial *rem = spec;
  std::vector<size_t> indices;
  bool enlarged = 0;
  for (size_t i = 0; i < rem->len(); i++) {
    Monomial *m = rem->get_mon(i);
    Term *t = m->get_term();
    bool flag = 0;
    if (t->degree() == 1) {
      push_mstack(m->copy());
      flag = 1;
    } else if (t->get_ref() > 1) {
      Gate *g = gate(t->get_var_num());

      for (auto &sub : g->get_parents()) {
        if (sub->get_elim())
          continue;
        Polynomial *sub_gc = sub->get_gate_constraint();
        if (sub_gc->len() != 2)
          continue;

        Monomial *inner_m = sub_gc->get_mon(1);
        Term *inner_t = inner_m->get_term();
        if (t != inner_t)
          continue;
        Term *lt = sub_gc->get_lt()->copy();
        Monomial *tmp = new Monomial(m->coeff, lt);

        if (proof_logging) {
          mpz_t neg;
          mpz_init(neg);
          mpz_neg(neg, m->coeff);
          Monomial *sub_mon = new Monomial(neg, term_x->copy());
          Monomial **intern_mstack = new Monomial *[1];
          intern_mstack[0] = sub_mon;
          Polynomial *sub_poly = new Polynomial(intern_mstack, 1, 2);
          factor_idx.push_back(sub_gc->get_idx());
          factor_polys.push_back(sub_poly);
        }

        push_mstack(tmp);
        flag = 1;
        break;
      }
    }
    if (!flag) {
      if (!enlarged)  // TODO this should not be here. this should be automatically handled when adding new gates
      {
        enlarge_gates(rem->len());
        enlarged = 1;
      }
      Term *rep_t = extend_var_gates(t);
      Monomial *tmp = new Monomial(m->coeff, rep_t->copy());

      if (proof_logging) {
        mpz_t neg;
        mpz_init(neg);
        mpz_neg(neg, m->coeff);
        Monomial *sub_mon = new Monomial(neg, term_x->copy());

        Monomial **intern_mstack = new Monomial *[1];
        intern_mstack[0] = sub_mon;
        Polynomial *sub_poly = new Polynomial(intern_mstack, 1, 2);
        Polynomial *ext_repl = gate(rep_t->get_var_num())->get_gate_constraint();
        factor_idx.push_back(ext_repl->get_idx());
        factor_polys.push_back(sub_poly);
      }

      push_mstack(tmp);
    }
  }

  Polynomial *tmp = build_poly();
  if (verbose > 2) {
    msg_nl("linearized spec ");
    tmp->print(stdout);
  }

  adjust_level_of_extended_gates();

  if (proof_logging) {
    factor_idx.push_back(x_spec->get_idx());
    factor_polys.push_back(new Polynomial());
    Polynomial *rem_x_tmp = multiply_poly_with_term(tmp, term_x);
    push_mstack(new Monomial(minus_one, 0));
    Polynomial *min_one = build_poly();
    delete (x_spec);
    x_spec = add_poly(rem_x_tmp, min_one);

    print_pac_vector_combi_rule(proof_file, factor_idx, factor_polys, x_spec);
  }

  return tmp;
}
/*------------------------------------------------------------------------*/

bool verify(const char *inp_f, Polynomial *spec, const char *out_f1, const char *out_f2, const char *out_f3) {
  assert(!proof_logging || inp_f);
  assert(!proof_logging || out_f1);
  assert(!proof_logging || out_f2);
  assert(!proof_logging || out_f3);

  FILE *f1 = 0, *f2 = 0, *f3 = 0;
  if (proof_logging) {
    if (!(f1 = fopen(out_f1, "w")))
      die(err_writing, "can not write output to '%s'", out_f1);

    if (!(f2 = fopen(out_f2, "w")))
      die(err_writing, "can not write output to '%s'", out_f2);

    if (!(f3 = fopen(out_f3, "w")))
      die(err_writing, "can not write output to '%s'", out_f3);
  }

  if (proof_logging) {
    polys_file = f1;
    print_circuit_poly(polys_file);
    print_dual_constraints(polys_file);
    proof_file = f2;
    init_proof_logging(-num_gates);
  }

  
  identify_final_stage_adder();
  
  if(!force_vanishing_off) find_vanishing_constraints_light();

  // Preprocessing of Circuit
  if (do_preprocessing && !force_guessing)  // guessing requires AIG nodes, otherwise they cannot be encoded to clauses
    preprocessing();

  // If needed spec will be linearized
  Polynomial *rem = spec->copy();
  

  if (proof_logging) {
    Polynomial *rem_x_tmp = multiply_poly_with_term(rem, term_x);
    push_mstack(new Monomial(minus_one, 0));
    Polynomial *min_one = build_poly();
    x_spec = add_poly(rem_x_tmp, min_one);

    pac_add_circuit_poly(polys_file, x_spec);
    print_refutation_spec(f3);
  }

  if (rem->degree() > 1) {
    Polynomial *tmp = linearize_spec(rem);
    delete (rem);
    rem = tmp;
  }
  assert(rem->degree() == 1);

 
  
  rem = reduce(rem);

  // Taking care of results
  bool res;
  print_hline();
  if (rem && !rem->is_constant_zero_poly()) {
    if (!check_inputs_only(rem)) {
      msg("REMAINDER IS");
      msg_nl(" ");
      rem->print(stdout);
      msg("");
      die(err_rem_poly, "internal sorting error - remainder polynomial contains non-inputs");
    }

    res = 0;
    msg("RESULT: INCORRECT MULTIPLIER");
    msg("");

    if (inp_f && gen_witness) {
      msg("REMAINDER IS");
      msg_nl(" ");
      rem->print(stdout);
      msg("");
      msg("GENERATING WITNESSES IS UNDER CONSTRUCTION");
      // generate_witness(rem, inp_f); TODO
    }
  } else {
    res = 1;
    msg("");
    msg("RESULT: CORRECT MULTIPLIER");
    if (proof_logging) {
      Polynomial *neg = multiply_poly_with_constant(x_spec, minus_one);
      print_pac_mul_const_rule(proof_file, x_spec, -1, neg);

      msg("");
      msg("writing gate constraints to '%s' ", out_f1);
      msg("writing proof certificate to '%s'", out_f2);
      msg("writing specification to '%s'    ", out_f3);
    }
  }

  if (rem) delete (rem);

  if (proof_logging) {
    fclose(f1);
    fclose(f2);
    fclose(f3);
  }
  return res;
}
