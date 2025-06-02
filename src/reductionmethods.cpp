/*------------------------------------------------------------------------*/
/*! \file reductionmethods.cpp
    \brief contains all reduction methods of polynomials

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#include "reductionmethods.h"

#include <ranges>

Polynomial* x_spec = 0;
std::map<Term*, Polynomial*> van_poly;
std::map<Term*, Polynomial*> dual_van_poly;
/*------------------------------------------------------------------------*/
Polynomial*
reduce_by_one_poly(Polynomial* p1, Polynomial* p2, bool non_lin_rewriting) {
  Polynomial* negfactor = divide_poly_by_term(p1, p2->get_lt());
  if (!negfactor || negfactor->is_constant_zero_poly())
    return p1->copy();

  if (mpz_cmp_si(p2->get_lm()->coeff, 0) == 1) {
    Polynomial* tmp = multiply_poly_with_constant(negfactor, minus_one);
    delete (negfactor);
    negfactor = tmp;
  }

  Polynomial* mult = multiply_poly(negfactor, p2);
  Polynomial* rem = add_poly(p1, mult);

  if (!proof_logging) {
    delete (mult);
    delete (negfactor);

  } else if (!non_lin_rewriting) {
    assert(proof_file);
   
    if (!negfactor->is_constant_one_poly())
      print_pac_combi_rule(proof_file, p2, negfactor, p1, 0, rem);
    else
      print_pac_add_rule(proof_file, p1, p2, rem);

    // print_pac_del_rule(proof_file,p1);
  } else {
    Polynomial * x_p1 = x_spec;
    Polynomial* neg_x_tmp = multiply_poly_with_term(negfactor, term_x);
    Polynomial* rem_x_tmp = multiply_poly_with_term(rem, term_x);
    push_mstack(new Monomial(minus_one, 0));
    Polynomial* min_one = build_poly();
    x_spec = add_poly(rem_x_tmp, min_one);
    
    print_pac_combi_rule(proof_file, p2, neg_x_tmp, x_p1, 0, x_spec);
  }

  return rem;
}
/*------------------------------------------------------------------------*/

Polynomial*
substitute_linear_poly(Polynomial* p1, Polynomial* p2) {
  assert(p1->degree() == 1 && p2->degree() == 1);
  Monomial* p1_m = 0;
  for (size_t i = 0; i < p1->len(); i++) {
    Monomial* p1_iter_mon = p1->get_mon(i);

    if (!p1_iter_mon->get_term())
      continue;

    int cmp = cmp_term(p1_iter_mon->get_term(), p2->get_lt());
    if (!cmp) {
      p1_m = p1_iter_mon;
      break;
    } else if (cmp == -1)
      break;
  }
  if (!p1_m)
    return p1->copy();

  Monomial* p2_m = p2->get_mon(0);

  if (!mpz_cmp(p1_m->coeff, p2_m->coeff)) {
    Polynomial* rem = sub_poly(p1, p2);

    if (proof_logging) {
      push_mstack(new Monomial(minus_one, term_x->copy()));
      Polynomial* min_one = build_poly();

      Polynomial* rem_x_tmp = multiply_poly_with_term(rem, term_x);
      push_mstack(new Monomial(minus_one, 0));
      Polynomial* min_one_rem = build_poly();
      Polynomial* x_spec_new = add_poly(rem_x_tmp, min_one_rem);

      print_pac_combi_rule(proof_file, p2, min_one, x_spec, 0, x_spec_new);
      delete (x_spec);
      x_spec = x_spec_new;
    }

    return rem;
  }

  Polynomial* res = 0;
  mpz_t coeff;
  mpz_init(coeff);
  mpz_tdiv_r(coeff, p1_m->coeff, p2_m->coeff);
  if (mpz_sgn(coeff) != 0) {
    die(1, "cannot use p2 to reduce p1");
  } else {
    mpz_tdiv_q(coeff, p1_m->coeff, p2_m->coeff);
    Polynomial* p2_lift = multiply_poly_with_constant(p2, coeff);
    res = sub_poly(p1, p2_lift);

    if (proof_logging) {
      mpz_t neg;
      mpz_init(neg);
      mpz_neg(neg, coeff);
      push_mstack(new Monomial(neg, term_x->copy()));
      Polynomial* min_one = build_poly();
      Polynomial* res_x_tmp = multiply_poly_with_term(res, term_x);
      push_mstack(new Monomial(minus_one, 0));
      Polynomial* min_one_rem = build_poly();
      Polynomial* x_spec_new = add_poly(res_x_tmp, min_one_rem);
      print_pac_combi_rule(proof_file, p2, min_one, x_spec, 0, x_spec_new);
      delete (x_spec);
      x_spec = x_spec_new;
    }

    delete (p2_lift);
  }
  mpz_clear(coeff);
  return res;
}
/*------------------------------------------------------------------------*/
int flipcount = 0;
Polynomial*
flip_var_in_poly(Polynomial* p1, Var* v, bool rem_van) {
  if (!proof_logging) {
    for (size_t i = 0; i < p1->len(); i++) {
      Monomial* m = p1->get_mon(i);
      if (!m->get_term()) {
        push_mstack(m->copy());
        continue;
      }
      if (m->get_term()->get_var_level() < v->get_level()) {
        for (size_t j = i; j < p1->len(); j++) {
          Monomial* next_m = p1->get_mon(j);
          push_mstack(next_m->copy());
        }
        break;
      }
      if (m->get_term()->contains(v)) {
        bool flag = 0;
        Term* t = divide_by_var(m->get_term(), v);

        Gate* g = gate(v->get_num());
        if (rem_van && t && v->is_dual() && g->get_van_twins().size() != 0) {
          for (auto& g_van : g->get_van_twins()) {
            if (t->contains(g_van->get_var())) {
              flag = 1;
              break;
            }
          }
        }
        if (!flag) {
          Term* vt = multiply_term_by_var(t, v->get_dual());
          mpz_t neg;
          mpz_init(neg);
          mpz_neg(neg, m->coeff);
          push_mstack(new Monomial(neg, vt));
          mpz_clear(neg);
        }
        push_mstack(new Monomial(m->coeff, t));

      } else {
        push_mstack(m->copy());
      }
    }

    Polynomial* rem = build_poly();

    return rem;
  } else {
    Polynomial* flip = 0;

    if (!v->is_dual()) {
      flip = gen_dual_constraint(v);
      if (proof_logging)
        flip->set_idx(gate(v->get_num())->get_dual_constraint()->get_idx());
    } else {
      flip = gate(v->get_num())->get_dual_constraint();
    }

    Polynomial* negfactor = divide_by_var(p1, flip->get_lt());
    Polynomial* rem;

    if (negfactor->is_constant_zero_poly()) {
      rem = p1->copy();

    } else {
      Polynomial* mult = multiply_poly(negfactor, flip);
      rem = add_poly(p1, mult);

      if (proof_logging) {
        if (!negfactor->is_constant_one_poly())
          print_pac_combi_rule(proof_file, flip, negfactor, p1, 0, rem);
        else
          print_pac_add_rule(proof_file, p1, flip, rem);
      }

      delete (mult);
    }
    if (!v->is_dual())
      delete (flip);

    delete (negfactor);
    // msg("flipped %i times", ++flipcount);
    return rem;
  }
}
/*------------------------------------------------------------------------*/
Polynomial*
mod_poly(Polynomial* p1, int exp) {
  if (!p1)
    return 0;

  mpz_t coeff;
  mpz_init(coeff);

  for (size_t i = 0; i < p1->len(); i++) {
    Monomial* m = p1->get_mon(i);
    mpz_tdiv_r_2exp(coeff, m->coeff, exp);
    if (mpz_sgn(coeff) != 0) {
      Monomial* tmp;
      if (m->get_term())
        tmp = new Monomial(coeff, m->get_term_copy());
      else
        tmp = new Monomial(coeff, 0);
      push_mstack_end(tmp);
    }
  }
  mpz_clear(coeff);
  Polynomial* out = build_poly();

  if (proof_logging) {
    mpz_t quot;
    mpz_init(quot);
    for (size_t i = 0; i < p1->len(); i++) {
      Monomial* m = p1->get_mon(i);

      mpz_tdiv_q_2exp(quot, m->coeff, exp);
      if (mpz_sgn(quot) != 0) {
        mpz_neg(quot, quot);
        Monomial* tmp;
        if (m->get_term())
          tmp = new Monomial(quot, m->get_term_copy());
        else
          tmp = new Monomial(quot, 0);
        push_mstack_end(tmp);
      }
    }
    mpz_clear(quot);

    Polynomial* p = build_poly();
    if (p) {
      push_mstack(new Monomial(one, term_x->copy()));
      Polynomial* mod = build_poly();
      Polynomial* px = multiply_poly(p, mod);
      Polynomial* mod_mul = multiply_poly_with_constant(px, mod_coeff);
      print_pac_mod_rule(proof_file, px, mod_mul);

      Polynomial* res_x_tmp = multiply_poly_with_term(out, term_x);
      push_mstack(new Monomial(minus_one, 0));
      Polynomial* min_one_rem = build_poly();
      Polynomial* x_spec_new = add_poly(res_x_tmp, min_one_rem);
      print_pac_add_rule(proof_file, x_spec, mod_mul, x_spec_new);
      delete (x_spec);
      x_spec = x_spec_new;

      delete (p);
      delete (mod);
    }
  }

  return out;
}
/*------------------------------------------------------------------------*/
Polynomial*
unflip_poly(Polynomial* p) {
  Polynomial* res = p->copy();
  Var* v = res->contains_dual_var();
  while (v) {
    Polynomial* tmp = flip_var_in_poly(res, v, 0);
    delete (res);
    res = tmp;
    v = res->contains_dual_var();
  }

  return res;
}

/*------------------------------------------------------------------------*/
Polynomial*
remove_vanishing_monomials(Polynomial* p,
                           std::vector<Polynomial*>* used_van_poly) {

  if (!proof_logging) {
    for (unsigned i = 0; i < p->len(); i++) {
      Monomial* m = p->get_mon(i);
      if (!m->get_term() || m->get_term()->degree() < 2) {
        push_mstack(m->copy());
        continue;
      }
      Term* t = m->get_term();

      bool flag = 0;
      while (t) {
        if (!t->get_var()->is_dual()) {
          Gate* g = gate(t->get_var_num());
          if (g->get_van_twins().size() != 0) {
            for (auto& g_van : g->get_van_twins()) {
              if (t->contains(g_van->get_var())) {
                van_mon_used_count++;
                flag = 1;
                break;
              }
            }
          }
          if (flag)
            break;
        }
        t = t->get_rest();
      }

      Term* shrunk = 0;
      if (!flag) {
        t = m->get_term();
        while (t) {
          if (!t->get_var()->is_dual()) {
            Gate* g = gate(t->get_var_num());
            for (auto& g_dv : g->get_dual_twins()) {
              if (m->get_term()->contains(g_dv->get_var())) {
                shrunk = divide_by_var(m->get_term(), g_dv->get_var());
                van_mon_used_count++;
                flag = 1;
              } else if (m->get_term()->contains(g_dv->get_var()->get_dual())){
                shrunk = 0;
                van_mon_used_count++;
                flag = 1;
              }
            }
          } 
          if (flag)
            break;
          t = t->get_rest();
        } 
      }

      if (!flag)
        push_mstack(m->copy());
      else if (shrunk)
        push_mstack(new Monomial(m->coeff, shrunk->copy()));
    }
    Polynomial* rem = build_poly();
    // rem->print(stdout);

    return rem;
  } else {  // with proof logging

    size_t i = 0, plen = p->len();
    Polynomial* rest = p->copy();
    for (; i < plen; i++) {
      Monomial* m = rest->get_mon(i);
      if (!m->get_term() || m->get_term()->degree() < 2) {
        continue;
      }
      Term* t = m->get_term();
      bool flag = 0;
      while (t) {
        if (!t->get_var()->is_dual()) {
          Gate* g = gate(t->get_var_num());
          if (g->get_van_twins().size() != 0) {
            for (auto& g_van : g->get_van_twins()) {
              if (t->contains(g_van->get_var())) {
                van_mon_used_count++;

                Term* t1 = new_quadratic_term(t->get_var(), g_van->get_var());

                Polynomial* p1 = van_poly[t1];
                bool find_flag = 0;
                if (used_van_poly) {
                  for (auto& vp : *used_van_poly) {
                    if (equal_poly(p1, vp)) {
                      find_flag = 1;
                      break;
                    }
                  }

                  if (!find_flag) {
                    if (proof_logging) {
                      fprintf(proof_file, "in0 %lu ", p1->get_idx());
                      p1->print(proof_file);
                    }
                    used_van_poly->push_back(p1);
                  }
                }
                assert(p1);
                Polynomial* p2 = reduce_by_one_poly(rest, p1);

                delete (rest);
                rest = p2;
                assert(rest);
                plen = rest->len();
                flag = 1;
                break;
              }
            }
          }
          if (flag)
            break;
        }
        t = t->get_rest();
      }

      if (!flag) {
        t = m->get_term();

        while (t) {
          if (!t->get_var()->is_dual()) {
            Gate* g = gate(t->get_var_num());
            for (auto& g_dv : g->get_dual_twins()) {
              if (m->get_term()->contains(g_dv->get_var())) {
                van_mon_used_count++;

                Term* t1 = new_quadratic_term(t->get_var(), g_dv->get_var());
                // t1->print(stdout);
                // t1->print_orig(stdout);
                Polynomial* p1 = dual_van_poly[t1];
                if (!p1) die(1, "did not find poly");

                if (used_van_poly) {
                  bool find_flag = 0;

                  for (auto& vp : *used_van_poly) {
                    if (equal_poly(p1, vp)) {
                      find_flag = 1;
                      break;
                    }
                  }
                  if (!find_flag) {
                    if (proof_logging) {
                      fprintf(proof_file, "in0 %lu ", p1->get_idx());
                      p1->print(proof_file);
                    }
                    used_van_poly->push_back(p1);
                  }
                }
                Polynomial* p2 = reduce_by_one_poly(rest, p1);

                delete (rest);
                rest = p2;
                assert(rest);
                plen = rest->len();
                flag = 1;
                break;
              }
            }
          }
          if (flag)
            break;
          t = t->get_rest();
        }
      }
    }
    
    return rest;
  }
}
/*------------------------------------------------------------------------*/

Polynomial*
unflip_poly_and_remove_van_mon(Polynomial* p) {
  Var* v = p->contains_dual_var();
  if (!v) {
    return remove_vanishing_monomials(p);
  }
  Polynomial * res = p->copy();
  
  while (v) {
    Polynomial* tmp = flip_var_in_poly(res, v, 1);
    if(!tmp) return 0;
   
    if (tmp->degree() > 1) {
      Polynomial* tmp2 = remove_vanishing_monomials(tmp);
      delete (tmp);
      tmp = tmp2;
    }
    delete (res);
    res = tmp;
    if(res) v = res->contains_dual_var();
    else return 0;
  }
  

  return res;
}
/*------------------------------------------------------------------------*/
struct {
  bool operator()(Gate* a, Gate* b) const {
    return a->get_var_level() > b->get_var_level();
  }
} cmpGateLvl;
/*------------------------------------------------------------------------*/
static Gate*
get_largest_node(Polynomial* p, std::set<Polynomial*> G) {
  std::list<Gate*> poly_var = get_var_of_poly(p, 0);

  std::list<Gate*> gate_var;
  for (Polynomial* g : G) {
    gate_var.push_back(gate(g->get_lt()->get_var_num()));
  }

  poly_var.sort(cmpGateLvl);

  for (auto it = poly_var.begin(); it != poly_var.end(); ++it) {
    if (std::find(gate_var.begin(), gate_var.end(), *it) != gate_var.end()) {
      Gate* res = *it;
      return res;
    }
  }

  return nullptr;
}

/*------------------------------------------------------------------------*/
static Polynomial * clean_phases(Polynomial * p1, Polynomial * p2){
  if (!p1) return 0;
  if (!p2) return p1->copy();

  if(p2->len() != 2 ) return p1->copy();

  Term * t = p2->get_tail_term();

  Polynomial * res = p1->copy();
  while (t && res){
    Polynomial * tmp1 = flip_var_in_poly(res, t->get_var(),0);
    if(!tmp1) return 0;
    Polynomial * tmp2 = flip_var_in_poly(res, t->get_var()->get_dual(),0);
    if (!tmp2) return 0;
    delete(res);
    if(tmp1->len() <= tmp2->len()){
      delete(tmp2);
      res = tmp1; 
    } else {
      delete(tmp1);
      res = tmp2; 
    }

    t = t->get_rest();
  }

  return res; 

}
/*------------------------------------------------------------------------*/
std::list<Gate*> dyn_red_guesses(Polynomial * rem, std::set<Polynomial*> G){
  std::list<Gate*> res;
  if (rem->degree() > 1 && rem->len() > 1 && rem->get_tail_poly()->degree()== 1){
    Term * t = rem->get_lt();
    while(t){
      Gate * g = gate(t->get_var_num());
  
      if(!g->get_input() && g->get_nf() && g->get_nf()->degree() > 1 ) res.push_back(g);

      t = t->get_rest();
    }
  }
  
  if(res.size() == 0) res.push_back(get_largest_node(rem, G));

  return res;
}

/*------------------------------------------------------------------------*/
bool reduce_to_zero(Polynomial* p, std::set<Polynomial*> G) {
  Polynomial* rem = p->degree() > 1 ? remove_vanishing_monomials(p) : p->copy();
  // if(rem)rem->print(stdout);

  std::list<Gate*> next_reduction;
  while (rem && !rem->is_constant_zero_poly()) {
    if(verbose > 2  && rem && rem-> len() <100) rem->print(stdout);
   
    if(next_reduction.size() == 0){
      next_reduction = dyn_red_guesses(rem, G);
    }
  
    Gate* v = next_reduction.front();
    next_reduction.pop_front();
 
    if (!v) {
      
      Polynomial* final = unflip_poly_and_remove_van_mon(rem);
    
      if (!final)
        return true;
      else if (final->is_constant_zero_poly()) {
        delete (final);
        return true;
      }
      delete (final);

      return false;
    }
   
    if (rem->len() > 8000) {
      std::cout << "Possibly wrong\n";
      return false;
    }

  // all of this should probably go outside the function?
   Polynomial* red = v->get_nf() ? v->get_nf() : v->get_gate_constraint();

  

    Polynomial* rem_unf = flip_var_in_poly(rem, red->get_lt()->get_var()->get_dual(), 1);
    delete (rem);
    if (!rem_unf || rem_unf->is_constant_zero_poly())
      return true;
    Polynomial * rem_unf1 = remove_vanishing_monomials(rem_unf);
    delete (rem_unf);
    if (!rem_unf1 || rem_unf1->is_constant_zero_poly())
      return true;

    

    Polynomial* res = reduce_by_one_poly(rem_unf1, red);
    if(!res) return true;
    Polynomial* res1 = remove_vanishing_monomials(res);
    if(!res1) return true;
    Polynomial * res_cleaned = clean_phases(res1, red);
    delete(res);
    if (res_cleaned) {
      msg_nl("d %2.i  lvl %5.i,  %6.i, reduced by ", gate(red->get_lt()->get_var_num())->get_dist(), red->get_lt()->get_var()->get_level(), res_cleaned->len());
      // red->print(stdout);
      if(verbose > 3 && res_cleaned-> len() <100) res_cleaned->print(stdout);
    }
    delete (rem_unf1);
    rem = res_cleaned;
  }

  return true;
}
