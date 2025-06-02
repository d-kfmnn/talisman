/*------------------------------------------------------------------------*/
/*! \file vanishing_constraints.cpp
    \brief contains functions used in the polynomial solver

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#include<list>
#include "vanishing_constraints.h"
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
static std::list<Gate *> collect_vanishing_pairs(Gate *g1, Gate *child, std::list<Gate *> prop) {
  for (unsigned g1_p_val : g1->get_aig_parents()) {
    if (g1_p_val & 1) continue;
    Gate *g1_p = gate(g1_p_val);
    //msg("propagated vanishing vanishing pair in collect %s and %s", g1_p->get_var_name(), child->get_var_name());
    g1_p->van_twins_push_back(child);
    child->van_twins_push_back(g1_p);
    van_mon_prop_count++;
    //msg("collected %s", g1_p->get_var_name());
    prop.push_back(g1_p);
    prop = collect_vanishing_pairs(g1_p, child, prop);
  }
  return prop;
}

/*------------------------------------------------------------------------*/
static void propagate_vanishing_pairs(Gate *g1, Gate *child, std::list<Gate *> prop) {
  for (unsigned g1_p_val : g1->get_aig_parents()) {
    if (g1_p_val & 1) continue;
    Gate *g1_p = gate(g1_p_val);

    //msg("propagated vanishing vanishing pair %s and %s", g1_p->get_var_name(), child->get_var_name());
    g1_p->van_twins_push_back(child);
    child->van_twins_push_back(g1_p);
    van_mon_prop_count++;
    for (auto &prop_elem : prop) {
      //msg("propagated vanishing vanishing prop elem pair %s and %s", g1_p->get_var_name(), prop_elem->get_var_name());
      prop_elem->van_twins_push_back(g1_p);
      g1_p->van_twins_push_back(prop_elem);
      van_mon_prop_count++;
    }

    propagate_vanishing_pairs(g1_p, child, prop);
  }
}
/*------------------------------------------------------------------------*/
static void gen_dual_van_constraint(Gate * triangle, Gate * n, Gate *nc1, Gate*nc2){
  Term * t1 = new_quadratic_term(triangle->get_var(), nc1->get_var());

  std::vector<int> indices;
    std::vector<const Polynomial*> co_factors;

  Polynomial *p1 =van_poly[t1]; 
  

  Term * t2 = new_quadratic_term(triangle->get_var(), nc2->get_var());
  Polynomial *p2 =van_poly[t2]; 
 

  Polynomial * p3 = n->get_aig_poly();
  Polynomial * p3_unf = unflip_poly(p3);
 
  indices.push_back(p3_unf->get_idx());

  push_mstack(new Monomial(one, triangle->get_aig_poly()->get_lt()));
  Polynomial *fac0 = build_poly();
  co_factors.push_back(fac0);


  Monomial * m = new Monomial(minus_one, nc1->get_aig_poly()->get_lt()->copy());
  push_mstack(m);
  Monomial * m1 = new Monomial(one, 0);
  push_mstack(m1);
  Polynomial * fac1 = build_poly();
 
  indices.push_back(p2->get_idx());
  co_factors.push_back(fac1);


  push_mstack(m1->copy());
  Polynomial * fac2 = build_poly();

  indices.push_back(p1->get_idx());
  co_factors.push_back(fac2);

  Monomial** intern_mstack = new Monomial*[2];
  Term* tt1 = new_quadratic_term(triangle->get_var(), n->get_var());
  Monomial* mm1 = new Monomial(minus_one, tt1);
  Term* tt2 = new_term(triangle->get_var());
  Monomial* mm2 = new Monomial(one, tt2);
  intern_mstack[0] = mm1;
  intern_mstack[1] = mm2;
  Polynomial* pp1 = new Polynomial(intern_mstack, 2, 2);
  
  print_pac_vector_combi_rule(proof_file, indices, co_factors, pp1);
  dual_van_poly.insert({tt1, pp1});


}

/*------------------------------------------------------------------------*/
static void find_vanishing_triangles() {
  for (unsigned i = 0; i < M; i++) {
    Gate *n = gates[i];
    if (n->get_input() > 0) continue;
    if (n->get_children().size() != 2) continue;
    
    Gate *ch1 = n->get_children().front();
    Gate *ch2 = n->get_children().back();

    if (!ch2->is_van_twin(ch1)) continue;
    std::vector<Gate *> triangle_candidates = {};
   // msg("n is %s", n->get_var_name());

    for (auto &ch1_van : ch1->get_van_twins()) {  // TODO: generalize to more than one case
      if (ch2->is_van_twin(ch1_van)) {
        triangle_candidates.push_back(ch1_van);
      }
    }
    if (triangle_candidates.size() > 0) {
      for (auto &triangle : triangle_candidates) {
        if(triangle == n) continue;
       if(verbose > 3) msg("found %s for %s", triangle->get_var_name(), n->get_var_name());

        triangle->dual_twins_push_back(n);
        //msg("dual twins push back (triangle) %s %s", triangle->get_var_name(), n->get_var_name());
        if(proof_logging) gen_dual_van_constraint(triangle, n, ch1, ch2);

        for (unsigned n_parents : n->get_aig_parents()) {
          if (!(n_parents & 1)) continue;
          Gate *np = gate(n_parents);
          if(verbose > 3) msg("found vanishing pair through triangle %s %s", np->get_var_name(), triangle->get_var_name());
          np->van_twins_push_back(triangle);
          triangle->van_twins_push_back(np);
          van_mon_prop_count++;
          std::list<Gate *> pos_parents_gp = {};
          pos_parents_gp = collect_vanishing_pairs(np, triangle, pos_parents_gp);
          propagate_vanishing_pairs(triangle, np, pos_parents_gp);
          
        }
        //propagate_dual_vanishing_pairs(triangle,n);
        if(verbose > 3) msg("found dual vanishing pair through triangle %s (1-%s)", triangle->get_var_name(), n->get_var_name());
        
      }
    }
  }
}

/*------------------------------------------------------------------------*/
static void propagate_xor_and(Gate * gp_gate, Gate *g, Gate * andg, Polynomial * p){
  
  for(auto &gpp_idx : gp_gate->get_aig_parents()){
    if(gpp_idx & 1) continue;
    Gate * gpp =gate(gpp_idx);

    aiger_and *and1 = is_model_and(gpp->get_var_num());
    if (!and1) continue;
    unsigned l = and1->rhs0, r = and1->rhs1;
    if (static_cast<int>(aiger_strip(l))== gp_gate->get_var_num() && aiger_sign(l)) continue;
    if (static_cast<int>(aiger_strip(r))== gp_gate->get_var_num() && aiger_sign(r)) continue;
   
    push_mstack(andg->get_aig_poly()->get_lm()->copy());
    Polynomial* f1 = build_poly();
    Polynomial * tmp0 = unflip_poly(gpp->get_aig_poly());
    Polynomial * left = multiply_poly(tmp0, f1);
    Polynomial * right = multiply_poly(p, tmp0->get_tail_poly());
    Polynomial * res = add_poly(left, right);

    std::vector<int> indices;
    std::vector<const Polynomial*> co_factors;
    indices.push_back(gpp->get_aig_poly()->get_idx());
    co_factors.push_back(f1);
    indices.push_back(p->get_idx());
    co_factors.push_back(gpp->get_aig_poly()->get_tail_poly());
    print_pac_vector_combi_rule(proof_file, indices,co_factors, res);

    if(g->is_in_pos_parents(gpp->get_var_num())) {
      van_poly.insert({res->get_lt(), res});
      //msg_nl("inserted van poly ");
      propagate_xor_and(andg, andg, gpp, res);
    }
  

    delete(tmp0);
    delete(f1);
    delete(left);
    delete(right);
    propagate_xor_and(gpp, g, andg, res);
    
  }
  
}

/*------------------------------------------------------------------------*/

static void gen_xor_and_van_constraint_and_propagate(Gate * gp_gate, Gate * xor_gate, Gate * andg){
 
  //Gate * xor_gate = gp_gate->children_front()->get_xor_gate() == 2 ? gp_gate->children_front() : gp_gate->children_back();

  Gate * l = xor_gate->children_front();
  Gate * r = xor_gate->children_back();

  Polynomial * l_unfl = unflip_poly(l->get_aig_poly());
  Polynomial * r_unfl = unflip_poly(r->get_aig_poly());
  Polynomial * xor_unfl = unflip_poly(xor_gate->get_aig_poly());

  Polynomial * xor_tmp = reduce_by_one_poly(xor_unfl, l_unfl);
  Polynomial * xor_red = reduce_by_one_poly(xor_tmp, r_unfl);

  Polynomial * tmp0 = unflip_poly(gp_gate->get_aig_poly());
  Polynomial * res = reduce_by_one_poly(tmp0, xor_red);

  Polynomial * and_tmp = unflip_poly(andg->get_aig_poly());

  
  std::vector<int> indices;
  std::vector<const Polynomial*> co_factors;

  // Step 1: 
  Polynomial * left0 = multiply_poly_with_term(and_tmp, xor_red->get_lt());
  push_mstack(new Monomial(one, xor_red->get_lt()));
  Polynomial * resp0 = build_poly();
  indices.push_back(and_tmp->get_idx());
  co_factors.push_back(resp0);

  
  Polynomial * right0 = multiply_poly(xor_red, and_tmp->get_tail_poly());
  indices.push_back(xor_red->get_idx());
  co_factors.push_back(and_tmp->get_tail_poly());

  indices.push_back(and_tmp->get_idx());
  push_mstack(new Monomial(minus_one, 0));
  Polynomial * fac = build_poly();
  co_factors.push_back(fac);

  Polynomial * result0 = add_poly(left0,right0);
  Polynomial * result1 = sub_poly(result0, and_tmp);
  print_pac_vector_combi_rule(proof_file, indices,co_factors, result1);
  dual_van_poly.insert({result1->get_lt(), result1});


  // Step 2:
  indices.clear();
  co_factors.clear();

  Polynomial * left = multiply_poly_with_monomial(and_tmp, res->get_lm());
  push_mstack(res->get_lm());
  Polynomial * resp = build_poly();
  Polynomial* tmp = multiply_poly_with_constant(and_tmp->get_tail_poly(), minus_one);
  Polynomial * right = multiply_poly(res, tmp);

  indices.push_back(and_tmp->get_idx());
  co_factors.push_back(resp);
  indices.push_back(res->get_idx());
  co_factors.push_back(tmp);

  Polynomial * result = add_poly(left,right);
  print_pac_vector_combi_rule(proof_file, indices,co_factors, result);
  van_poly.insert({result->get_lt(), result});
 

  propagate_xor_and(gp_gate, gp_gate, andg, result);
  
}
/*------------------------------------------------------------------------*/

static void gen_xor_child_van_constraints(Gate * l, Gate * r){
 
  //Gate * xor_gate = gp_gate->children_front()->get_xor_gate() == 2 ? gp_gate->children_front() : gp_gate->children_back();

  Polynomial * l_unfl = unflip_poly(l->get_aig_poly());
  Polynomial * r_unfl = unflip_poly(r->get_aig_poly());
  Polynomial * result = multiply_poly(l_unfl, r_unfl);
  print_pac_mul_rule(proof_file, l_unfl, r_unfl, result);

  
  van_poly.insert({result->get_lt(), result});
 

  
}
/*------------------------------------------------------------------------*/
static void identify_vanishing_pairs(Gate *g) {
 
  Gate * lg = g->children_front();
  Gate * rg = g->children_back();


  if(verbose > 3) msg("found vanishing xor child pair %s and %s", lg->get_var_name(), rg->get_var_name());
  lg->van_twins_push_back(rg);
  rg->van_twins_push_back(lg);
  if(proof_logging) gen_xor_child_van_constraints(lg, rg);

  if (g->get_aig_output()) return;

  Gate * llg = g->children_front()->children_front();
  Gate * lrg = g->children_front()->children_back();

  std::vector<Gate *> ands;
  for (auto &llg_p : llg->get_parents()) {
    if (g->is_child(llg_p)) continue;     // check whether internal xor node
    if (!llg_p->is_child(lrg)) continue;  // check whether both llg and lrg are parents
    if(llg_p->children_size() > 2) continue;
    ands.push_back(llg_p);
  }
  if (ands.size() == 0) return;
   
  if(do_vanishing_constraints){
  for (unsigned gp_negp : g->get_neg_parents()) {
    Gate *gp_gate = gate(gp_negp);
    for (auto &andg : ands) {
      if(verbose > 3) msg("found vanishing pair %s and %s", gp_gate->get_var_name(), andg->get_var_name());
      gp_gate->van_twins_push_back(andg);
      andg->van_twins_push_back(gp_gate);
      van_mon_prop_count++;
      
      if(proof_logging){
        gen_xor_and_van_constraint_and_propagate(gp_gate, g, andg);
      }

      for (unsigned gp_posp : gp_gate->get_pos_parents()) {
        Gate *gp_pos_gate = gate(gp_posp);
        if (gp_pos_gate->is_van_twin(andg)) continue;
        if(verbose > 3) msg("found vanishing pair2 %s and %s", gp_pos_gate->get_var_name(), andg->get_var_name());
        gp_pos_gate->van_twins_push_back(andg);
        andg->van_twins_push_back(gp_pos_gate);
        van_mon_prop_count++;

        
        for (unsigned andg_posp : andg->get_pos_parents()) {
          Gate *and_posp_gate = gate(andg_posp);
          if(verbose > 3) msg("found vanishing pair3 %s and %s", gp_pos_gate->get_var_name(), and_posp_gate->get_var_name());
          gp_pos_gate->van_twins_push_back(and_posp_gate);
          and_posp_gate->van_twins_push_back(gp_pos_gate);
          van_mon_prop_count++;
        }
      }
      
    }
  }
}

  if (ands.size() == 1) {
    Gate *and1 = *ands.begin();
    and1->set_xor_and(g);
    g->set_xor_and(and1);
    and1->dual_twins_push_back(g);
    if(verbose > 3) msg("dual twins push back (identify) %s %s", and1->get_var_name(), g->get_var_name());
  }
}
/*------------------------------------------------------------------------*/
static void find_and_propagate_xor_and(){
  for (unsigned i = 0; i < M; i++) {
    Gate *n = gates[i];
    //if(!n->get_fsa()) continue;
   
    if (n->get_xor_gate() != 1) continue;
    if(n->children_size() != 2) continue;
   
   
    identify_vanishing_pairs(n);
     
    
  }
}
/*------------------------------------------------------------------------*/
static void gen_xor_and_van_constraint(Gate * xor_gate, Gate * andg){
 
  //Gate * xor_gate = gp_gate->children_front()->get_xor_gate() == 2 ? gp_gate->children_front() : gp_gate->children_back();

  Gate * l = xor_gate->children_front();
  Gate * r = xor_gate->children_back();

  Polynomial * l_unfl = unflip_poly(l->get_aig_poly());
  Polynomial * r_unfl = unflip_poly(r->get_aig_poly());
  Polynomial * xor_unfl = unflip_poly(xor_gate->get_aig_poly());

  Polynomial * xor_tmp = reduce_by_one_poly(xor_unfl, l_unfl);
  Polynomial * xor_red = reduce_by_one_poly(xor_tmp, r_unfl);

  Polynomial * and_tmp = unflip_poly(andg->get_aig_poly());

  
  std::vector<int> indices;
  std::vector<const Polynomial*> co_factors;

  // Step 1: 
  Polynomial * left0 = multiply_poly_with_term(and_tmp, xor_red->get_lt());
  push_mstack(new Monomial(one, xor_red->get_lt()));
  Polynomial * resp0 = build_poly();
  indices.push_back(and_tmp->get_idx());
  co_factors.push_back(resp0);

  
  Polynomial * right0 = multiply_poly(xor_red, and_tmp->get_tail_poly());
  indices.push_back(xor_red->get_idx());
  co_factors.push_back(and_tmp->get_tail_poly());

  indices.push_back(and_tmp->get_idx());
  push_mstack(new Monomial(minus_one, 0));
  Polynomial * fac = build_poly();
  co_factors.push_back(fac);

  Polynomial * result0 = add_poly(left0,right0);
  Polynomial * result1 = sub_poly(result0, and_tmp);
  print_pac_vector_combi_rule(proof_file, indices,co_factors, result1);
  dual_van_poly.insert({result1->get_lt(), result1});
}
/*------------------------------------------------------------------------*/


static void find_xor_and(){
  for (unsigned i = 0; i < M; i++) {
    Gate *g = gates[i];

    if (g->get_xor_gate() != 1) continue;
    
    Gate * and1 = g->get_xor_and_gate();
    if(!and1) continue;

    if(!force_vanishing_off)  {
      and1->dual_twins_push_back(g);
      if(proof_logging) gen_xor_and_van_constraint(g, and1);
      van_mon_poly_count++;
      if (verbose > 1) msg("dual twins push back (identify) %s %s", and1->get_var_name(), g->get_var_name());
    }
  }
}

/*------------------------------------------------------------------------*/
void find_vanishing_constraints(){
  find_and_propagate_xor_and();
  //find_vanishing_triangles();
}

/*------------------------------------------------------------------------*/
void find_vanishing_constraints_light(){
  find_xor_and();
}

