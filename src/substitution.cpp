/*------------------------------------------------------------------------*/
/*! \file substitution.cpp
    \brief Used to identify special subcircuits

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#include <list>

#include "substitution.h"
/*------------------------------------------------------------------------*/
// ERROR CODES:

/*------------------------------------------------------------------------*/
// Local variables

unsigned aig_idx;

bool no_cin;
bool single_gen_gate;

static Gate * carry_out;
static Gate * carry_in;
static std::vector<Gate*> outputs;
static std::vector<unsigned> original_outputs;
static std::vector<unsigned> rewritten_outputs;
static std::vector<Gate*> inputs;
static std::list<unsigned>plain_inputs;
static std::vector<Gate*> c_ins;
/*------------------------------------------------------------------------*/

static bool all_single_output() {
  for (unsigned i = 0; i < NN-1; i++) {
    Gate * n = gate(slit(i));
    if (n->parents_size() > 1) return 0;
  }
  return 1;
}

/*------------------------------------------------------------------------*/

static bool all_outputs_are_xor() {
  unsigned lit = slit(NN-1);
  if(lit < 2) return 0;

  for (unsigned i = 1; i < NN-1; i++) {
    unsigned lit = slit(i);
    if(lit < 2) return 0;
    Gate * n = gate(slit(i));
    if (!n->get_xor_gate()) return 0;
  }
  return 1;
}

/*------------------------------------------------------------------------*/

static bool slice_two_needs_carry_in_slice_zero() {
  unsigned lit2 = slit(2);
  if(lit2 < 2) return 0;
  Gate * out2 = gate(slit(2));
  unsigned lit = slit(0);
  if(lit < 2) return 0;
  Gate * out0 = gate(slit(0));
  if (out2->parents_size() > 3 && out0->parents_size() == 1) return 0;
  return 1;
}

/*------------------------------------------------------------------------*/

static bool cin_in_slice_0() {
  unsigned lit = slit(0);
  if(lit < 2) return 0;
  Gate * n = gate(slit(0));
  if (n->parents_size() > 1) return 1;
  return 0;
}

/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/

static void push_to_inputs(Gate * n) {
  inputs.push_back(n);
  n->inc_fsa_inp();
  n->mark_fsa();
}

/*------------------------------------------------------------------------*/

static void push_to_outputs(Gate * n, int i) {
  outputs.push_back(n);
  if (verbose >= 2)  msg("found output %i %s", i, n->get_var_name());
}

/*------------------------------------------------------------------------*/

static void push_to_cins(Gate * n, int i) {
  c_ins.push_back(n);
  carry_in = n;
  n->mark_fsa();
  if (verbose >= 2)  msg("found cin of slice %i %s", i, n->get_var_name());
}

/*------------------------------------------------------------------------*/

static void set_carry_in(Gate * n) {
  carry_in = n;
  n->mark_fsa();
  if (verbose >= 3)  msg("identified carry in %s ", n->get_var_name());
}

/*------------------------------------------------------------------------*/

static void identify_carry_out() {
  Gate * largest_aig_output = gate(slit(NN-1));

  if (largest_aig_output->get_xor_gate() != 1) {
    carry_out = largest_aig_output;
    push_to_outputs(carry_out, NN-1);
  } else {
    Gate * l = xor_left_child(largest_aig_output);
    Gate * r = xor_right_child(largest_aig_output);
    if (r->get_var()->get_level() > l->get_var()->get_level()) carry_out = r;
    else
      carry_out = l;
    push_to_outputs(carry_out, -1);
  }
  if (verbose >= 3)  msg("identified carry out %s", carry_out->get_var_name());
}

/*------------------------------------------------------------------------*/

static bool identify_propagate_and_generate_gates() {
// slice NN-1 contains carry, slice 0 is no XOR
  for (int i = NN-2; i > 0; i--) {
    Gate * n = gate(slit(i));

    if (i == 2 && n->parents_size() > 3) {
      assert(gate(slit(0))->parents_size() > 1);

      push_to_outputs(n, 2);
      push_to_outputs(gate(slit(1)), 1);
      push_to_outputs(gate(slit(0)), 0);

      push_to_inputs(n);
      push_to_inputs(gate(slit(1)));
      set_carry_in(gate(slit(0)));
      return 1;
    }

    Gate * internal_xor, *l = 0, *r = 0;
    if (i == 1 && n->parents_size() > 1) { internal_xor = n;
    } else {
      l = xor_left_child(n);
      r = xor_right_child(n);
      internal_xor = l->get_xor_gate() ? l : r;
    }


    int cmp = NN-1;
    if (internal_xor->parents_size() < 3) break;

    if (internal_xor->parents_size() == 3 && i < 3*cmp/4
        && !cin_in_slice_0()) {
      if (all_single_output()) break;
      else if (!booth) break;  // sp-ba-csv
      
    }


    internal_xor->mark_prop_gen_gate();
    if (verbose >= 2)
      msg("found propagate gate %s", internal_xor->get_var_name());
    Gate *g_0 = 0, *g_1 = 0;

    if (internal_xor->get_xor_gate()==1 &&
        xor_left_child(internal_xor)->parents_size() != 2 &&
        xor_right_child(internal_xor)->parents_size() != 2 &&
        (i != 1 || n->parents_size() == 1 || booth)) {
          
      Gate * internal_and = internal_xor->get_xor_and_gate();
      
      internal_and->mark_prop_gen_gate();
      if (verbose >= 2)
        msg("found generate gate %s", internal_and->get_var_name());

      aiger_and * par = is_model_and(internal_and->get_var_num());
      g_0 = gate(par->rhs0);
      g_1 = gate(par->rhs1);
      g_0->set_neg(aiger_sign(par->rhs0));
      g_1->set_neg(aiger_sign(par->rhs1));
      push_to_inputs(g_0);
      push_to_inputs(g_1);
    } else if (booth) {
      push_to_inputs(internal_xor);

      if (verbose >= 3)  msg("pushed xor %s", internal_xor->get_var_name());
      single_gen_gate = 1;
    }

    push_to_outputs(n, i);
    if (i != 1 || n->parents_size() == 1) {
      if (l->get_xor_gate()) push_to_cins(r, i);
      else
        push_to_cins(l, i);
    } else  {
      Gate * c = gate(slit(0));
      if (c->parents_size() > 1) {
         push_to_cins(c, i);
         push_to_outputs(c, 0);
      } else if (booth &&(g_0->get_xor_gate() || g_1->get_xor_gate())) {
        // bp-ba-csv
        Gate * not_xor_cin = g_0->get_xor_gate() ? g_1 : g_0;
        push_to_cins(not_xor_cin, i);
        no_cin = 1;
      } 
    }
  }
  return 1;
}

/*------------------------------------------------------------------------*/

static void fix_inputs() {
  if (!cin_in_slice_0()) return;
  

  static std::vector<Gate*> inputs_cpy;
  for (std::vector<Gate*>::const_iterator it = inputs.begin();
      it != inputs.end(); ++it) {
    Gate * n = *it;
    if (!n->get_prop_gen_gate()) { inputs_cpy.push_back(n);
    } else {
      aiger_and * and1 = is_model_and(n->get_var_num());
      if (aiger_sign(and1->rhs0) != aiger_sign(and1->rhs1)) {
        if (aiger_sign(and1->rhs0)) inputs_cpy.push_back(gate(and1->rhs0));
        if (aiger_sign(and1->rhs1)) inputs_cpy.push_back(gate(and1->rhs1));
      }  
    }
  }
  inputs = inputs_cpy;
}

/*------------------------------------------------------------------------*/

bool follow_path_and_mark_gates(Gate * n, bool init) {
  if (n->get_input() && !n->get_fsa_inp()) return 0;

  n->mark_fsa();
  if(verbose > 3) msg("marked %s", n->get_var_name());

  if (n == carry_in) return 1;
  if (n->get_fsa_inp()) return 1;

  aiger_and * and1 = is_model_and(n->get_var_num());
  Gate * l = gate(and1->rhs0);
  Gate * r = gate(and1->rhs1);

  if (!r->get_prop_gen_gate() && carry_in == r && init  && !r->get_neg()) {
     r->set_neg(aiger_sign(and1->rhs1));
  }
  if (!follow_path_and_mark_gates(r, init)) return 0;

  if (!l->get_prop_gen_gate() && carry_in == l && init && !l->get_neg()) {
     l->set_neg(aiger_sign(and1->rhs0));
  }
  if (!follow_path_and_mark_gates(l, init)) return 0;


  return 1;
}

/*------------------------------------------------------------------------*/

bool follow_all_output_paths_and_mark_gates() {
  msg("checking last stage adder");
  for (std::vector<Gate*>::const_iterator it = outputs.begin();
      it != outputs.end(); ++it) {
    Gate * n = *it;
    if (verbose >= 3)  msg("follow path starting with %s", n->get_var_name());
    bool init =(it == outputs.begin()) ? 1 : 0;
    if (!follow_path_and_mark_gates(n, init)) return 0;
  }
  return 1;
}

/*------------------------------------------------------------------------*/

static void correctly_mark_inputs() {
    for (unsigned i = 0; i< inputs.size(); i++) {
      if (inputs[i]->get_prop_gen_gate()) continue;
      if (!inputs[i]->get_aig_output()) inputs[i]->reset_fsa_inp();
    }

    for (unsigned i = M-1; i >= 1; i--) {
      Gate * n = gates[i];
      if (!n->get_prop_gen_gate()) continue;
      if (single_gen_gate && n->get_fsa_inp()) continue;

      n->reset_fsa_inp();

      aiger_and * and1 = is_model_and(n->get_var_num());
      unsigned l = aiger_strip(and1->rhs0), r = aiger_strip(and1->rhs1);

      if (!n->get_xor_gate()) {
        gate(l)->inc_fsa_inp();
        gate(r)->inc_fsa_inp();
      }
    }

    carry_in->inc_fsa_inp();

    if (single_gen_gate) {
      for (unsigned i = 0; i< inputs.size(); i++) {
        if (!inputs[i]->get_fsa_inp()) {
          inputs[i]->inc_fsa_inp();
        }
      }
    }
}

/*------------------------------------------------------------------------*/
void unmark_fsa(){
  for (unsigned i = 0; i < M; i++) {
    Gate *n = gates[i];
    n->remove_fsa();
  }
}
/*------------------------------------------------------------------------*/
bool follow_all_output_paths_cs(Gate * n, bool final) {
  if (n->get_input()) {msg("%s", n->get_var_name()); return 0;}
  if (n->get_fsa()) return 1;
  n->mark_fsa();
  if(verbose > 3) msg("marked %s %i", n->get_var_name(), n->get_xor_gate());


  if (n->get_xor_and_inp()) return 1;
  if (final && n->get_xor_gate() ==1 ) return 1;

  aiger_and * and1 = is_model_and(n->get_var_num());
  Gate * l = gate(and1->rhs0);
  Gate * r = gate(and1->rhs1);

  if (!follow_all_output_paths_cs(r, final)) return 0;
  if (!follow_all_output_paths_cs(l,final)) return 0;


  return 1;
}

/*------------------------------------------------------------------------*/

static bool try_carry_save(){
  unsigned lit = slit(NN-1);
  if(lit < 2) return 0;

  bool flag = 0;

  for (unsigned i = NN-2; i > 0; i--) {
    unsigned lit = slit(i);
    if(lit < 2) return 0;
    Gate*n = gate(slit(i));

    bool res = follow_all_output_paths_cs(n,0);
    if(!res) return 0;
    if(i > 1 && n->get_xor_gate()==1 && gate(slit(i-1))->get_xor_gate()==1 && flag) break;
    if(i > 1 && n->get_xor_gate()==1 && gate(slit(i-1))->get_xor_gate()==1 && !flag) flag = 1;
    
    
    
  }

  
  Gate * n = gate(slit(NN-1));
  if (n->get_xor_gate() == 1){
    Gate * n1 = xor_left_child(n);
    Gate * n2 = xor_right_child(n);
    if((n1->get_xor_gate() == 1) != (n2->get_xor_gate() == 1)) return 0;
    n = (n1->get_xor_gate() == 1) ? n2 : n1;
  } 
  bool res = follow_all_output_paths_cs(n,1);
  if(!res) return 0;
  




  return 1;
}



/*------------------------------------------------------------------------*/
bool identify_final_stage_adder() {
 
  if (!all_outputs_are_xor()) {
    if(try_carry_save()) return 1;
    msg("substitution not possible - not all outputs are XOR");
    unmark_fsa();
    return 0;
  } 
  if (!slice_two_needs_carry_in_slice_zero()) {
    msg("substitution not possible - carry in slice 0 not found");
    unmark_fsa();
    return 0;
  }
 
  identify_carry_out(); 
 
  if (!identify_propagate_and_generate_gates()) {
    msg("substitution not possible - propagate and generate gates not found");
    unmark_fsa();
    return 0;
  } 
  fix_inputs();

  if (!follow_all_output_paths_and_mark_gates()) {
    msg("substitution not possible - no clear boundaries");
    unmark_fsa();
    return 0;
  }

  correctly_mark_inputs();
  
  return 1;
}
/*----------------------------------------------------------------------------*/
void mark_bottom_of_circuit(Gate *g){
  g->mark_fsa();
  for (unsigned i = 0; i < num_gates; i++) {
    Gate *n = gates[i];
    if(n->get_var_level() <= g->get_var_level()) n->mark_fsa();
  }
}
