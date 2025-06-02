/*------------------------------------------------------------------------*/
/*! \file extensions.cpp
    \brief contains functions to add extension variables

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#include "extensions.h"
/*------------------------------------------------------------------------*/
unsigned extended_gates = 0;

static void set_parents_and_children_of_extension_var(Term *t, Gate *g) {
  while (t) {
    Gate *tmp = gate(t->get_var_num());
    g->children_push_back(tmp);
    tmp->parents_push_back(g);
    t = t->get_rest();
  }
}

/*------------------------------------------------------------------------*/


Term *extend_var_gates(Term *t) {
  if (num_gates == size_gates)
    die(2, "gates too small");

  int level = -2 - 2 * extended_gates;
  std::string name = std::string("t") + std::to_string(extended_gates++);
  gates[num_gates] = new Gate(M - num_gates - 1, name, level, 0, 0);
  num_gates++;
  Gate *g = gates[num_gates - 1];
  g->set_ext();

  Monomial **intern_mstack = new Monomial *[2];

  intern_mstack[0] = new Monomial(minus_one, new_term(g->get_var()));
  intern_mstack[1] = new Monomial(one, t);
  Polynomial *p = new Polynomial(intern_mstack, 2, t->degree());
  g->set_gate_constraint(p);
  set_parents_and_children_of_extension_var(t, g);

  if (verbose >= 2) {
    msg("added extension var: %s", g->get_var_name());
    msg_nl("extension poly: ");
    p->print(stdout);
  }

  if(proof_logging){
    print_pac_extension_rule_for_mon(proof_file, g, t, p);
  }

 

  return p->get_lt();
}
/*------------------------------------------------------------------------*/
void adjust_level_of_extended_gates() {

  if (verbose > 3)
    msg("adding extension variables at correct level");
  for (unsigned i = M - 1 + MM; i < num_gates; i++) {
    Gate *g = gates[i];
    if (verbose > 3)
      msg("old level %s %i", g->get_var_name(), g->get_var_level());
    g->set_var_level(g->get_var_level() + 2 * NN + 2);
  }

  for (unsigned i = 0; i < NN; i++) {
    Gate *g = gates[i];
    if (verbose > 3)
      msg("old level %s %i", g->get_var_name(), g->get_var_level());
    g->set_var_level(g->get_var_level() - 2 * extended_gates - 2);
  }

  if (verbose > 3) {
    for (unsigned i = M - 1 + MM; i < num_gates; i++) {
      Gate *g = gates[i];
      msg("new level %s %i", g->get_var_name(), g->get_var_level());
    }

    for (unsigned i = 0; i < NN; i++) {
      Gate *g = gates[i];
      msg("new level %s %i", g->get_var_name(), g->get_var_level());
    }
  }
 
}

/*------------------------------------------------------------------------*/