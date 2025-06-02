/*------------------------------------------------------------------------*/
/*! \file subcircuit.cpp
    \brief contains functions to identify a subcircuit

  This file contains all functions to cut out a subcircuit

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#include "subcircuit.h"

#include <flint/fmpq.h>

#include <algorithm>
#include <array>
#include <cstdio>
#include <flint/fmpq_mat.h>
#include <iostream>
#include <memory>
#include <random>
#include <tuple>
#include <unordered_set>

#include "matrix.h"
#include "term.h"
/*------------------------------------------------------------------------*/
std::set<Gate*, LargerGate> var;
std::set<Gate*, LargerGate> sc_inputs;  // check whether sets are actually the
                                        // most appropriate data structures
std::set<Gate*, SmallerGate> gate_poly;

size_t fanout_size_last_call = 0;



/*------------------------------------------------------------------------*/
static void
print_subcircuit_guessing(Gate* root, int depth, size_t fanout_size) {
  msg("");
  msg("subcircuit with root %s at dist %i, depth %i, fanout size %i:",
      root->get_var_name(),
      root->get_dist(),
      depth,
      fanout_size);
  msg("%i gates:", gate_poly.size());
  for (auto& g : gate_poly) {
    msg_nl("  %s, dist %i, parentsize %i   ",
           g->get_var_name(),
           g->get_dist(),
           g->aig_parents_size());
    g->print_gate_constraint(stdout);
  }
  msg("");
  msg("%i inputs:", sc_inputs.size());
  for (auto& g : sc_inputs) {
    if (g->get_dist()) {
      msg_nl("  %s, dist %i, parentsize %i   ",
             g->get_var_name(),
             g->get_dist(),
             g->aig_parents_size());
      g->print_gate_constraint(stdout);
    } else {
      msg("  %s, dist %i, parentsize %i   ",
          g->get_var_name(),
          g->get_dist(),
          g->aig_parents_size());
    }
  }
  msg("");
  msg("");
}
/*------------------------------------------------------------------------*/
static void
add_children_guessing(Gate* g, Gate* root, int depth, size_t fanout_size, bool init) {
  if (g->get_input())
    return;
  if (!depth)
    return;
  if (fanout_size && !init && g->aig_parents_size() >= fanout_size && !root->is_aig_child(g)) {
    fanout_size_last_call = g->aig_parents_size();
    return;
  }

  if (!gate_poly.contains(g)) {
    gate_poly.insert(g);
  }
  sc_inputs.erase(g);
  var.insert(g);
  if (verbose > 3)
    msg("added child %s, parentsize %i", g->get_var_name(), g->parents_size());

  for (auto& gc : g->get_aig_children()) {
    if (gc->get_elim())
      continue;
    var.insert(gc);
    if (!gate_poly.contains(gc)) {
      if (verbose > 2)
        msg("inserted %s to sc_inputs", gc->get_var_name());
      sc_inputs.insert(gc);
    }
  }

  for (auto& gc : g->get_children()) {
    add_children_guessing(gc, root, depth - 1, fanout_size, 0);
  }
}
/*------------------------------------------------------------------------*/
static void
push_inputs_guessing(size_t fanout_size) {
  std::vector<Gate*> placeholders;
  for (auto& g : sc_inputs) {
    if (g->aig_parents_size() == 1 && g->aig_parents_size() < fanout_size && !g->get_xor_and_gate() && !g->get_input()) {
      gate_poly.insert(g);
      placeholders.push_back(g);
      
      if (verbose > 1)
        msg("pushed single input %s, parentsize %i",
            g->get_var_name(),
            g->aig_parents_size());

      for (auto& gc : g->get_aig_children()) {
        if (gc->get_elim())
          continue;
        var.insert(gc);
        if (!gate_poly.contains(gc)) {
          //     msg("inserted %s to sc_inputs", gc->get_var_name());
          sc_inputs.insert(gc);
        }
      }
    }

    if (g->get_dist() > 0) {
      bool flag = 1;
      for (auto& gc : g->get_aig_children()) {
        if (!sc_inputs.contains(gc)) {
          flag = 0;
          break;
        }
      }
      if (flag) {
        gate_poly.insert(g);
        placeholders.push_back(g);
        if (verbose > 1)
          msg("pushed input %s whose inputs are inputs, parentsize %i",
              g->get_var_name(),
              g->aig_parents_size());
      }
    }
  }

  for (auto& g : placeholders) {
    sc_inputs.erase(g);
  }
}
/*------------------------------------------------------------------------*/

static void
push_pp_guessing() {
  std::vector<Gate*> placeholders;
  for (auto& g : sc_inputs) {
    if (g->get_pp()) {
      gate_poly.insert(g);
      placeholders.push_back(g);
      if (verbose > 1)
        msg(
            "pushed pp %s, parentsize %i", g->get_var_name(), g->aig_parents_size());

      for (auto& gc : g->get_aig_children()) {
        var.insert(gc);
        sc_inputs.insert(gc);
      }
    }
  }

  for (auto& g : placeholders) {
    sc_inputs.erase(g);
  }
}
/*------------------------------------------------------------------------*/
static void
add_spouses_guessing(Gate* g) {
  for (auto& gc : g->get_aig_children()) {
    if (gc->get_input())
      continue;
    for (auto& g_sib_val : gc->get_aig_parents()) {
      Gate * g_sib = gate(g_sib_val);
      if (g == g_sib)
        continue;
      if (g_sib->get_var_level() > g->get_var_level())
        continue;
      if (g_sib->get_elim())
        continue;

      if (!gate_poly.contains(g_sib) && !g_sib->get_input()) {
        gate_poly.insert(g_sib);
      }
      sc_inputs.erase(g_sib);
      var.insert(g_sib);
      if (verbose > 1)
        msg(
            "added spouse %s, dist %i", g_sib->get_var_name(), g_sib->get_dist());

      for (auto& gsib_c : g_sib->get_aig_children()) {
        if (gsib_c->get_elim())
          continue;
        var.insert(gsib_c);
        if (!gate_poly.contains(gsib_c))
          sc_inputs.insert(gsib_c);
      }
    }
  }
}

/*------------------------------------------------------------------------*/
static void
add_parents_guessing(Gate* node, Gate* g) {
  for (auto& node_p_val : node->get_aig_parents()) {
    Gate * node_p = gate(node_p_val);
    if (node_p->get_var_level() > g->get_var_level())
      continue;
    if (node_p->get_output())
      continue;
    if (node_p->get_elim())
      continue;

    bool flag = 0;
    for (auto& node_p_c : node_p->get_aig_children()) {
      if (node == node_p_c)
        continue;
      if (!var.contains(node_p_c)) {
        flag = 1;
        break;
      }
    }
    if (flag)
      continue;

    if (!gate_poly.contains(node_p) && !node_p->get_input()) {
      gate_poly.insert(node_p);

      sc_inputs.erase(node_p);
      if (verbose > 1)
        msg("added parent %s, dist %i", node_p->get_var_name(), g->get_dist());
      var.insert(node_p);
      add_parents_guessing(node_p, g);
    }
  }
}

/*------------------------------------------------------------------------*/
static void
add_common_ancestors_guessing(Gate* g) {
  for (auto& node : var) {
    for (auto& node_p_val : node->get_aig_parents()) {
      Gate * node_p = gate(node_p_val);
      if (node_p == g)
        continue;
      if (node_p->get_var_level() > g->get_var_level())
        continue;

      if (node_p->get_output())
        continue;
      if (node_p->get_elim())
        continue;

      bool flag = 0;

      for (auto& node_p_c : node_p->get_aig_children()) {
        if (node == node_p_c)
          continue;
        if (!var.contains(node_p_c)) {
          flag = 1;
          break;
        }
      }

      if (flag)
        continue;

      if (!gate_poly.contains(node_p) && !node_p->get_input()) {
        if (verbose > 1)
          msg("added common ancestor %s, dist %i",
              node_p->get_var_name(),
              g->get_dist());
        gate_poly.insert(node_p);
        sc_inputs.erase(node_p);
        var.insert(node_p);
        add_parents_guessing(node_p, g);
      }
    }
  }
}
/*------------------------------------------------------------------------*/
static void
add_ancestors_with_same_dist_guessing(Gate* g) {
  for (auto& node : var) {
    for (auto& node_p_val : node->get_aig_parents()) {
      Gate * node_p = gate(node_p_val);
      if (node_p == g)
        continue;
      if (node_p->get_dist() != g->get_dist())
        continue;
      if (node_p->get_var_level() > g->get_var_level())
        continue;

      if (node_p->get_output())
        continue;
      if (node_p->get_elim())
        continue;

      for (auto& node_p_c : node_p->get_aig_children()) {
        if (node == node_p_c)
          continue;
        if (!var.contains(node_p_c)) {
          var.insert(node_p_c);
          sc_inputs.insert(node_p_c);
        }
      }

      if (!gate_poly.contains(node_p) && !node_p->get_input()) {
        if (verbose > 1)
          msg("added same dist ancestor %s, dist %i",
              node_p->get_var_name(),
              g->get_dist());
        gate_poly.insert(node_p);
        sc_inputs.erase(node_p);
        var.insert(node_p);
        add_parents_guessing(node_p, g);
      }
    }
  }
}
/*------------------------------------------------------------------------*/
static bool
expand_inputs_guessing(Gate* inp_g, int depth, size_t fanout_size) {
  bool flag_exit = 1;
  for (auto& g : sc_inputs) {
    if (!g->get_input()) {
      flag_exit = 0;
    }
  }
  if (flag_exit) return 1;

  // check whether there is a suitable candidate for expansion
  Gate* exp = 0;
  std::vector<Gate*> placeholders;
  for (auto& g : sc_inputs) {
    if(g->get_input()) continue;
    if (g->aig_parents_size() < fanout_size) {
      exp = g;
      break;
    }
    if (g->get_dist() + depth > inp_g->get_dist()) {
      exp = g;
      break;
    }
  }

  bool flag = 1;
  int i = 1;
  while (!exp && flag) {
    flag = 0;
    for (auto& g : sc_inputs) {
      if (g->get_dist() > 1 && g->aig_parents_size() < fanout_size + i) {
        exp = g;
        break;
      }
      if (g->get_dist() > 1) flag = 1;
    }
    if (flag) i++;
  }
  if (!exp) return 0;
  // count how many suitable candidates for expansions are there:

  for (auto& g : sc_inputs) {
    if (g->get_dist() > 0 && g->aig_parents_size() <= fanout_size) {
      placeholders.push_back(g);
    }
  }
  // If only 3 suitable nodes can be expanded, we expand all of them
  if (placeholders.size() > 0 && placeholders.size() < 4) {
    for (auto& exp : placeholders) {
      gate_poly.insert(exp);
      if(verbose > 1) msg("expanded by %s", exp->get_var_name());
      sc_inputs.erase(exp);

      for (auto& gc : exp->get_aig_children()) {
        if (gc->get_elim())
          continue;
        var.insert(gc);
        if (!gate_poly.contains(gc))
          sc_inputs.insert(gc);
      }
    }

    return 1;
  }

  // expand single gate with max distance
  for (auto& g : sc_inputs) {
    if (g->get_dist() > exp->get_dist() && g->aig_parents_size() < fanout_size)
      exp = g;
  }

  gate_poly.insert(exp);
  sc_inputs.erase(exp);

  if (verbose > 1)
    msg("expand input %s, parentsize %i",
        exp->get_var_name(),
        exp->aig_parents_size());

  for (auto& gc : exp->get_aig_children()) {
    if (gc->get_elim())
      continue;
    var.insert(gc);
    if (!gate_poly.contains(gc))
      sc_inputs.insert(gc);
  }
  return 1;
}
/*------------------------------------------------------------------------*/



static bool subcircuit_for_guessing(Gate* g,
  int depth,
  size_t fanout_size,
  int init,
  bool single_expand){
bool expand = true;
    if (!single_expand)  // collect all children for certain depth
      add_children_guessing(g, g, depth, fanout_size, 1);
    else
      expand = expand_inputs_guessing(g, depth, fanout_size);  // expand single inputs
    if (!expand)
      return false;

    // expand based on inputs and root node
    add_ancestors_with_same_dist_guessing(g);
    add_spouses_guessing(g);
    push_inputs_guessing(fanout_size);
    push_pp_guessing();
    add_common_ancestors_guessing(g);

    if (verbose > 1)
      print_subcircuit_guessing(g, depth, fanout_size);

    return true;
  }


/*------------------------------------------------------------------------*/
// Identify Sub-Circuit no force guess
/*------------------------------------------------------------------------*/
static void
add_children(Gate* g, Gate* root, int depth, size_t fanout_size, bool init) {
  if (g->get_input())
    return;
  if (!depth)
    return;
  if (fanout_size && !init && g->parents_size() >= fanout_size && !root->is_child(g)) {
    fanout_size_last_call = g->parents_size();
    return;
  }

  if (!gate_poly.contains(g)) {
    gate_poly.insert(g);
  }
  sc_inputs.erase(g);
  var.insert(g);
  if (verbose > 3)
    msg("added child %s, parentsize %i", g->get_var_name(), g->parents_size());

  for (auto& gc : g->get_children()) {
    if (gc->get_elim())
      continue;
    var.insert(gc);
    if (!gate_poly.contains(gc)) {
      if (verbose > 2)
        msg("inserted %s to sc_inputs", gc->get_var_name());
      sc_inputs.insert(gc);
    }
  }

  for (auto& gc : g->get_children()) {
    add_children(gc, root, depth - 1, fanout_size, 0);
  }
}
/*------------------------------------------------------------------------*/
static void
push_inputs(size_t fanout_size) {
  std::vector<Gate*> placeholders;
  for (auto& g : sc_inputs) {
    if (g->parents_size() == 1 && g->parents_size() < fanout_size && !g->get_xor_and_gate() && !g->get_input()) {
      gate_poly.insert(g);
      placeholders.push_back(g);
      
      if (verbose > 1)
        msg("pushed single input %s, parentsize %i",
            g->get_var_name(),
            g->parents_size());

      for (auto& gc : g->get_children()) {
        if (gc->get_elim())
          continue;
        var.insert(gc);
        if (!gate_poly.contains(gc)) {
          //     msg("inserted %s to sc_inputs", gc->get_var_name());
          sc_inputs.insert(gc);
        }
      }
    }

    if (g->get_dist() > 0) {
      bool flag = 1;
      for (auto& gc : g->get_children()) {
        if (!sc_inputs.contains(gc)) {
          flag = 0;
          break;
        }
      }
      if (flag) {
        gate_poly.insert(g);
        placeholders.push_back(g);
        if (verbose > 1)
          msg("pushed input %s whose inputs are inputs, parentsize %i",
              g->get_var_name(),
              g->parents_size());
      }
    }
  }

  for (auto& g : placeholders) {
    sc_inputs.erase(g);
  }
}
/*------------------------------------------------------------------------*/

static void
push_pp() {
  std::vector<Gate*> placeholders;
  for (auto& g : sc_inputs) {
    if (g->get_pp()) {
      gate_poly.insert(g);
      placeholders.push_back(g);
      if (verbose > 1)
        msg(
            "pushed pp %s, parentsize %i", g->get_var_name(), g->parents_size());

      for (auto& gc : g->get_children()) {
        var.insert(gc);
        sc_inputs.insert(gc);
      }
    }
  }

  for (auto& g : placeholders) {
    sc_inputs.erase(g);
  }
}
/*------------------------------------------------------------------------*/
static void
add_spouses(Gate* g) {
  for (auto& gc : g->get_children()) {
    if (gc->get_input())
      continue;
    for (auto& g_sib : gc->get_parents()) {
      if (g == g_sib)
        continue;
      if (g_sib->get_var_level() > g->get_var_level())
        continue;
      if (g_sib->get_elim())
        continue;

      if (!gate_poly.contains(g_sib) && !g_sib->get_input()) {
        gate_poly.insert(g_sib);
      }
      sc_inputs.erase(g_sib);
      var.insert(g_sib);
      if (verbose > 1)
        msg(
            "added spouse %s, dist %i", g_sib->get_var_name(), g_sib->get_dist());

      for (auto& gsib_c : g_sib->get_children()) {
        if (gsib_c->get_elim())
          continue;
        var.insert(gsib_c);
        if (!gate_poly.contains(gsib_c))
          sc_inputs.insert(gsib_c);
      }
    }
  }
}

/*------------------------------------------------------------------------*/
static void
add_parents(Gate* node, Gate* g) {
  for (auto& node_p : node->get_parents()) {
    if (node_p->get_var_level() > g->get_var_level())
      continue;
    if (node_p->get_output())
      continue;
    if (node_p->get_elim())
      continue;

    bool flag = 0;
    for (auto& node_p_c : node_p->get_children()) {
      if (node == node_p_c)
        continue;
      if (!var.contains(node_p_c)) {
        flag = 1;
        break;
      }
    }
    if (flag)
      continue;

    if (!gate_poly.contains(node_p) && !node_p->get_input()) {
      gate_poly.insert(node_p);

      sc_inputs.erase(node_p);
      if (verbose > 1)
        msg("added parent %s, dist %i", node_p->get_var_name(), g->get_dist());
      var.insert(node_p);
      add_parents(node_p, g);
    }
  }
}

/*------------------------------------------------------------------------*/
static void
add_common_ancestors(Gate* g) {
  for (auto& node : var) {
    for (auto& node_p : node->get_parents()) {
      if (node_p == g)
        continue;
      if (node_p->get_var_level() > g->get_var_level())
        continue;

      if (node_p->get_output())
        continue;
      if (node_p->get_elim())
        continue;

      bool flag = 0;

      for (auto& node_p_c : node_p->get_children()) {
        if (node == node_p_c)
          continue;
        if (!var.contains(node_p_c)) {
          flag = 1;
          break;
        }
      }

      if (flag)
        continue;

      if (!gate_poly.contains(node_p) && !node_p->get_input()) {
        if (verbose > 1)
          msg("added common ancestor %s, dist %i",
              node_p->get_var_name(),
              g->get_dist());
        gate_poly.insert(node_p);
        sc_inputs.erase(node_p);
        var.insert(node_p);
        add_parents(node_p, g);
      }
    }
  }
}
/*------------------------------------------------------------------------*/
static void
add_ancestors_with_same_dist(Gate* g) {
  for (auto& node : var) {
    for (auto& node_p : node->get_parents()) {
      if (node_p == g)
        continue;
      if (node_p->get_dist() != g->get_dist())
        continue;
      if (node_p->get_var_level() > g->get_var_level())
        continue;

      if (node_p->get_output())
        continue;
      if (node_p->get_elim())
        continue;

      for (auto& node_p_c : node_p->get_children()) {
        if (node == node_p_c)
          continue;
        if (!var.contains(node_p_c)) {
          var.insert(node_p_c);
          sc_inputs.insert(node_p_c);
        }
      }

      if (!gate_poly.contains(node_p) && !node_p->get_input()) {
        if (verbose > 1)
          msg("added same dist ancestor %s, dist %i",
              node_p->get_var_name(),
              g->get_dist());
        gate_poly.insert(node_p);
        sc_inputs.erase(node_p);
        var.insert(node_p);
        add_parents(node_p, g);
      }
    }
  }
}
/*------------------------------------------------------------------------*/
static bool
expand_inputs(Gate* inp_g, int depth, size_t fanout_size) {
  bool flag_exit = 1;
  for (auto& g : sc_inputs) {
    if (!g->get_input()) {
      flag_exit = 0;
    }
  }
  if (flag_exit) return 1;

  // check whether there is a suitable candidate for expansion
  Gate* exp = 0;
  std::vector<Gate*> placeholders;
  for (auto& g : sc_inputs) {
    if(g->get_input()) continue;
    if (g->parents_size() < fanout_size) {
      exp = g;
      break;
    }
    if (g->get_dist() + depth > inp_g->get_dist()) {
      exp = g;
      break;
    }
  }

  bool flag = 1;
  int i = 1;
  while (!exp && flag) {
    flag = 0;
    for (auto& g : sc_inputs) {
      if (g->get_dist() > 1 && g->parents_size() < fanout_size + i) {
        exp = g;
        break;
      }
      if (g->get_dist() > 1) flag = 1;
    }
    if (flag) i++;
  }
  if (!exp) return 0;
  // count how many suitable candidates for expansions are there:

  for (auto& g : sc_inputs) {
    if (g->get_dist() > 0 && g->parents_size() <= fanout_size) {
      placeholders.push_back(g);
    }
  }
  // If only 3 suitable nodes can be expanded, we expand all of them
  if (placeholders.size() > 0 && placeholders.size() < 4) {
    for (auto& exp : placeholders) {
      gate_poly.insert(exp);
      if(verbose > 1) msg("expanded by %s", exp->get_var_name());
      sc_inputs.erase(exp);

      for (auto& gc : exp->get_children()) {
        if (gc->get_elim())
          continue;
        var.insert(gc);
        if (!gate_poly.contains(gc))
          sc_inputs.insert(gc);
      }
    }

    return 1;
  }

  // expand single gate with max distance
  for (auto& g : sc_inputs) {
    if (g->get_dist() > exp->get_dist() && g->parents_size() < fanout_size)
      exp = g;
  }

  gate_poly.insert(exp);
  sc_inputs.erase(exp);

  if (verbose > 1)
    msg("expand input %s, parentsize %i",
        exp->get_var_name(),
        exp->parents_size());

  for (auto& gc : exp->get_children()) {
    if (gc->get_elim())
      continue;
    var.insert(gc);
    if (!gate_poly.contains(gc))
      sc_inputs.insert(gc);
  }
  return 1;
}
/*------------------------------------------------------------------------*/

static void
print_subcircuit(Gate* root, int depth, size_t fanout_size) {
  msg("");
  msg("subcircuit with root %s at dist %i, depth %i, fanout size %i:",
      root->get_var_name(),
      root->get_dist(),
      depth,
      fanout_size);
  msg("%i gates:", gate_poly.size());
  for (auto& g : gate_poly) {
    msg_nl("  %s, dist %i, parentsize %i   ",
           g->get_var_name(),
           g->get_dist(),
           g->parents_size());
    g->print_gate_constraint(stdout);
  }
  msg("");
  msg("%i inputs:", sc_inputs.size());
  for (auto& g : sc_inputs) {
    if (g->get_dist()) {
      msg_nl("  %s, dist %i, parentsize %i   ",
             g->get_var_name(),
             g->get_dist(),
             g->parents_size());
      g->print_gate_constraint(stdout);
    } else {
      msg("  %s, dist %i, parentsize %i   ",
          g->get_var_name(),
          g->get_dist(),
          g->parents_size());
    }
  }
  msg("");
  msg("");
}
/*------------------------------------------------------------------------*/
static void
gen_fsa_subcircuit(Gate* g) {
  for (unsigned i = num_gates-1; i > 0; i--) {
    Gate* n = gates[i];
    if (n->get_elim() == 1)
      continue;
    if (n->get_var_level() > g->get_var_level())
      continue;
    if (!n->get_fsa())
      continue;

    if (n->get_input()) {
      sc_inputs.insert(n);
      var.insert(n);
    } else {
      bool flag = 0;
      for (auto& nc : n->get_children()) {
        if (!nc->get_fsa()) {
          flag = 1;
          break;
        }
      }
      if (flag) {
        sc_inputs.insert(n);
        var.insert(n);
      } else {
        gate_poly.insert(n);
        var.insert(n);
      }
    }
  }
}
/*------------------------------------------------------------------------*/
bool is_internal_fsa(Gate* g) {
  if (!g->get_fsa())
    return false;
  if(g->get_input()) return false;
  for (auto& gc : g->get_children()) {
    if (!gc->get_fsa())
      return false;
  }
  return true;
}
/*------------------------------------------------------------------------*/
static bool
get_subcircuit(Gate* g,
               int depth,
               size_t fanout_size,
               int init,
               bool single_expand) {
  if (init == 1) {
    var.clear();
    gate_poly.clear();
    sc_inputs.clear();
  }

  // if g belongs to specially marked circuit collect all nodes with same marking
  if (is_internal_fsa(g)) {
    gen_fsa_subcircuit(g);
    if (verbose > 1)
      print_subcircuit(g, depth, fanout_size);
 
    return true;
  } else if (force_guessing){
    return subcircuit_for_guessing(g,depth,fanout_size,init,single_expand);

  } else {
    bool expand = true;
    if (!single_expand)  // collect all children for certain depth
      add_children(g, g, depth, fanout_size, 1);
    else
      expand = expand_inputs(g, depth, fanout_size);  // expand single inputs
    if (!expand)
      return false;

    // expand based on inputs and root node
    add_ancestors_with_same_dist(g);
    add_spouses(g);
    push_inputs(fanout_size);
    push_pp();
    add_common_ancestors(g);

 
    if (verbose > 1)
      print_subcircuit(g, depth, fanout_size);

    return true;
  }
}
/*------------------------------------------------------------------------*/
static inline Normalized_poly
normalize(const Polynomial* p, std::map<Var*, size_t>& var_to_id) {
  Normalized_poly np;

  for (size_t i = 0; i < p->len(); i++) {
    Monomial* m = p->get_mon(i);
    mpz_class c(m->coeff);
    Term* t = m->get_term();
    std::vector<size_t> nt;
    if (!t)
      nt.push_back(0);
    while (t) {
      Var* v = t->get_var();
      nt.push_back(var_to_id[v]);
      t = t->get_rest();
    }
    np.emplace_back(c, nt);
  }

  return np;
}
/*------------------------------------------------------------------------*/
static void
compress_subcircuit(std::set<Gate*, SmallerGate>& subcircuit,
                    std::vector<Normalized_poly>& res,
                    std::map<Var*, size_t>& var_to_id) {
  res.clear();
  var_to_id.clear();
  // save id = 0 for constant coefficient
  size_t id = 1;

  for (const auto& gate : subcircuit) {
    Polynomial* g;
    if (gate->get_nf())
      g = gate->get_nf();
    else {
      g = gate->get_gate_constraint();
      g = unflip_poly_and_remove_van_mon(g);
      gate->set_nf(g);
    }

    for (size_t i = 0; i < g->len(); i++) {
      Term* t = g->get_mon(i)->get_term();
      while (t) {
        Var* v = t->get_var();
        if (var_to_id.find(v) == var_to_id.end()) {
          var_to_id[v] = id++;
        }
        t = t->get_rest();
      }
    }
    Normalized_poly h = normalize(g, var_to_id);
    res.push_back(h);
    if (!gate->get_nf())
      delete (g);
  }
}
/*------------------------------------------------------------------------*/
bool get_and_compress_subcircuit(Gate* g,
                                 int depth,
                                 size_t fanout_size,
                                 int init,
                                 bool single_expand,
                                 std::vector<Normalized_poly>& normalized,
                                 std::map<Var*, size_t>& var_to_id) {
  double pre_circuit_time = process_time();
  if (!get_subcircuit(g, depth, fanout_size, init, single_expand)){
    find_circuit_time += (process_time() - pre_circuit_time);
    return false;
  }
    

  // we do not compress internal-fsa as it will not be cached
  if (!is_internal_fsa(g) || force_fglm) {
    compress_subcircuit(gate_poly, normalized, var_to_id);
  }

  find_circuit_time += (process_time() - pre_circuit_time);
  return true;
}
/*------------------------------------------------------------------------*/
// Computing normalforms
/*------------------------------------------------------------------------*/
std::vector<Polynomial*>
compute_normalforms(std::vector<Polynomial*>* used_van_poly, std::vector<Polynomial*>* new_nf_poly) {
  if (gate_poly.size() == 0)
    return {};

  if (verbose > 2)
    msg("starting computing normal forms top down");

  if (verbose > 2) {
    msg("input:");
    for (const auto& gatep : gate_poly) {
      if (gatep->get_nf()) {
        msg_nl("recycled nf ");
        gatep->print_nf(stdout);
      } else {
        msg_nl("");
        gatep->print_gate_constraint(stdout);
      }
    }

    msg("");
  }

  std::vector<Polynomial*> input_poly;
  for (auto rit = gate_poly.rbegin(); rit != gate_poly.rend(); rit++) {
    Gate* gatep = *rit;
    Polynomial* gpol_raw = 0;
    if (gatep->get_nf()) {
      gpol_raw = gatep->get_nf();
      input_poly.push_back(gpol_raw);
    } else {
      die(123,
          "mismatch - gate %s has no normal form from compress circuit",
          gatep->get_var_name());
    }
  }

  std::vector<Polynomial*> rewritten;
  for (auto it = input_poly.begin(); it != input_poly.end(); it++) {
    Polynomial* gpol = *it;

    for (auto rit_inner = it; rit_inner != input_poly.end(); ++rit_inner) {
      Polynomial* gatep_inner = *rit_inner;

      if (gpol->get_lt() == gatep_inner->get_lt())
        continue;
      if (gpol->len() == 1)
        continue;

      if (verbose > 2) {
        msg_nl("reducing by:");
        gatep_inner->print(stdout);
      }
      Polynomial* tmp = reduce_by_one_poly(gpol, gatep_inner);
      if (tmp->degree() > 1) {
        Polynomial* tmp2 = remove_vanishing_monomials(tmp, used_van_poly);
        delete (tmp);
        tmp = tmp2;
      }

      if(!proof_logging) check_if_propagate(tmp);

      if (verbose > 2) {
        msg_nl("result:");
        tmp->print(stdout);
      }
      delete (gpol);
      gpol = tmp;
    }

    if(!proof_logging) check_if_propagate(gpol);

    rewritten.push_back(gpol);

    gate(gpol->get_lt()->get_var_num())->set_nf(gpol);
    if(proof_logging) new_nf_poly->push_back(gpol);
  }

  if (verbose > 2) {
    msg("Output of normal forms");
    for (const auto& gatep : rewritten) {
      gatep->print(stdout);
    }
  }
  assert(rewritten.size() > 0);

  return rewritten;
}

/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
std::deque<std::map<Gate*, bool>> collected_assignments;
/*------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
static bool call_kissat(std::vector<std::vector<int>> cnf_clauses, std::map<int, Gate*> ids) {
  count_kissat_call++;

  int max_lit = ids.size();
  kissat* solver = kissat_init();

  // Add CNF clauses to the solver
  for (const auto& clause : cnf_clauses) {
    for (int lit : clause) {
      kissat_add(solver, lit);  // Add literal
    }
    kissat_add(solver, 0);  // Terminate clause
  }

  // Solve the CNF
  kissat_set_option(solver, "quiet", 1);
  int result = kissat_solve(solver);

  // Print result
  if (result == 10) {  // SAT
    if (verbose > 2) std::cout << "SATISFIABLE\n";
    std::map<Gate*, bool> assignment;

    // Retrieve and print the satisfying assignment
    for (int var = 1; var <= max_lit; var++) {
      int value = kissat_value(solver, var);
      if (var <= max_lit) assignment.insert({ids[var], value > 0 ? 1 : 0});
    }
    collected_assignments.push_back(assignment);

  } else if (result == 20) {  // UNSAT
    if (verbose > 2) std::cout << "UNSATISFIABLE\n";
  } else {
    if (verbose > 2) std::cout << "UNKNOWN result\n";
  }

  // Clean up
  kissat_release(solver);
  return result == 20;
}
/*------------------------------------------------------------------------*/
static std::tuple<std::map<Gate*, int>, std::map<int, Gate*>>var_cnf_mapping(std::vector<Var*> vars){
    // map literals to names
    std::map<Gate*, int> lit_id;
    std::map<int, Gate*> inverse_lit_id;
    lit_id.clear();
    int id = 1;
  
    // Generate the mappings
    unsigned count = 0;
    id = 1;
    for(auto & v: vars){
      Gate * g = gate(v->get_num());
  
      // check that vars are unique
      if(lit_id.find(g) != lit_id.end()) die(4, "Gate %s already has id %i", lit_id[g]);
  
      inverse_lit_id[id] = g;
      lit_id[g] = id++;
    }

  return {lit_id, inverse_lit_id}; 
}


/*------------------------------------------------------------------------*/
static std::vector<std::vector<int>> translate_aig_part_to_cnf(std::map<Gate*, int> lit_id){
  PB2CNF pb2cnf;
  std::vector<std::vector<int>> cnf_clauses;
  int firstFreshVariable = lit_id.size() + 1;
  int rhs = 1;
  std::vector<int64_t> weights2 = {1, 1};
  std::vector<int64_t> weights3 = {1, 1, 1};

  // Encoding the AIG
  for(auto&g: gate_poly){
    if(!g->is_extension()){
      aiger_and* and1 = is_model_and(g->get_var_num());
      int lit_id_lhs = lit_id[g];
      int lit_id_rhs0 = aiger_sign(and1->rhs0) ? -lit_id[gate(and1->rhs0)] : lit_id[gate(and1->rhs0)];
      int lit_id_rhs1 = aiger_sign(and1->rhs1) ? -lit_id[gate(and1->rhs1)] : lit_id[gate(and1->rhs1)];

      std::vector<int> literals = {-lit_id_lhs, lit_id_rhs0};
      firstFreshVariable = pb2cnf.encodeGeq(weights2, literals, rhs, cnf_clauses, firstFreshVariable) + 1;

      std::vector<int> literals2 = {-lit_id_lhs, lit_id_rhs1};
      firstFreshVariable = pb2cnf.encodeGeq(weights2, literals2, rhs, cnf_clauses, firstFreshVariable) + 1;

      std::vector<int> literals3 = {lit_id_lhs, -lit_id_rhs0, -lit_id_rhs1};
      firstFreshVariable = pb2cnf.encodeGeq(weights3, literals3, rhs, cnf_clauses, firstFreshVariable) + 1;

    } else {
      int lit_id_lhs = lit_id[g];
      int lit_id_rhs0 = lit_id[g->children_front()];
      int lit_id_rhs1 = lit_id[g->children_back()];

      std::vector<int> literals = {-lit_id_lhs, lit_id_rhs0};
      firstFreshVariable = pb2cnf.encodeGeq(weights2, literals, rhs, cnf_clauses, firstFreshVariable) + 1;

      std::vector<int> literals2 = {-lit_id_lhs, lit_id_rhs1};
      firstFreshVariable = pb2cnf.encodeGeq(weights2, literals2, rhs, cnf_clauses, firstFreshVariable) + 1;

      std::vector<int> literals3 = {lit_id_lhs, -lit_id_rhs0, -lit_id_rhs1};
      firstFreshVariable = pb2cnf.encodeGeq(weights3, literals3, rhs, cnf_clauses, firstFreshVariable) + 1;
    }
  }  
  return cnf_clauses;
}
/*------------------------------------------------------------------------*/

static std::vector<std::vector<int>> translate_poly_to_cnf(Polynomial *p ,std::map<Gate*, int> lit_id, std::vector<std::vector<int>> cnf_clauses, bool second_call){
  PB2CNF pb2cnf;
    // Stores generated CNF clauses
  int firstFreshVariable = lit_id.size() + 1;
  
  // Encoding the target polynomial
  Polynomial* p_print = p;
  if (second_call) p_print = multiply_poly_with_constant(p, minus_one);

  std::vector<int64_t> poly_weights;
  std::vector<int> poly_ids;

  for (size_t i = 0; i < p_print->len(); i++) {
    Monomial* m = p_print->get_mon(i);
    if (m->get_term()) {
      poly_weights.push_back(mpz_get_si(m->coeff));
      poly_ids.push_back(lit_id[gate(m->get_term()->get_var_num())]);
    }
  }

  if (!p_print->get_mon((p_print->len()) - 1)->get_term()) {  // Last term is not constant

    int rhs_p = -1 * mpz_get_si(p_print->get_mon((p_print->len()) - 1)->coeff) + 1;
    firstFreshVariable = pb2cnf.encodeGeq(poly_weights, poly_ids, rhs_p, cnf_clauses, firstFreshVariable) + 1;

  } else {
    firstFreshVariable = pb2cnf.encodeGeq(poly_weights, poly_ids, 1, cnf_clauses, firstFreshVariable) + 1;
  }
  if (second_call) delete (p_print);
  return cnf_clauses;
}


/*------------------------------------------------------------------------*/
std::random_device seed;
std::mt19937 generator(seed());
std::uniform_int_distribution<uint32_t> uniform;
std::unordered_map<Var*,int> var_to_col;

static void sample_subcircuit(fmpq_mat_t mat, int row_idx) {
  int i = 0;
  uint32_t rand = 0;

  // constant term
  fmpq_set_si(fmpq_mat_entry(mat, row_idx, fmpq_mat_ncols(mat) - 1), 1, 1);

  // set all inputs
  for (auto& g : sc_inputs) {
    if (i++ % 32 == 0)
      rand = uniform(generator);
    
    Var* v = g->get_var();
    int val = (int)(rand & 1U);
    rand >>= 1U;
    
    v->set_value(val);
    v->get_dual()->set_value(1 - val);
    
    fmpq_set_si(fmpq_mat_entry(mat, row_idx, var_to_col[v]), val, 1);
  }

  // compute outputs
  for (auto& gate : gate_poly) {
    Polynomial* g = gate->get_aig_poly();

    int val = g->evaluate();
    Var* v = gate->get_var();

    v->set_value(val);
    v->get_dual()->set_value(1 - val);

    fmpq_set_si(fmpq_mat_entry(mat, row_idx, var_to_col[v]), val, 1);
  }
}
/*------------------------------------------------------------------------*/
static void sample_trivial(fmpq_mat_t mat) {
  for(int val = 0; val < 2; val++) {

    // constant term
    fmpq_set_si(fmpq_mat_entry(mat, val, fmpq_mat_ncols(mat) - 1), 1, 1);
    
    // set all inputs
    for (auto& g : sc_inputs) {
      Var* v = g->get_var();
      v->set_value(val);
      v->get_dual()->set_value(1 - val);
      fmpq_set_si(fmpq_mat_entry(mat, val, var_to_col[v]), val, 1);
    }

    // compute outputs
    for (auto& gate : gate_poly) {
      Polynomial* g = gate->get_aig_poly();

      int val_g = g->evaluate();
      Var* v = gate->get_var();

      v->set_value(val_g);
      v->get_dual()->set_value(1 - val_g);

      fmpq_set_si(fmpq_mat_entry(mat, val, var_to_col[v]), val_g, 1);
    }
  }
}

/*------------------------------------------------------------------------*/
static void sample_dual(fmpq_mat_t mat, int row_idx) {

  // constant term
  fmpq_set_si(fmpq_mat_entry(mat, row_idx,  fmpq_mat_ncols(mat) - 1), 1, 1);
  
  for (auto& g : sc_inputs) {
    Var* v = g->get_var();
    int val = v->get_value();
    v->set_value(val);
    v->get_dual()->set_value(1 - val);
    
    fmpq_set_si(fmpq_mat_entry(mat, row_idx, var_to_col[v]), val, 1);
  }

  // compute outputs
  for (auto& gate : gate_poly) {
    Polynomial* g = gate->get_aig_poly();

    int val = g->evaluate();
    Var* v = gate->get_var();

    v->set_value(val);
    v->get_dual()->set_value(1 - val);

    fmpq_set_si(fmpq_mat_entry(mat, row_idx, var_to_col[v]), val, 1);
  }
}
/*------------------------------------------------------------------------*/
Polynomial *
verify_guess(Polynomial* p, std::set<Polynomial*>& gb, std::vector<std::vector<int>> aig_clauses, int& eval_count, int& sat_count, std::map<Gate*, int> lit_id, std::map<int, Gate*> inverse_lit_id) {
  evaluated_guess_count++;
  eval_count++;
  if (use_algebra_reduction) { // use ideal membership
    if (reduce_to_zero(p, gb)) {
      if (verbose > 3) std::cout << "===== CORRECT =====\n";
      if (verbose > 3) p->print(stdout);
      Gate* p_lt = gate(p->get_lt()->get_var_num());
      p_lt->set_nf(p->copy());
      return p;
    } else {
      sat_count++;
      if (verbose > 3) std::cout << "===== WRONG =====" << std::endl;
      if (verbose > 3) p->print(stdout);
      delete (p);
      return nullptr;
    }
  } else {  // use KISSAT

  
    auto cnf_clauses = translate_poly_to_cnf(p, lit_id, aig_clauses, 0);
    bool run1 = call_kissat(cnf_clauses, inverse_lit_id);

    bool run2 = false;
    if(run1){
      cnf_clauses = translate_poly_to_cnf(p, lit_id, aig_clauses, 1);
      run2 = call_kissat(cnf_clauses, inverse_lit_id);
    }


    if (run1 && run2) {
      correct_guess_count++;
      if(proof_logging){
        pac_add_circuit_poly(polys_file, p);
      }
      
      Gate* p_lt = gate(p->get_lt()->get_var_num());
      p_lt->set_nf(p->copy());
      p_lt->update_gate_poly(p->copy());



      if (verbose > 1) std::cout << "===== CORRECT =====" << std::endl;
      if (verbose > 1) p->print(stdout);
      return p;
    } else {
      sat_count++;
      if (verbose > 1) std::cout << "===== WRONG =====" << std::endl;
      if (verbose > 1) p->print(stdout);
      delete (p);
      return nullptr;
    }
  }
}
/*------------------------------------------------------------------------*/
void append_collected_assignments(fmpq_mat_t mat) {

  if(collected_assignments.size() == 0)
    return;

  // count number of nonzero rows of mat
  int n = 0;
  while(n < fmpq_mat_nrows(mat) && !row_is_zero(mat, n)) n++;

  fmpq_mat_t extended;
  fmpq_mat_init(extended, n+collected_assignments.size(), fmpq_mat_ncols(mat));

  // copy nonzero rows
  int i = 0;
  for(; i < n; i++)
    for(int j = 0; j < fmpq_mat_ncols(extended); j++)
      fmpq_set(fmpq_mat_entry(extended,i,j), fmpq_mat_entry(mat,i,j));

  // append samples
  for(; i < fmpq_mat_nrows(extended); i++) {
    auto sample = collected_assignments.front();
    collected_assignments.pop_front();

    // constant term
    fmpq_set_si(fmpq_mat_entry(extended, i, fmpq_mat_ncols(mat) - 1), 1, 1);
    
    for (auto& g : sc_inputs) {
      Var* v = g->get_var();
      int val = sample[g];
      fmpq_set_si(fmpq_mat_entry(extended, i, var_to_col[v]), val, 1);
    }

    // compute outputs
    for (auto& g : gate_poly) {
      Var* v = g->get_var();
      int val = sample[g];
      fmpq_set_si(fmpq_mat_entry(extended, i, var_to_col[v]), val, 1);
    }
  }

  fmpq_mat_swap(mat, extended);
  fmpq_mat_clear(extended);
}

/*------------------------------------------------------------------------*/
std::vector<Polynomial*>
guess_linear() {
  count_guess_call++;

  std::vector<Polynomial*> result;

  std::vector<Var*> vars;
  for (auto& g : sc_inputs)
    vars.push_back(g->get_var());

  for (auto& g : gate_poly)
    vars.push_back(g->get_var());

  // sort vars in decreasing order
  std::vector<Var*> vars_sorted(vars.begin(), vars.end());
  std::sort(vars_sorted.begin(), vars_sorted.end(), [](auto& v1, auto& v2) {
    return v1->get_level() > v2->get_level();
  });
  var_to_col.clear();
  for(int j = 0; j < vars_sorted.size(); j++)
    var_to_col[vars_sorted[j]] = j; 

  double pre_guess_time = process_time();

  // Need one additional column for constant term
  int n = vars.size() + 1;
  int N = std::min(10 * n, 10'000) + 2;

  fmpq_mat_t mat;
  fmpq_mat_t K;
  fmpq_mat_init(mat, N, n);
  fmpq_mat_init(K, 1, 1);
 
  sample_trivial(mat);
 
  for (int i = 2; i < N; i+=2) {
    sample_subcircuit(mat, i);
    sample_dual(mat, i+1);
  }

  std::vector<Term*> terms;
  for (const auto& v : vars_sorted) {
    Term* t = new_term(v, 0);
    terms.push_back(t);
  }
  terms.push_back(0);

  std::set<Polynomial*> gb;
  for (const auto& gate : gate_poly) {
    if(!gate->get_nf())
      gate->set_nf(gate->get_gate_constraint()->copy());
  }
  mpz_t c;
  mpz_init(c);
  bool found_root = false;
  guess_time += (process_time() - pre_guess_time);


  // Prove Part
  int iteration_count = 0;
  int eval_count = 0, sat_count = 0;

  // Initialize CNF Translation by generating mapping and translating aig part to cnf
  auto mappings = var_cnf_mapping(vars_sorted);
  auto lit_id = std::get<0>(mappings);
  auto inverse_lit_id = std::get<1>(mappings);
  auto aig_clauses = translate_aig_part_to_cnf(lit_id);


  while(!found_root) {
    eval_count = 0, sat_count = 0;
    iteration_count++; total_iterations_count++;
    pre_guess_time = process_time();
    int nr_assignments = collected_assignments.size();
    append_collected_assignments(mat);
    // msg("M dim = %li, %li  using %i collected assignments", fmpq_mat_nrows(mat), fmpq_mat_ncols(mat), nr_assignments);

    // clear result from previous iteration (if existent)
    for(auto& p : result)
      delete(p);      
    result.clear();
    
    fmpq_mat_clear(K);
    kernel(mat, K);

    if(fmpq_mat_nrows(K) == 0) {
      // msg("NO LINEAR POLIES IN SUBCIRCUIT");
      break;
    }
    
    if (fmpq_is_zero(fmpq_mat_entry(K,0,0))) {
      // msg("NO LINEAR POLY FOR ROOT");
      break;
    }
      
    int nr_lin_polies = fmpq_mat_nrows(K);
    // msg("done guessing - building %i polies", nr_lin_polies);
    guess_time += (process_time() - pre_guess_time);
    total_guesses_count += nr_lin_polies;
    if(nr_lin_polies>max_guesses_count) max_guesses_count = nr_lin_polies;

 
    result.reserve(nr_lin_polies);

    bool all_already_linear = true;
    for (int i = 0; i < nr_lin_polies; i++) {
      if (!normalize_row(K,i))
        continue;

      for (int j = 0; j < n; j++) {
        if (fmpq_is_zero(fmpq_mat_entry(K, i, j)))
          continue;        
        Term* t = terms[j];
        fmpz_get_mpz(c, fmpq_mat_entry_num(K, i, j));
        Monomial* m = new Monomial(c, t ? t->copy() : 0);
        push_mstack(m);
      }
      Polynomial* p = build_poly();

      // if we already have a linear polynomial, continue
      Polynomial * nf = gate(p->get_lt()->get_var_num())->get_nf();
      if(nf and nf->degree() <= 1) {
        delete(p);
        continue;
      }
      
      all_already_linear = false;
      // verify correctness
      double pre_proof_time = process_time();
      p = verify_guess(p, gb, aig_clauses, eval_count, sat_count, lit_id, inverse_lit_id);
      double after_proof_time = process_time();
      proof_time += (after_proof_time - pre_proof_time);
      if(p) {
        found_root = found_root or (i==0);
        result.push_back(p);
      } 
    }
     
    // msg("evaluated %i polies, %i remain sat", eval_count, sat_count);
    accuracy[iteration_count-1]+= (static_cast<double>(eval_count-sat_count)/eval_count*100);
    iteration_on_level[iteration_count-1]+=1;
    
    if(all_already_linear) break;
  }
  if(iteration_count > max_iterations_count) max_iterations_count = iteration_count;
  
  collected_assignments.clear();
  mpz_clear(c);
  fmpq_mat_clear(K);
  fmpq_mat_clear(mat);

  return result;
}

/*------------------------------------------------------------------------*/
