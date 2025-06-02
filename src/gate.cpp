/*------------------------------------------------------------------------*/
/*! \file gate.cpp
    \brief contains the class Gate and further functions to
    organize the gate structure, such as initializing the gate constraints

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#include "gate.h"

#include <list>
#include <string>
/*------------------------------------------------------------------------*/
// Global variables
int add_var = 0;
int max_dist = 0;
/*------------------------------------------------------------------------*/
// ERROR CODES:
static int err_allocate = 91;  // failed to allocate gates

/*------------------------------------------------------------------------*/
Gate::Gate(int n_, std::string name_, int level_, bool input_, bool output_)
    : v(new Var(name_, level_, n_, 0)), input(input_), output(output_) {
  std::string negname = name_;
  if (proof_logging && negname.size() > 0) {
    negname.insert(1, 1, '_');
  } else if (negname.size() > 0) {
    negname.insert(0, "(1-");
    negname.insert(negname.size(), ")");
  }

  Var *d = new Var(negname, level_ + 1, n_, 1);
  v->set_dual_var(d);
  d->set_dual_var(v);
}
/*------------------------------------------------------------------------*/
void Gate::set_elim() {
  if (elim) return;

  for (auto &gc : get_children()) {
    gc->parents_remove(this);
  }
  elim = 1;
  if (verbose > 3)
    msg("eliminated %s", get_var_name());
  if (gate_constraint) {
    delete (gate_constraint);
    gate_constraint = 0;
  }
}
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
Gate *xor_left_child(const Gate *n) {
  if (!n->get_xor_gate()) return 0;

  aiger_and *and1 = is_model_and(n->get_var_num());
  if (!and1) return 0;
  unsigned l = and1->rhs0;
  if (!aiger_sign(l)) return 0;
  l = aiger_strip(l);

  aiger_and *land = is_model_and(l);
  if (!land) return 0;

  unsigned ll = land->rhs0;
  return gate(ll);
}

/*------------------------------------------------------------------------*/
bool upper_half_xor_output() {
  for (unsigned i = num_gates - 2; i > M - 1; i--) {
    Gate *n = gates[i];
    if (!n->children_size()) return 0;
    n = n->children_front();
    if (!n->get_xor_gate()) return 0;
  }

  if (!gates[M - 1]->children_size()) return 0;

  return 1;
}

/*------------------------------------------------------------------------*/

Gate *xor_right_child(const Gate *n) {
  if (!n->get_xor_gate()) return 0;

  aiger_and *and1 = is_model_and(n->get_var_num());
  if (!and1) return 0;
  unsigned l = and1->rhs0;
  if (!aiger_sign(l)) return 0;
  l = aiger_strip(l);

  aiger_and *land = is_model_and(l);
  if (!land) return 0;

  unsigned lr = land->rhs1;
  return gate(lr);
}

/*------------------------------------------------------------------------*/

static Polynomial *positive_poly(Var *v) {
  Term *t = new_term(v);
  Monomial *m = new Monomial(one, t);
  push_mstack_end(m);

  return build_poly();
}

/*------------------------------------------------------------------------*/
static Polynomial *negative_poly(Var *v) {
  Term *t1 = new_term(v);
  Monomial *m1 = new Monomial(minus_one, t1);
  push_mstack(m1);

  Monomial *m2 = new Monomial(one, 0);
  push_mstack(m2);

  return build_poly();
}

/*------------------------------------------------------------------------*/

static Polynomial *get_node_constraint(Gate *g, unsigned sign, bool flipped = 1) {
  if (g) {
    Var *v1 = g->get_var();
    if (sign && flipped)
      return positive_poly(v1->get_dual());
    else if (sign)
      return negative_poly(v1);
    else
      return positive_poly(v1);
  } else {
    if (sign) {
      push_mstack_end(new Monomial(one, 0));
      return build_poly();
    } else
      return 0;
  }
}

/*------------------------------------------------------------------------*/
Polynomial *gen_gate_constraint(unsigned i) {
  assert(i >= NN && i < M + MM - 1);
  Polynomial *p = 0;
  Gate *n = gates[i];
  // gate constraint
  if (i < M - 1 && (n->get_xor_gate() != 1 || !do_local_lin)) {
    assert(!n->get_input());

    aiger_and *and1 = is_model_and(n->get_var_num());
    assert(and1);

    unsigned l = and1->rhs0, r = and1->rhs1;
    Gate *l_gate = gate(l), *r_gate = gate(r);

    Var *v = n->get_var();
    Term *t1 = new_term(v);
    Monomial *m1 = new Monomial(minus_one, t1);

    push_mstack_end(m1);

    Var *v1 = aiger_sign(and1->rhs0) ? l_gate->get_var()->get_dual() : l_gate->get_var();
    Var *v2 = aiger_sign(and1->rhs1) ? r_gate->get_var()->get_dual() : r_gate->get_var();
    Monomial *m2 = new Monomial(one, new_quadratic_term(v1, v2));
    push_mstack_end(m2);

    p = build_poly();
  } else if (i < M - 1) {
    p = gen_xor_constraint(n);
    lin_xor_constraint_count++;

  }

  else {  // output

    assert(n->get_output());

    unsigned lit = slit(i - M + 1);
    Var *v = n->get_var();
    Term *t1 = new_term(v);
    Monomial *m1 = new Monomial(minus_one, t1);
    push_mstack_end(m1);

    if (lit == 1) {
      push_mstack_end(new Monomial(one, 0));

    } else if (aiger_sign(slit(i - M + 1))) {
      Term *t1 = new_term(gate(slit(i - M + 1))->get_var());
      Monomial *m2 = new Monomial(minus_one, t1);
      push_mstack_end(m2);
      Monomial *m3 = new Monomial(one, 0);
      push_mstack_end(m3);
    } else {
      Term *t = new_term(gate(slit(i - M + 1))->get_var());
      Monomial *m = new Monomial(one, t);
      push_mstack_end(m);
    }

    p = build_poly();
  }
  // p->set_idx(2+i-NN);
  return p;
}

//-------------------------------------------------------------
Polynomial *gen_xor_constraint(Gate *n) {
  assert(!n->get_input());

  aiger_and *and1 = is_model_and(n->get_var_num());
  assert(and1);

  unsigned l = and1->rhs0, r = and1->rhs1;
  Gate *l_gate = gate(l), *r_gate = gate(r);

  Gate *smaller = l_gate->get_var_level() < r_gate->get_var_level() ? l_gate : r_gate;
  aiger_and *smand1 = is_model_and(smaller->get_var_num());
  assert(smand1);
  unsigned ll = smand1->rhs0, rr = smand1->rhs1;
  Gate *ll_gate = gate(ll), *rr_gate = gate(rr);

  Var *v = n->get_var();
  Term *t1 = new_term(v);
  Monomial *m1 = new Monomial(minus_one, t1);
  push_mstack_end(m1);

  Term *t2 = new_term(smaller->get_var());
  Monomial *m2 = new Monomial(minus_two, t2);
  push_mstack_end(m2);

  Polynomial *p_h = build_poly();
  Polynomial *p1 = get_node_constraint(ll_gate, aiger_sign(smand1->rhs0), 0);
  Polynomial *p2 = get_node_constraint(rr_gate, aiger_sign(smand1->rhs1), 0);
  Polynomial *p_tl = add_poly(p1, p2);

  Polynomial *p = add_poly(p_h, p_tl);
  delete (p_h);
  delete (p_tl);

  delete (p1);
  delete (p2);

  return p;
}
/*------------------------------------------------------------------------*/

static void init_gate_constraint(unsigned i) {
  assert(i >= NN && i < M + MM - 1);
  Gate *n = gates[i];

  Polynomial *p = gen_gate_constraint(i);

  n->set_gate_constraint(p);
  n->set_aig_poly(p->copy());
  
}
/*------------------------------------------------------------------------*/

Polynomial *Gate::get_gate_constraint() const {
  if (!gate_constraint) {
    // output aig are 0, -1, ...-NN+2
    if (output)
      init_gate_constraint(-1 * get_var_num() + M - 1);
    // gates are numbered 2,4,6,8,..
    else
      init_gate_constraint(get_var_num() / 2 - 1);
  }
  return gate_constraint;
}

/*------------------------------------------------------------------------*/
Polynomial *Gate::get_dual_constraint() {
  if (!dual_constraint) {
    Var *v = get_var();
    if (!v->is_dual()) v = v->get_dual();
    dual_constraint = gen_dual_constraint(v);
  }
  return dual_constraint;
}

/*------------------------------------------------------------------------*/
std::list<Gate *> get_var_of_poly(Polynomial *p, bool tail) {
  std::list<Gate *> res;
  size_t j = tail ? 1 : 0;
  for (size_t i = j; i < p->len(); i++) {
    Monomial *m = p->get_mon(i);
    Term *t = m->get_term();
    while (t) {
      Gate *tmp = gate(t->get_var_num());
      auto it = std::find(res.begin(), res.end(), tmp);
      if (it == res.end()) {
        res.push_back(tmp);
      }
      t = t->get_rest();
    }
  }

  return res;
}

/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/

void Gate::update_gate_poly(Polynomial *p, bool rec) {
  if (gate_constraint)
    delete (this->gate_constraint);

  this->set_gate_constraint(p);

  auto g_orig_children = this->get_children();
  for (auto &gc : g_orig_children) {
    gc->parents_remove(this);
  }

  auto g_update_children = get_var_of_poly(this->get_gate_constraint());
  this->set_children(g_update_children);

  for (auto &gp : g_update_children) {
    gp->parents_push_back(this);
  }
  
}
/*------------------------------------------------------------------------*/
bool Gate::is_dual_twin(const Gate *n) const {
  for (auto it = dual_twins.begin(); it != dual_twins.end(); ++it) {
    Gate *g = *it;
    if (g == n)
      return 1;
  }
  return 0;
}
/*------------------------------------------------------------------------*/
bool Gate::is_van_twin(const Gate *n) const {
  for (auto it = van_twins.begin(); it != van_twins.end(); ++it) {
    Gate *g = *it;
    if (g == n)
      return 1;
  }
  return 0;
}
/*------------------------------------------------------------------------*/
bool Gate::is_child(const Gate *n) const {
  for (auto it = children.begin(); it != children.end(); ++it) {
    Gate *child = *it;
    if (child == n)
      return 1;
  }
  return 0;
}
/*------------------------------------------------------------------------*/
bool Gate::is_aig_child(const Gate *n) const {
  for (auto it = aig_children.begin(); it != aig_children.end(); ++it) {
    Gate *child = *it;
    if (child == n)
      return 1;
  }
  return 0;
}
/*------------------------------------------------------------------------*/
bool Gate::is_in_parents(const Gate *n) const {
  for (auto it = parents.begin(); it != parents.end(); ++it) {
    Gate *parents = *it;
    if (parents == n)
      return 1;
  }
  return 0;
}
/*------------------------------------------------------------------------*/
bool Gate::is_in_aig_parents(unsigned n) const {
  for (auto it = aig_parents.begin(); it != aig_parents.end(); ++it) {
    unsigned parents = *it;
    if (parents == n)
      return 1;
  }
  return 0;
}

/*------------------------------------------------------------------------*/
bool Gate::is_in_neg_parents(unsigned n) const {
  for (auto it = neg_parents.begin(); it != neg_parents.end(); ++it) {
    unsigned parents = *it;
    if (parents == n)
      return 1;
  }
  return 0;
}
/*------------------------------------------------------------------------*/
bool Gate::is_in_pos_parents(unsigned n) const {
  for (auto it = pos_parents.begin(); it != pos_parents.end(); ++it) {
    unsigned parents = *it;
    if (parents == n)
      return 1;
  }
  return 0;
}
/*------------------------------------------------------------------------*/
Gate::~Gate() {
  delete (v->get_dual());
  delete (v);
  if (gate_constraint)
    delete (gate_constraint);
  if (normal_form)
    delete (normal_form);
  if (dual_constraint)
    delete (dual_constraint);
  if (aig_poly)
    delete (aig_poly);
}

/*------------------------------------------------------------------------*/
Gate **gates;
unsigned num_gates = 0;
unsigned size_gates = 0;

/*------------------------------------------------------------------------*/

static void mark_aig_outputs() {
  for (unsigned i = 0; i < MM; i++) {
    unsigned lit = slit(i);
    if (lit < 2)
      continue;
    Gate *n = gate(lit);
    n->mark_aig_output();
  }
}
/*------------------------------------------------------------------------*/
static void allocate_gates() {
  unsigned aiger;
  num_gates = M + MM - 1;

  msg("allocating %i gates", num_gates);
  gates = new Gate *[num_gates + MM];
  size_gates = num_gates + MM;

  if (!gates)
    die(err_allocate, "failed to allocate gates");
  int level = 0;
  if (mult_spec) {
    // inputs a
    for (unsigned i = a0; i <= al; i += ainc) {
      aiger = 2 * (i + 1);
      assert(is_model_input(aiger));

      std::string name = std::string("a") + std::to_string((i - a0) / ainc);
      gates[i] = new Gate(aiger, name, level += 2, 1);
    }

    // inputs b
    for (unsigned i = b0; i <= bl; i += binc) {
      aiger = 2 * (i + 1);
      assert(is_model_input(aiger));

      std::string name = std::string("b") + std::to_string((i - b0) / binc);
      gates[i] = new Gate(aiger, name, level += 2, 1);
    }
  } else {
    for (unsigned i = 0; i < NN; i++) {
      aiger = 2 * (i + 1);
      assert(is_model_input(aiger));

      std::string name = std::string("i") + std::to_string(i);
      gates[i] = new Gate(aiger, name, level += 2, 1);
      if (verbose > 3) msg("allocated inp %s", gates[i]->get_var_name());
    }
  }

  // internal gates
  for (unsigned i = NN; i < M - 1; i++) {
    aiger = 2 * (i + 1);
    assert(is_model_and(aiger));

    std::string name = std::string("l") + std::to_string(aiger);
    gates[i] = new Gate(aiger, name, 0);
    if (verbose > 3) msg("allocated gate %s", gates[i]->get_var_name());
  }

  for (unsigned i = NN; i < M - 1; i++) {
    Gate *n = gates[i];

    assert(is_model_and(n->get_var_num()));
    aiger_and *and1 = is_model_and(n->get_var_num());

    if (!and1)
      continue;
    unsigned l = and1->rhs0, r = and1->rhs1;

    if (l < 2 || r < 2) {
      n->set_dist(1);
      if (verbose > 3) msg("gate %s has distance %i", n->get_var_name(), 1);
    } else {
      Gate *l_gate = gate(l), *r_gate = gate(r);
      int dist_l = l_gate->get_dist();
      int dist_r = r_gate->get_dist();
      int dist_n = dist_l > dist_r ? dist_l + 1 : dist_r + 1;

      n->set_dist(dist_n);
      if (verbose > 3) msg("gate %s has distance %i", n->get_var_name(), dist_n);
      if (dist_n > max_dist)
        max_dist++;
    }
  }
  msg("max dist is %i", max_dist);

  mark_aig_outputs();

  for (int dist_it = 1; dist_it <= max_dist; dist_it++) {
    for (unsigned i = NN; i < M - 1; i++) {
      Gate *n = gates[i];
      if (n->get_dist() == dist_it) {
        level += 2;
        n->set_var_level(level);
      }
    }
  }

  // output s
  for (unsigned i = M - 1; i < M - 1 + MM; i++) {
    aiger = i - M + 1;
    std::string name = std::string("s") + std::to_string(aiger);
    gates[i] = new Gate(M - i - 1, name, 2 * (i + 1), 0, 1);
    if (verbose > 3) msg("allocated outp %s", gates[i]->get_var_name());
  }
}

/*------------------------------------------------------------------------*/

static void init_gate_constraints() {
  for (unsigned i = NN; i < M - 1; i++) {
    init_gate_constraint(i);
  }

  for (unsigned i = 0; i < MM; i++) {
    init_gate_constraint(i + M - 1);
  }
}

/*------------------------------------------------------------------------*/
static void set_parents_and_children() {
  unsigned pp = 0;

  // for (unsigned i = NN; i < M; i++) {
  for (unsigned i = M - 1; i >= NN; i--) {
    Gate *n = gates[i];
    assert(!n->get_input());

    aiger_and *and1 = is_model_and(n->get_var_num());

    if (!and1)
      continue;
    unsigned l = and1->rhs0, r = and1->rhs1;
    Gate *l_gate = gate(l), *r_gate = gate(r);
    n->children_push_back(l_gate);
    n->children_push_back(r_gate);
    n->aig_children_push_back(l_gate);
    n->aig_children_push_back(r_gate);
    if (verbose >= 4)
      msg("node %s has children %s, %s", n->get_var_name(), l_gate->get_var_name(), r_gate->get_var_name());

    if (l_gate && r_gate) {
      if (l_gate->get_input() && r_gate->get_input() && !aiger_sign(l) && !aiger_sign(r)) {
        n->mark_pp();
        pp++;
        if (verbose >= 4)
          msg("partial product %s", n->get_var_name());
      }
    }
    if (l_gate) {
      l_gate->parents_push_back(n);
      if (aiger_sign(l)) {
        l_gate->aig_parents_push_back(n->get_var_num() + 1);
        l_gate->neg_parents_push_back(n->get_var_num() + 1);
      } else {
        l_gate->aig_parents_push_back(n->get_var_num());
        if (n->neg_parents_size() > 0) {
          l_gate->pos_parents_push_back(n->get_var_num());
          //       msg("pushed back1 %s, %s", l_gate->get_var_name(), n->get_var_name());
        }
        for (unsigned n_pos : n->get_pos_parents()) {
          l_gate->pos_parents_push_back(n_pos);
          //      msg("pushed back2 %s, %s", l_gate->get_var_name(), gate(n_pos)->get_var_name());
        }
      }
    }

    if (r_gate) {
      r_gate->parents_push_back(n);
      if (aiger_sign(r)) {
        r_gate->aig_parents_push_back(n->get_var_num() + 1);
        r_gate->neg_parents_push_back(n->get_var_num() + 1);
      } else {
        r_gate->aig_parents_push_back(n->get_var_num());
        if (n->neg_parents_size() > 0) {
          r_gate->pos_parents_push_back(n->get_var_num());
          //      msg("pushed back3 %s, %s", r_gate->get_var_name(), n->get_var_name());
        }
        for (unsigned n_pos : n->get_pos_parents()) {
          r_gate->pos_parents_push_back(n_pos);
          //         msg("pushed back4 %s, %s", l_gate->get_var_name(), gate(n_pos)->get_var_name());
        }
      }
    }
  }

  // set children for extra outputs
  for (unsigned i = 0; i < MM; i++) {
    Gate *n = gates[i + M - 1];
    assert(n->get_output());
    unsigned lit = slit(i);
    if (lit < 2)
      continue;
    Gate *model_output_gate = gate(lit);
    n->children_push_back(model_output_gate);
    if (verbose >= 4)
      msg("node %s has child %s", n->get_var_name(), model_output_gate->get_var_name());
    model_output_gate->parents_push_back(n);
  }

  if (verbose >= 1)
    msg("found %i partial products", pp);

  if(pp != NN/2*NN/2)   booth = 1;
}
/*------------------------------------------------------------------------*/


static void mark_xor_and() {
  for (unsigned i = 0; i < M; i++) {
    Gate *g = gates[i];

    if (g->get_xor_gate() != 1) continue;
    if (g->children_size() != 2) continue;

    //if (g->get_aig_output()) continue;

    Gate *llg = g->children_front()->children_front();
    Gate *lrg = g->children_front()->children_back();

    std::vector<Gate *> ands;
    for (auto &llg_p : llg->get_parents()) {
      if (g->is_child(llg_p)) continue;     // check whether internal xor node
      if (!llg_p->is_child(lrg)) continue;  // check whether both llg and lrg are parents

      ands.push_back(llg_p);
    }
    if (ands.size() == 0) continue;

    if (ands.size() == 1) {
      Gate *and1 = *ands.begin();
      and1->set_xor_and(g);
      g->set_xor_and(and1);    
      llg->mark_xor_and_inp();
      lrg->mark_xor_and_inp();  
    }
  }
}
/*------------------------------------------------------------------------*/
void set_xor() {
  int found_xor = 0;
  for (unsigned i = 0; i < M; i++) {
    Gate *n = gates[i];
    if (n->get_input() > 0) continue;
    if (n->get_xor_gate() > 0) continue;

    aiger_and *and1 = is_model_and(n->get_var_num());
    if (!and1) continue;
    unsigned l = and1->rhs0, r = and1->rhs1;
    if (!aiger_sign(l)) continue;
    if (!aiger_sign(r)) continue;
    if (l == r || l == aiger_not(r)) continue;
    l = aiger_strip(l);
    r = aiger_strip(r);
    aiger_and *land = is_model_and(l);
    if (!land) continue;
    aiger_and *rand = is_model_and(r);
    if (!rand) continue;

    unsigned ll = land->rhs0, lr = land->rhs1;
    unsigned rl = rand->rhs0, rr = rand->rhs1;
    if ((ll == aiger_not(rl) && lr == aiger_not(rr)) ||
        (ll == aiger_not(rr) && lr == aiger_not(rl))) {
      gate(l)->set_xor_gate(2);
      gate(r)->set_xor_gate(2);
      n->set_xor_gate(1);
      found_xor++;

      if (verbose >= 4) msg("xor-gate %s", n->get_var_name());
    }
  }
  if (verbose >= 1) msg("found %i xor-gates", found_xor);
}

/*------------------------------------------------------------------------*/

void init_gates() {
  allocate_gates();
  set_parents_and_children();
  init_gate_constraints();
  set_xor();
  mark_xor_and();
}
/*------------------------------------------------------------------------*/
void enlarge_gates(int added_size) {
  uint64_t new_size_gates = size_gates + added_size;
  Gate **new_gates_table = new Gate *[new_size_gates];
  if (size_gates > 0) memcpy(new_gates_table, gates, size_gates * sizeof(Gate *));
  delete[] gates;
  gates = new_gates_table;
  size_gates = new_size_gates;
}
/*------------------------------------------------------------------------*/

Gate *gate(int lit) {
  if (lit <= 0)
    return gates[M - lit - 1];
  if (lit < 2)
    return 0;
  return gates[lit / 2 - 1];
}

/*------------------------------------------------------------------------*/

void delete_gates() {
  for (unsigned i = 0; i < num_gates; i++) {
    if (gates[i]) msg("delete %s", gates[i]->get_var_name());
    if (gates[i]) delete (gates[i]);
  }
  delete[] gates;
}
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
Gate *search_for_parent(Term *t, Gate *exclude) {
  assert(t);

  Gate *g = gate(t->get_var_num());
  for (auto &parent : g->get_parents()) {
    if (parent == exclude)
      continue;
    Polynomial *gc = parent->get_gate_constraint();
    if (gc->len() != 2)
      continue;

    if (t == gc->get_tail_term())
      return parent;
  }
  return 0;
}
/*------------------------------------------------------------------------*/
Gate *search_for_parent_dual(Term *t) {
  assert(t);
  Gate *g = gate(t->get_var_num());

  for (auto &parent : g->get_parents()) {
    Polynomial *gc = parent->get_gate_constraint();

    if (gc->len() != 2)
      continue;
    if (equal_up_to_duality(t, gc->get_tail_term()))
      return parent;
  }
  return 0;
}
/*------------------------------------------------------------------------*/
bool equal_children(const Gate *g1, const Gate *g2) {
  if (g1->children_size() != g2->children_size())
    return 0;

  for (auto &child : g1->get_children()) {
    if (!g2->is_child(child))
      return 0;
  }
  return 1;
}
