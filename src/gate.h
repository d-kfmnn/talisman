/*------------------------------------------------------------------------*/
/*! \file gate.h
    \brief contains the class Gate and further functions to
    organize the gate structure, such as initializing the gate constraints

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_GATE_H_
#define TALISMAN_SRC_GATE_H_
/*------------------------------------------------------------------------*/
#include <list>
#include <map>
#include <queue>
#include <string>
#include <map>

#include "aig.h"
#include "polynomial.h"
/*------------------------------------------------------------------------*/
// / set to true when a signed or unsigned multiplier is verified
extern int add_var;
extern int max_dist;
/*------------------------------------------------------------------------*/
extern unsigned size_gates;
extern std::map<Polynomial*, size_t> van_constr;


/** \class Gate
  Internal structure to represent the AIG graph.
*/
class Gate
{
  // / Variable of the gate, as used in the polynomials
  Var *v;

  Var *replace_var = 0;

  // / True if gate is an input
  bool input;

  // / True if the gate is an output s_i
  bool output;

  // / True if the gate is an output in the aig
  bool aig_output = 0;

  // / True if gate is identified as a partial product
  bool partial_product = 0;

  // / is set to 1 for root node, 2 for internal nodes of XORs
  int xor_gate = 0;
  Gate * xor_and = 0;

  int distance = 0;

   // / True if circuit is a prop_gen_gate(-substitute)
   bool prop_gen_gate = 0;

   // / True if gate is identified to belong to complex fsa(-substitute)
   bool fsa = 0;
 
   // / True if gate is input of complex fsa(-substitute)
   int fsa_inp = 0;
 
   // / True if gate occurs negative(-substitute)
   bool neg = 0;

  // / True if gate is eliminated during preprocessing
  bool elim = 0;

  bool xor_and_inp = 0;

  bool extension = 0;

  // / Polynomial implied by the aig gate
  Polynomial *gate_constraint = 0;
  Polynomial * dual_constraint = 0;
  Polynomial * normal_form = 0;
  Polynomial * aig_poly = 0;

  // / list of gates that create vanishing monomials
  std::vector<Gate *> van_twins;

    // / list of gates that create vanishing monomials
  std::vector<Gate *> dual_twins;

  // / list of gates that are parents
  std::list<Gate *> parents;

  std::vector<unsigned> aig_parents;
  
  std::vector<unsigned> pos_parents;
  std::vector<unsigned> neg_parents;

  // / list of gates that are children
  std::list<Gate *> children;
  std::list<Gate*> aig_children;

public:
  /**
      Constructor
      Calls constructor of Var

      @param n_ value, corresponding to aiger value
      @param name_ string name of the variable
      @param level position in order of the variable
  */
  Gate(int n_, std::string name_, int level_, bool input_ = 0, bool output_ = 0);

  // Var
  Var *get_var() const
  {
    return v;
  }
  int get_var_num() const
  {
    return v->get_num();
  }
  int get_var_level() const
  {
    return v->get_level();
  }
  const char *get_var_name() const
  {
    return v->get_name();
  }

  void set_var_level(int l)
  {
    v->set_level(l);
    v->get_dual()->set_level(l + 1);
  }

  // Replace_var
  Var *get_rep_var() const
  {
    return replace_var;
  }
  void set_rep_var(Var *v)
  {
    replace_var = v;
  }

  // Input & output
  bool get_input() const
  {
    return input;
  }
  bool get_output() const
  {
    return output;
  }



    /**
      Getter for prop_gen_gate

      @return member prop_gen_gate
  */
 bool get_prop_gen_gate() const {return prop_gen_gate;}

 /**
     Sets prop_gen_gate to true
 */
 void mark_prop_gen_gate() {prop_gen_gate = 1;}

 /**
     Sets prop_gen_gate to false
 */
 void unmark_prop_gen_gate() {prop_gen_gate = 0;}

 /**
     Getter for fsa

     @return member fsa
 */
 bool get_fsa() const {return fsa;}

 /**
     Sets fsa to true
 */
 void mark_fsa() {fsa = 1;}
 void remove_fsa() {fsa = 0;}



bool get_xor_and_inp() const {return xor_and_inp;}


void mark_xor_and_inp() {xor_and_inp = 1;}

 void set_ext() {extension = 1;}
 bool is_extension() {return extension;}

 /**
     Getter for fsa_inp

     @return member fsa_inp
 */
 int get_fsa_inp() const {return fsa_inp;}

 /**
     Increases fsa_inp
 */
 void inc_fsa_inp() {fsa_inp++;}

 /**
     Sets fsa_inp to 0
 */
 void reset_fsa_inp() {fsa_inp = 0;}

 /**
     Getter for neg

     @return member neg
 */
 bool get_neg() const {return neg;}

 /**
     Setter for neg

     @param val Boolean
 */
 void set_neg(bool val) {neg = val;}

  // AIG output
  bool get_aig_output() const
  {
    return aig_output;
  }
  void mark_aig_output()
  {
    aig_output = 1;
  }

  // Partial products
  bool get_pp() const
  {
    return partial_product;
  }
  void mark_pp()
  {
    partial_product = 1;
  }

  // xor gate
  int get_xor_gate() const
  {
    return xor_gate;
  }

  void set_xor_gate(int val)
  {
    xor_gate = val;
  }

    // xor gate
  Gate* get_xor_and_gate() const
  {
    return xor_and;
  }

  void set_xor_and(Gate *g)
  {
    xor_and = g;
  }

  // Distance to inputs
  int get_dist() const
  {
    return distance;
  }
  void set_dist(int l)
  {
    distance = l;
  }

  // Red
  bool get_elim() const
  {
    return elim;
  }
  void set_elim();

  // Gate constraint
  Polynomial *get_gate_constraint() const;
  Polynomial * get_dual_constraint();
  void update_gate_poly(Polynomial *p, bool rec = 1);
  void set_gate_constraint(Polynomial *p)
  {
    gate_constraint = p;
  }
  void print_gate_constraint(FILE *file) const
  {
    gate_constraint->print(file);
  }

    // Normal form 
  Polynomial *get_nf() const {
    return normal_form;
  }

  void set_nf(Polynomial *p)
  {
    normal_form = p;
  }

  Polynomial *get_aig_poly() const {
    return aig_poly;
  }

  void set_aig_poly(Polynomial *p)
  {
    aig_poly = p;
  }
  
  void print_nf(FILE *file) const
  {
    normal_form->print(file);
  }

  // Vanishing twins
  bool is_van_twin(const Gate *n) const;
  size_t van_twins_size() const
  {
    return van_twins.size();
  }
  void van_twins_push_back(Gate *n)
  {
   
      this->van_twins.push_back(n);
      //msg("inserted %s to van twins of %s", n->get_var_name(), this->get_var_name());
  
  }
 std::vector<Gate *> get_van_twins() const
  {
    return van_twins;
  }

    // DUal twins
  bool is_dual_twin(const Gate *n) const;
  size_t dual_twins_size() const
  {
    return dual_twins.size();
  }
  void dual_twins_push_back(Gate *n)
  {
    if(!this->is_dual_twin(n)) {
      this->dual_twins.push_back(n);
      //msg("inserted %s to van twins of %s", n->get_var_name(), this->get_var_name());
    }
  }
 std::vector<Gate *> get_dual_twins() const
  {
    return dual_twins;
  }

  // Children
  std::list<Gate *> get_children() const
  {
    return children;
  }
  size_t children_size() const
  {
    return children.size();
  }

  bool is_child(const Gate *n) const;
 

  void set_children(std::list<Gate *> c)
  {
    children = c;
  }
  void children_push_back(Gate *n)
  {
    if (n)
      children.push_back(n);
  }

  void children_remove(Gate *n)
  {
    children.remove(n);
  }
  void delete_children()
  {
    children = std::list<Gate *>();
  }

  Gate *children_front() const
  {
    return children.front();
  }
  Gate *children_back() const
  {
    return children.back();
  }

  // Parents
  std::list<Gate *> get_parents()
  {
    return parents;
  }
  size_t parents_size() const
  {
    return parents.size();
  }

  bool is_in_parents(const Gate *n) const;
  void parents_push_back(Gate *n)
  {
    parents.push_back(n);
  }
  void parents_remove(Gate *n)
  {
    parents.remove(n);
  }


  

  std::list<Gate *>::const_iterator parents_begin() const
  {
    return parents.begin();
  }
  std::list<Gate *>::const_iterator parents_end() const
  {
    return parents.end();
  }

  Gate *parents_front() const
  {
    return parents.front();
  }
  Gate *parents_back() const
  {
    return parents.back();
  }



  // aig Parents
  std::vector<unsigned> get_aig_parents()
  {
    return aig_parents;
  }
  size_t aig_parents_size() const
  {
    return aig_parents.size();
  }
   unsigned aig_parents_front() const
  {
    return aig_parents.front();
  }

  bool is_in_aig_parents(unsigned n) const;
  void aig_parents_push_back(unsigned n)
  {
    aig_parents.push_back(n);
  }

  // aig Children
std::list<Gate*> get_aig_children()
{
    return aig_children;
}

size_t aig_children_size() const
{
    return aig_children.size();
}

Gate* aig_children_front() const
{
    return aig_children.front();
}
bool is_aig_child(const Gate* n) const;


void aig_children_push_back(Gate* n)
{
    aig_children.push_back(n);
}


  // pos Parents
  std::vector<unsigned> get_pos_parents()
  {
    return pos_parents;
  }
  size_t pos_parents_size() const
  {
    return pos_parents.size();
  }
   unsigned pos_parents_front() const
  {
    return pos_parents.front();
  }

  bool is_in_pos_parents(unsigned n) const;
  void pos_parents_push_back(unsigned n)
  {
    pos_parents.push_back(n);
  }

 // neg Parents
  std::vector<unsigned> get_neg_parents()
  {
    return neg_parents;
  }
  size_t neg_parents_size() const
  {
    return neg_parents.size();
  }
   unsigned neg_parents_front() const
  {
    return neg_parents.front();
  }

  bool is_in_neg_parents(unsigned n) const;
  void neg_parents_push_back(unsigned n)
  {
    neg_parents.push_back(n);
  }

  ~Gate();
};

/*------------------------------------------------------------------------*/
// / Gate ** where all gates are stored
extern Gate **gates;

// / Counts the number of gates
extern unsigned num_gates;

/*------------------------------------------------------------------------*/
void init_gates();
void enlarge_gates(int added_size);

/**
    Returns the gate with aiger value 'lit'

    @param lit unsigned integer

    @returns Gate*
*/
Gate *gate(int lit);



Polynomial *gen_gate_constraint(unsigned i);
Polynomial *gen_xor_constraint(Gate * n);

/**
    Deletes Gate** gates by calling destructur of Gate
*/
void delete_gates();

/*------------------------------------------------------------------------*/

/**
    Identifies whether the output gates of slice N until NN-2 are
    XOR gates

    @return True when output gates of bigger slices are all xor gates
*/
bool upper_half_xor_output();

/**
    Returns the 'left' child of an XOR

    @param n Gate that is root of XOR

    @return Gate*
*/
Gate * xor_left_child(const Gate * n);

/**
    Returns the 'right' child of an XOR

    @param n Gate that is root of XOR

    @return Gate*
*/
Gate * xor_right_child(const Gate * n);

Gate *search_for_parent(Term *t, Gate *exclude = 0);
Gate *search_for_parent_dual(Term *t);
bool equal_children(const Gate *g1, const Gate *g2);

std::list<Gate *> get_var_of_poly(Polynomial * p, bool tail = 1);

struct LargerGate {
  bool operator()(const Gate* a, const Gate* b) const {
    return a->get_var_level() > b->get_var_level();
  }
};

struct SmallerGate {
  bool operator()(const Gate* a, const Gate* b) const {
    return a->get_var_level() < b->get_var_level();
  }
};

#endif // TALISMAN_SRC_GATE_H_
