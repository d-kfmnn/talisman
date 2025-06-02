/*------------------------------------------------------------------------*/
/*! \file subcircuit.h
    \brief contains functions to identify a subcircuit

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_SUBCIRCUIT_H_
#define TALISMAN_SRC_SUBCIRCUIT_H_
/*------------------------------------------------------------------------*/
#include <gmpxx.h>
#include <set>
#include <string.h>

#include <pblib/pb2cnf.h>
extern "C" {
  #include "kissat.h"  // Include Kissat header
}

#include "gate.h"
#include "pac.h"
#include "propagate.h"
/*------------------------------------------------------------------------*/
extern size_t fanout_size_last_call;
/*------------------------------------------------------------------------*/
extern std::set<Gate*, SmallerGate> gate_poly;
extern std::set<Gate*, LargerGate> var;
extern std::set<Gate*, LargerGate> sc_inputs;

struct Normalized_poly {
  std::vector<mpz_class> coeffs;
  std::vector<std::vector<size_t>> terms;

  Normalized_poly()
    : coeffs()
    , terms() {};

  void inline emplace_back(mpz_class& c, std::vector<size_t>& t) {
    coeffs.push_back(std::move(c));
    terms.push_back(std::move(t));
  }

  size_t inline size() const { return terms.size(); }

  void inline print() const {
    for(size_t i = 0; i < terms.size(); i++) {
      std::cout << coeffs[i] << ": (";
      for(auto v : terms[i])
        std::cout << v << ", ";
      std::cout << "), ";
    }
    std::cout << "\n";
  }
};

struct circuit_hash {
  size_t operator()(const std::vector<Normalized_poly>& circuit) const {
    size_t seed = 0;
    std::hash<std::string> h1;
    std::hash<size_t> h2;
    for(const auto& p : circuit) {
      for(size_t i = 0; i < p.size(); i++) {
        seed
          ^= h1(p.coeffs[i].get_str()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        for(auto t : p.terms[i])
          seed ^= h2(t) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      } 
      // return seed;
    }
    return seed;
  }
};

bool is_internal_fsa(Gate *g);

bool
get_and_compress_subcircuit(Gate* g,
                            int depth,
                            size_t fanout_size,
                            int init, bool single_expand,
                            std::vector<Normalized_poly>& normalized,
                            std::map<Var*, size_t>& var_to_id);

std::vector<Polynomial*>
compute_normalforms(std::vector<Polynomial*> *used_van_poly, std::vector<Polynomial*> *new_nf_poly);


std::vector<Polynomial*> guess_linear();


#endif// TALISMAN_SRC_SUBCIRCUIT_H_
