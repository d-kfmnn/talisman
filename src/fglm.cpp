/*------------------------------------------------------------------------*/
/*! \file fglm.cpp
    \brief contains functions used in the polynomial solver

  This file contains all functions used for preprocessing the Gröbner basis
  and for reducing the specification by the slices.

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#include "fglm.h"

#include <gmpxx.h>
#include <stdio.h>

#include <ranges>
#include <unordered_map>
#include <unordered_set>

#include "gmp.h"
#include "matrix.h"
#include "polynomial.h"
#include "reductionmethods.h"
#include "signal_statistics.h"
#include "subcircuit.h"
#include "substitution.h"

#include <iostream>

typedef std::vector<std::pair<mpz_class, size_t>> compressed_polynomial;
std::vector<std::vector<int>> indices;

compressed_polynomial
compress_linear(Polynomial* g, std::map<Var*, size_t>& var_to_id) {
  assert(g->degree() <= 1);

  compressed_polynomial p;
  for(size_t j = 0; j < g->len(); j++) {
    Monomial* m = g->get_mon(j);
    mpz_class c(m->coeff);
    size_t id = m->get_term() ? var_to_id[m->get_term()->get_var()] : 0;
    p.emplace_back(c, id);
  }
  return p;
}
/*------------------------------------------------------------------------*/
/*
 * Possible optimizations:
 *
 *   Store normal forms
 */
static std::vector<compressed_polynomial>
run_fglm(std::vector<Polynomial*>& normal_forms,
         std::map<Var*, size_t>& var_to_id) {
  count_fglm_call++;
  if(verbose > 2)
    msg("========= running run_fglm =========");
  std::vector<compressed_polynomial> compressed_res;
 
  // extract row and column terms
  std::unordered_set<Term*> term_set;
  std::vector<Term*> rows;
  std::vector<std::pair<Term*, int>> cols;
  for(size_t i = 0; i < normal_forms.size(); i++) {
    Polynomial* g = normal_forms[i];
    cols.emplace_back(g->get_lt(), i);
    term_set.insert(g->get_lt());
    for(size_t j = 1; j < g->len(); j++) {
      Term* t = g->get_mon(j)->get_term();
      if(term_set.insert(t).second)
        rows.push_back(t);
    }
  }

  // get linear terms for matrix and make
  // term-to-index map

  std::unordered_map<Term*, size_t> term_to_id;
  size_t i = 0;
  for(const auto& t : rows) {
    term_to_id[t] = i++;
    if(!t or t->degree() == 1)
      cols.emplace_back(t, -1);
  }

  // sort columns by decreasing leading terms
  std::sort(cols.begin(), cols.end(), [](auto& t1, auto& t2) {
    return cmp_term(t1.first, t2.first) == 1;
  });

  // set up matrix
  size_t n_rows = rows.size();
  size_t n_cols = cols.size();
  fmpq_mat_t mat;
  fmpq_mat_init(mat, n_rows, n_cols);
  long j = 0;

  for(const auto& [t, id] : cols) {
    // those are the linear terms
    if(id < 0) {
      fmpq_set_si(fmpq_mat_entry(mat, term_to_id[t], j), 1, 1);
    }
    // those are the normal forms
    else {
      Polynomial* g = normal_forms[id];
      // we assume that leading coefficient is +- 1
      assert(mpz_cmpabs_ui(g->get_lm()->coeff, 1) == 0);
      int sign = mpz_sgn(g->get_lm()->coeff);
      for(size_t k = 1; k < g->len(); k++) {
        Monomial* m = g->get_mon(k);
        size_t i = term_to_id[m->get_term()];
        mpz_class c(m->coeff);
        c = sign > 0 ? -c : c;
        fmpq_set_si(fmpq_mat_entry(mat, i, j), c.get_si(), 1);
      }
    }
    j++;
  }

  fmpq_mat_t K;
  kernel(mat, K);
  
  bool is_zero = true;
  for(long i = 0; i < fmpq_mat_nrows(K); i++)
    if(!row_is_zero(K,i)) {
      is_zero = false;
      break;
    }
  // kernel is zero --> we found no linear polynomials
  if(is_zero) {
    return compressed_res;
  }

  mpz_t tmp;
  mpz_init(tmp);
  // construct compressed form of linear polynomials
  compressed_res.reserve(fmpq_mat_nrows(K));
  indices.clear();
  for(long i = 0; i < fmpq_mat_nrows(K); i++) {
    //  check if row contains a denominator -> that would be bad
    if(!is_denom_free(K,i))
      continue;
    j = 0;
    while(fmpq_is_zero(fmpq_mat_entry(K,i,j))) j++;
   
    compressed_polynomial p;
    std::vector<int> indices_p;
    for(; j < fmpq_mat_ncols(K); j++) {
      if(fmpq_is_zero(fmpq_mat_entry(K,i,j)))
        continue;
      size_t id = 0;
      // constant coefficient gets id = 0
      Term* t = cols[j].first;
      if(t) {
        Var* v = t->get_var();
        assert(var_to_id.find(v) != var_to_id.end());
        id = var_to_id[v];
      }
      mpz_class c;
      fmpq_get_mpz_frac(c.get_mpz_t(), tmp, fmpq_mat_entry(K,i,j));
      p.emplace_back(c, id);

      // log coeff * nf = new_poly
      // cols[j].second >= 0 ensures that those come from normal_form
      if(proof_logging and cols[j].second >= 0) {
        Polynomial* nf = normal_forms[cols[j].second];
        c = -c;
        if(static_cast<int>(c.get_si()) != 1) {
          Polynomial* pp = multiply_poly_with_constant(nf, c.get_mpz_t());
          print_pac_mul_const_rule(
            proof_file, nf, static_cast<int>(c.get_si()), pp);
          indices_p.push_back(pp->get_idx());
          delete(pp);

        } else {
          indices_p.push_back(nf->get_idx());
        }
      }
    }
    compressed_res.push_back(p);
    indices.push_back(indices_p);
  }
 
  return compressed_res;
}

/*------------------------------------------------------------------------*/
std::vector<Polynomial*> linear_polies;
static void
construct_linear_polynomials(
                             std::vector<compressed_polynomial>& compressed,
                             std::map<Var*, size_t>& var_to_id) {
  linear_polies.clear();
  std::vector<Term*> id_to_term(var_to_id.size() + 1);
  for(auto [var, id] : var_to_id) {
    Term* t = new_term(var, 0);
    id_to_term[id] = t;
  }
  id_to_term[0] = 0;

  linear_polies.reserve(compressed.size());
  for(const auto& f : compressed) {
    for(const auto& [c, id] : f) {
      Term* t = id_to_term[id];
      Monomial* m
        = new Monomial(const_cast<mpz_ptr>(c.get_mpz_t()), t ? t->copy() : 0);
      push_mstack(m);
    }
    linear_polies.push_back(build_poly());
  }
}

static bool
update_gates(Gate* g) {
  bool flag = 0;
  
  if(verbose > 2)
    msg("Found the following linear polynomials:");
  for(size_t i = 0; i < linear_polies.size(); i++) {
    const auto& p = linear_polies[i];
    if(verbose > 2)
      p->print(stdout);
    if(proof_logging && indices.size() > 0) {
      print_pac_vector_add_rule(proof_file, indices[i], p);
    }
    //check_if_propagate(p);
  }

  for(const auto& p : linear_polies) {
    Gate* p_g = gate(p->get_lt()->get_var_num());
    if(p_g == g) {
      flag = 1;
      break;
    }
  }

  if(flag) {
    for(const auto& p : linear_polies) {
      Gate* p_g = gate(p->get_lt()->get_var_num());
      Polynomial* p_g_constraint = p_g->get_gate_constraint();

      if(p_g_constraint->degree() > 1 or cmp_poly(p_g_constraint, p) == 1) {
        bool case_one = p_g_constraint->degree() > 1;
        p_g->update_gate_poly(p, 0);

        if(verbose > 3) {
          if(case_one) {
            msg_nl("updated gate poly of %s to ", p_g->get_var_name());
          } else {
            msg_nl("updated linear gate poly of %s to ", p_g->get_var_name());
          }
          p_g->get_gate_constraint()->print(stdout);
        }
      } else if(!proof_logging) {
        delete p;
      }
    }
  }
  indices.clear();
  if(verbose > 2)
    msg("========= finished run_fglm =========");
  return flag;
}

/*------------------------------------------------------------------------*/
bool
cmp_np(const Normalized_poly& p1, const Normalized_poly& p2) {
  return p1.coeffs == p2.coeffs and p1.terms == p2.terms;
}

auto circuit_cmp = [](const std::vector<Normalized_poly>& c1,
                      const std::vector<Normalized_poly>& c2) -> bool {
  if(c1.size() != c2.size())
    return false;
  for(size_t i = 0; i < c1.size(); i++) {
    if(c1[i].coeffs != c2[i].coeffs or c1[i].terms != c2[i].terms)
      return false;
  }
  return true;
};

std::unordered_map<std::vector<Normalized_poly>,
                   std::vector<compressed_polynomial>,
                   circuit_hash,
                   decltype(circuit_cmp)>
  cached_circuits;

std::map<Var*, size_t> var_to_id;
std::vector<Normalized_poly> circuit;
std::vector<compressed_polynomial> cache;
std::vector<Polynomial*> normal_forms;
std::map<size_t, std::vector<Polynomial*>> used_van_mon;

/*------------------------------------------------------------------------*/
static bool
linearize_via_msolve(Gate* g) {
  count_msolve_call++;

  std::string vars;
  std::string poly;
  std::string red;
  std::vector<Gate*> gate_vec;

  int randomNumber = std::rand();

  // Construct the string
  std::string output = "./tmp/" + std::to_string(randomNumber) + ".ms";
  std::string respath = "./tmp/" + std::to_string(randomNumber) + ".out";
  std::string polypath = "./tmp/" + std::to_string(randomNumber) + ".line";

  FILE* f = fopen(output.c_str(), "w");
  if(!f)
    die(2, "cannot open file %s", output.c_str());

  for(const auto& tag : var) {
    fprintf(f, "%s", tag->get_var_name());
    if(tag != *var.rbegin())
      fprintf(f, ",");
  }

  fprintf(f, "\n");
  fprintf(f, "1073741827\n");

  for(const auto& gatep : gate_poly) {
    Polynomial* gpol = gatep->get_gate_constraint();
    Polynomial* tmp = unflip_poly(gpol);
    tmp->print(f, 0);
    fprintf(f, ",\n");
    delete(tmp);
  }

  for(const auto& gatep : var) {
    fprintf(f, "-%s^2+%s", gatep->get_var_name(), gatep->get_var_name());
    if(gatep != *var.rbegin())
      fprintf(f, ",\n");
  }
  fclose(f);

  std::string msolvecall = "msolve -f " + output + " -g 2 ";
  std::string add_grep
    = msolvecall + "| grep -m2 " + g->get_var_name() + " | tail -n1 | ";
  std::string remove_ones = add_grep
                            + "sed 's/\\(\\^\\)1\\b//g; s/\\[//g; "
                              "s/+1073741826/-1/g ; s/+1073741825/-2/g' > "
                            + respath;
  std::string call = remove_ones;

  system(call.c_str());
  Polynomial* target = parse_specification_polynomial(respath.c_str());
  std::string delcall = "rm " + output + " " + respath;
  system(delcall.c_str());
  if(target->degree() > 1)
    return 0;

  g->update_gate_poly(target);

  return 1;
}

/*------------------------------------------------------------------------*/
static int
internal_linearize(Gate* g,
                            int depth,
                            size_t fanout_size,
                            int init,
                            bool single_expand) {
  // statistics
  total_circuit_lin_count++;

  double call_init_time = process_time();

  bool res = 0;

  var_to_id.clear();
  circuit.clear();
  cache.clear();
  normal_forms.clear();

  // get subcircuit
  
  if(!get_and_compress_subcircuit(
       g, depth, fanout_size, init, single_expand, circuit, var_to_id))
    return -1;
 


  circuit_hash hasher;
  std::size_t hash_value = hasher(circuit);// Get hash value
  std::vector<size_t> indices_input_new_pattern;

  for(const auto& gatep : gate_poly) {
    if(gatep->get_nf())
      indices_input_new_pattern.push_back(gatep->get_nf()->get_idx());
  }

  std::vector<Polynomial*> guessed;
  std::vector<Polynomial*> new_nf_poly;
  // check cache
  bool found_cache = false;
  if(cached_circuits.find(circuit) != cached_circuits.end()) {
    found_cache = true;
    cache = cached_circuits[circuit];

    if(verbose > 1)
      msg("found a cached circuit at dist %i", g->get_dist());
    circut_cached_count++;

  } else if(!msolve) {
    
  
    if(is_internal_fsa(g) && !force_fglm) {   
      double pre_gap_time = process_time(); 
      linear_polies = guess_linear();

      if(linear_polies.size() == 0) {
        unmark_fsa();
        gate_poly.clear();
        sc_inputs.clear();
        var.clear();
      }
      
      
      gap_time += (process_time() - pre_gap_time);
      linearization_time += process_time() - call_init_time;
      return update_gates(g);
    } else if (force_guessing) {
      double pre_gap_time = process_time(); 
      linear_polies = guess_linear();
      gap_time += (process_time() - pre_gap_time);
      std::vector<compressed_polynomial> cache;
      for(auto& poly: linear_polies){
        cache.push_back(compress_linear(poly, var_to_id));
      }

      if(do_caching) cached_circuits[circuit] = cache;
      linearization_time += (process_time() - call_init_time);
      return update_gates(g);
   
    } else {
      double pre_fglm_time = process_time(); 
      // reduce gates
      std::vector<Polynomial*> used_van_poly;
      
   
      int i = 0;
      if(proof_logging && do_caching) {
        fprintf(proof_file, "pattern_new %lu {\n", hash_value);

        for(auto& v : var_to_id) {
          v.first->set_id(v.second);
        }

        for(const auto& gatep : gate_poly) {
          fprintf(proof_file, "in%i %lu ", i++, gatep->get_nf()->get_idx());

          gatep->print_nf(proof_file);
        }
      }


      double pre_nf_time = process_time();
      normal_forms = compute_normalforms(&used_van_poly, &new_nf_poly);
      nf_time += (process_time() - pre_nf_time);

      if(proof_logging && do_caching) {
        used_van_mon.insert({ hash_value, used_van_poly });
      }

    
      

      assert(normal_forms.size() > 0);

      // exit if linearisation was found during normal form computation
      if(g->get_gate_constraint()->degree() == 1) {
        msg("found desired linear poly during computing normal forms");
        res = 1;
        fglm_time += (process_time() - pre_fglm_time);
        goto clean_up;
      } else {
        // run fglm
        double pre_matrix_time = process_time();
        cache = run_fglm(normal_forms, var_to_id);
        matrix_time += (process_time() - pre_matrix_time);

        // cache result
        if(do_caching)
          cached_circuits[circuit] = cache;
      }
      fglm_time += (process_time() - pre_fglm_time);
    }
  } else {
    res = linearize_via_msolve(g);
    if(res) {
      compressed_polynomial g_compr
        = compress_linear(g->get_gate_constraint(), var_to_id);
      cache.push_back(g_compr);
      cached_circuits[circuit] = cache;
    }
    linearization_time += (process_time() - call_init_time);
    return res;
  }

  // construct polynomials
  construct_linear_polynomials(cache, var_to_id);

  res = update_gates(g);


  if(proof_logging && do_caching) {
    if(!found_cache) {
      int i = 0;
      for(auto& p : linear_polies) {
        fprintf(proof_file, "out%i %lu;\n", i++, p->get_idx());
      }
      for(auto& p : new_nf_poly) {
        fprintf(proof_file, "out%i %lu;\n", i++, p->get_idx());
      }
      fprintf(proof_file, "};\n");
    }
   

    for(auto& v : var_to_id) {
      v.first->set_id(0);
    }


    fprintf(proof_file, "pattern_apply %lu {\n", hash_value);
    for(auto& v : var_to_id) {
      fprintf(proof_file, "v%lu  %s;\n", v.second, v.first->get_name());
    }


    int i = 0;
 

    for(const auto& gatep : indices_input_new_pattern) {
      fprintf(proof_file, "in%i %lu;\n", i++, gatep);
    }

    std::vector<Polynomial*> van_p_pattern = used_van_mon[hash_value];
    for(const auto& gatep : van_p_pattern) {
      fprintf(proof_file, "in%i %lu;\n", i++, gatep->get_idx());
    }

    int j = print_pac_pattern_out_rules(proof_file, linear_polies, 0);

    print_pac_pattern_out_rules(proof_file, new_nf_poly, j);


    fprintf(proof_file, "};\n");
  }

clean_up:


  if(res) {
    for(size_t i = 0; i < normal_forms.size(); i++) {
      Gate* tmp = gate(normal_forms[i]->get_lt()->get_var_num());
      delete(tmp->get_nf());
      tmp->set_nf(0);
    }
  }
  linearization_time += (process_time() - call_init_time);
  return res;
}

/*------------------------------------------------------------------------*/
// Function from outside
bool
linearize_via_fglm_or_gap(Gate* g) {
  count_unique_gb_call++;
  int max_depth = g->get_dist();
  int depth = sc_depth;
  size_t fanout_size = sc_fanout;
  fanout_size_last_call = 0;

  int count = 1;
  int res = internal_linearize(g, depth, fanout_size, count++, 0);

  while(!res && depth < max_depth) {
    circuit_enlarged_count++;
    if(count % 15 == 0) {
      res = internal_linearize(g, depth, fanout_size_last_call + 1, 1, 0);
      count++;
      if(!res) {
       res = internal_linearize(g, ++depth, fanout_size, 1, 0);
       count ++;
      }
    }

    if(!res) res = internal_linearize(g, depth, fanout_size, count++, 1);


    if(res == -1 && max_depth <= 6) {
      if(count-2 > max_depth_count)
        max_depth_count = count-2;
      return 0;
    }

 
  }
  if(count-2 > max_depth_count)
    max_depth_count = count-2;

  return res;
}
