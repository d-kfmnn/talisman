/*------------------------------------------------------------------------*/
/*! \file signal_statistics.h
    \brief used to handle signals, messages and statistics

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_SIGNAL_STATISTICS_H_
#define TALISMAN_SRC_SIGNAL_STATISTICS_H_
/*------------------------------------------------------------------------*/
#include <signal.h>
#include <stdarg.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <vector>

#include <iostream>
/*------------------------------------------------------------------------*/
extern void (*original_SIGINT_handler)(int);
extern void (*original_SIGSEGV_handler)(int);
extern void (*original_SIGABRT_handler)(int);
extern void (*original_SIGTERM_handler)(int);

// Specifications
extern bool miter_spec;  // single output, spec: s = 1
extern bool mult_spec;   // spec: sum(2^is_i) = sum(2^ia_i)*sum(2^ib_i)
extern bool assert_spec; // multiple outputs, sum(outputs) - num(outputs) = 0
extern bool no_spec;     // spec is provided in extra file

// Ablation studies
extern bool do_preprocessing;
extern bool do_vanishing_constraints;
extern bool do_caching;
extern bool do_local_lin;
extern bool normal_form_top_down;
extern bool msolve;
extern bool use_algebra_reduction;
extern bool force_fglm;
extern bool force_guessing;
extern bool force_vanishing_off;

// Statistic counters
extern int van_mon_depth_count;
extern int l_f_count;
extern int node_has_only_two_grandchildren_count;
extern int children_share_l_f_count;
extern int count_guess_call;
extern int count_kissat_call;
extern int count_fglm_call;
extern int count_msolve_call;
extern int total_circuit_lin_count;
extern int equiv_gate_count;
extern int lin_xor_constraint_count;
extern int non_linear_count;
extern int linear_count;
extern int circut_cached_count;
extern int van_mon_poly_count;
extern int van_mon_prop_count;
extern int van_mon_used_count;
extern int count_unique_gb_call;
extern int circuit_enlarged_count;
extern int max_depth_count;
extern int correct_guess_count;
extern int evaluated_guess_count;
extern int total_guesses_count;
extern int max_guesses_count;
extern int max_iterations_count;
extern int total_iterations_count;
extern std::vector<double>accuracy;
extern std::vector<int>iteration_on_level;

// / Enables proof logging
extern bool proof_logging;

extern FILE *proof_file;
extern FILE *polys_file;

// Subcircuit sizes
extern size_t sc_depth;
extern size_t sc_fanout;

extern bool booth;



extern double linearization_time;
extern double fglm_time;
extern double nf_time;
extern double matrix_time;
extern double find_circuit_time;
extern double gap_time;
extern double proof_time;
extern double guess_time;

extern struct timeval start_tv;

// / Level of output verbosity, ranges from 0 to 4
extern int verbose;

/**
    Initialize all signals
*/
void init_all_signal_handers();

void init_time();

/**
    Resets all signal handlers
*/
void reset_all_signal_handlers();

/*------------------------------------------------------------------------*/

/**
    Prints an error message to stderr and exits the program

    @param char* fmt message
    @param int error_code
*/
void die(int error_code, const char *fmt, ...);

/**
    Prints a message to stdout

    @param char* fmt message
*/
void msg(const char *fmt, ...);
void msg_nl(const char *fmt, ...);

void print_hline();
/*------------------------------------------------------------------------*/

// / Time measures used for verify/certify modus

/**
    Determines the used process time
*/
double process_time();

/**
    Print statistics of maximum memory and used process time depending on
    selected modus

    @param modus integer
*/
void print_statistics();

#endif // TALISMAN_SRC_SIGNAL_STATISTICS_H_
