/*------------------------------------------------------------------------*/
/*! \file signal_statistics.cpp
    \brief used to handle signals, messages and statistics

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/

#include "signal_statistics.h"

/*------------------------------------------------------------------------*/
// Global variable
void (*original_SIGINT_handler)(int);
void (*original_SIGSEGV_handler)(int);
void (*original_SIGABRT_handler)(int);
void (*original_SIGTERM_handler)(int);
/*------------------------------------------------------------------------*/

// Specifications
bool miter_spec = 0;  // single output, spec: s = 1
bool mult_spec = 0;   // spec: sum(2^is_i) = sum(2^ia_i)*sum(2^ib_i)
bool assert_spec = 0; // multiple outputs, sum(outputs) - num(outputs) = 0
bool no_spec = 1;     // spec is provided in extra file

// Subcircuit sizes
size_t sc_depth = 2;
size_t sc_fanout = 4;

// Ablation studies
bool do_preprocessing = 1;
bool do_vanishing_constraints = 0;
bool do_caching = 1;
bool do_local_lin = 0;
bool normal_form_top_down = 1;
bool msolve = 0;
bool use_algebra_reduction = 0;
bool force_fglm = 0;
bool force_guessing = 0;
bool proof_logging = 0;
bool force_vanishing_off = 0;

// Statistics
int van_mon_depth_count = 0;
int van_mon_prop_count = 0;
int l_f_count = 0;
int node_has_only_two_grandchildren_count = 0;
int children_share_l_f_count = 0;
int lin_xor_constraint_count = 0;
int total_circuit_lin_count = 0;
int equiv_gate_count = 0;
int circut_cached_count = 0;
int count_fglm_call = 0;
int count_unique_gb_call = 0;
int count_msolve_call = 0;
int circuit_enlarged_count = 0;
int max_depth_count = 0;
bool booth = 0;

int correct_guess_count = 0;
int count_guess_call = 0;
int count_kissat_call = 0;
int evaluated_guess_count= 0;
int total_guesses_count = 0;
int max_guesses_count = 0;
int max_iterations_count = 0;
int total_iterations_count = 0;
std::vector<double>accuracy (100, 0.0);;
std::vector<int>iteration_on_level (100, 0);;


// Time
double linearization_time = 0;
double fglm_time = 0;
double nf_time = 0;
double matrix_time = 0;
double find_circuit_time = 0;
double gap_time = 0;
double proof_time = 0;
double guess_time = 0;

int van_mon_poly_count = 0;
int van_mon_used_count = 0;

// Reductions
int non_linear_count = 0;
int linear_count = 0;

FILE *proof_file = NULL;
FILE *polys_file = NULL;

// Global variable
int verbose = 1;

/*------------------------------------------------------------------------*/

static const char *signal_name(int sig)
{
  switch (sig)
  {
  case SIGINT:
    return "SIGINT";
  case SIGSEGV:
    return "SIGSEGV";
  case SIGABRT:
    return "SIGABRT";
  case SIGTERM:
    return "SIGTERM";
  default:
    return "SIGUNKNOWN";
  }
}

/*------------------------------------------------------------------------*/

static void catch_signal(int sig)
{
  printf("c\nc caught signal '%s'(%d)\nc\n", signal_name(sig), sig);
  printf("c\nc raising signal '%s'(%d) again\n", signal_name(sig), sig);
  reset_all_signal_handlers();
  fflush(stdout);
  raise(sig);
}

/*------------------------------------------------------------------------*/

void init_all_signal_handers()
{
  original_SIGINT_handler = signal(SIGINT, catch_signal);
  original_SIGSEGV_handler = signal(SIGSEGV, catch_signal);
  original_SIGABRT_handler = signal(SIGABRT, catch_signal);
  original_SIGTERM_handler = signal(SIGTERM, catch_signal);
}

/*------------------------------------------------------------------------*/

void reset_all_signal_handlers()
{
  (void)signal(SIGINT, original_SIGINT_handler);
  (void)signal(SIGSEGV, original_SIGSEGV_handler);
  (void)signal(SIGABRT, original_SIGABRT_handler);
  (void)signal(SIGTERM, original_SIGTERM_handler);
}

/*------------------------------------------------------------------------*/

void msg_nl(const char *fmt, ...)
{
  va_list ap;
#ifdef HAVEUNLOCKEDIO
  fputs_unlocked("[talisman] ", stdout);
#else
  fputs("[talisman] ", stdout);
#endif
  va_start(ap, fmt);
  vfprintf(stdout, fmt, ap);
  va_end(ap);
  fflush(stdout);
}
/*------------------------------------------------------------------------*/

void msg(const char *fmt, ...)
{
  va_list ap;
#ifdef HAVEUNLOCKEDIO
  fputs_unlocked("[talisman] ", stdout);
#else
  fputs("[talisman] ", stdout);
#endif
  va_start(ap, fmt);
  vfprintf(stdout, fmt, ap);
  va_end(ap);
#ifdef HAVEUNLOCKEDIO
  fputc_unlocked('\n', stdout);
#else
  fputc('\n', stdout);
#endif
  fflush(stdout);
}

/*------------------------------------------------------------------------*/
void print_hline()
{
#ifdef HAVEUNLOCKEDIO
  fputs_unlocked("[talisman] ", stdout);
  fputs_unlocked("-------------------------------------------------------", stdout);
  fputc_unlocked('\n', stdout);
#else
  fputs("[talisman] ", stdout);
  fputs("-------------------------------------------------------", stdout);
  fputc('\n', stdout);
#endif
  fflush(stdout);
}

/*------------------------------------------------------------------------*/

void die(int error_code, const char *fmt, ...)
{
  fflush(stdout);
  va_list ap;
  fprintf(stderr, "*** [talisman] error code %i \n", error_code);
#ifdef HAVEUNLOCKEDIO
  fputs_unlocked("*** [talisman] ", stderr);
#else
  fputs("*** [talisman] ", stderr);
#endif
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  fputc('\n', stderr);
  fflush(stderr);
  exit(error_code);
}

/*------------------------------------------------------------------------*/

static size_t maximum_resident_set_size()
{
  struct rusage u;
  if (getrusage(RUSAGE_SELF, &u))
    return 0;
  return ((size_t)u.ru_maxrss) << 10;
}

/*------------------------------------------------------------------------*/
struct timeval start_tv;
void init_time()
{
  gettimeofday(&start_tv, NULL);
}

/*------------------------------------------------------------------------*/

double process_time()
{
  struct timeval tv;
  double elapsed = 0.0;

  gettimeofday(&tv, NULL);
  elapsed = (tv.tv_sec - start_tv.tv_sec) + (tv.tv_usec - start_tv.tv_usec) / 1000000.0;
  return elapsed;

  /* struct rusage u;
   if (getrusage(RUSAGE_SELF, &u)) return 0;
   double res = u.ru_utime.tv_sec + 1e-6 * u.ru_utime.tv_usec;
   res += u.ru_stime.tv_sec + 1e-6 * u.ru_stime.tv_usec;
   return res;*/
}

static double percent(double a, double b)
{
  return b ? 100.0 * a / b : 0;
}
static double average(double a, double b)
{
  return b ? a / b : 0;
}
/*------------------------------------------------------------------------*/

void print_statistics()
{
  msg("");
  print_hline();
  msg("STATISTICS:");
  msg("");
  msg("LINEARIZATION");
  msg("total linearization calls: %13i", total_circuit_lin_count);
  msg("unique linearization calls:%13i", count_unique_gb_call);
  msg("sub-circuits enlarged:     %13i (max: %i times)", circuit_enlarged_count, max_depth_count); //buggy
  msg("");
  msg("cached circuits found:     %13i (%6.2f%% of total linearizations)", circut_cached_count, percent(circut_cached_count, total_circuit_lin_count));
  msg("new computations:          %13i (%6.2f%% of total linearizations)", total_circuit_lin_count-circut_cached_count, percent(total_circuit_lin_count-circut_cached_count, total_circuit_lin_count));
  int unique = total_circuit_lin_count-circut_cached_count;
  msg("  guess and prove calls:   %13i (%6.2f%% of new computations)", count_guess_call, percent(count_guess_call, unique));
  msg("    kissat calls:          %13i", count_kissat_call);
  msg("    guessed poly:          %13i (max: %2i, avg: %3.1f)", total_guesses_count, max_guesses_count, average(total_guesses_count, total_iterations_count));
  msg("    evaluated guessed poly:%13i (%6.2f%% of total guesses)", evaluated_guess_count, percent(evaluated_guess_count, total_guesses_count));
  msg("    correct guessed poly:  %13i (%6.2f%% of evaluated guesses)", correct_guess_count, percent(correct_guess_count, evaluated_guess_count));
  msg("    iterations:            %13i (max: %2i, avg: %3.1f)", total_iterations_count, max_iterations_count, average(total_iterations_count, count_guess_call));
  msg_nl("    average accuracies:               ");
   
    for(int i = 0; i<max_iterations_count;i++){
      fprintf(stdout, "%6.2f%% ", average(accuracy[i], iteration_on_level[i]));
      if(i< max_iterations_count-1) fprintf(stdout, "- ");
    }  
  fprintf(stdout, "\n");
  msg("");
  msg("  fglm calls:              %13i (%6.2f%% of new computations)", count_fglm_call, percent(count_fglm_call, unique));
  msg("  msolve calls:            %13i (%6.2f%% of new computations)", count_msolve_call, percent(count_msolve_call, unique));
  
  msg("");
  msg("VANISHING MONOMIALS: ");
  msg("total number:              %13i", van_mon_poly_count);
  msg("propagations:              %13i", van_mon_prop_count);
  msg("applications:              %13i", van_mon_used_count);
  
  msg("");
  msg("REDUCTIONS: ");
  msg("total reductions:          %13i", linear_count + non_linear_count);
  msg("linear reductions:         %13i (%6.2f %)", linear_count, percent(linear_count, linear_count + non_linear_count));
  msg("non-linear reductions:     %13i (%6.2f %)", non_linear_count,
      percent(non_linear_count, linear_count + non_linear_count));
  msg("");
  msg("TIME AND MEMORY: ");
  msg("maximum resident set size:     %12.2f MB", maximum_resident_set_size() / static_cast<double>((1 << 20)));
  double end_time = process_time();
  msg("total process time:            %13.3f seconds", end_time);
  msg("");
  msg("linearization time:            %13.3f seconds (%2.2f %% of total time)", linearization_time, percent(linearization_time, end_time));
  msg("  getting circuits time:       %13.3f seconds (%2.2f %% of linearization time)", find_circuit_time, percent(find_circuit_time, linearization_time));
  msg("  fglm time:                   %13.3f seconds (%2.2f %% of linearization time)", fglm_time, percent(fglm_time, linearization_time));
  msg("  guess-and-prove time:        %13.3f seconds (%2.2f %% of linearization time)", gap_time, percent(gap_time, linearization_time));  
  msg("");
  msg("  fglm time:                   %13.3f seconds (%2.2f %% of linearization time)", fglm_time, percent(fglm_time, linearization_time));
  msg("   used time for normal forms :        %4.3f seconds (%2.2f %% of fglm time)", nf_time, percent(nf_time, fglm_time));
  msg("   used time for linear combinations : %4.3f seconds (%2.2f %% of fglm time)", matrix_time, percent(matrix_time, fglm_time));
  msg("");
  msg("  guess-and-prove time:         %13.3f seconds (%2.2f %% of linearization time)", gap_time, percent(gap_time, linearization_time));
  msg("   used time for guessing :             %4.3f seconds (%2.2f %% of g&p time)", guess_time, percent(guess_time, gap_time));
  msg("   used time for proving :             %4.3f seconds (%2.2f %% of g&p time)", proof_time, percent(proof_time, gap_time));
  print_hline();
}
