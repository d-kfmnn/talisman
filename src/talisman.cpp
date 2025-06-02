/*------------------------------------------------------------------------*/
/*! \file talisman.cpp
    \brief main file of our tool TalisMan

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#define VERSION "1.0"
/*------------------------------------------------------------------------*/
// / Manual of TalisMan, will be printed with command line '-h'
static const char *USAGE =
    "\n"
    "### USAGE ###\n"
    "usage : talisman <input file> <spec-mode> [proof-logging] [options] \n"
    "\n"
    "General Help\n"
    "------------\n"
    "    -h | --help           Displays usage information and exits.\n"
    "\n"
    "\n"
    "Input File <input file>\n"
    "--------------\n"
    "    <input file>          File containing the AIG graph"
    "\n"
    "\n"
    "Specifications <spec-mode>\n"
    "--------------\n"
    "    Specify a type of verification spec:\n"
    "    -miter-spec           Use the miter specification.\n"
    "    -mult-spec            Use the multiplier specification.\n"
    "    -assert-spec          Use the assertion specification.\n"
    "    <spec_file>           If none of the above pre-determined specs is used, a spec file has to be provided.\n"
    "\n"
    "    Note: Only one specification can be selected. Attempting to select multiple specs will result in an error.\n"
    "\n"
    "\n"
    "Proof Logging [proof-logging]\n"
    "-------------\n"
    "    -proofs [level] <poly> <steps> <spec>    Enable proof logging mode.\n"
    "      Note: Only one proof format can be selected. If multiple formats are specified, an error will occur.\n"
    "\n"
    "    If -proofs is set three output files have to provided in the following order:\n"
    "      <poly>              File for proof axioms.\n"
    "      <steps>             File for proof steps.\n"
    "      <spec>              File for proof spec.\n"
    "\n"
    "\n"
    "ADDITIONAL OPTIONS"
    "\n"
    "Counter-Example Generation\n"
    "--------------------------\n"
    "  -nce | --no-counter-examples     Disables counter-example generation in case of incorrect circuit.\n"
    "\n"
    "Sub-Circuit size\n"
    "--------------------------\n"
    "  -f <int>                         Non-negative value for fanout size, 0 turns fanout limit off (default value: 4).\n"
    "  -d <int>                         Positive value for depth (default: 2).\n"
    "\n"
    "Ablation\n"
    "--------------------------\n"
    "  -npp  | --no-preprocessing        Disables the preprocessing phase. (no rewriting of AIG).\n"
    "  -nvc  | --no-vanishing            Turns on vanishing constraints \n"
    "  -nch  | --no-caching              Turns off caching of circuits \n"
    "  -dll  | --do-local-linearization  Enables the local linearization and only uses FGLM to linearize.\n"
    "  -alg  | --algebraic-reduction     Use algebraic reductions instead of SAT in guess and proof\n"
    "  -gap  | --force-guessing          Forces the linearization to only use guess-and-proof\n"
    "  -fglm | --force-fglm              Forces the linearization to only use fglm\n"
    "\n"
    "\n"
    "Verbosity Levels\n"
    "----------------\n"
    "    Control the level of output detail:\n"
    "    -v0                   Minimal output (silent mode).\n"
    "    -v1                   Low verbosity (default).\n"
    "    -v2                   Medium verbosity.\n"
    "    -v3                   High verbosity.\n"
    "    -v4                   Maximum verbosity (debug-level).\n"
    "\n"

    "Example Usages\n"
    "-------------\n"
    "    talisman input.aig -mult-spec\n"
    "    talisman input.aig spec.txt\n"
    "    talisman -v3 -proofs -p2 -miter-spec input.aig output1.txt output2.txt output3.txt\n";

/*------------------------------------------------------------------------*/
#include <cctype>  // for std::isdigit
#include <cstdlib> // for std::atoi
#include <ctime>

#include "gate.h"
#include "parser.h"
#include "polynomial_solver.h"
#include "specpoly.h"
#include "substitution.h"
/*------------------------------------------------------------------------*/
// / Name of the input file
static const char *input_name = 0;
static const char *spec_name = 0;
static bool spec_selected = 0;

// / \brief
// / Name of first output file, which stores the CNF miter in '-substitute', and
// / the gate constraints in '-certify'.
static const char *output_name1 = 0;

// / \brief
// / Name of second output file, storing the rewritten AIG in '-substitute',
// / and the core proof in '-certify'.
static const char *output_name2 = 0;

// / Name of third output file. Stores the specification in '-certify'.
static const char *output_name3 = 0;


/*------------------------------------------------------------------------*/
// ERROR CODES:

static int err_no_file = 10;    // no input file given
static int err_spec_sel = 11;   // mode has already been selected/not selected
static int err_wrong_arg = 12;  // wrong number of arguments given
static int err_proof_form = 13; // too many proof formats selected

/*------------------------------------------------------------------------*/
/**
    Calls the deallocaters of the involved data types
    @see reset_all_signal_handlers()
    @see delete_gates()
    @see deallocate_terms()
    @see deallocate_mstack()
    @see clear_mpz()
*/
static void reset_all()
{
  reset_all_signal_handlers();
  //delete_gates(); 
  deallocate_terms();
  deallocate_mstack();
  clear_mpz();
}
/*------------------------------------------------------------------------*/
bool is_number(const std::string &s)
{
  return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
}
//----------------------------------------------------------------------

/**
    Main Function of TalisMan.
    Reads the given AIG and depending on the selected mode, either
    calls the substution engine or the polynomial solver.

    Prints statistics to stdout after finishing.
*/
int main(int argc, char **argv)
{
  unsigned int seed= time (0);
  std::srand (seed);
  init_time();

  msg("TalisMan Version " VERSION);
  msg("");
  msg("Copyright(C) 2025 Daniela Kaufmann, TU Wien, Austria");
  msg("                  Clemens Hofstadler, JKU Linz, Austria");
  msg("");

  for (int i = 1; i < argc; i++)
  {
    if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help"))
    {
      fputs(USAGE, stdout);
      fflush(stdout);
      exit(0);
    }
    else if (!strcmp(argv[i], "-v0"))
    {
      verbose = 0;
    }
    else if (!strcmp(argv[i], "-v1"))
    {
      verbose = 1;
    }
    else if (!strcmp(argv[i], "-v2"))
    {
      verbose = 2;
    }
    else if (!strcmp(argv[i], "-v3"))
    {
      verbose = 3;
    }
    else if (!strcmp(argv[i], "-v4"))
    {
      verbose = 4;
    }
    else if (!strcmp(argv[i], "-proofs"))
    {
      proof_logging = true;
    }
    else if (!strcmp(argv[i], "-f") && i + 1 < argc)
    {

      std::string arg_value = argv[i + 1];

      if (is_number(arg_value))
      {
        int value = std::stoi(arg_value); // Convert to integer
        if (value < 0)
          die(123, "-f needs to be followed by a non-negative integer");
        else
          sc_fanout = value;
      }
      else
     
        die(123, "-f needs to be followed by a non-negative integer");
      i++;
    }
    else if (!strcmp(argv[i], "-d") && i + 1 < argc)
    {
      std::string arg_value = argv[i + 1];

      if (is_number(arg_value))
      {
        int value = std::stoi(arg_value); // Convert to integer
        if (value <= 0)
          die(123, "-d needs to be followed by a positive integer");
        else
          sc_depth = value;
      }
      else
     
        die(123, "-d needs to be followed by a positive integer");
      i++;
    }

    else if (!strcmp(argv[i], "-miter-spec"))
    {
      if (spec_selected)
        die(err_spec_sel, "only one specification is allowed (try '-h')");

      miter_spec = 1;
      no_spec = 0;
      spec_selected = 1;
    }
    else if (!strcmp(argv[i], "-mult-spec"))
    {
      if (spec_selected)
        die(err_spec_sel, "only one specification is allowed (try '-h')");

      mult_spec = 1;
      no_spec = 0;
      spec_selected = 1;
    }
    else if (!strcmp(argv[i], "-assert-spec"))
    {
      if (spec_selected)
        die(err_spec_sel, "only one specification is allowed (try '-h')");

      assert_spec = 1;
      no_spec = 0;
      spec_selected = 1;
    }
    else if (!strcmp(argv[i], "--no-counter-examples") || (!strcmp(argv[i], "-nce")))
    {
      gen_witness = 0;
    }
    else if (!strcmp(argv[i], "--no-vanishing") || (!strcmp(argv[i], "-nvc")))
    {
      force_vanishing_off = 1;
    }
    else if (!strcmp(argv[i], "--no-caching") || (!strcmp(argv[i], "-nch")))
    {
      do_caching = 0;
    }
    else if (!strcmp(argv[i], "--algebraic-reduction") || (!strcmp(argv[i], "-alg")))
    {
      use_algebra_reduction = 1;
    }
    else if (!strcmp(argv[i], "--no-preprocessing") || (!strcmp(argv[i], "-npp")))
    {
      do_preprocessing = 0;
    }
    else if (!strcmp(argv[i], "--do-local-linearization") || (!strcmp(argv[i], "-dll")))
    {
      do_local_lin = 1;
    }
    else if (!strcmp(argv[i], "--force-fglm") || (!strcmp(argv[i], "-fglm")))
    {
      force_fglm = 1;
    }
    else if (!strcmp(argv[i], "--force-guessing") || (!strcmp(argv[i], "-gap")))
    {
      force_guessing = 1;
      sc_depth = 4;
    }
    else if (!strcmp(argv[i], "--msolve") || (!strcmp(argv[i], "-m")))
    {
      msolve = 1;
    }
    else if (output_name3)
    {
      die(err_wrong_arg, "too many arguments '%s', '%s', '%s', '%s' and '%s'(try '-h')", input_name, output_name1,
          output_name2, output_name3, argv[i]);
    }
    else if (output_name2)
    {
      output_name3 = argv[i];
    }
    else if (output_name1)
    {
      output_name2 = argv[i];
    }
    else if (spec_name || (input_name && !no_spec))
    {
      output_name1 = argv[i];
    }
    else if (input_name && no_spec)
    {
      spec_name = argv[i];
    }
    else
    {
      input_name = argv[i];
    }
  }

  if (!input_name)
    die(err_no_file, "no input file given (try '-h')");
  if (no_spec && !spec_name)
  {

    die(err_no_file, "no spec file in modus 'no_spec' given (try '-h')");
  }

  if(proof_logging && (!output_name1 || !output_name2 || !output_name3))
  {
    die(err_proof_form, "proof logging requires three output files (try '-h')");
  }

  if(!proof_logging && output_name1){
    die(err_proof_form, "invalid option '%s' (try '-h')", output_name1);
  }

  if(proof_logging && msolve){
    die(123, "invalid combination of options: proof logging is not supported by msolve (try '-h')");
  }

  if(force_fglm && force_guessing){
    die(123, "invalid combination of options: fglm and guessing cannot be forced at the same time (try '-h')");
  }

  if(force_fglm && use_algebra_reduction){
    die(123, "invalid combination of options: algebra reduction can only be used in guessing (try '-h')");
  }



  print_hline();
  msg("SETTINGS");

  msg("seed: %i", seed);
  msg("preprocessing: %s", do_preprocessing ? "enabled" : "disabled");
  msg("vanishing constraints: %s", do_vanishing_constraints ? "enabled" : (force_vanishing_off ? "disabled" : "partially enabled"));
  msg("local linearization: %s", do_local_lin ? "enabled" : "disabled");
  msg("caching: %s", do_caching ? "enabled" : "disabled");
  msg("");
  msg("fanout limitation: %s", sc_fanout ? "enabled" : "disabled");
  if(sc_fanout) msg("subcircuit fanout: %i", sc_fanout);
  msg("subcircuit depth: %i", sc_depth);
  msg("");
  msg("linearization: %s", msolve ? "Groebner basis using msolve" : "Matrix-based using normal forms");
  msg("reduction: %s", use_algebra_reduction ? "Ideal membership" : "Kissat");
  msg("");

  if (no_spec)
    msg("spec from file %s will be used", spec_name);
  else if (miter_spec)
    msg("specification: miter");
  else if (mult_spec)
    msg("specification: unsigned multiplier");
  else if (assert_spec)
    msg("specification: assertion");

  if (proof_logging)
  {
   msg("proof logging: enabled");
  }

  print_hline();

  // Initialization Phase
  init_all_signal_handers();
  init_nonces();

  parse_aig(input_name);
  bool res;
  init_mpz(NN);
  init_gates();
  print_hline();


  // Generate or Parse Spec
  Polynomial *spec;
  if (mult_spec)
    spec = mult_spec_poly();
  else if (miter_spec)
    spec = miter_spec_poly();
  else if (assert_spec)
    spec = assertion_spec_poly();
  else
    {
      msg("reading specification polynomial from '%s'", spec_name);
      spec = parse_specification_polynomial(spec_name);
    }
  if(verbose > 1){
    msg("generated spec poly: ");
    msg_nl("");
    spec->print(stdout);
    print_hline();
  }

  res = verify(input_name, spec, output_name1, output_name2, output_name3); 

  // Resetting
  if (spec)
    delete (spec);

  reset_aig_parsing();
  reset_all();

  print_statistics();
  //msg("called build_poly %i", running_idx);
  return res;
}
