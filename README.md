
TalisMan - Trusted Algebraic LInearization of Sub-Circuits with Matrix-based Algorithms using Normalforms
================================================================================

Our tool TalisMan 1.0 is able to parse and verify AIGs.

For further information we refer to the paper

Clemens Hofstadler, Daniela Kaufmann. 
`Guess and Prove: A Hybrid Approach to Linear Polynomial Recovery in Circuit Verification`
In Proc. 31st Intl. Conference on Principles and Practice of Constraint Programming (CP) 2025.


  
----------------------------------------------------------------  
  
# Dependencies:  
`libgmp` (https://gmplib.org/)  
              `lflint` (https://flintlib.org/)  
              `lpb` (https://github.com/master-keying/pblib)  
              `lkissat` (https://github.com/arminbiere/kissat)  

# Installation 

Use `./configure.sh && make` to configure and build `TalisMan 1.0`.


# Usage 
    talisman <input file> <spec-mode> [options] 

 General Help
 ------------
    -h | --help           Displays usage information and exits.


Input File <input file>
--------------
    <input file>          File containing the AIG graph

Specifications <spec-mode>
--------------
    Specify a type of verification spec:
    -miter-spec           Use the miter specification.
    -mult-spec            Use the multiplier specification.
    -assert-spec          Use the assertion specification.
    <spec_file>           If none of the above pre-determined specs is used, a spec file has to be provided.

    Note: Only one specification can be selected. Attempting to select multiple specs will result in an error.


Counter-Example Generation
--------------------------
    -nce | --no-counter-examples     Disables counter-example generation in case of incorrect circuit.

Sub-Circuit size
--------------------------
    -f <int>                         Non-negative value for fanout size, 0 turns fanout limit off (default value: 4).
    -d <int>                         Positive value for depth (default: 2).

Ablation
--------------------------
    -npp  | --no-preprocessing        Disables the preprocessing phase. (no rewriting of AIG).
    -nvc  | --no-vanishing            Turns on vanishing constraints 
    -nch  | --no-caching              Turns off caching of circuits 
    -dll  | --do-local-linearization  Enables the local linearization and only uses FGLM to linearize.
    -alg  | --algebraic-reduction     Use algebraic reductions instead of SAT in guess and proof
    -gap  | --force-guessing          Forces the linearization to only use guess-and-proof
    -fglm | --force-fglm              Forces the linearization to only use fglm


Verbosity Levels
----------------
    Control the level of output detail:
    -v0                   Minimal output (silent mode).
    -v1                   Low verbosity (default).
    -v2                   Medium verbosity.
    -v3                   High verbosity.
    -v4                   Maximum verbosity (debug-level).

Example Usages
-------------
    talisman examples/abc4.aig -mult-spec
    talisman examples/abc4.aig -mult-spec -nch -fglm

