// This file exists only to assert some assumed properties on solver types
// that were used to optimize Minisat::vec<T> for them.
// If it fails to compile, it is not neccessarily wrong, just weird.

#include "minisat/core/SolverTypes.h"

#include <type_traits>

namespace Minisat {
    static_assert(std::is_trivially_destructible<Lit>::value, "Minisat::Lit is expected to be trivially destructible");
    static_assert(std::is_trivially_destructible<Var>::value, "Minisat::Var is expected to be trivially destructible");
    static_assert(std::is_trivially_destructible<lbool>::value, "Minisat::lbool is expected to be trivially destructible");
    static_assert(std::is_trivially_destructible<CRef>::value, "Minisat::CRef is expected to be trivially destructible");
}
