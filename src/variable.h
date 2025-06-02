/*------------------------------------------------------------------------*/
/*! \file variable.h
    \brief contains the class Var

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_VARIABLE_H_
#define TALISMAN_SRC_VARIABLE_H_
/*------------------------------------------------------------------------*/
#include <cstdlib>
#include <string>

#include "hash_val.h"
/*------------------------------------------------------------------------*/

/** \class Var
    represents a variable is assigned to a gate(see <gate.h>) and is
    used to represent variables in terms, see <term.h>
*/

class Var
{
    // / name of variable
    const std::string name;

    // / Hash value of variables, used for storing terms
    int hash;

    // / Increasing value that indicates the order of the variable
    int level = 0;

    int id = 0;

    // / corresponding value used to relate AIG gates to Gate class
    const int num;

    // / links dual variable
    Var *dual;

    // / True if Variable is a dual var
  const bool d = 0;

  int value = -1;

  public:
    /** Constructor

       @param name_ name
       @param level_ level
       @param num_ num, default is 0

    */
    Var(std::string name_, int level_, int num_ = 0, bool dual_ = 0) : name(name_), level(level_), num(num_), d(dual_)
    {
        hash = hash_string(name_);
    }

    /** Getter for member name, and converts string to char*

        @return const char *
    */
    const char *get_name() const
    {
        return name.c_str();
    }

    /** Getter for member hash

        @return integer
    */
    int get_hash() const
    {
        return hash;
    }

    int get_id() const
    {
        return id;
    }

    void set_id(int i)
    {
        id = i;
    }

    /** Getter for member d

        @return bool
    */
    bool is_dual() const
    {
        return d;
    }

    /** Getter for member level

        @return integer
    */
    int get_level() const
    {
        return level;
    }

    /** Setter for member level
        @param l int
    */
    void set_level(int l)
    {
        level = l;
    }

    /** Getter for member num

        @return integer
    */
    int get_num() const
    {
        return num;
    }

    /** Getter for member dual

      @return integer
    */
    Var *get_dual() const
    {
        return dual;
    }

    /** Setter for dual variable
        @param d Variable *
    */
    void set_dual_var(Var *d)
    {
        dual = d;
    }

  bool has_value() { return value >= 0; }

  int get_value() { return value; }

  void set_value(int v) { value = v; }
};

#endif // TALISMAN_SRC_VARIABLE_H_
