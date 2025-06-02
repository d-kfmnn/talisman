/*------------------------------------------------------------------------*/
/*! \file polynomial.cpp
    \brief contains arithmetic operations for polynomials

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#include "polynomial.h"
/*------------------------------------------------------------------------*/


Polynomial::Polynomial() {}
Polynomial::Polynomial(Monomial** m, size_t len, int d)
  : mon(m)
  , num_mon(len)
  , deg(d) {}

Monomial*
Polynomial::get_mon(size_t i) const {
  if(i < num_mon)
    return mon[i];
  else
    return 0;
}

/*------------------------------------------------------------------------*/

Polynomial*
Polynomial::copy() {

 for (size_t i = 0 ; i < len(); i++) {
    Monomial * m = get_mon(i);
    push_mstack_end(m->copy());
  }
  Polynomial * out = build_poly();
  out->set_idx(idx);
  return out;
}
/*------------------------------------------------------------------------*/


void
Polynomial::print(FILE* file, bool end) const {
  if(len() == 0) {
    #ifdef HAVEUNLOCKEDIO
      fputs_unlocked("0", file);
    #else
      fputs("0", file);
    #endif
    } else {
    for(size_t i = 0; i < len(); i++) {
      Monomial* m = get_mon(i);
      if(i == 0)
        m->print(file, 1);
      else
        m->print(file, 0);
    }
  }
  if(end)
    fputs(";\n", file);
}

/*------------------------------------------------------------------------*/
Polynomial::~Polynomial() {

  for(size_t i = 0; i < len(); i++) {
    Monomial* m = get_mon(i);
    deallocate_monomial(m);
  }
  delete[] mon;
}

/*------------------------------------------------------------------------*/

bool
Polynomial::is_constant_zero_poly() const {
  return len() == 0;
}

/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
bool
Polynomial::is_constant_one_poly() const {
  if(len() != 1)
    return 0;

  Monomial* m = get_mon(0);
  if(m->get_term())
    return 0;
  if(mpz_cmp_si(m->coeff, 1) != 0)
    return 0;

  return 1;
}
/*------------------------------------------------------------------------*/
Var*
Polynomial::contains_dual_var() const {
  for(size_t i = 0; i < len(); i++) {
    Term* t = get_mon(i)->get_term();
    if(!t)
      return 0;
     Var* v = t->extract_first_dual_var();
    if(v)
      return v;
  }
  return 0;
}
/*------------------------------------------------------------------------*/
Monomial*
Polynomial::get_largest_mon() const {
  for(size_t i = 0; i < len(); i++) {
    Term* t = get_mon(i)->get_term();
    if(t->degree() == degree())
      return get_mon(i);
  }
  return 0;
}
/*------------------------------------------------------------------------*/
Polynomial*
Polynomial::get_tail_poly() const {
  for(size_t i = 1; i < len(); i++) {
    push_mstack_end(get_mon(i));
  }
  return build_poly();
}

/*------------------------------------------------------------------------*/
Term*
Polynomial::get_largest_term() const {
  for(size_t i = 0; i < len(); i++) {
    Term* t = get_mon(i)->get_term();
    if(t->degree() == degree())
      return t;
  }
  return 0;
}

/*------------------------------------------------------------------------*/
int
Polynomial::min_term_size() const {
  int length = INT_MAX;

  for(size_t i = 0; i < len(); i++) {
    Monomial* m = get_mon(i);

    int tlen = 0;
    if(m->get_term())
      tlen = m->get_term_size();
    if(tlen < length)
      length = tlen;
  }
  return length;
}

/*------------------------------------------------------------------------*/
// Local variables
static size_t size_mstack;   // /< size of mstack
static size_t num_mstack = 0;// /< number of elements in mstack
static Monomial** mstack;    // /< Monomial** used for building poly

/*------------------------------------------------------------------------*/

static void
enlarge_mstack() {
  size_t new_size_mstack = size_mstack ? 2 * size_mstack : 1;

  Monomial** newArr = new Monomial*[new_size_mstack];
  if(size_mstack > 0)  memcpy(newArr, mstack, size_mstack * sizeof(Monomial*));
  delete[] mstack;
  mstack = newArr;
  size_mstack = new_size_mstack;
}

/*------------------------------------------------------------------------*/

static void
clear_mstack() {
  num_mstack = 0;
  size_mstack = 0;
  mstack = 0;
}

/*------------------------------------------------------------------------*/

void
deallocate_mstack() {
  delete[](mstack);
}

/*------------------------------------------------------------------------*/
void
push_mstack_end(Monomial* m) {
  if(size_mstack == num_mstack)
    enlarge_mstack();

  assert(m);
  if(mpz_sgn(m->coeff) == 0) {
    deallocate_monomial(m);
    return;
  }

  mstack[num_mstack++] = m;
}

/*------------------------------------------------------------------------*/

void
push_mstack(Monomial* m) {
  assert(m);
  if(mpz_sgn(m->coeff) == 0) {
    deallocate_monomial(m);
    return;
  }

  if(size_mstack == num_mstack)
    enlarge_mstack();
  if(num_mstack == 0) {
    mstack[num_mstack++] = m;
    return;
  }

  if(!m->get_term()) {
    Monomial* tmp = mstack[num_mstack - 1];

    if(tmp->get_term()) {
      mstack[num_mstack++] = m;
    } else {
      mpz_t coeff;
      mpz_init(coeff);
      mpz_add(coeff, tmp->coeff, m->coeff);
      deallocate_monomial(m);
      deallocate_monomial(tmp);

      if(mpz_sgn(coeff) != 0)
        mstack[num_mstack - 1] = new Monomial(coeff, 0);
      else
        --num_mstack;

      mpz_clear(coeff);
    }
  } else {
    assert(num_mstack > 0);
    int i = num_mstack - 1;
    int cmp = -1;
    Monomial* tmp = 0;

    while(i >= 0) {
      tmp = mstack[i];
      cmp = cmp_term(tmp->get_term(), m->get_term());

      if(cmp >= 0)
        break;
      i--;
    }

    if(cmp == 0) {
      mpz_t coeff;
      mpz_init(coeff);
      mpz_add(coeff, tmp->coeff, m->coeff);

      if(mpz_sgn(coeff) == 0) {
        for(unsigned j = i; j < num_mstack - 1; j++)
          mstack[j] = mstack[j + 1];
        num_mstack--;
      } else {
        mstack[i] = new Monomial(coeff, m->get_term_copy());
      }
      deallocate_monomial(m);
      deallocate_monomial(tmp);
      mpz_clear(coeff);
    } else {
      for(int j = num_mstack; j > i + 1; j--)
        mstack[j] = mstack[j - 1];
      mstack[i + 1] = m;
      num_mstack++;
    }
  }
}

/*------------------------------------------------------------------------*/
size_t running_idx = 1;

Polynomial*
build_poly() {

  if(num_mstack == 0)
    return 0;
  size_t deg = 0, tmp_deg = 0;

  for(unsigned j = 0; j < num_mstack; j++) {
    Term* tmp = mstack[j]->get_term();
    if(tmp)
      tmp_deg = tmp->degree();
    if(deg < tmp_deg)
      deg = tmp_deg;
  }

  Polynomial* res = new Polynomial(mstack, num_mstack, deg);
  running_idx++;
  res->set_idx(running_idx);
  clear_mstack();
  return res;
}

/*------------------------------------------------------------------------*/
Polynomial*
gen_dual_constraint( Var* v) {
  Var* d = v->get_dual();

  push_mstack_end(new Monomial(minus_one, new_term(v)));
  push_mstack_end(new Monomial(minus_one, new_term(d)));
  push_mstack_end(new Monomial(one, 0));

  Polynomial* p = build_poly();
  return p;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
bool
equal_poly(Polynomial* p1, Polynomial* p2) {
  if(p1->len() != p2->len())
    return 0;
  if(p1->degree() != p2->degree())
    return 0;

  for(size_t i = 0; i < p1->len(); i++) {
    if(p1->get_mon(i)->get_term() != p2->get_mon(i)->get_term())
      return 0;
    if(mpz_cmp(p1->get_mon(i)->coeff, p2->get_mon(i)->coeff) != 0)
      return 0;
  }
  return 1;
}

/*------------------------------------------------------------------------*/
int
cmp_poly(Polynomial* p1, Polynomial* p2) {

  if(equal_poly(p1, p2))
    return 0;

  size_t i = 0;

  Monomial* m1 = p1->get_mon(i);
  Monomial* m2 = p2->get_mon(i);
  mpz_t coeff;
  mpz_init(coeff);

  while(i < p1->len() && i < p2->len()) {
    if(!m1->get_term())
      return -1;
    if(!m2->get_term())
      return 1;

    int cmp = cmp_term(m1->get_term(), m2->get_term());
    if(cmp)
      return cmp;
    else {
      i++;
      m1 = p1->get_mon(i);
      m2 = p2->get_mon(i);
    }
  }

  if(i < p1->len())
    return 1;
  return -1;
}
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/

Polynomial*
add_poly(Polynomial* p1, Polynomial* p2) {
  if(!p1)
    return p2->copy();
  if(!p2)
    return p1->copy();
  assert(p1);
  assert(p2);

  size_t i = 0, j = 0;

  Monomial* m1 = p1->get_mon(i);
  Monomial* m2 = p2->get_mon(j);
  mpz_t coeff;
  mpz_init(coeff);

  while(i < p1->len() && j < p2->len()) {
    if(!m1->get_term() || !m2->get_term()) {
      if(!m1->get_term() && !m2->get_term()) {
        mpz_add(coeff, m1->coeff, m2->coeff);
        if(mpz_sgn(coeff) != 0) {
          Monomial* m = new Monomial(coeff, 0);
          push_mstack_end(m);
        }
        m1 = p1->get_mon(++i);
        m2 = p2->get_mon(++j);
      } else if(!m1->get_term()) {
        push_mstack_end(m2->copy());
        m2 = p2->get_mon(++j);
      } else {
        push_mstack_end(m1->copy());
        m1 = p1->get_mon(++i);
      }
    } else {
      int cmp = cmp_term(m1->get_term(), m2->get_term());
      if(cmp == 1) {
        push_mstack_end(m1->copy());
        m1 = p1->get_mon(++i);
      } else if(cmp == -1) {
        push_mstack_end(m2->copy());
        m2 = p2->get_mon(++j);
      } else {
        mpz_add(coeff, m1->coeff, m2->coeff);
        if(mpz_sgn(coeff) != 0) {
          Monomial* m = new Monomial(coeff, m1->get_term_copy());
          push_mstack_end(m);
        }
        m1 = p1->get_mon(++i);
        m2 = p2->get_mon(++j);
      }
    }
  }
  mpz_clear(coeff);

  while(i < p1->len()) {
    push_mstack_end(m1->copy());
    m1 = p1->get_mon(++i);
  }
  while(j < p2->len()) {
    push_mstack_end(m2->copy());
    m2 = p2->get_mon(++j);
  }

  Polynomial* p = build_poly();
  return p;
}

/*------------------------------------------------------------------------*/
Polynomial*
sub_poly(Polynomial* p1, Polynomial* p2) {
  Polynomial* tmp = multiply_poly_with_constant(p2, minus_one);
  Polynomial* add = add_poly(tmp, p1);
  delete(tmp);
  return add;
}

/*------------------------------------------------------------------------*/

Polynomial*
multiply_poly(Polynomial* p1, Polynomial* p2) {
  if(!p1 || !p2)
    return 0;
  assert(p1);
  assert(p2);

  mpz_t coeff;
  mpz_init(coeff);

  for(size_t i = 0; i < p1->len(); i++) {
    Monomial* m1 = p1->get_mon(i);
    for(size_t j = 0; j < p2->len(); j++) {
      Monomial* m2 = p2->get_mon(j);

      push_mstack(multiply_monomial(m1, m2));
    }
  }
  Polynomial* p = build_poly();
  mpz_clear(coeff);
  return p;
}

/*------------------------------------------------------------------------*/

Polynomial*
multiply_poly_with_constant(Polynomial* p1, mpz_t c) {
  if(mpz_sgn(c) == 0)
    return 0;
  mpz_t coeff;
  mpz_init(coeff);

  for(size_t i = 0; i < p1->len(); i++) {
    Monomial* m = p1->get_mon(i);
    mpz_mul(coeff, m->coeff, c);
    if(m->get_term())
      push_mstack_end(new Monomial(coeff, m->get_term_copy()));
    else
      push_mstack_end(new Monomial(coeff, 0));
  }
  Polynomial* tmp = build_poly();
  mpz_clear(coeff);
  return tmp;
}

/*------------------------------------------------------------------------*/
Polynomial*
multiply_poly_with_term(Polynomial* p1, Term* t) {
  if(!t)
    return p1->copy();
  if(!p1) return 0;

  for(size_t i = 0; i < p1->len(); i++) {
    Term* t1 = p1->get_mon(i)->get_term();
    Term* res = multiply_term(t1, t);
    if(!t1) {
      res = t->copy();

    }    
    push_mstack_end(new Monomial(p1->get_mon(i)->coeff, res));
  }
  Polynomial* tmp = build_poly();
  return tmp;
}

/*------------------------------------------------------------------------*/
Polynomial*
multiply_poly_with_monomial(Polynomial* p1, Monomial * m) {
  if(!m)
    return p1->copy();
  if(!p1) return 0;

  mpz_t coeff;
  mpz_init(coeff);

  for(size_t i = 0; i < p1->len(); i++) {
    Monomial* m1 = p1->get_mon(i);
    Term * t1 = m1->get_term();
    Term* res = multiply_term(t1, m->get_term());
    mpz_mul(coeff, m->coeff, m1->coeff);
    if(!t1) {
      res = m->get_term()->copy();

    }    
    push_mstack_end(new Monomial(coeff, res));
  }
  Polynomial* tmp = build_poly();
  mpz_clear(coeff);
  return tmp;
}

/*------------------------------------------------------------------------*/
Polynomial*
divide_by_var(Polynomial* p1, const Term* t) {
  assert(t->degree() == 1);
   Var* v = t->get_var();

  for(size_t i = 0; i < p1->len(); i++) {
    Monomial* lm_tmp = p1->get_mon(i);

    if(!lm_tmp->get_term())
      break;
    if(cmp_term(lm_tmp->get_term(), t) == -1)
      break;

    if(lm_tmp->get_term()->contains(v)) {
      Term* t_rem = divide_by_var(lm_tmp->get_term(), v);
      if(t_rem) {
        push_mstack_end(new Monomial(lm_tmp->coeff, t_rem->copy()));
      } else {
        push_mstack_end(new Monomial(lm_tmp->coeff, 0));
        break;
      }
    }
  }
  Polynomial* p = build_poly();
  return p;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/

Polynomial*
divide_poly_by_term(Polynomial* p1, const Term* t) {
  if(t->degree() == 1)
    return divide_by_var(p1, t);
  else {
    for(size_t i = 0; i < p1->len(); i++) {

      Monomial* lm_tmp = p1->get_mon(i);

      if(!lm_tmp->get_term())
        continue;

      if(lm_tmp->get_term()->contains_subterm(t)) {
        Term* t_rem = divide_by_term(lm_tmp->get_term(), t);
        push_mstack_end(new Monomial(lm_tmp->coeff, t_rem));
      }
    }
    return build_poly();
  }
}
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
mpz_t one;
mpz_t minus_one;
mpz_t minus_two;
mpz_t base;
mpz_t mod_coeff;

/*------------------------------------------------------------------------*/

void
init_mpz(unsigned exp) {
  mpz_init_set_ui(one, 1);
  mpz_init_set_si(minus_one, -1);
  mpz_init_set_si(minus_two, -2);
  mpz_init_set_ui(base, 2);

  mpz_init(mod_coeff);
  mpz_pow_ui(mod_coeff, base, exp);
}

/*------------------------------------------------------------------------*/

void
clear_mpz() {
  mpz_clear(one);
  mpz_clear(minus_one);
  mpz_clear(minus_two);
  mpz_clear(base);
  mpz_clear(mod_coeff);
}

/*------------------------------------------------------------------------*/
