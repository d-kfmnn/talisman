/*------------------------------------------------------------------------*/
/*! \file specpoly.cpp
    \brief Used to generate specification polynomials

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#include "specpoly.h"
/*------------------------------------------------------------------------*/
// Global var

FILE *parse_file;
unsigned lineno;
unsigned charno;
unsigned lineno_at_start_of_last_token;
/*------------------------------------------------------------------------*/
typedef const char *Token;

/**
    Determines whether char is either character, number of '_'

    @param ch char

    @return bool

*/
static bool is_valid_variable_letter(char ch)
{
    // if (ch == '_') return 1;
    if ('0' <= ch && ch <= '9')
        return 1;
    //  if ('a' <= ch && ch <= 'z') return 1;
    //  if ('A' <= ch && ch <= 'Z') return 1;
    return 0;
}

/*------------------------------------------------------------------------*/
/**
    Determines whether char is a character

    @param ch char

    @return bool

*/
static bool is_valid_variable_first_letter(char ch)
{
    if ('a' <= ch && ch <= 'z')
        return 1;
    if ('A' <= ch && ch <= 'Z')
        return 1;
    return 0;
}

/*------------------------------------------------------------------------*/
#ifndef NDEBUG
/**
    Determines whether name is a valid variable name

    @param name const char *

    @return bool

*/
static bool is_valid_variable_name(const char *name)
{
    if (!name)
        return 0;
    if (!is_valid_variable_first_letter(name[0]))
        return 0;
    for (const char *p = name + 1; *p; p++)
        if (!is_valid_variable_letter(*p))
            return 0;
    return 1;
}

#endif

static Var *var_from_string_via_gate(const char *name)
{
    if (!name)
        return 0;
    std::string str(name);

    for (unsigned i = 0; i < num_gates; i++)
    {
        Gate *g = gates[i];
        std::string g_name(g->get_var_name());
        if (g_name == str)
            return g->get_var();
    }

    die(1,
        "Variable %s from specification not contained in input AIG! \n"
        "   Following formats are assumed for output variables 's<num>', internal variables 'l<num>', and primary "
        "inputs 'i<num>', \n"
        "   where <num> corresponds to the number in the AIG",
        str.c_str());
    return 0;
}

/*------------------------------------------------------------------------*/
/// size of the buffer
static size_t size_buffer;
/// number of elements in the buffer
static size_t num_buffer;
/// buffer to store char
static char *buffer;

/*------------------------------------------------------------------------*/

/**
    Enlarges the allocated buffer
*/
static void enlarge_buffer()
{
    size_t new_size_buffer = size_buffer ? 2 * size_buffer : 1;
    buffer = reinterpret_cast<char *>(realloc(buffer, new_size_buffer));
    size_buffer = new_size_buffer;
}

/*------------------------------------------------------------------------*/
/**
    Pushes a character to the buffer
*/
static void push_buffer(char ch)
{
    if (size_buffer == num_buffer)
        enlarge_buffer();
    buffer[num_buffer++] = ch;
}

/*------------------------------------------------------------------------*/
/**
   Clears the buffer, i.e. sets the number to 0
*/
static void clear_buffer()
{
    num_buffer = 0;
}

/*------------------------------------------------------------------------*/

void deallocate_buffer()
{
    free(buffer);
}

/*------------------------------------------------------------------------*/
/// currently saved char
static int saved_char;
/// determines whether we currently saved a char
static bool char_saved;
/*------------------------------------------------------------------------*/
/**
    stores the next character

    @return int
*/
static int next_char()
{
    int res;
    if (char_saved)
    {
        res = saved_char;
        char_saved = 0;
    }
    else
    {
#ifdef HAVEUNLOCKEDIO
        res = getc_unlocked(parse_file);
#else
        res = getc(parse_file);
#endif
    }
    if (res == '\n')
        lineno++;
    if (res != EOF)
        charno++;
    return res;
}

/*------------------------------------------------------------------------*/
/**
    In case we need to undo reading a character, we store the current character

    @param ch int
*/
static void prev_char(int ch)
{
    assert(!char_saved);
    if (ch == '\n')
    {
        assert(lineno > 0);
        lineno--;
    }
    else if (ch != EOF)
    {
        assert(charno > 0);
        charno--;
    }
    saved_char = ch;
    char_saved = 1;
}

/*------------------------------------------------------------------------*/

/// We use Token to identify components in the proof
typedef const char *Token;

/// Token, used to store the components in the proof
static Token token;

static Token END_OF_FILE_TOKEN = "end-of-file";
static Token MINUS_TOKEN = "minus operator";
static Token PERCENT_TOKEN = "linear combination operator";
static Token PLUS_TOKEN = "addition operator";
static Token MULTIPLY_TOKEN = "multiplication operator";
static Token COMMA_TOKEN = "comma separator";
static Token SEMICOLON_TOKEN = "semicolon separator";
static Token NUMBER_TOKEN = "number";
static Token VARIABLE_TOKEN = "variable";
static Token EQUALITY_TOKEN = "equal";
static Token EXPONENT_TOKEN = "exponent";
static Token L_PARENTHESIS_TOKEN = "open parenthesis";
static Token R_PARENTHESIS_TOKEN = "close parenthesis";
static Token L_SQUARE_BRACKET_TOKEN = "open bracket";
static Token R_SQUARE_BRACKET_TOKEN = "close bracket";
/*------------------------------------------------------------------------*/

bool is_delete_token()
{
    if (!buffer)
        return 0;
    if (buffer[0] != 'd')
        return 0;
    if (buffer[1])
        return 0;
    return 1;
}
/*------------------------------------------------------------------------*/

static bool is_separator_token()
{
    if (token == COMMA_TOKEN)
        return 1;
    if (token == SEMICOLON_TOKEN)
        return 1;
    return 0;
}

static void parse_error(const char *msg, ...);
/*------------------------------------------------------------------------*/
Token get_token()
{
    return token;
}
/*------------------------------------------------------------------------*/
/**
    stores a new token

    @param t token
*/
static Token new_token(Token t)
{
    push_buffer(0);
    token = t;
    return token;
}

/*------------------------------------------------------------------------*/

Token next_token()
{
    clear_buffer();
    for (;;)
    {
        int ch = next_char();
        if (ch == ' ' || ch == '\t' || ch == '\r' || ch == '\n')
            continue;
        lineno_at_start_of_last_token = lineno;
        if (ch == EOF)
            return new_token(END_OF_FILE_TOKEN);
        push_buffer(ch);

        if ('0' <= ch && ch <= '9')
        {
            while ('0' <= (ch = next_char()) && ch <= '9')
                push_buffer(ch);
            prev_char(ch);
            return new_token(NUMBER_TOKEN);
        }
        if (is_valid_variable_first_letter(ch))
        {
            while (is_valid_variable_letter(ch = next_char()))
                push_buffer(ch);
            prev_char(ch);
            return new_token(VARIABLE_TOKEN);
        }
        if (ch == '-')
            return new_token(MINUS_TOKEN);
        if (ch == '+')
            return new_token(PLUS_TOKEN);
        if (ch == '*')
            return new_token(MULTIPLY_TOKEN);
        if (ch == '%')
            return new_token(PERCENT_TOKEN);
        if (ch == '^')
            return new_token(EXPONENT_TOKEN);
        if (ch == ',')
            return new_token(COMMA_TOKEN);
        if (ch == ';')
            return new_token(SEMICOLON_TOKEN);
        if (ch == '=')
            return new_token(EQUALITY_TOKEN);
        if (ch == '(')
            return new_token(L_PARENTHESIS_TOKEN);
        if (ch == ')')
            return new_token(R_PARENTHESIS_TOKEN);
        if (ch == '[')
            return new_token(L_SQUARE_BRACKET_TOKEN);
        if (ch == ']')
            return new_token(R_SQUARE_BRACKET_TOKEN);
        if (isprint(ch))
            parse_error("invalid character");
        else
            parse_error("invalid character code 0x%02x", ch);
    }
}

/*------------------------------------------------------------------------*/
bool is_semicolon_token()
{
    return token == SEMICOLON_TOKEN;
}
bool is_comma_token()
{
    return token == COMMA_TOKEN;
}
bool is_plus_token()
{
    return token == PLUS_TOKEN;
}
bool is_multiply_token()
{
    return token == MULTIPLY_TOKEN;
}
bool is_open_parenthesis_token()
{
    return token == L_PARENTHESIS_TOKEN;
}
bool is_close_parenthesis_token()
{
    return token == R_PARENTHESIS_TOKEN;
}
bool is_equality_token()
{
    return token == EQUALITY_TOKEN;
}
bool is_lin_combi_token()
{
    return token == PERCENT_TOKEN;
}
bool is_eof_token()
{
    return token == END_OF_FILE_TOKEN;
}
bool following_token_is_EOF()
{
    return next_token() == END_OF_FILE_TOKEN;
}

/*------------------------------------------------------------------------*/

void parse_error(const char *msg, ...)
{
    fflush(stdout);
    fprintf(stderr, "*** parse error in line %i", lineno_at_start_of_last_token);
    if (buffer[0] && isprint(buffer[0]))
        fprintf(stderr, " at '%s'", buffer);
    else if (token == END_OF_FILE_TOKEN)
        fputs(" at end-of-file", stderr);
    fputs(": ", stderr);
    va_list ap;
    va_start(ap, msg);
    vfprintf(stderr, msg, ap);
    va_end(ap);
    fputc('\n', stderr);
    fflush(stderr);
    exit(1);
}

/*------------------------------------------------------------------------*/

/**
    parses a variable, a new variable is only allocated if new_var_allowed = 1

    @param new_var_allowed bool
*/
static Var *parse_variable()
{
    assert(is_valid_variable_name(buffer));

    Var *res = var_from_string_via_gate(buffer);

    return res;
}

/*------------------------------------------------------------------------*/

/**
    parses a term, a new variable is only allocated if new_var_allowedd = 1

    @param new_var_allowed bool
*/
static Term *parse_term()
{

    while (token == VARIABLE_TOKEN)
    {
        Var *variable = parse_variable();
        assert(variable);

        add_to_vstack(variable);
        next_token();
        if (token == MULTIPLY_TOKEN)
            next_token();
        if (token == EXPONENT_TOKEN)
            die(1, "exponents currently not supported");
    }
    Term *res = build_term_from_stack(1);
    return res;
}

/*------------------------------------------------------------------------*/

/**
    parses a monomial, a new variable is only allocated if new_var_allowed = 1

    @param sign int
    @param new_var_allowed bool
*/
static Monomial *parse_monomial(bool sign)
{
    mpz_t tmp_gmp;
    mpz_init(tmp_gmp);
    if (token == NUMBER_TOKEN)
    {
        mpz_set_str(tmp_gmp, buffer, 10);
        next_token();
    }
    else if (token == VARIABLE_TOKEN)
    {
        mpz_set_ui(tmp_gmp, 1);
    }
    else
    {
        parse_error("expected monomial");
    }
    if (sign)
        mpz_neg(tmp_gmp, tmp_gmp);
    if (token == MULTIPLY_TOKEN)
        next_token();
    Term *term = parse_term();

    Monomial *res = new Monomial(tmp_gmp, term);
    mpz_clear(tmp_gmp);
    return res;
}

/*------------------------------------------------------------------------*/

static Polynomial *parse_polynomial()
{
    next_token();
    bool sign;
    if (token == MINUS_TOKEN)
    {
        next_token();
        if (token == NUMBER_TOKEN && buffer[0] == '0')
            parse_error("unexpected '0' after '-'");
        sign = 1;
    }
    else
    {
        sign = 0;
    }
    for (;;)
    {
        Monomial *monomial = parse_monomial(sign);
        push_mstack(monomial);
        if (is_separator_token() || is_eof_token())
            break;
        if (token == MINUS_TOKEN)
        {
            sign = 1;
            next_token();
        }
        else if (token == PLUS_TOKEN)
        {
            sign = 0;
            next_token();
        }
        else
        {
            parse_error("unexpected %s", token);
        }
    }

    Polynomial *res = build_poly();
    return res;
}

/*------------------------------------------------------------------------*/

Polynomial *parse_specification_polynomial(const char *file_name)
{

    parse_file = fopen(file_name, "r");
    if (!parse_file)
        die(1, "can not open '%s' for reading", file_name);

    lineno = 1;
    charno = 0;

    // msg("reading specification polynomial from '%s'", file_name);

    Polynomial *spec = parse_polynomial();
    assert(is_separator_token() || is_eof_token());

    if (fclose(parse_file))
        die(1, "failed to close '%s'", file_name);
    if (verbose >= 3)
        msg("read %i bytes from '%s'", charno, file_name);

    return spec;
}

/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
Polynomial *mult_spec_poly()
{
    mpz_t coeff;
    mpz_init(coeff);

    // outputs
    for (int i = NN - 1; i >= 0; i--)
    {
        Var *v = gates[i + M - 1]->get_var();

        mpz_pow_ui(coeff, base, i);
        mpz_neg(coeff, coeff);

        Term *t1 = new_term(v);
        Monomial *m1 = new Monomial(coeff, t1);
        push_mstack(m1);
    }

    // partial products
    for (int i = NN / 2 - 1; i >= 0; i--)
    {
        Var *a = gates[a0 + i * ainc]->get_var();

        for (int j = NN / 2 - 1; j >= 0; j--)
        {
            Var *b = gates[b0 + j * binc]->get_var();
            mpz_pow_ui(coeff, base, i + j);
            add_to_vstack(b);
            add_to_vstack(a);
            Term *t1 = build_term_from_stack();
            Monomial *m1 = new Monomial(coeff, t1);
            push_mstack(m1);
        }
    }
    mpz_clear(coeff);

    Polynomial *p = build_poly();
    return p;
}
/*------------------------------------------------------------------------*/
Polynomial *miter_spec_poly()
{

    assert(MM == 1);

    // outputs
    Var *v = gates[M - 1]->get_var();

    Term *t1 = new_term(v);
    Monomial *m1 = new Monomial(one, t1);
    push_mstack(m1);

    Polynomial *p = build_poly();
    return p;
}
/*------------------------------------------------------------------------*/
Polynomial *assertion_spec_poly()
{

    // outputs
    for (int i = MM - 1; i >= 0; i--)
    {
        Var *v = gates[i + M - 1]->get_var();

        Monomial *m = new Monomial(one, new_term(v));
        push_mstack(m);
    }

    mpz_t coeff;
    mpz_init(coeff);
    int output_num = -1 * MM;
    mpz_set_si(coeff, output_num);
    Monomial *m = new Monomial(coeff, 0);
    push_mstack(m);
    mpz_clear(coeff);

    Polynomial *p = build_poly();
    return p;
}