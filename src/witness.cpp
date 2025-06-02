/*------------------------------------------------------------------------*/
/*! \file witness.cpp
    \brief contains functions used in the polynomial solver

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/

#include "witness.h"
/*------------------------------------------------------------------------*/

// ERROR CODES:
static int err_writing = 61; // cannot write to
static int err_witness = 62; // cannot write to
/*------------------------------------------------------------------------*/

bool check_inputs_only(const Polynomial *p)
{
    for (size_t i = 0; i < p->len(); i++)
    {
        Monomial *m = p->get_mon(i);
        if (!m->get_term())
        {
            continue;
        }
        else
        {
            const Var *v = m->get_term()->get_var();
            Gate *n = gate(v->get_num());
            if (!n->get_input())
                return 0;
        }
    }
    return 1;
}

/*------------------------------------------------------------------------*/

static void write_witness_vector(const Term *t, FILE *file)
{
    msg_nl("");

    if (ainc == 2)
    {
        for (unsigned i = 0; i <= NN / 2 - 1; i++)
        {
            Var *v = gates[a0 + i * ainc]->get_var();

            if (t->contains(v))
            {
                fprintf(file, "1");
                fprintf(stdout, "%s = ", v->get_name());
            }
            else
            {
                fprintf(file, "0");
            }

            Var *w = gates[b0 + i * binc]->get_var();

            if (t->contains(w))
            {
                fprintf(file, "1");
                fprintf(stdout, "%s = ", w->get_name());
            }
            else
            {
                fprintf(file, "0");
            }
        }
    }
    else if (ainc == 1)
    {
        for (unsigned i = 0; i <= NN / 2 - 1; i++)
        {
            Var *v = gates[a0 + i * ainc]->get_var();

            if (t->contains(v))
            {
                fprintf(file, "1");
                fprintf(stdout, "%s = ", v->get_name());
            }
            else
            {
                fprintf(file, "0");
            }
        }

        for (unsigned i = 0; i <= NN / 2 - 1; i++)
        {
            Var *w = gates[b0 + i * binc]->get_var();

            if (t->contains(w))
            {
                fprintf(file, "1");
                fprintf(stdout, "%s = ", w->get_name());
            }
            else
            {
                fprintf(file, "0");
            }
        }
    }

    fprintf(stdout, "1, all other inputs = 0;\n");
    fprintf(file, "\n");
}

/*------------------------------------------------------------------------*/

void write_witnesses(const Polynomial *p, FILE *file)
{
    assert(check_inputs_only(p));

    int len = p->min_term_size();
    if (len == 0)
    {
        msg("  all inputs = 0;\n");
        for (unsigned i = 0; i <= NN / 2 - 1; i++)
            fprintf(file, "00");

        fprintf(file, "\n");
    }
    else
    {
        for (size_t i = 0; i < p->len(); i++)
        {
            Monomial *m = p->get_mon(i);
            if (m->get_term())
            {
                int tlen = m->get_term_size();
                if (tlen == len)
                    write_witness_vector(m->get_term(), file);
            }
        }
    }
    fprintf(file, ".");
}

/*------------------------------------------------------------------------*/

void generate_witness(const Polynomial *p, const char *name)
{
    if (!check_inputs_only(p))
        die(err_witness, "cannot generate witness, as remainder polynomial contains non-inputs");

#define WITNESSLEN 100
    char witness_name[WITNESSLEN];
    memset(witness_name, '\0', sizeof(witness_name));
    for (int i = 0; name[i] != '.'; i++)
    {
        witness_name[i] = name[i];
    }
    snprintf(witness_name + strlen(witness_name), WITNESSLEN - strlen(witness_name), "%s", ".cex");

    FILE *witness_file;
    if (!(witness_file = fopen(witness_name, "w")))
        die(err_writing, "cannot write output to '%s'", witness_name);

    msg("");
    msg("COUNTER EXAMPLES ARE: ");

    write_witnesses(p, witness_file);

    msg("");
    msg("");
    msg("Counter examples are written to %s", witness_name);
    msg("You can run 'aigsim' from the AIGER library (http://fmv.jku.at/aiger/)");
    msg("to simulate the provided counter example(s).");
    msg("");
    msg("Note: 'aiger/aigsim %s %s' produces output in the form:", name, witness_name);
    if (ainc == 2)
    {
        msg_nl(" ");
        if (NN == 2)
            fprintf(stdout, "  a[0]b[0]  s[0]\n");
        else if (NN == 4)
            fprintf(stdout, "  a[0]b[0]a[1]b[1]  s[0]s[1]s[2]s[3]\n");
        else
            fprintf(stdout, "  a[0]b[0]a[1]b[1]...a[%u]b[%u]  s[0]s[1]s[2]...s[%u]\n", NN / 2 - 1, NN / 2 - 1, NN - 1);
    }
    else
    {
        msg_nl(" ");
        if (NN == 2)
            fprintf(stdout, "  a[0]b[0]  s[0]\n");
        else if (NN == 4)
            fprintf(stdout, "  a[0]a[1]b[0]b[1]  s[0]s[1]s[2]s[3]\n");
        else
            fprintf(stdout, "  a[0]a[1]...a[%u]b[0]b[1]...b[%u]  s[0]s[1]s[2]...s[%u]\n", NN / 2 - 1, NN / 2 - 1,
                    NN - 1);
    }

    fclose(witness_file);
}
