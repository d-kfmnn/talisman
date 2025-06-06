/**************************************************************************************[Options.cc]
Copyright (c) 2008-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#include "minisat/mtl/Sort.h"
#include "minisat/utils/Options.h"
#include "minisat/utils/ParseUtils.h"

using namespace Minisat;

void Minisat::parseOptions(int& argc, char** argv, bool strict)
{
    int i, j;
    for (i = j = 1; i < argc; i++){
        const char* str = argv[i];
        if (match(str, "--") && match(str, Option::getHelpPrefixString()) && match(str, "help")){
            if (*str == '\0')
                printUsageAndExit(argc, argv);
            else if (match(str, "-verb"))
                printUsageAndExit(argc, argv, true);
        } else {
            bool parsed_ok = false;

            for (int k = 0; !parsed_ok && k < Option::getOptionList().size(); k++){
                parsed_ok = Option::getOptionList()[k]->parse(argv[i]);

                // fprintf(stderr, "checking %d: %s against flag <%s> (%s)\n", i, argv[i], Option::getOptionList()[k]->name, parsed_ok ? "ok" : "skip");
            }

            if (!parsed_ok) {
                if (strict && match(argv[i], "-"))
                    fprintf(stderr, "ERROR! Unknown flag \"%s\". Use '--%shelp' for help.\n", argv[i], Option::getHelpPrefixString()), exit(1);
                else
                    argv[j++] = argv[i];
            }
        }
    }

    argc -= (i - j);
}


void Minisat::setUsageHelp      (const char* str){ Option::getUsageString() = str; }
void Minisat::setHelpPrefixStr  (const char* str){ Option::getHelpPrefixString() = str; }

void Minisat::printUsageAndExit (int, char** argv, bool verbose)
{
    const char* usage = Option::getUsageString();
    if (usage != NULL)
        fprintf(stderr, usage, argv[0]);

    sort(Option::getOptionList(), Option::OptionLt());

    const char* prev_cat  = NULL;
    const char* prev_type = NULL;

    for (int i = 0; i < Option::getOptionList().size(); i++){
        const char* cat  = Option::getOptionList()[i]->category;
        const char* type = Option::getOptionList()[i]->type_name;

        if (cat != prev_cat)
            fprintf(stderr, "\n%s OPTIONS:\n\n", cat);
        else if (type != prev_type)
            fprintf(stderr, "\n");

        Option::getOptionList()[i]->help(verbose);

        prev_cat  = Option::getOptionList()[i]->category;
        prev_type = Option::getOptionList()[i]->type_name;
    }

    fprintf(stderr, "\nHELP OPTIONS:\n\n");
    fprintf(stderr, "  --%shelp        Print help message.\n", Option::getHelpPrefixString());
    fprintf(stderr, "  --%shelp-verb   Print verbose help message.\n", Option::getHelpPrefixString());
    fprintf(stderr, "\n");
    exit(0);
}

Minisat::Option::~Option() = default;

bool Minisat::DoubleOption::parse(const char * str) {
    const char* span = str;

    if (!match(span, "-") || !match(span, name) || !match(span, "="))
        return false;

    char*  end;
    double tmp = strtod(span, &end);

    if (end == NULL)
        return false;
    else if (tmp >= range.end && (!range.end_inclusive || tmp != range.end)) {
        fprintf(stderr, "ERROR! value <%s> is too large for option \"%s\".\n", span, name);
        exit(1);
    } else if (tmp <= range.begin && (!range.begin_inclusive || tmp != range.begin)) {
        fprintf(stderr, "ERROR! value <%s> is too small for option \"%s\".\n", span, name);
        exit(1);
    }

    value = tmp;
    // fprintf(stderr, "READ VALUE: %g\n", value);

    return true;
}

void Minisat::DoubleOption::help(bool verbose) {
    fprintf(stderr, "  -%-12s = %-8s %c%4.2g .. %4.2g%c (default: %g)\n",
            name, type_name,
            range.begin_inclusive ? '[' : '(',
            range.begin,
            range.end,
            range.end_inclusive ? ']' : ')',
            value);
    if (verbose) {
        fprintf(stderr, "\n        %s\n", description);
        fprintf(stderr, "\n");
    }
}

bool Minisat::IntOption::parse(const char * str) {
    const char* span = str;

    if (!match(span, "-") || !match(span, name) || !match(span, "="))
        return false;

    char*   end;
    int32_t tmp = strtol(span, &end, 10);

    if (end == NULL)
        return false;
    else if (tmp > range.end) {
        fprintf(stderr, "ERROR! value <%s> is too large for option \"%s\".\n", span, name);
        exit(1);
    } else if (tmp < range.begin) {
        fprintf(stderr, "ERROR! value <%s> is too small for option \"%s\".\n", span, name);
        exit(1);
    }

    value = tmp;

    return true;
}

void Minisat::IntOption::help(bool verbose) {
    fprintf(stderr, "  -%-12s = %-8s [", name, type_name);
    if (range.begin == INT32_MIN)
        fprintf(stderr, "imin");
    else
        fprintf(stderr, "%4d", range.begin);

    fprintf(stderr, " .. ");
    if (range.end == INT32_MAX)
        fprintf(stderr, "imax");
    else
        fprintf(stderr, "%4d", range.end);

    fprintf(stderr, "] (default: %d)\n", value);
    if (verbose) {
        fprintf(stderr, "\n        %s\n", description);
        fprintf(stderr, "\n");
    }
}

bool Minisat::StringOption::parse(const char * str) {
    const char* span = str;

    if (!match(span, "-") || !match(span, name) || !match(span, "="))
        return false;

    value = span;
    return true;
}

void Minisat::StringOption::help(bool verbose) {
    fprintf(stderr, "  -%-10s = %8s\n", name, type_name);
    if (verbose) {
        fprintf(stderr, "\n        %s\n", description);
        fprintf(stderr, "\n");
    }
}

bool Minisat::BoolOption::parse(const char * str) {
    const char* span = str;

    if (match(span, "-")) {
        bool b = !match(span, "no-");

        if (strcmp(span, name) == 0) {
            value = b;
            return true;
        }
    }

    return false;
}

void Minisat::BoolOption::help(bool verbose) {

    fprintf(stderr, "  -%s, -no-%s", name, name);

    for (uint32_t i = 0; i < 32 - strlen(name) * 2; i++)
        fprintf(stderr, " ");

    fprintf(stderr, " ");
    fprintf(stderr, "(default: %s)\n", value ? "on" : "off");
    if (verbose) {
        fprintf(stderr, "\n        %s\n", description);
        fprintf(stderr, "\n");
    }
}

 #ifndef _MSC_VER

bool Minisat::Int64Option::parse(const char * str) {
    const char* span = str;

    if (!match(span, "-") || !match(span, name) || !match(span, "="))
        return false;

    char*   end;
    int64_t tmp = strtoll(span, &end, 10);

    if (end == NULL)
        return false;
    else if (tmp > range.end) {
        fprintf(stderr, "ERROR! value <%s> is too large for option \"%s\".\n", span, name);
        exit(1);
    } else if (tmp < range.begin) {
        fprintf(stderr, "ERROR! value <%s> is too small for option \"%s\".\n", span, name);
        exit(1);
    }

    value = tmp;

    return true;
}

void Minisat::Int64Option::help(bool verbose) {

#       ifndef PRIi64
#           ifdef __MINGW32__
#                   define PRIi64 "I64i"
#               else
#                   error PRIi64 not defined
#           endif
#       endif

    fprintf(stderr, "  -%-12s = %-8s [", name, type_name);
    if (range.begin == INT64_MIN)
        fprintf(stderr, "imin");
    else
        fprintf(stderr, "%4" PRIi64, range.begin);

    fprintf(stderr, " .. ");
    if (range.end == INT64_MAX)
        fprintf(stderr, "imax");
    else
        fprintf(stderr, "%4" PRIi64, range.end);

    fprintf(stderr, "] (default: %" PRIi64 ")\n", value);
    if (verbose) {
        fprintf(stderr, "\n        %s\n", description);
        fprintf(stderr, "\n");
    }
}

 #endif
