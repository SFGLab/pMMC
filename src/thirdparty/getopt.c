/*
 * getopt.c - Portable getopt for Windows
 *
 * Public domain implementation compatible with POSIX getopt().
 * Based on the classic AT&T implementation.
 *
 * Only compiled on Windows (_WIN32). On Linux, the system libc
 * provides getopt via <unistd.h>.
 */

#ifdef _WIN32

#include <stdio.h>
#include <string.h>
#include "getopt.h"

char *optarg = NULL;
int optind = 1;
int opterr = 1;
int optopt = '?';

int getopt(int argc, char *const argv[], const char *optstring) {
    static int sp = 1;
    int c;
    const char *cp;

    if (sp == 1) {
        /* Check if we've run out of arguments or hit a non-option */
        if (optind >= argc ||
            argv[optind][0] != '-' ||
            argv[optind][1] == '\0') {
            return -1;
        }
        /* Handle "--" end-of-options marker */
        if (argv[optind][0] == '-' && argv[optind][1] == '-' &&
            argv[optind][2] == '\0') {
            optind++;
            return -1;
        }
    }

    c = argv[optind][sp];
    optopt = c;

    /* Look up the option character in the option string */
    cp = strchr(optstring, c);
    if (c == ':' || cp == NULL) {
        if (opterr)
            fprintf(stderr, "%s: unknown option '-%c'\n", argv[0], c);
        if (argv[optind][++sp] == '\0') {
            optind++;
            sp = 1;
        }
        return '?';
    }

    /* Check if this option requires an argument */
    if (*(cp + 1) == ':') {
        if (argv[optind][sp + 1] != '\0') {
            /* Argument is the rest of this argv element */
            optarg = &argv[optind][sp + 1];
        } else if (optind + 1 < argc) {
            /* Argument is the next argv element */
            optarg = argv[++optind];
        } else {
            if (opterr)
                fprintf(stderr, "%s: option '-%c' requires an argument\n",
                        argv[0], c);
            optind++;
            sp = 1;
            return (optstring[0] == ':') ? ':' : '?';
        }
        optind++;
        sp = 1;
    } else {
        optarg = NULL;
        if (argv[optind][++sp] == '\0') {
            optind++;
            sp = 1;
        }
    }

    return c;
}

#endif /* _WIN32 */
