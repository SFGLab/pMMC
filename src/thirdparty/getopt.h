/*
 * getopt.h - Portable getopt for Windows
 *
 * Public domain implementation compatible with POSIX getopt().
 * On Linux/POSIX systems, use <unistd.h> instead of this header.
 */

#ifndef GETOPT_H_
#define GETOPT_H_

#ifdef _WIN32

#ifdef __cplusplus
extern "C" {
#endif

extern char *optarg;
extern int optind;
extern int opterr;
extern int optopt;

int getopt(int argc, char *const argv[], const char *optstring);

#ifdef __cplusplus
}
#endif

#endif /* _WIN32 */

#endif /* GETOPT_H_ */
