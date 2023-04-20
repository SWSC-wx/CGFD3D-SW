#ifndef H_PRINT_EXECINFO
#define H_PRINT_EXECINFO
#include <execinfo.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
/* Obtain a backtrace and print it to stdout. */
void
print_trace (void)
{
	void *array[10];
	size_t size;
	char **strings;
	size_t i;

	size = backtrace (array, 10);
	strings = backtrace_symbols (array, size);

	printf ("Obtained %zd stack frames.\n", size);

	for (i = 0; i < size; i++)
	printf ("%s\n", strings[i]);

	free (strings);
}

/* A dummy function to make the backtrace more interesting. */
void
sigtrace (int sig)
{
	printf("Printing trace on signal %d:\n", sig);
	print_trace ();
	exit(1);
}
void
sigtrace_noquit (int sig)
{
	printf("Printing trace on signal %d:\n", sig);
	print_trace ();
}


int signal_bind(void)
{
	signal(SIGSEGV, sigtrace);
	signal(SIGUSR1, sigtrace_noquit);
	return 0;
}
#endif
