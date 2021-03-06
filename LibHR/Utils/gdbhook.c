#include "gdbhook.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

// To be able to run GDB on MPI jobs, as desribed in
// https://www.open-mpi.org/faq/?category=debugging#serial-debuggers
void gdb_hook() {

  char *debugging_status = getenv("GDBMPI");
  if (0 == strcmp("1", debugging_status)) {

    volatile int i = 0;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("PID %d on %s ready for attach\n", getpid(), hostname);
    fflush(stdout);
    while (0 == i)
      sleep(5);
  }
}
