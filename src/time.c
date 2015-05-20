#include "time.h"
#include "mem.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

TIME *CreateClock(clock_t t){
  TIME *Time = (TIME *) Calloc(1, sizeof(TIME));
  Time->cpu_start = t;
  return Time;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void StopTimeNDRM(TIME *Time, clock_t t){
  Time->cpu_ndrm = t;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void StopCalcAll(TIME *Time, clock_t t){
  Time->cpu_total = t - Time->cpu_start;
  fprintf(stdout, "Total cpu time: %.2g minutes.\n", ((double) Time->cpu_total 
  / CLOCKS_PER_SEC) / 60.);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void RemoveClock(TIME *Time){
  Free(Time);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
