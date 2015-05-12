#ifndef LEVELS_H_INCLUDED
#define LEVELS_H_INCLUDED

#include "defs.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// COMPRESSION LEVELS
//
// LEVEL_0 IS RESERVED TO THE ALGORITHM
#define LEVEL_1 " 1: -m 13:100:0:0/0 "
#define LEVEL_2 " 2: -m 13:100:1:0/0 "
#define LEVEL_3 " 3: -m 14:200:0:0/0 "
#define LEVEL_4 " 4: -m 14:200:1:0/0 "
#define LEVEL_5 " 5: -m 18:500:1:0/0 -c 30 "
#define LEVEL_6 " 6: -m 18:500:1:1/100 -c 35 -g 0.85 "
#define LEVEL_7 " 7: -m 18:500:1:1/100 -m 13:100:0:0/0 -c 35 -g 0.85 "
#define LEVEL_8 " 8: -m 20:1000:1:3/100 -m 13:100:0:0/0 -m 8:1:0:0/0 -c 40 -g 0.85 "
#define LEVEL_9 " 9: -m 20:1000:1:3/100 -m 14:200:1:0/0 -m 10:10:0:0/0 -m 4:1:0:0/0 -c 40 -g 0.85 "
#define LEVEL_10 " 10: -m 20:1000:1:3/100 -m 14:200:1:0/0 -m 12:20:0:0/0 -m 4:1:0:0/0 -c 50 -g 0.85 "

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

char    *GetLevels  (uint8_t);
void    PrintLevels (void);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif

