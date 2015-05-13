#ifndef PARSER_H_INCLUDED
#define PARSER_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include "defs.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

typedef struct{
  uint8_t  sym;
  uint8_t  type;
  uint8_t  header;
  uint8_t  line;
  uint8_t  dna;
  }
PARSER;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void     FileType      (PARSER *, FILE *);
PARSER   *CreateParser (void);
int32_t  ParseSym      (PARSER *, uint8_t);
void     RemoveParser  (PARSER *);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
