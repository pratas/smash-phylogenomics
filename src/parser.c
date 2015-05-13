#include <string.h>
#include "parser.h"
#include "mem.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// CREATE PARSER
//
PARSER *CreateParser(void){
  PARSER *PA = (PARSER *) Calloc(1, sizeof(PARSER));
  PA->sym    = 0;
  PA->type   = 0;
  PA->header = 0;
  return PA;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// FILE TYPE
//
void FileType(PARSER *PA, FILE *IN){
  rewind(IN);
  PA->sym = fgetc(IN);
    switch(PA->sym){
    case '>': PA->type = 1; break;
    case '@': PA->type = 2; break;
    default : PA->type = 0;
    }
  rewind(IN);
  } 

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// PARSE SYM
//
int32_t ParseSym(PARSER *PA, uint8_t sym){

  switch(PA->type){
    // IS A FASTA FILE
    case 1:
      switch(sym){
        case '>':  PA->header = 1; return -1;
        case '\n': PA->header = 0; return -1;
        default:   if(PA->header==1) return -1;
        }
    break;

    // IS A FASTQ FILE
    case 2:
      /*
      switch(line){
        case 0: if(sym == '\n'){ line = 1; dna = 1; begin = 0; } break;
        case 1: if(sym == '\n'){ line = 2; dna = 0; }            break;
        case 2: if(sym == '\n'){ line = 3; dna = 0; }            break;
        case 3: if(sym == '\n'){ line = 0; dna = 0; }            break;
        }
      if(dna == 0 || sym == '\n') continue;
      */
    break;

    // OTHER (SUCH AS DNA SEQ)
    default: ;
    }

  // NUCLEOTIDE PARSE
  if(sym != 'A' && sym != 'C' && sym != 'G' && sym != 'T')
    return -1;

  return sym;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// REMOVE PARSER
//
void RemoveParser(PARSER *PA){
  Free(PA);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
