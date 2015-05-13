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
  // IS A FASTA FILE
  if(PA->type == 1){
    switch(sym){
      case '>':  PA->header = 1; return -1;
      case '\n': PA->header = 0; return -1;
      default:   if(PA->header==1) return -1;
      }
    }

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
