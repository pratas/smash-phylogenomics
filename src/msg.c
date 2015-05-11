#include "msg.h"
#include <stdio.h>
#include <stdlib.h>

void ModelsExplanation(void){
  fprintf(stderr,
  "                                                                         \n"
  "  -m <c>:<d>:<i>:<m/e>  reference context model (ex:-m 13:100:0:0/0),    \n"
  "  -m <c>:<d>:<i>:<m/e>  reference context model (ex:-m 18:1000:0:1/1000),\n"
  "  ...                                                                    \n"
  "                        templates use <c> for context-order size, <d> for\n"
  "                        alpha (1/<d>), <i> (0 or 1) to set the usage of  \n"
  "                        inverted repeats (set 1 to use) and <m> to the   \n"
  "                        maximum allowed mutation on the context without  \n"
  "                        being discarded (usefull in deep contexts), under\n"
  "                        the estimator <e>.                               \n"
  "                                                                       \n");
  }

void PrintMenu(void){
  fprintf(stderr,
  "Usage: smash-phylog [OPTION]... [FILE1]:[FILE2]:[FILE3]:...               \n"
  "A similarity matrix builder for Phylogenomics computation                \n"
  "                                                                         \n"
  "Non-mandatory arguments:                                                 \n"
  "                                                                         \n"
  "  -h                          give this help,                            \n"
  "  -s                          show compression levels,                   \n"
  "  -v                          verbose mode (more information),           \n"
  "  -V                          display version number,                    \n"
  "  -f                          force overwrite of output,                 \n"
  "  -l <level>                  level of compression [1;9],                \n"
  "  -c <cache>                  maximum collisions for hash cache. Memory  \n"
  "                              values are higly dependent of the parameter\n"
  "                              specification,                             \n"
  "                                                                         \n"
  "Mandatory arguments:                                                     \n"
  "                                                                         \n"
  "  <FILE>:<FILE>:<FILE>:<...>  input files (last argument). Use \":\" for \n"
  "                              file name splitting.                       \n"
  "                                                                         \n"
  "Report bugs to <{pratas,ap,pjf}@ua.pt>.                                \n");
  }

void PrintVersion(void){
  fprintf(stderr,
  "                                                                         \n"
  "                         ====================                            \n"
  "                         | smash-phylog %u.%u |                          \n"
  "                         ====================                            \n"
  "                                                                         \n"
  "smash-phylog: similarity matrix builder for phylogenomics computation.   \n"
  "Copyright (C) 2014-2015 University of Aveiro. This is a Free software.   \n"
  "You may redistribute copies of it under the terms of the GNU - General   \n"
  "Public License v2 <http://www.gnu.org/licenses/gpl.html>. There is NOT   \n"
  "ANY WARRANTY, to the extent permitted by law. Developed and Written by   \n"
  "Diogo Pratas, Armando J. Pinho and Paulo J. S. G. Ferreira.\n\n", VERSION,
  RELEASE);
  }

