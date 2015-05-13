#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <time.h>
#include <pthread.h>
#include "mem.h"
#include "defs.h"
#include "param.h"
#include "msg.h"
#include "parser.h"
#include "buffer.h"
#include "levels.h"
#include "common.h"
#include "context.h"

CModel **Models;

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - C O M P R E S S O R - - - - - - - - - - - - - -
/*
void Compress(Parameters *P, CModel **cModels, uint8_t id, uint32_t 
refNModels, INF *I){
  FILE        *Reader  = Fopen(P->tar[id], "r");
  char        *name    = concatenate(P->tar[id], ".co");
  FILE        *Writter = Fopen(name, "w");
  uint32_t    n, x, k, cModel, totModels, idxPos;
  int32_t     idx = 0;
  uint64_t    compressed = 0, nSymbols = 0, nBases = 0;
  double      *cModelWeight, cModelTotalWeight = 0;
  uint8_t     *readerBuffer, *symbolBuffer, sym, irSym, *pos, type = 0, 
              header = 1, line = 0, dna = 0;
  PModel      **pModel, *MX;
  FloatPModel *PT;
  #ifdef PROGRESS
  uint64_t    i = 0;
  #endif

  if(P->verbose)
    fprintf(stdout, "Analyzing data and creating models ...\n");

  FILE *IAE = NULL;
  char *IAEName = NULL;
  IAEName = concatenate(P->tar[id], ".iae");
  IAE = Fopen(IAEName, "w");
  
  sym = fgetc(Reader);
  switch(sym){
    case '>': type = 1; break;
    case '@': type = 2; break;
    default : type = 0;
    }
  rewind(Reader);

  switch(type){
    case 1:  nBases = NDNASymInFasta(Reader); break;
    case 2:  nBases = NDNASymInFastq(Reader); break;
    default: nBases = NDNASyminFile (Reader); break;
    }
 
  nSymbols      = NBytesInFile(Reader);
  
  // EXTRA MODELS DERIVED FROM EDITS
  totModels = P->nModels;
  for(n = 0 ; n < P->nModels ; ++n) 
    if(P->model[n].edits != 0)
      totModels += 1;

  pModel        = (PModel  **) Calloc(totModels, sizeof(PModel *));
  for(n = 0 ; n < totModels ; ++n)
    pModel[n]   = CreatePModel(ALPHABET_SIZE);
  MX            = CreatePModel(ALPHABET_SIZE);
  PT            = CreateFloatPModel(ALPHABET_SIZE);
  readerBuffer  = (uint8_t *) Calloc(BUFFER_SIZE, sizeof(uint8_t));
  symbolBuffer  = (uint8_t *) Calloc(BUFFER_SIZE + BGUARD, sizeof(uint8_t));
  symbolBuffer += BGUARD;
  cModelWeight  = (double   *) Calloc(totModels, sizeof(double));

  for(n = 0 ; n < totModels ; ++n)
    cModelWeight[n] = 1.0 / totModels;

  for(n = 0 ; n < P->nModels ; ++n){
    if(P->model[n].type == TARGET){
      cModels[n] = CreateCModel(P->model[n].ctx, P->model[n].den, 
      P->model[n].ir, TARGET, P->col, P->model[n].edits, P->model[n].eDen);
      }
    }

  if(P->verbose){
    fprintf(stdout, "Done!\n");
    fprintf(stdout, "Compressing target sequence %d [bases: %"PRIu64"] ...\n", 
    id + 1, nBases);
    }

  while((k = fread(readerBuffer, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos){
      #ifdef PROGRESS
      CalcProgress(nSymbols, ++i);
      #endif

      sym = readerBuffer[idxPos];
      if(type == 1){  // IS A FAST[A] FILE
        if(sym == '>'){ header = 1; continue; }
        if(sym == '\n' && header == 1){ header = 0; continue; }
        if(sym == '\n') continue;
        if(header == 1) continue;
        }
      else if(type == 2){ // IS A FAST[Q] FILE
        switch(line){
          case 0: if(sym == '\n'){ line = 1; dna = 1; } break;
          case 1: if(sym == '\n'){ line = 2; dna = 0; } break;
          case 2: if(sym == '\n'){ line = 3; dna = 0; } break;
          case 3: if(sym == '\n'){ line = 0; dna = 0; } break;
          }
        if(dna == 0 || sym == '\n') continue;
        }

      // REMOVE SPECIAL SYMBOLS [WINDOWS TXT ISSUES]
      if(sym < 65 || sym > 122) continue; 

      // FINAL FILTERING DNA CONTENT
      if(sym != 'A' && sym != 'C' && sym != 'G' && sym != 'T'){
        #ifdef ESTIMATE
        if(P->estim != 0)
          fprintf(IAE, "0\n");
        #endif
        continue;
        }

      symbolBuffer[idx] = sym = DNASymToNum(sym);
      memset((void *)PT->freqs, 0, ALPHABET_SIZE * sizeof(double));

      n = 0;
      pos = &symbolBuffer[idx-1];
      for(cModel = 0 ; cModel < P->nModels ; ++cModel){
        GetPModelIdx(pos, cModels[cModel]);
        ComputePModel(cModels[cModel], pModel[n], cModels[cModel]->pModelIdx,
        cModels[cModel]->alphaDen);
        ComputeWeightedFreqs(cModelWeight[n], pModel[n], PT);

        if(cModels[cModel]->edits != 0){
          ++n;
          cModels[cModel]->SUBS.seq->buf[cModels[cModel]->SUBS.seq->idx] = sym;
          cModels[cModel]->SUBS.idx = GetPModelIdxCorr(cModels[cModel]->SUBS.
          seq->buf+cModels[cModel]->SUBS.seq->idx-1, cModels[cModel], cModels
          [cModel]->SUBS.idx);
          ComputePModel(cModels[cModel], pModel[n], cModels[cModel]->SUBS.idx, 
          cModels[cModel]->SUBS.eDen);
          ComputeWeightedFreqs(cModelWeight[n], pModel[n], PT);
          }

        ++n;
        }

      MX->sum  = MX->freqs[0] = 1 + (unsigned) (PT->freqs[0] * MX_PMODEL);
      MX->sum += MX->freqs[1] = 1 + (unsigned) (PT->freqs[1] * MX_PMODEL);
      MX->sum += MX->freqs[2] = 1 + (unsigned) (PT->freqs[2] * MX_PMODEL);
      MX->sum += MX->freqs[3] = 1 + (unsigned) (PT->freqs[3] * MX_PMODEL);

      fprintf(IAE, "%.3g\n", 2-(PModelSymbolNats(MX, sym) / M_LN2));

      cModelTotalWeight = 0;
      for(n = 0 ; n < totModels ; ++n){
        cModelWeight[n] = Power(cModelWeight[n], P->gamma) * (double) 
        pModel[n]->freqs[sym] / pModel[n]->sum;
        cModelTotalWeight += cModelWeight[n];
        }

      for(n = 0 ; n < P->nModels ; ++n){
        if(cModels[n]->ref == TARGET){
          UpdateCModelCounter(cModels[n], sym, cModels[n]->pModelIdx);
          if(cModels[n]->ir != 0){                // REVERSE COMPLEMENTS
            irSym = GetPModelIdxIR(symbolBuffer+idx, cModels[n]);
            UpdateCModelCounter(cModels[n], irSym, cModels[n]->pModelIdxIR);
            }
          }
        }

      for(n = 0 ; n < totModels ; ++n)
        cModelWeight[n] /= cModelTotalWeight; // RENORMALIZE THE WEIGHTS

      n = 0;
      for(cModel = 0 ; cModel < P->nModels ; ++cModel){
        if(cModels[cModel]->edits != 0){
          CorrectCModelSUBS(cModels[cModel], pModel[++n], sym);
          }
        ++n;
        }

      if(++idx == BUFFER_SIZE){
        memcpy(symbolBuffer-BGUARD, symbolBuffer+idx-BGUARD, BGUARD);
        idx = 0;
        }

      ++compressed;
      }

  fclose(Writter);
  fclose(IAE);
  Free(IAEName);

  Free(MX);
  Free(name);
  Free(cModelWeight);
  for(n = 0 ; n < P->nModels ; ++n)
    if(P->model[n].type == REFERENCE)
      ResetCModelIdx(cModels[n]);
    else
      FreeCModel(cModels[n]);
  for(n = 0 ; n < totModels ; ++n){
    Free(pModel[n]->freqs);
    Free(pModel[n]);
    }
  Free(pModel);
  Free(PT);
  Free(readerBuffer);
  Free(symbolBuffer-BGUARD);
  fclose(Reader);

  if(P->verbose == 1)
    fprintf(stdout, "Done!                          \n");  // SPACES ARE VALID 

  I[id].bytes = 2222222;
  I[id].size  = compressed;
  }
*/


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - F I L T E R I N G - - - - - - - - - - - - -

void FilterTarget(Threads T){
  //XXX: REF , TAR , wTAR ?



  return;
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - T H R E A D I N G - - - - - - - - - - - - - - -

void *FilterThread(void *Thr){
  Threads *T = (Threads *) Thr;
  FilterTarget(T[0]);
  pthread_exit(NULL);
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - R E F E R E N C E - - - - - - - - - - - - -

void LoadReference(Threads T){
  FILE     *Reader = Fopen(P->files[T.id], "r");
  uint32_t n, k, idxPos;
  PARSER   *PA = CreateParser();
  CBUF     *symBuf = CreateCBuffer(BUFFER_SIZE, BGUARD);
  uint8_t  *readBuf = (uint8_t *) Calloc(BUFFER_SIZE, sizeof(uint8_t));
  uint8_t  sym, irSym;

  FileType(PA, Reader);

  while((k = fread(readBuf, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos){

      if(ParseSym(PA, (sym = readBuf[idxPos])) == -1) continue;
      symBuf->buf[symBuf->idx] = sym = DNASymToNum(sym);

      for(n = 0 ; n < P->nModels ; ++n){
        GetPModelIdx(symBuf->buf+symBuf->idx-1, Models[n]);
        UpdateCModelCounter(Models[n], sym, Models[n]->pModelIdx);
        if(Models[n]->ir == 1){                         // INVERTED REPEATS
          irSym = GetPModelIdxIR(symBuf->buf+symBuf->idx, Models[n]);
          UpdateCModelCounter(Models[n], irSym, Models[n]->pModelIdxIR);
          }
        }

      UpdateCBuffer(symBuf);
      }
 
  for(n = 0 ; n < P->nModels ; ++n)
    ResetCModelIdx(Models[n]);
  Free(readBuf);
  RemoveCBuffer(symBuf);
  RemoveParser(PA);
  fclose(Reader);
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - C O M P R E S S O R - - - - - - - - - - - - - -

void CompressAction(Threads *T, uint32_t ref){
  uint32_t n;

  Models = (CModel **) Malloc(P->nModels * sizeof(CModel *));
  for(n = 0 ; n < P->nModels ; ++n)
    Models[n] = CreateCModel(T[ref].model[n].ctx, T[ref].model[n].den, 
    T[ref].model[n].ir, REFERENCE, P->col, T[ref].model[n].edits, 
    T[ref].model[n].eDen);

  LoadReference(T[ref]);

  pthread_t t[P->nThreads];

  //XXX: if n = ref ??? (DIAGONAL)
  ref = 0;
  do{
    for(n = 0 ; n < P->nThreads ; ++n)
      pthread_create(&(t[n+1]), NULL, FilterThread, (void *) &(T[ref+n]));
    for(n = 0 ; n < P->nThreads ; ++n) // DO NOT JOIN FORS!
      pthread_join(t[n+1], NULL);
    }
  while((ref += P->nThreads) < P->nFiles && ref + P->nThreads <= P->nFiles);

  if(ref < P->nFiles){ // EXTRA - OUT OF THE MAIN LOOP
    for(n = ref ; n < P->nFiles ; ++n)
      pthread_create(&(t[n+1]), NULL, FilterThread, (void *) &(T[n]));
    for(n = ref ; n < P->nFiles ; ++n) // DO NOT JOIN FORS!
      pthread_join(t[n+1], NULL);
    }

/*// HERE, WITH THREADING, MEMORY IS QUADRATIC!
  for(n = 0 ; n < P->nFiles ; ++n)
    if(T.id != n){
      //LoadReference(T, Models, IN); //NEW REFERENCE
      //CompressTarget(T, Models, T.id, n);
      }
*/
 
  for(n = 0 ; n < P->nModels ; ++n)
    FreeCModel(Models[n]);
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int32_t main(int argc, char *argv[]){
  char        **p = *&argv, **xargv, *xpl = NULL;
  int32_t     xargc = 0;
  uint32_t    n, k, col, ref;
  clock_t     start = clock();
  double      gamma;
  Threads     *T;
  
  P = (Parameters *) Malloc(1 * sizeof(Parameters));
  if((P->help = ArgsState(DEFAULT_HELP, p, argc, "-h")) == 1 || argc < 2){
    PrintMenu();
    return EXIT_SUCCESS;
    }

  if(ArgsState(DEF_VERSION, p, argc, "-V")){
    PrintVersion();
    return EXIT_SUCCESS;
    }

  if(ArgsState(0, p, argc, "-s")){
    PrintLevels(); 
    return EXIT_SUCCESS;
    }

  P->verbose  = ArgsState  (DEFAULT_VERBOSE, p, argc, "-v" );
  P->force    = ArgsState  (DEFAULT_FORCE,   p, argc, "-f" );
  P->level    = ArgsNum    (0, p, argc, "-l", MIN_LEVEL, MAX_LEVEL);
  P->nThreads = ArgsNum    (DEFAULT_THREADS, p, argc, "-t", MIN_THREADS, 
  MAX_THREADS);

  gamma = DEFAULT_GAMMA;
  for(n = 1 ; n < xargc ; ++n) 
    if(strcmp(xargv[n], "-g") == 0) 
      gamma = atof(xargv[n+1]);

  col = MAX_COLLISIONS;
  for(n = 1 ; n < xargc ; ++n) 
    if(strcmp(xargv[n], "-c") == 0) 
      col = atoi(xargv[n+1]);

  P->col      = ArgsNum    (col,   p, argc, "-c", 1, 200);
  P->gamma    = ArgsDouble (gamma, p, argc, "-g");
  P->gamma    = ((int)(P->gamma * 65536)) / 65536.0;
  P->nFiles   = ReadFNames (P, argv[argc-1]);

  P->nModels  = 0;
  for(n = 1 ; n < argc ; ++n)
    if(strcmp(argv[n], "-m") == 0)
      P->nModels += 1;

  if(P->nModels == 0 && P->level == 0)
    P->level = DEFAULT_LEVEL;

  if(P->level != 0){
    xpl = GetLevels(P->level);
    xargc = StrToArgv(xpl, &xargv);
    for(n = 1 ; n < xargc ; ++n)
      if(strcmp(xargv[n], "-m") == 0)
        P->nModels += 1;
    }

  if(P->nModels == 0){
    fprintf(stderr, "Error: at least you need to use a context model!\n");
    return EXIT_FAILURE;
    }

  // READ MODEL PARAMETERS FROM XARGS & ARGS
  T = (Threads *) Calloc(P->nFiles, sizeof(Threads));
  for(ref = 0 ; ref < P->nFiles ; ++ref){
    T[ref].model = (ModelPar *) Calloc(P->nModels, sizeof(ModelPar));
    T[ref].id = ref;
    k = 0;
    for(n = 1 ; n < argc ; ++n)
      if(strcmp(argv[n], "-m") == 0)
        T[ref].model[k++] = ArgsUniqModel(argv[n+1], 0);
    if(P->level != 0){
      for(n = 1 ; n < xargc ; ++n)
        if(strcmp(xargv[n], "-m") == 0)
          T[ref].model[k++] = ArgsUniqModel(xargv[n+1], 0);
      }
    }

  if(P->verbose) PrintArgs(P, T[0]);

  P->size   = (uint64_t *) Calloc(P->nFiles, sizeof(uint64_t));
  P->matrix = (double  **) Calloc(P->nFiles, sizeof(double *));
  for(n = 0 ; n < P->nFiles ; ++n){
    P->matrix[n] = (double *) Calloc(P->nFiles, sizeof(double));
    }

  if(P->nThreads > P->nFiles){
    fprintf(stderr, "Error: the number of threads must not be higher than the "
    "number of files\n");
    exit(1);
    }

  for(n = 0 ; n < P->nFiles ; ++n){
    fprintf(stderr, "Running reference %u ...\n", n+1);
    CompressAction(T, n);
    fprintf(stderr, "Done!\n");
    }

  fprintf(stdout, "Final matrix:\n");
  for(n = 0 ; n < P->nFiles ; ++n){
    for(k = 0 ; k < P->nFiles ; ++k)
      fprintf(stdout, "%7.5g ", P->matrix[n][k]);
    fprintf(stdout, "\n");
    }

  //TODO: human readable & min protection
  fprintf(stdout, "Total cpu time: %.2g minutes.\n", (((double) (clock() - 
  start)) / CLOCKS_PER_SEC) / 60.);

  return EXIT_SUCCESS;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
