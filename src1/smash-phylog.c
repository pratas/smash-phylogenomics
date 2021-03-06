#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <time.h>
#include <pthread.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/uio.h>
#include <sys/mman.h>
#include "mem.h"
#include "time.h"
#include "defs.h"
#include "param.h"
#include "msg.h"
#include "parser.h"
#include "buffer.h"
#include "levels.h"
#include "common.h"
#include "cmodel.h"
#include "filters.h"
#include "paint.h"

CModel **Models;  // MEMORY SHARED BY THREADING

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - - P A I N T - - - - - - - - - - - - - - - -

void PaintMatrix(void){
  FILE *Plot = Fopen("plot.svg", "w");
  Painter *Paint;
  uint32_t ref, tar;

  Paint = CreatePainter(DEFAULT_CX * 2 + ((DEFAULT_WIDTH + DEFAULT_SPACE) * 
  P->nFiles) + 5);

  PrintHead(Plot, (2 * DEFAULT_CX) + (((Paint->width + DEFAULT_SPACE) *
  P->nFiles) - DEFAULT_SPACE), Paint->size + EXTRA);

  Rect(Plot, (2 * DEFAULT_CX) + (((Paint->width + DEFAULT_SPACE) *
  P->nFiles) - DEFAULT_SPACE), Paint->size + EXTRA, 0, 0, "#ffffff");

  // PRINT HEATMAP SCALE
  uint32_t size = (Paint->width + DEFAULT_SPACE) * P->nFiles - DEFAULT_SPACE;
  for(ref = 0 ; ref < size ; ++ref){
    char color[12];
    Rect(Plot, Paint->width, 1, DEFAULT_CX - (DEFAULT_WIDTH*2), Paint->cy + ref, 
    HeatMapColor(((double) ref / size), color));
    }
  if(P->nFiles > 4)
    Text90d(Plot, -DEFAULT_CX - ((size/2)+46), DEFAULT_CX-(DEFAULT_WIDTH*2 + 2), 
    "SIMILARITY");
  Text   (Plot, DEFAULT_CX-(DEFAULT_WIDTH*2 + 14), Paint->cy+13, "1");
  Text   (Plot, DEFAULT_CX-(DEFAULT_WIDTH*2 + 14), Paint->cy+size, "0");

  for(ref = 0 ; ref < P->nFiles ; ++ref){
    for(tar = P->nFiles ; tar-- ; ){ // INVERT LOOP: INCREASE SPEED OF LEARNING
      char color[12];
      Rect(Plot, Paint->width, Paint->width, Paint->cx, Paint->cy, 
      HeatMapColor(BoundDouble(0.0, P->matrix[ref][tar], 1.0), color));
      Paint->cx += Paint->width + DEFAULT_SPACE;
      }
    // TEXT HAS 16 PX -> CALCULATE AVERAGE POSITION
    Text   (Plot, Paint->cx + 4, (Paint->cy+Paint->width/2)+6, P->files[ref]);
    Text90d(Plot, 4-DEFAULT_CX, (Paint->cy+Paint->width/2)+10, 
    P->files[P->nFiles-1-ref]);
    Paint->cx =  DEFAULT_CX;
    Paint->cy += Paint->width + DEFAULT_SPACE;
    }
  
  PrintFinal(Plot);
  RemovePainter(Paint);
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - C O M P R E S S I N G - - - - - - - - - - - - - 

void CompressTarget(Threads T){
  FILE        *Reader  = Fopen(P->files[T.id], "r");
  double      *cModelWeight, cModelTotalWeight = 0, bits = 0, instance = 0;
  uint64_t    nBase = 0;
  uint32_t    n, k, idxPos, totModels, cModel;
  PARSER      *PA = CreateParser();
  CBUF        *symBuf = CreateCBuffer(BUFFER_SIZE, BGUARD);
  uint8_t     *readBuf = (uint8_t *) Calloc(BUFFER_SIZE, sizeof(uint8_t));
  uint8_t     sym, *pos;
  PModel      **pModel, *MX;
  CModel      **Shadow;
  FloatPModel *PT;

  totModels = P->nModels; // EXTRA MODELS DERIVED FROM EDITS
  for(n = 0 ; n < P->nModels ; ++n) 
    if(T.model[n].edits != 0)
      totModels += 1;

  Shadow = (CModel **) Calloc(P->nModels, sizeof(CModel *));
  for(n = 0 ; n < P->nModels ; ++n)
    Shadow[n] = CreateShadowModel(Models[n]); 
  pModel        = (PModel  **) Calloc(totModels, sizeof(PModel *));
  for(n = 0 ; n < totModels ; ++n)
    pModel[n]   = CreatePModel(ALPHABET_SIZE);
  MX            = CreatePModel(ALPHABET_SIZE);
  PT            = CreateFloatPModel(ALPHABET_SIZE);
  cModelWeight  = (double   *) Calloc(totModels, sizeof(double));

  for(n = 0 ; n < totModels ; ++n)
    cModelWeight[n] = 1.0 / totModels;

  FileType(PA, Reader);

  nBase = 0;
  while((k = fread(readBuf, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos){
      if(ParseSym(PA, (sym = readBuf[idxPos])) == -1) continue;
      symBuf->buf[symBuf->idx] = sym = DNASymToNum(sym);

      memset((void *)PT->freqs, 0, ALPHABET_SIZE * sizeof(double));

      n = 0;
      pos = &symBuf->buf[symBuf->idx-1];
      for(cModel = 0 ; cModel < P->nModels ; ++cModel){
        GetPModelIdx(pos, Shadow[cModel]);
        ComputePModel(Models[cModel], pModel[n], Shadow[cModel]->pModelIdx,
        Shadow[cModel]->alphaDen);
        ComputeWeightedFreqs(cModelWeight[n], pModel[n], PT);
        if(Shadow[cModel]->edits != 0){
          ++n;
          Shadow[cModel]->SUBS.seq->buf[Shadow[cModel]->SUBS.seq->idx] = sym;
          Shadow[cModel]->SUBS.idx = GetPModelIdxCorr(Shadow[cModel]->SUBS.
          seq->buf+Shadow[cModel]->SUBS.seq->idx-1, Shadow[cModel], Shadow
          [cModel]->SUBS.idx);
          ComputePModel(Models[cModel], pModel[n], Shadow[cModel]->SUBS.idx, 
          Shadow[cModel]->SUBS.eDen);
          ComputeWeightedFreqs(cModelWeight[n], pModel[n], PT);
          }
        ++n;
        }

      MX->sum  = (MX->freqs[0] = 1 + (unsigned) (PT->freqs[0] * MX_PMODEL));
      MX->sum += (MX->freqs[1] = 1 + (unsigned) (PT->freqs[1] * MX_PMODEL));
      MX->sum += (MX->freqs[2] = 1 + (unsigned) (PT->freqs[2] * MX_PMODEL));
      MX->sum += (MX->freqs[3] = 1 + (unsigned) (PT->freqs[3] * MX_PMODEL));
      bits += (instance = PModelSymbolLog(MX, sym));
      nBase++;

      cModelTotalWeight = 0;
      for(n = 0 ; n < totModels ; ++n){
        cModelWeight[n] = Power(cModelWeight[n], P->gamma) * (double) 
        pModel[n]->freqs[sym] / pModel[n]->sum;
        cModelTotalWeight += cModelWeight[n];
        }

      for(n = 0 ; n < totModels ; ++n)
        cModelWeight[n] /= cModelTotalWeight; // RENORMALIZE THE WEIGHTS

      n = 0;
      for(cModel = 0 ; cModel < P->nModels ; ++cModel){
        if(Shadow[cModel]->edits != 0){
          CorrectCModelSUBS(Shadow[cModel], pModel[++n], sym);
          }
        ++n;
        }

      UpdateCBuffer(symBuf);
      }

  Free(cModelWeight);
  for(n = 0 ; n < totModels ; ++n)
    RemovePModel(pModel[n]);
  Free(pModel);
  RemovePModel(MX);
  RemoveFPModel(PT);
  for(n = 0 ; n < P->nModels ; ++n)
    FreeShadow(Shadow[n]);
  Free(Shadow);
  Free(readBuf);
  RemoveCBuffer(symBuf);
  RemoveParser(PA);
  fclose(Reader);

  P->matrix[P->ref][T.id] = nBase == 0 ? 101 : bits / 2 / nBase; // 101 -> nan
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - F   T H R E A D I N G - - - - - - - - - - - - - - -

void *CompressThread(void *Thr){
  Threads *T = (Threads *) Thr;
  CompressTarget(T[0]);
  pthread_exit(NULL);
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - R E F E R E N C E - - - - - - - - - - - - -

void LoadReference(Threads T){
  FILE     *Reader = Fopen(P->files[T.id], "r");
  uint32_t n;
  PARSER   *PA = CreateParser();
  CBUF     *symBuf = CreateCBuffer(BUFFER_SIZE, BGUARD);
  uint8_t  sym, irSym, *readBuf;
  FileType(PA, Reader);
  fclose(Reader);
  struct   stat s;
  size_t   size, k;
  long     fd = open(P->files[T.id], O_RDONLY);

  fstat (fd, & s);
  size = s.st_size;
  readBuf = (uint8_t *) mmap(0, size, PROT_READ, MAP_PRIVATE, fd, 0);
  for(k = 0 ; k < size ; ++k){
    if(ParseSym(PA, (sym = *readBuf++)) == -1) continue;
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
  RemoveCBuffer(symBuf);
  RemoveParser(PA);
  close(fd);
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - C O M P R E S S O R - - - - - - - - - - - - - -

void CompressAction(Threads *T, uint32_t ref){
  uint32_t n, k;
  pthread_t t[P->nThreads];
  P->ref = ref;

  Models = (CModel **) Malloc(P->nModels * sizeof(CModel *));
  for(n = 0 ; n < P->nModels ; ++n)
    Models[n] = CreateCModel(T[ref].model[n].ctx, T[ref].model[n].den, 
    T[ref].model[n].ir, REFERENCE, P->col, T[ref].model[n].edits, 
    T[ref].model[n].eDen);

  fprintf(stderr, "  [+] Loading reference %u ... ", ref+1);
  LoadReference(T[ref]);
  fprintf(stderr, "Done!\n");
  
  fprintf(stderr, "      [+] Compressing %u targets ... ", P->nFiles);
  ref = 0;
  do{
    for(n = 0 ; n < P->nThreads ; ++n)
      pthread_create(&(t[n+1]), NULL, CompressThread, (void *) &(T[ref+n]));
    for(n = 0 ; n < P->nThreads ; ++n) // DO NOT JOIN FORS!
      pthread_join(t[n+1], NULL);
    }
  while((ref += P->nThreads) < P->nFiles && ref + P->nThreads <= P->nFiles);

  if(ref < P->nFiles){ // EXTRA - OUT OF THE MAIN LOOP
    for(n = ref, k = 0 ; n < P->nFiles ; ++n)
      pthread_create(&(t[++k]), NULL, CompressThread, (void *) &(T[n]));
    for(n = ref, k = 0 ; n < P->nFiles ; ++n) // DO NOT JOIN FORS!
      pthread_join(t[++k], NULL);
    }
  fprintf(stderr, "Done!\n");

  for(n = 0 ; n < P->nModels ; ++n)
    FreeCModel(Models[n]);
  Free(Models);
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int32_t main(int argc, char *argv[]){
  char        **p = *&argv, **xargv, *xpl = NULL;
  int32_t     xargc = 0;
  uint32_t    n, k, col, ref, index;
  double      gamma, threshold;
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
  P->force    = ArgsState  (DEFAULT_FORCE,   p, argc, "-F" );
  P->level    = ArgsNum    (0, p, argc, "-l", MIN_LEVEL, MAX_LEVEL);
  P->nThreads = ArgsNum    (DEFAULT_THREADS, p, argc, "-n", MIN_THREADS, 
  MAX_THREADS);

  P->nModels = 0;
  for(n = 1 ; n < argc ; ++n)
    if(strcmp(argv[n], "-m") == 0)
      P->nModels += 1;

  if(P->nModels == 0 && P->level == 0)
    P->level = DEFAULT_LEVEL;

  P->nModels  = 0;
  if(P->level != 0){
    xpl = GetLevels(P->level);
    xargc = StrToArgv(xpl, &xargv);
    for(n = 1 ; n < xargc ; ++n)
      if(strcmp(xargv[n], "-m") == 0)
        P->nModels += 1;
    }

  gamma = DEFAULT_GAMMA;
  for(n = 1 ; n < xargc ; ++n) 
    if(strcmp(xargv[n], "-g") == 0) 
      gamma = atof(xargv[n+1]);

  threshold = DEFAULT_THRESHOLD;
  for(n = 1 ; n < xargc ; ++n)
    if(strcmp(xargv[n], "-t") == 0)
      threshold = atof(xargv[n+1]);

  index = DEFAULT_INDEX;
  for(n = 1 ; n < xargc ; ++n)
    if(strcmp(xargv[n], "-i") == 0)
      index = atoi(xargv[n+1]);

  col = MAX_COLLISIONS;
  for(n = 1 ; n < xargc ; ++n) 
    if(strcmp(xargv[n], "-c") == 0) 
      col = atoi(xargv[n+1]);

  P->col        = ArgsNum    (col,   p, argc, "-c", 1, 200);
  P->gamma      = ArgsDouble (gamma, p, argc, "-g");
  P->gamma      = ((int)(P->gamma * 65536)) / 65536.0;
  P->nFiles     = ReadFNames (P, argv[argc-1]);
  P->index      = ArgsNum    (index, p, argc, "-i", 1, P->nFiles);

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

  fprintf(stderr, "\n");
  if(P->verbose) PrintArgs(P, T[0]);

  P->matrix = (double  **) Calloc(P->nFiles, sizeof(double *));
  P->size   = (uint64_t *) Calloc(P->nFiles, sizeof(uint64_t));
  for(n = 0 ; n < P->nFiles ; ++n){
    P->size[n]   = FopenBytesInFile(P->files[n]);
    P->matrix[n] = (double *) Calloc(P->nFiles, sizeof(double));
    }

  if(P->nThreads > P->nFiles){
    fprintf(stderr, "Error: the number of threads must not be higher than the "
    "number of files\n");
    exit(1);
    }

  fprintf(stderr, "==[ PROCESSING ]====================\n");
  TIME *Time = CreateClock(clock());

  for(n = 0 ; n < P->nFiles ; ++n)
    CompressAction(T, n);
  StopTimeNDRM(Time, clock());

  fprintf(stderr, "  [+] Painting heatmap ... ");
  PaintMatrix();
  fprintf(stderr, "Done! \n\n");

  fprintf(stderr, "==[ RESULTS ]=======================\n");
  fprintf(stderr, "Normalized Dissimilarity Rate matrix:\n");
  for(n = 0 ; n < P->nFiles ; ++n){
    for(k = 0 ; k < P->nFiles ; ++k)
      fprintf(stderr, "%.4lf\t", P->matrix[n][k]);
    fprintf(stderr, "\n");
    }
  fprintf(stderr, "\n");

  fprintf(stderr, "==[ STATISTICS ]====================\n");
  StopCalcAll(Time, clock());
  fprintf(stderr, "\n");

  RemoveClock(Time);

  for(ref = 0 ; ref < P->nFiles ; ++ref)
    Free(T[ref].model);
  Free(T);

  return EXIT_SUCCESS;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
