#include <float.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

#include "../lib/klib/khash.h"
#include "../lib/klib/kseq.h"
#include "dude.h"
#include "seq.h"
#include <version.h>

#define VERSION "0.4"

// default values of the options
#define DEFAULT_K 5
#define DEFAULT_I 2
#define DEFAULT_L 200
#define DEFAULT_R 20
#define DEFAULT_C 0

// context table: the key is the context without the centre symbol, the value
// is the index of the count vector of the key in counts
KHASH_MAP_INIT_INT64(ctx, int)

typedef struct {
  khash_t(ctx) * index;
  // dim counts for each key: how many times each centre symbol was seen
  int *counts;
  int dim;
  // count vectors used and allocated
  int n, m;
} ctx_table_t;

// rolling window of the 2k + 1 symbols around the centre
// substitutions: 2-bit base codes, homopolymers: 6-bit symbols (hom_symbol)
typedef struct {
  int k;
  int bits;
  uint64_t side_mask;
  // the k symbols on each side of the centre, the oldest in the high bits
  uint64_t left, right;
  int centre;
  // symbols since the start or since the last base that is not A, C, G or T
  int filled;
} window_t;

// kseq init
__KS_GETC(gzread, BUFFER_SIZE)

__KS_GETUNTIL(gzread, BUFFER_SIZE)

__KSEQ_READ(static)

// log file
FILE *LOG_FILE;

// denoise
void denoise_substitutions(dude_t *dude, context_helper_t *helper, char *output,
                           char *reads);
void denoise_homopolymers(dude_t *dude, context_helper_t *helper, char *output,
                          char *reads);

void print_help();
void print_version();

// dude
ctx_table_t *create_table(int dim);

int *get_counts(ctx_table_t *table, uint64_t key, int add);

void destroy_table(ctx_table_t *table);

void window_init(window_t *w, dude_t *dude);

int window_push(window_t *w, int symbol);

uint64_t window_key(window_t *w);

int window_slot(window_t *w);

int get_symbol(dude_t *dude, char c);

void count_contexts(dude_t *dude, const char *seq, ctx_table_t *table);

char get_min_base(dude_t *dude, context_helper_t *helper, ctx_table_t *table,
                  window_t *w, int position, int quality);

kseq_t *denoise_seq(dude_t *dude, context_helper_t *helper, kseq_t *seq,
                    ctx_table_t *table);

// utilities
void write_seq(FILE *file, kseq_t *sequence);

void write_denoised_seq(FILE *file, kseq_t *denoised, kseq_t *sequence);

FILE *open_output(char *path);

gzFile open_reads(char *path);

int run_tests();

// ================================
//
//           Program Main
//
// ================================
int main(int argc, char *argv[]) {

  int k = DEFAULT_K, i = DEFAULT_I, l = DEFAULT_L, r = DEFAULT_R, c = DEFAULT_C;
  char *reads = NULL, *matrix = NULL, *output = NULL, *type = "read";

  init_base_codes();

  int p;
  while ((p = getopt(argc, argv, "q:m:k:i:o:c:l:r:M:tvh")) != -1)
    switch (p) {
    case 'q':
      reads = optarg;
      break;
    case 'm':
      matrix = optarg;
      break;
    case 'k':
      k = atoi(optarg);
      break;
    case 'i':
      i = atoi(optarg);
      break;
    case 'o':
      output = optarg;
      break;
    case 'l':
      l = atoi(optarg);
      break;
    case 'r':
      r = atoi(optarg);
      break;
    case 'c':
      c = atoi(optarg);
      // Mode { 0 - subs, 1 hom, 2 both {SH}, 3 {HS} };
      if (c > 3 || c < 0) {
        fprintf(stderr, "Invalid correction type: -c must be 0, 1, 2 or 3.\n");
        print_help();
        exit(EXIT_FAILURE);
      }
      break;
    case 'M':
      type = optarg;
      if (strncmp(type, "qual", 5) != 0 && strncmp(type, "read", 5) != 0) {
        fprintf(stderr, "Invalid matrix type: -M must be read or qual.\n");
        print_help();
        exit(EXIT_FAILURE);
      }
      break;
    case 't':
      return run_tests() ? EXIT_FAILURE : 0;
    case 'v':
      print_version();
      return 0;
    case 'h':
      print_help();
      return 0;
    case ':':
    case '?':
      print_help();
      exit(EXIT_FAILURE);
    default:
      print_help();
      exit(EXIT_FAILURE);
    }

  if (reads == NULL || matrix == NULL || output == NULL) {
    fprintf(stderr, "Missing arguments: -q, -m and -o are required.\n");
    print_help();
    exit(EXIT_FAILURE);
  }

  if (k < 1 || i < 1 || l < 1 || r < 1) {
    fprintf(stderr,
            "Invalid arguments: -k, -i, -l and -r must be greater than 0.\n");
    print_help();
    exit(EXIT_FAILURE);
  }

  if (k > MAX_K || i > MAX_I) {
    fprintf(stderr,
            "Invalid arguments: -k must be %d or less and -i must be %d or "
            "less.\n",
            MAX_K, MAX_I);
    print_help();
    exit(EXIT_FAILURE);
  }

  char mode_str[3] = "S-";
  if (c == 1) {
    mode_str[0] = 'H';
  } else if (c == 2) {
    mode_str[0] = 'S';
    mode_str[1] = 'H';
  } else if (c == 3) {
    mode_str[0] = 'H';
    mode_str[1] = 'S';
  }

  print_version();
  printf("\nArguments\n");
  printf("# reads                 : %s\n", reads);
  printf("# channel matrix        : %s\n", matrix);
  printf("# context length subs   : %d\n", k);
  printf("# context length indels : %d\n", i);
  printf("# max read length       : %d\n", l);
  printf("# matrix ranges size    : %d\n", r);
  printf("# correction type       : %s\n", mode_str);
  printf("# denoise type          : %s\n", type);

  printf("\nOutput\n");
  printf("# denoise substitutions : %s_subs.fastq\n", output);
  printf("# denoise indels        : %s_hom.fastq\n", output);

  // create the DUDE params and context helper
  // TODO-QUAL
  dude_t *dude = create_dude_params(k, i, l, r, type);

  load_matrices(matrix, dude);

  context_helper_t *helper = create_helper();

  char *subs_output_name = append_str(output, "_subs.fastq");
  char *indels_output_name = append_str(output, "_hom.fastq");

  // ========== Log file ============== //
  /*
  time_t rawtime;
  struct tm * timeinfo;
  time (&rawtime);
  timeinfo = localtime (&rawtime);

  char log_name[30];

  int ln = strftime(log_name, 80, "info-%Y-%m-%d_%H%M.log", timeinfo);
  if (ln < 0) {
       error_exit();
  }
  printf("Creating log file. File: %s\n", log_name);

  LOG_FILE = fopen(log_name, "w+");

  */
  // ========== Denoise =============== //
  printf("Starting the denoise\n");

  switch (c) {
  case 0:
    denoise_substitutions(dude, helper, subs_output_name, reads);
    break;
  case 1:
    denoise_homopolymers(dude, helper, indels_output_name, reads);
    break;
  case 2:
    // Substitutions and then Homopolymers
    denoise_substitutions(dude, helper, subs_output_name, reads);
    denoise_homopolymers(dude, helper, indels_output_name, subs_output_name);
    break;
  case 3:
    // Homopolymers and then Substitutions
    denoise_homopolymers(dude, helper, indels_output_name, reads);
    denoise_substitutions(dude, helper, subs_output_name, indels_output_name);
    break;
  }

  // TODO clean up dude
  destroy_helper(helper);

  // ========== Log file ============== //
  // fclose(LOG_FILE);

  fprintf(stdout, "Completed\n");

  return 0;
}

void print_help() {
  print_version();
  printf("\nDUDE denoiser for Ion Torrent reads.\n\n");
  printf("Usage: iondude -q reads.fastq -m matrices.txt -o denoised "
         "[options]\n\n");
  printf("Required:\n");
  printf("  -q <file>    reads, FASTQ (plain or gzip)\n");
  printf("  -m <file>    channel matrices made by the estimator\n");
  printf("  -o <prefix>  output prefix: <prefix>_subs.fastq and "
         "<prefix>_hom.fastq\n\n");
  printf("Options:\n");
  printf(
      "  -k <int>     context length for substitutions                  [%d]\n",
      DEFAULT_K);
  printf(
      "  -i <int>     context length for homopolymers                   [%d]\n",
      DEFAULT_I);
  printf(
      "  -l <int>     max read length, as used by the estimator         [%d]\n",
      DEFAULT_L);
  printf(
      "  -r <int>     width of the position ranges, as in the estimator [%d]\n",
      DEFAULT_R);
  printf(
      "  -c <int>     correction type                                   [%d]\n",
      DEFAULT_C);
  printf("                 0  substitutions only\n");
  printf("                 1  homopolymers only\n");
  printf("                 2  substitutions, then homopolymers\n");
  printf("                 3  homopolymers, then substitutions\n");
  printf("  -M <str>     matrices by read position (read) or by base quality "
         "(qual) [read]\n");
  printf("  -t           run the self-tests\n");
  printf("  -v           print the version\n");
  printf("  -h           print this help\n");
}

void print_version() {
  printf("iondude %s (branch %s, commit %s)\n", VERSION, GIT_BRANCH,
         GIT_COMMIT_HASH);
}

// ====================================================
//                  SUBSTITUTIONS
//=====================================================
void denoise_substitutions(dude_t *dude, context_helper_t *helper, char *output,
                           char *reads) {

  dude->mode = SUBSTITUTIONS;

  // allocate a new hash table
  ctx_table_t *table_subs = create_table(SUBS_MATRIX_DIM);

  // set auxiliary counters
  helper->counter = (float *)calloc(SUBS_MATRIX_DIM, sizeof(float));
  if (!helper->counter) {
    error_exit();
  }

  helper->counters_inv = (float *)calloc(SUBS_MATRIX_DIM, sizeof(float));
  if (!helper->counters_inv) {
    error_exit();
  }

  FILE *subs_output = open_output(output);

  // Initialize Input Sequence
  gzFile fp = open_reads(reads);
  kseq_t *seq = kseq_init(fp);
  int l;

  // count the contexts for substitutions
  // windows with a base that is not A, C, G or T are not counted
  while ((l = kseq_read(seq)) >= 0) {
    count_contexts(dude, seq->seq.s, table_subs);
  }

  // rewind the file and start the denoise process
  gzrewind(fp);
  kseq_rewind(seq);

  // Denoise substitutions
  while ((l = kseq_read(seq)) >= 0) {
    // iondude does not correct reads with uncalled, ambiguous or lowercase
    // bases
    if (is_acgt_seq(seq->seq.s)) {
      kseq_t *denoised = denoise_seq(dude, helper, seq, table_subs);
      if (denoised != NULL) {
        write_denoised_seq(subs_output, denoised, seq);
        kseq_free(denoised);
      } else {
        fprintf(stderr,
                KRED
                "[subs denoise] sequence '%s' shorter than threshold\n" KNRM,
                seq->name.s);
        // save as it is
        write_seq(subs_output, seq);
      }
    } else {
      // save as it is
      write_seq(subs_output, seq);
    }
  }

  // close
  fclose(subs_output);

  // clean up the table
  destroy_table(table_subs);

  // close the original file
  kseq_destroy(seq);
  gzclose(fp);

  // reset helper
  free(helper->counter);
  free(helper->counters_inv);
}

// =================================================
//                  HOMOPOLYMERS
//==================================================
void denoise_homopolymers(dude_t *dude, context_helper_t *helper, char *output,
                          char *reads) {

  dude->mode = HOMOPOLYMERS;

  ctx_table_t *table_hom = create_table(HOM_MATRIX_DIM);

  gzFile fp = open_reads(reads);
  kseq_t *seq = kseq_init(fp);

  FILE *hom_output = open_output(output);

  helper->counter = (float *)calloc(HOM_MATRIX_DIM, sizeof(float));
  helper->counters_inv = (float *)calloc(HOM_MATRIX_DIM, sizeof(float));

  int l = 0;

  // count the contexts for homopolymers
  while ((l = kseq_read(seq)) >= 0) {
    // iondude does not correct reads with uncalled, ambiguous or lowercase
    // bases
    if (is_acgt_seq(seq->seq.s)) {
      // the sequence is encoded [homopolymers are encoded]
      kseq_t *encoded = encode_homopolymers_seq(seq);

      count_contexts(dude, encoded->seq.s, table_hom);

      kseq_free(encoded);
    }
  }

  gzrewind(fp);
  kseq_rewind(seq);

  // denoise homopolymers
  while ((l = kseq_read(seq)) >= 0) {
    if (is_acgt_seq(seq->seq.s)) {
      kseq_t *encoded_hom = encode_homopolymers_seq(seq);
      kseq_t *denoised = denoise_seq(dude, helper, encoded_hom, table_hom);

      if (denoised == NULL) {
        kseq_free(encoded_hom);
        fprintf(stderr,
                "[homopolymers denoise] sequence '%s' denoise not performed, "
                "sequence too short\n",
                seq->name.s);
        // save as it is
        write_seq(hom_output, seq);
        continue;
      }

      kseq_t *decoded = decode_homopolymers_seq(denoised);

      if (decoded != NULL) {
        write_denoised_seq(hom_output, decoded, seq);
        kseq_free(denoised);
        kseq_free(decoded);
        kseq_free(encoded_hom);
      } else {
        kseq_free(denoised);
        kseq_free(encoded_hom);
        fprintf(stderr, "[homopolymers denoise] sequence '%s' denoise failed\n",
                seq->name.s);
        // save as it is
        write_seq(hom_output, seq);
      }
    } else {
      // save as it is
      write_seq(hom_output, seq);
    }
  }

  // close
  fclose(hom_output);
  kseq_destroy(seq);
  gzclose(fp);

  // clean up the table
  destroy_table(table_hom);

  // reset helper
  free(helper->counter);
  free(helper->counters_inv);
}

// create an empty context table, dim counts for each key
ctx_table_t *create_table(int dim) {
  ctx_table_t *table = (ctx_table_t *)calloc(1, sizeof(ctx_table_t));
  if (table == NULL) {
    error_exit();
  }
  table->index = kh_init(ctx);
  table->dim = dim;
  return table;
}

// get the count vector of a key
// if the key is not in the table: add it when add is 1, else return NULL
// the pointer is valid until the next call that adds a key
int *get_counts(ctx_table_t *table, uint64_t key, int add) {
  khint_t it = kh_get(ctx, table->index, key);
  if (it != kh_end(table->index)) {
    return table->counts + (size_t)kh_val(table->index, it) * table->dim;
  }
  if (!add) {
    return NULL;
  }

  // grow the count vectors array
  if (table->n == table->m) {
    table->m = table->m ? table->m * 2 : 1024;
    int *grown = (int *)realloc(table->counts,
                                (size_t)table->m * table->dim * sizeof(int));
    if (grown == NULL) {
      error_exit();
    }
    table->counts = grown;
  }

  int ret;
  it = kh_put(ctx, table->index, key, &ret);
  if (ret < 0) {
    error_exit();
  }
  kh_val(table->index, it) = table->n;

  int *counts = table->counts + (size_t)table->n * table->dim;
  memset(counts, 0, table->dim * sizeof(int));
  table->n++;
  return counts;
}

// free the memory of a context table
void destroy_table(ctx_table_t *table) {
  kh_destroy(ctx, table->index);
  free(table->counts);
  free(table);
}

// start an empty window for the current mode
void window_init(window_t *w, dude_t *dude) {
  int subs = (dude->mode == SUBSTITUTIONS);
  w->k = subs ? dude->k : dude->h;
  w->bits = subs ? 2 : 6;
  // one side has at most 32 bits (MAX_K and MAX_I)
  w->side_mask = (1ULL << (w->bits * w->k)) - 1;
  w->left = 0;
  w->right = 0;
  w->centre = 0;
  w->filled = 0;
}

// shift the next symbol into the right side of the window
// the oldest symbol of the right side becomes the centre, and the old centre
// goes into the left side
// a symbol -1 (not A, C, G or T) empties the window
// returns 1 when the window has 2k + 1 symbols
int window_push(window_t *w, int symbol) {
  if (symbol < 0) {
    w->filled = 0;
    return 0;
  }
  int b = w->bits, k = w->k;
  w->left = ((w->left << b) | (uint64_t)w->centre) & w->side_mask;
  w->centre = (int)(w->right >> (b * (k - 1))) & ((1 << b) - 1);
  w->right = ((w->right << b) | (uint64_t)symbol) & w->side_mask;
  w->filled++;
  return w->filled > 2 * k;
}

// the key is the left and the right sides, without the centre
// homopolymers: the base of the centre run is also in the key, because the
// denoise only changes the run length
uint64_t window_key(window_t *w) {
  uint64_t key = (w->left << (w->bits * w->k)) | w->right;
  if (w->bits == 2) {
    return key;
  }
  return (key << 2) | (uint64_t)(w->centre >> 4);
}

// the index of the centre in the count vector: the base code (substitutions)
// or the run length - 1 (homopolymers)
int window_slot(window_t *w) {
  return (w->bits == 2) ? w->centre : w->centre & 15;
}

// the symbol of a sequence char for the current mode, -1 if not valid
int get_symbol(dude_t *dude, char c) {
  if (dude->mode == SUBSTITUTIONS) {
    return base_code[(unsigned char)c];
  }
  return hom_symbol(c);
}

void count_contexts(dude_t *dude, const char *seq, ctx_table_t *table) {
  // go over each base of the sequence and increase the counter of the centre
  // symbol for the context. Windows with a base that is not A, C, G or T are
  // not counted
  window_t w;
  window_init(&w, dude);

  for (int i = 0; seq[i] != '\0'; i++) {
    if (window_push(&w, get_symbol(dude, seq[i]))) {
      get_counts(table, window_key(&w), 1)[window_slot(&w)]++;
    }
  }
}

// Denoise every base on the sequence
kseq_t *denoise_seq(dude_t *dude, context_helper_t *helper, kseq_t *seq,
                    ctx_table_t *table) {

  int k = (dude->mode == SUBSTITUTIONS) ? dude->k : dude->h;

  // do not touch if shorter
  if (seq->seq.l == 0 || (k * 2) + 1 > seq->seq.l) {
    return NULL;
  }

  kseq_t *ret = kseq_dup(seq);

  window_t w;
  window_init(&w, dude);

  // for every base on the sequence, obtain the base that minimize
  for (int i = 0; i < seq->seq.l; i++) {
    if (!window_push(&w, get_symbol(dude, seq->seq.s[i]))) {
      continue;
    }
    int centre = i - k;

    // TODO: - optional threshold
    //       - configurable fastq quality type
    // 1 is PHRED is used for iontorrent
    int quality = get_quality_num(seq->qual.s[centre], 1, centre);

    // extreme high quality bases aren't denoised
    if (quality >= MAX_QUAL) {
      continue;
    }

    // save the base that minimize
    ret->seq.s[centre] = get_min_base(dude, helper, table, &w, centre, quality);
  }

  return ret;
}

// obtain the base that minimize
char get_min_base(dude_t *dude, context_helper_t *helper, ctx_table_t *table,
                  window_t *w, int position, int quality) {

  int subs = (dude->mode == SUBSTITUTIONS);
  int upper_bound = table->dim;

  // the column vector with the context counts for each possible base
  int *counts = get_counts(table, window_key(w), 0);
  for (int i = 0; i < upper_bound; i++) {
    helper->counter[i] = counts ? counts[i] : 0;
  }

  // observed base, or homopolymer length - 1
  int encoded_base = window_slot(w);

  // The matrix index depends on the running type
  int position_index = 0;
  if (dude->type == POSITION) {
    position_index = get_index_for_matrix(position, dude);
  } else {
    position_index = quality;
  }

  // row vector m^T PI^-1, the counters times the inverted channel matrix
  for (int x = 0; x < upper_bound; x++) {
    float value = 0;
    for (int y = 0; y < upper_bound; y++) {
      if (subs) {
        value +=
            helper->counter[y] * dude->subs_channel_i[position_index][y][x];
      } else {
        value += helper->counter[y] * dude->hom_channel_i[position_index][y][x];
      }
    }
    helper->counters_inv[x] = value;
  }

  // find the base that minimize
  float min = FLT_MAX;
  int base_that_minimize = 0;
  for (int argmin_symbol = 0; argmin_symbol < upper_bound; argmin_symbol++) {

    // Has the cached product wise values of the channel transition matrix and
    // the loss matrix Multiply both columns to get and if lower than the
    // minimum then store
    float innerMin = 0;
    for (int n = 0; n < upper_bound; n++) {
      if (subs) {
        innerMin += helper->counters_inv[n] *
                    dude->subs_schur_product[position_index][encoded_base]
                                            [argmin_symbol][n];
      } else {
        innerMin += helper->counters_inv[n] *
                    dude->hom_schur_product[position_index][encoded_base]
                                           [argmin_symbol][n];
      }
    }

    if (innerMin < min) {
      base_that_minimize = argmin_symbol;
      min = innerMin;
    }
  }

  if (subs) {
    return "ACGT"[base_that_minimize];
  }
  // homopolymers: the same base with the new run length
  return encode_h_base("ACGT"[w->centre >> 4], base_that_minimize);
}

// open the output file, exit on error
FILE *open_output(char *path) {
  FILE *file = fopen(path, "w+");
  if (file == NULL) {
    fprintf(stderr, "Could not open the output file '%s'.\n", path);
    exit(EXIT_FAILURE);
  }
  return file;
}

// open the reads file (plain or gzip), exit on error
gzFile open_reads(char *path) {
  gzFile file = gzopen(path, "r");
  if (file == NULL) {
    fprintf(stderr, "Could not open the reads file '%s'.\n", path);
    exit(EXIT_FAILURE);
  }
  return file;
}

// write sequence
void write_seq(FILE *file, kseq_t *sequence) {
  if (sequence->name.s) {
    fprintf(file, "@%s\n", sequence->name.s);
  } else {
    fprintf(file, "@%s\n", "failed");
  }
  if (sequence->seq.s) {
    fprintf(file, "%s\n", sequence->seq.s);
  } else {
    fprintf(file, "%s\n", "failed");
  }

  if (sequence->comment.s) {
    fprintf(file, "+%s\n", sequence->comment.s);
  } else {
    fputs("+\n", file);
  }
  fprintf(file, "%s\n", sequence->qual.s);
}

// write sequence, use the metadata from the original sequence
void write_denoised_seq(FILE *file, kseq_t *denoised, kseq_t *sequence) {
  if (sequence->name.s) {
    fprintf(file, "@%s\n", sequence->name.s);
  } else {
    fprintf(file, "@%s\n", "failed");
  }
  if (denoised->seq.s) {
    fprintf(file, "%s\n", denoised->seq.s);
  } else {
    fprintf(file, "%s\n", "failed");
  }

  if (sequence->comment.s) {
    fprintf(file, "+%s\n", sequence->comment.s);
  } else {
    fputs("+\n", file);
  }
  fprintf(file, "%s\n", denoised->qual.s);
}

// ================================================ //
//         Encode/Decode homopolymers RLE Tests    //
// ============================================== //
// returns the number of failed reads
int test_rle() {

  gzFile fp;
  kseq_t *seq;

  fp = gzopen("unit_test/data.fastq", "r");
  seq = kseq_init(fp);

  int i = 0, s = 0, f = 0;
  while (kseq_read(seq) >= 0) {

    kseq_t *encoded = encode_homopolymers_seq(seq);
    kseq_t *decoded = decode_homopolymers_seq(encoded);

    int test_success = strlen(decoded->seq.s) > 0 &&
                       strcmp(decoded->seq.s, seq->seq.s) == 0 &&
                       strlen(decoded->qual.s) == strlen(decoded->seq.s);

    if (test_success) {
      s++;
    } else {
      f++;
    }

    i++;

    free(encoded->name.s);
    free(encoded->comment.s);
    free(encoded->seq.s);
    free(encoded->qual.s);
    free(encoded);
    free(decoded->name.s);
    free(decoded->comment.s);
    free(decoded->seq.s);
    free(decoded->qual.s);
    free(decoded);
  }

  if (f == 0 && s == i) {
    printf("RLE tests passed\n");
  }

  if (f != 0 && s != i) {
    fprintf(stderr, "some RLE test failed\n");
  }

  if (f != 0 && s == 0) {
    fprintf(stderr, "every RLE test failed\n");
  }

  kseq_destroy(seq);
  gzclose(fp);
  return f;
}

// ================================================ //
//              Context encoding Tests              //
// ============================================== //
// compare the keys and slots of all the full windows of a read with the
// expected values. HOMOPOLYMERS: the read is run-length encoded first
// returns 1 if the test failed
int test_keys(const char *name, denoise_m mode, int k, char *read, int n,
              const uint64_t *keys, const int *slots) {
  dude_t dude;
  memset(&dude, 0, sizeof dude);
  dude.mode = mode;
  dude.k = k;
  dude.h = k;

  char *seq = read;
  kseq_t *encoded = NULL;
  if (mode == HOMOPOLYMERS) {
    // the quality values are not used, the read is a valid placeholder
    kseq_t raw;
    memset(&raw, 0, sizeof raw);
    raw.seq.s = read;
    raw.seq.l = strlen(read);
    raw.qual = raw.seq;
    encoded = encode_homopolymers_seq(&raw);
    seq = encoded->seq.s;
  }

  window_t w;
  window_init(&w, &dude);

  int found = 0, failed = 0;
  for (int i = 0; seq[i] != '\0'; i++) {
    if (!window_push(&w, get_symbol(&dude, seq[i]))) {
      continue;
    }
    if (found >= n || window_key(&w) != keys[found] ||
        window_slot(&w) != slots[found]) {
      failed = 1;
    }
    found++;
  }

  if (failed || found != n) {
    fprintf(stderr, "test '%s' failed\n", name);
    failed = 1;
  }
  kseq_free(encoded);
  return failed;
}

// count the contexts of 3 reads and check the count vectors
// returns 1 if the test failed
int test_counts() {
  dude_t dude;
  memset(&dude, 0, sizeof dude);
  dude.mode = SUBSTITUTIONS;
  dude.k = 2;

  ctx_table_t *table = create_table(SUBS_MATRIX_DIM);
  count_contexts(&dude, "ACGTTGCA", table);
  count_contexts(&dude, "ACTTTGCA", table);
  // no full window without the N
  count_contexts(&dude, "ACNTT", table);

  // AC_TT (31) has G once and T once, TT_CA (244) has G twice
  int *c31 = get_counts(table, 31, 0);
  int *c244 = get_counts(table, 244, 0);
  int failed = table->n != 6 || get_counts(table, 0, 0) != NULL ||
               c31 == NULL || c31[0] != 0 || c31[1] != 0 || c31[2] != 1 ||
               c31[3] != 1 || c244 == NULL || c244[0] != 0 || c244[1] != 0 ||
               c244[2] != 2 || c244[3] != 0;

  if (failed) {
    fprintf(stderr, "test 'context counts' failed\n");
  }
  destroy_table(table);
  return failed;
}

// returns the number of failed tests
int test_encoding() {
  int f = 0;

  // the examples of the plan (encoding_demo.c)
  f += test_keys("substitutions k 2", SUBSTITUTIONS, 2, "ACGTTGCA", 4,
                 (uint64_t[]){31, 110, 185, 244}, (int[]){2, 3, 3, 2});
  // runs A3 C1 G2 T1: centre C1 (key 645, length 1) and G2 (key 4290, length 2)
  f += test_keys("homopolymers i 1", HOMOPOLYMERS, 1, "AAACGGT", 2,
                 (uint64_t[]){645, 4290}, (int[]){0, 1});

  // k 16: the key has 64 bits, the first base (C) is in the 2 high bits
  f += test_keys("substitutions k 16", SUBSTITUTIONS, 16,
                 "CAAAAAAAAAAAAAAAGAAAAAAAAAAAAAAAT", 1,
                 (uint64_t[]){(1ULL << 62) | 3}, (int[]){2});
  f += test_keys("substitutions k 16 all T", SUBSTITUTIONS, 16,
                 "TTTTTTTTTTTTTTTTGTTTTTTTTTTTTTTTT", 1,
                 (uint64_t[]){UINT64_MAX}, (int[]){2});
  // i 5: runs T10 G10 T10 G10 T10 [C3] T10 G10 T10 G10 T10, the key has 62 bits
  f += test_keys("homopolymers i 5", HOMOPOLYMERS, 5,
                 "TTTTTTTTTTGGGGGGGGGGTTTTTTTTTTGGGGGGGGGGTTTTTTTTTT"
                 "CCC"
                 "TTTTTTTTTTGGGGGGGGGGTTTTTTTTTTGGGGGGGGGGTTTTTTTTTT",
                 1, (uint64_t[]){4154459029623433701ULL}, (int[]){2});

  // bases that are not A, C, G or T empty the window
  f += test_keys("N before", SUBSTITUTIONS, 2, "ACNACGTTGCA", 4,
                 (uint64_t[]){31, 110, 185, 244}, (int[]){2, 3, 3, 2});
  f += test_keys("N inside", SUBSTITUTIONS, 2, "ACGTTNGCA", 1,
                 (uint64_t[]){31}, (int[]){2});
  f += test_keys("lowercase", SUBSTITUTIONS, 2, "ACGTtGCA", 0, NULL, NULL);
  f += test_keys("short read", SUBSTITUTIONS, 2, "ACGT", 0, NULL, NULL);

  f += test_counts();

  if (f == 0) {
    printf("context encoding tests passed\n");
  }
  return f;
}

// returns the number of failed tests
int run_tests() {
  int f = test_rle() + test_encoding();
  if (f == 0) {
    printf("every test passed, nice\n");
  }
  return f;
}
