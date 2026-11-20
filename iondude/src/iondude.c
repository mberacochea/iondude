#include <float.h>
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

// Hash table init
// name part of init must be unique for the key, value types.
// in this instance 34 is arbitrary symbolic name for a hashtable
// that contains string keys and int values.
const int khStrInt = 34;
KHASH_MAP_INIT_STR(khStrInt, int)

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
void count_contexts(dude_t *dude, context_helper_t *helper, char *seq,
                    khash_t(khStrInt) * table);

int get_context(dude_t *dude, context_helper_t *helper, int position,
                char *seq);

int load_modified_context(dude_t *dude, context_helper_t *helper, char base);

char get_min_base(dude_t *dude, context_helper_t *helper,
                  khash_t(khStrInt) * table, char base, int position,
                  int quality);

kseq_t *denoise_seq(dude_t *dude, context_helper_t *helper, kseq_t *seq,
                    khash_t(khStrInt) * table);

// utilities
void write_seq(FILE *file, kseq_t *sequence);

void write_denoised_seq(FILE *file, kseq_t *denoised, kseq_t *sequence);

void deallocate_table(khash_t(khStrInt) * table);

FILE *open_output(char *path);

gzFile open_reads(char *path);

void run_tests();

// ================================
//
//           Program Main
//
// ================================
int main(int argc, char *argv[]) {

  int k = DEFAULT_K, i = DEFAULT_I, l = DEFAULT_L, r = DEFAULT_R, c = DEFAULT_C;
  char *reads = NULL, *matrix = NULL, *output = NULL, *type = "read";

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
      run_tests();
      return 0;
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

  context_helper_t *helper = create_helper(k);

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
  khash_t(khStrInt) *table_subs = kh_init(khStrInt);

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
  while ((l = kseq_read(seq)) >= 0) {
    count_contexts(dude, helper, seq->seq.s, table_subs);
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
  deallocate_table(table_subs);

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

  khash_t(khStrInt) *table_hom = kh_init(khStrInt);

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

      count_contexts(dude, helper, encoded->seq.s, table_hom);

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
  deallocate_table(table_hom);

  // reset helper
  free(helper->counter);
  free(helper->counters_inv);
}

void count_contexts(dude_t *dude, context_helper_t *helper, char *seq,
                    khash_t(khStrInt) * table) {
  // This method will go over each base on the provided Sequence
  // and will encode the context and increase the counter for the context on the
  // table
  int len = strlen(seq), i = 0;

  khint_t hash_value;

  int k = (dude->mode == SUBSTITUTIONS) ? dude->k : dude->h;

  // sequences shorter than the provided K value are ignored
  if (len < k * 2) {
    return;
  }

  for (i = k; i < len - k; i++) {

    int ctx_get = get_context(dude, helper, i, seq);

    if (ctx_get == 0) {
      fprintf(stderr, "[context counts] sequence '%s' and position '%d' failed",
              seq, i);
      continue;
    }

    // query the hash table
    hash_value = kh_get(khStrInt, table, helper->context);

    if (hash_value != kh_end(table)) {
      int counter = kh_val(table, hash_value);
      kh_value(table, hash_value) = counter + 1;
    } else {
      int ret;
      hash_value = kh_put(khStrInt, table, strdup(helper->context), &ret);
      if (ret) {
        kh_value(table, hash_value) = 1;
      }
    }
  }
}

// change the base on a
int load_modified_context(dude_t *dude, context_helper_t *helper, char base) {
  // clean and change
  int k = (dude->mode == SUBSTITUTIONS) ? dude->k : dude->h;

  int i = 0;
  for (; i <= k * 2; i++) {
    helper->context_modified[i] = helper->context[i];
  }
  helper->context_modified[i] = '\0';
  helper->context_modified[k] = base;
  return 1;
}

// get the context of a sequence from a starting point
int get_context(dude_t *dude, context_helper_t *helper, int position,
                char *seq) {

  // k [Base] k + 1
  int k = (dude->mode == SUBSTITUTIONS) ? dude->k : dude->h;

  int start = position - k, end = position + k + 1;

  int new_ctx_len = end - start;

  if (k * 2 > new_ctx_len) {
    return 0;
  }

  int n = 0, j = 0;
  for (n = 0, j = start; j < end; n++, j++) {
    helper->context[n] = seq[j];
  }

  helper->context[n] = '\0';

  return 1;
}

// Denoise every base on the sequence
kseq_t *denoise_seq(dude_t *dude, context_helper_t *helper, kseq_t *seq,
                    khash_t(khStrInt) * table) {

  int k = (dude->mode == SUBSTITUTIONS) ? dude->k : dude->h;

  // do not touch if shorter
  if (seq->seq.l == 0 || (k * 2) + 1 > seq->seq.l) {
    return NULL;
  }

  kseq_t *ret = kseq_dup(seq);

  // for every base on the sequence, obtain the base that minimize
  for (int i = k; i < seq->seq.l - k; i++) {
    // TODO: - optional threshold
    //       - configurable fastq quality type
    // 1 is PHRED is used for iontorrent
    int quality = get_quality_num(seq->qual.s[i], 1, i);

    // extreme high quality bases aren't denoised
    if (quality >= MAX_QUAL) {
      ret->seq.s[i] = seq->seq.s[i];
      continue;
    }

    // get the context for this base
    if (get_context(dude, helper, i, seq->seq.s)) {
      // save the base that minimize
      ret->seq.s[i] =
          get_min_base(dude, helper, table, seq->seq.s[i], i, quality);
    } else {
      ret->seq.s[i] = seq->seq.s[i];
      fprintf(stderr,
              KRED "[denoise seq] context for seq '%s' and position '%d'",
              seq->name.s, i);
    }
  }

  return ret;
}

// obtain the base that minimize
char get_min_base(dude_t *dude, context_helper_t *helper,
                  khash_t(khStrInt) * table, char base, int position,
                  int quality) {

  // build the occurrences col vector
  int subs = (dude->mode == SUBSTITUTIONS);
  int upper_bound = (subs) ? SUBS_MATRIX_DIM : HOM_MATRIX_DIM;

  khint_t hash_value;

  // get the start and finish points
  int offset = 0;
  if (!subs) {
    if (base >= T_HOM) {
      offset = T_HOM;
    } else if (base >= G_HOM) {
      offset = G_HOM;
    } else if (base >= C_HOM) {
      offset = C_HOM;
    } else if (base >= A_HOM) {
      offset = A_HOM;
    }
  }

  // build the column vector with the context counts for each possible base
  int cwb = 0;
  for (int i = 0; i < upper_bound; i++) {
    if (subs) {
      cwb = load_modified_context(dude, helper, decode_s_base(i));
    } else {
      cwb = load_modified_context(dude, helper, (char)(offset + i));
    }

    int val = 0;
    if (cwb == 1) {
      // get counter value
      hash_value = kh_get(khStrInt, table, helper->context_modified);
      if (hash_value != kh_end(table)) {
        val = kh_val(table, hash_value);
      }
    }
    helper->counter[i] = val;
  }

  // find the base that minimize
  float min = FLT_MAX;
  int base_that_minimize = 0;
  int encoded_base = 0;
  if (subs) {
    // observed base
    encoded_base = encode_s_base(base);
  } else {
    // homopolymer length
    encoded_base = base - offset;
  }

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

  // argmin
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

    // printf("\n%.6f", innerMin);

    if (innerMin < min) {
      base_that_minimize = (subs) ? argmin_symbol : offset + argmin_symbol;
      min = innerMin;
    }
  }

  if (subs) {
    return decode_s_base(base_that_minimize);
  } else {
    return base_that_minimize;
  }
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

// free the memory of a hash table
void deallocate_table(khash_t(khStrInt) * table) {
  // clean up the table
  khint_t kv;
  for (kv = kh_begin(table); kv != kh_end(table); ++kv) {
    if (kh_exist(table, kv)) {
      free((char *)kh_key(table, kv));
    }
  }
  kh_destroy(khStrInt, table);
}

// ================================================ //
//         Encode/Decode homopolymers RLE Tests    //
// ============================================== //
void run_tests() {

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
    printf("every test passed, nice\n");
  }

  if (f != 0 && s != i) {
    fprintf(stderr, "some test failed\n");
  }

  if (f != 0 && s == 0) {
    fprintf(stderr, "every test failed");
  }

  kseq_destroy(seq);
  gzclose(fp);
}
