#ifndef DUDE
#define DUDE

#include "../lib/cjson/cJSON.h"

#define SUBS_MATRIX_DIM 4
#define HOM_MATRIX_DIM 10

// max context lengths: the context key must fit in 64 bits
// substitutions: 2k bases of 2 bits
#define MAX_K 16
// homopolymers: 2i runs of 6 bits and the base of the centre run (2 bits)
#define MAX_I 5

// POSITION type: bases with a Phred score of MAX_QUAL or more are not denoised
// QUAL type: the max_length of the matrices file gives the limit
#define MAX_QUAL 50

// the matrices file format (iondude/matrices/README.md)
#define MATRIX_FORMAT "iondude-matrices"
#define MATRIX_FORMAT_VERSION 1

// size of the error messages of the matrices file checks
#define MATRIX_ERROR_SIZE 256

// Colors
#define KNRM "\x1B[0m"
#define KRED "\x1B[31m"

typedef enum denoise_mode {
  SUBSTITUTIONS,
  HOMOPOLYMERS,
} denoise_m;

typedef enum denoise_type {
  POSITION,
  QUAL,
} denoise_t;

typedef float subs_matrix_type[SUBS_MATRIX_DIM][SUBS_MATRIX_DIM];
typedef float hom_matrix_type[HOM_MATRIX_DIM][HOM_MATRIX_DIM];

// this types are used to store the cached product wise multipicaton of each
// combination of loss fn matrix col and the transition matrix col
typedef float subs_schur_product[SUBS_MATRIX_DIM][SUBS_MATRIX_DIM]
                                [SUBS_MATRIX_DIM];
typedef float hom_schur_product[HOM_MATRIX_DIM][HOM_MATRIX_DIM][HOM_MATRIX_DIM];

typedef struct {
  // channel matrix values and the inverted matrix
  // loss function matrix

  // substitutions
  subs_matrix_type *subs_channel;
  subs_matrix_type *subs_channel_i;
  subs_matrix_type subs_loss;
  subs_schur_product *subs_schur_product;

  // indels on homopolymers
  hom_matrix_type *hom_channel;
  hom_matrix_type *hom_channel_i;
  hom_matrix_type hom_loss;
  hom_schur_product *hom_schur_product;

  // width of each slot: read positions (POSITION) or 1 (QUAL)
  int range;

  // number of matrices of each kind
  // POSITION: max read length / range, rounded up. QUAL: max quality score + 1
  int slots;

  // substitutions context length
  int k;
  // indels context length
  int h;

  // by position in read OR by quality score of the base
  enum denoise_type type;

  // substitutions or homopolymers
  enum denoise_mode mode;
} dude_t;

// the values of the matrices file that create_dude_params needs
typedef struct {
  enum denoise_type type;
  int range;
  int max_length;
  // the source fields, only to print them
  char tool[64];
  char version[64];
  char date[16];
  char input[1024];
} matrix_meta_t;

// context helper struct
typedef struct {
  // auxiliary counters
  float *counter;
  float *counters_inv;
} context_helper_t;

// read and parse a JSON file, exit on an error. Free it with cJSON_Delete
cJSON *parse_matrix_file(const char *path);

// read the metadata of a matrices file, exit on an error
matrix_meta_t read_matrix_meta(const char *path);

// read the matrices of a file into dude, exit on an error
void load_matrices(const char *path, dude_t *dude);

// check the metadata of a parsed matrices file
// returns 0, or 1 with the message in error (MATRIX_ERROR_SIZE)
int parse_matrix_meta(const cJSON *root, matrix_meta_t *meta, char *error);

// check the matrices of a parsed matrices file and copy them into dude
// returns 0, or 1 with the message in error (MATRIX_ERROR_SIZE)
int parse_matrices(const cJSON *root, dude_t *dude, char *error);

dude_t *create_dude_params(int k, int h, matrix_meta_t meta);

void destroy_dude_params(dude_t *dude);

context_helper_t *create_helper(void);

void destroy_helper(context_helper_t *helper);

int get_index_for_matrix(int position_in_read, dude_t *dude);

void error_exit();

#endif
