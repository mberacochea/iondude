#ifndef DUDE
#define DUDE

#define SUBS_MATRIX_DIM 4
#define HOM_MATRIX_DIM 10

// QUAL mode: one matrix per Phred score from 0 to MAX_QUAL - 1
// bases with a higher score are not denoised
#define MAX_QUAL 50

// flags for the matrices loaded on each slot
#define LOADED_SUBS 1
#define LOADED_SUBS_I 2
#define LOADED_HOM 4
#define LOADED_HOM_I 8

// get_matrix_slot value for a matrix without a range
#define ALL_SLOTS -2

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

  // matrices position ranges
  int range;
  int range_max;

  // number of matrices per type (range_max for POSITION, MAX_QUAL for QUAL)
  int slots;
  // LOADED_* flags for each slot
  int *loaded;

  // substitutions context length
  int k;
  // indels context length
  int h;

  // by position in read OR by quality score of the base
  enum denoise_type type;

  // substitutions or homopolymers
  enum denoise_mode mode;
} dude_t;

// context helper struct
typedef struct {
  // context
  char *context;
  // context modified with other base in the i position
  char *context_modified;
  // auxiliary counters
  float *counter;
  float *counters_inv;
} context_helper_t;

void load_matrices(char *estimation, dude_t *dude);

void load_substitution_matrix(dude_t *dude, FILE *file, char *header,
                              int isChannel);

void load_homopolymers_matrix(dude_t *dude, FILE *file, char *header,
                              int isChannel);

void complete_matrices(dude_t *dude);

void validate_matrices(dude_t *dude, int slot, float *channel, float *inverse,
                       int dim);

dude_t *create_dude_params(int k, int h, int read_len, int range, char *mode);

int get_matrix_slot(dude_t *dude, char *header);

context_helper_t *create_helper(int k);

void destroy_helper(context_helper_t *helper);

int get_index_for_matrix(int position_in_read, dude_t *dude);

void error_exit();

#endif
