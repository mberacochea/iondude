#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dude.h"

// max difference allowed on the matrices validation
#define MATRIX_TOLERANCE 1e-3

static const char *SUBS_SYMBOLS[SUBS_MATRIX_DIM] = {"A", "C", "G", "T"};

// write the message in error, returns 1
static int fail(char *error, const char *format, ...) {
  va_list args;
  va_start(args, format);
  vsnprintf(error, MATRIX_ERROR_SIZE, format, args);
  va_end(args);
  return 1;
}

static void matrix_error(const char *path, const char *error) {
  fprintf(stderr, "Invalid matrices file '%s': %s\n", path, error);
  exit(1);
}

cJSON *parse_matrix_file(const char *path) {
  FILE *file = fopen(path, "rb");
  if (!file) {
    fprintf(stderr, "Could not open the matrices file '%s'.\n", path);
    exit(1);
  }

  char *text = NULL;
  size_t size = 0, used = 0, n;
  do {
    if (used + 4096 + 1 > size) {
      size = 2 * size + 4096 + 1;
      text = realloc(text, size);
      if (text == NULL) {
        error_exit();
      }
    }
    n = fread(text + used, 1, size - used - 1, file);
    used += n;
  } while (n > 0);
  text[used] = '\0';
  fclose(file);

  cJSON *root = cJSON_Parse(text);
  free(text);
  if (root == NULL) {
    matrix_error(path, "not a valid JSON file");
  }
  return root;
}

// an integer value of the field, as a JSON number
static int get_int(const cJSON *object, const char *name, int *value) {
  const cJSON *item = cJSON_GetObjectItemCaseSensitive(object, name);
  if (!cJSON_IsNumber(item) || item->valuedouble < INT_MIN ||
      item->valuedouble > INT_MAX ||
      item->valuedouble != (int)item->valuedouble) {
    return 0;
  }
  *value = (int)item->valuedouble;
  return 1;
}

// copy a string field to out (size bytes)
static int get_string(const cJSON *object, const char *name, char *out,
                      size_t size) {
  const char *value =
      cJSON_GetStringValue(cJSON_GetObjectItemCaseSensitive(object, name));
  if (value == NULL) {
    return 0;
  }
  snprintf(out, size, "%s", value);
  return 1;
}

// YYYY-MM-DD
static int is_iso_date(const char *date) {
  if (strlen(date) != 10 || date[4] != '-' || date[7] != '-') {
    return 0;
  }
  for (int i = 0; i < 10; i++) {
    if (i != 4 && i != 7 && !isdigit((unsigned char)date[i])) {
      return 0;
    }
  }
  return 1;
}

int parse_matrix_meta(const cJSON *root, matrix_meta_t *meta, char *error) {
  memset(meta, 0, sizeof *meta);

  if (!cJSON_IsObject(root)) {
    return fail(error, "file: must be a JSON object");
  }

  const char *format =
      cJSON_GetStringValue(cJSON_GetObjectItemCaseSensitive(root, "format"));
  if (format == NULL || strcmp(format, MATRIX_FORMAT) != 0) {
    return fail(error, "format: must be \"%s\"", MATRIX_FORMAT);
  }

  int version;
  if (!get_int(root, "format_version", &version) ||
      version != MATRIX_FORMAT_VERSION) {
    return fail(error, "format_version: must be %d", MATRIX_FORMAT_VERSION);
  }

  const cJSON *source = cJSON_GetObjectItemCaseSensitive(root, "source");
  if (!cJSON_IsObject(source)) {
    return fail(error, "source: missing");
  }
  if (!get_string(source, "tool", meta->tool, sizeof meta->tool)) {
    return fail(error, "source.tool: missing");
  }
  if (!get_string(source, "version", meta->version, sizeof meta->version)) {
    return fail(error, "source.version: missing");
  }
  if (!get_string(source, "date", meta->date, sizeof meta->date)) {
    return fail(error, "source.date: missing");
  }
  if (!is_iso_date(meta->date)) {
    return fail(error, "source.date: must be an ISO 8601 date (YYYY-MM-DD)");
  }
  if (!get_string(source, "input", meta->input, sizeof meta->input)) {
    return fail(error, "source.input: missing");
  }

  if (!cJSON_IsObject(cJSON_GetObjectItemCaseSensitive(root, "estimation"))) {
    return fail(error, "estimation: missing");
  }

  const char *slot_type =
      cJSON_GetStringValue(cJSON_GetObjectItemCaseSensitive(root, "slot_type"));
  if (slot_type != NULL && strcmp(slot_type, "position") == 0) {
    meta->type = POSITION;
  } else if (slot_type != NULL && strcmp(slot_type, "quality") == 0) {
    meta->type = QUAL;
  } else {
    return fail(error, "slot_type: must be one of position, quality");
  }

  if (!get_int(root, "range", &meta->range) || meta->range < 1) {
    return fail(error, "range: must be an integer greater than 0");
  }
  if (meta->type == QUAL && meta->range != 1) {
    return fail(error, "range: must be 1 for slot_type quality");
  }
  if (!get_int(root, "max_length", &meta->max_length) ||
      meta->max_length < 1) {
    return fail(error, "max_length: must be an integer greater than 0");
  }

  return 0;
}

// copy a dim x dim JSON matrix to out (row-major), or check that it has only
// integers if out is NULL
static int read_matrix(const cJSON *json, float *out, int dim,
                       const char *field, char *error) {
  if (!cJSON_IsArray(json) || cJSON_GetArraySize(json) != dim) {
    return fail(error, "%s: must be %d x %d", field, dim, dim);
  }
  for (int i = 0; i < dim; i++) {
    const cJSON *row = cJSON_GetArrayItem(json, i);
    if (!cJSON_IsArray(row) || cJSON_GetArraySize(row) != dim) {
      return fail(error, "%s: must be %d x %d", field, dim, dim);
    }
    for (int j = 0; j < dim; j++) {
      const cJSON *value = cJSON_GetArrayItem(row, j);
      if (!cJSON_IsNumber(value)) {
        return fail(error, "%s: must have only numbers", field);
      }
      if (out == NULL && value->valuedouble != floor(value->valuedouble)) {
        return fail(error, "%s: must have only integers", field);
      }
      if (out != NULL) {
        out[i * dim + j] = (float)value->valuedouble;
      }
    }
  }
  return 0;
}

// the channel rows must sum to 1 and channel x inverse must be the identity
static int validate_matrices(const float *channel, const float *inverse,
                             int dim, const char *field, char *error) {
  for (int i = 0; i < dim; i++) {
    float sum = 0;
    for (int j = 0; j < dim; j++) {
      sum += channel[i * dim + j];

      float product = 0;
      for (int n = 0; n < dim; n++) {
        product += channel[i * dim + n] * inverse[n * dim + j];
      }
      if (fabsf(product - (i == j)) > MATRIX_TOLERANCE) {
        return fail(error, "%s: channel x inverse is not the identity matrix",
                    field);
      }
    }
    if (fabsf(sum - 1) > MATRIX_TOLERANCE) {
      return fail(error, "%s: channel row %d does not sum to 1", field, i);
    }
  }
  return 0;
}

// cache the Schur product of each channel column and each loss column:
// schur[observed][x][y] = channel[y][observed] * loss[y][x]
static void set_schur_product(const float *channel, const float *loss,
                              float *schur, int dim) {
  for (int observed = 0; observed < dim; observed++) {
    for (int x = 0; x < dim; x++) {
      for (int y = 0; y < dim; y++) {
        schur[(observed * dim + x) * dim + y] =
            channel[y * dim + observed] * loss[y * dim + x];
      }
    }
  }
}

// check one section (substitutions or homopolymers) and copy its slots
// symbols: the symbol names, NULL for the numbers 1 to dim
// channel, inverse and schur: the arrays of dude, slots x dim x dim (x dim)
static int parse_section(const cJSON *root, const char *kind, int dim,
                         const char **symbols, dude_t *dude, float *channel,
                         float *inverse, float *schur, const float *loss,
                         char *error) {
  const cJSON *section = cJSON_GetObjectItemCaseSensitive(root, kind);
  if (!cJSON_IsObject(section)) {
    return fail(error, "%s: missing", kind);
  }

  const cJSON *names = cJSON_GetObjectItemCaseSensitive(section, "symbols");
  int symbols_ok = cJSON_IsArray(names) && cJSON_GetArraySize(names) == dim;
  for (int i = 0; symbols_ok && i < dim; i++) {
    const cJSON *name = cJSON_GetArrayItem(names, i);
    symbols_ok = symbols ? cJSON_IsString(name) &&
                               strcmp(name->valuestring, symbols[i]) == 0
                         : cJSON_IsNumber(name) && name->valuedouble == i + 1;
  }
  if (!symbols_ok) {
    return fail(error, "%s.symbols: must be %s", kind,
                symbols ? "[\"A\", \"C\", \"G\", \"T\"]"
                        : "[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]");
  }

  const cJSON *slots = cJSON_GetObjectItemCaseSensitive(section, "slots");
  if (!cJSON_IsArray(slots) || cJSON_GetArraySize(slots) != dude->slots) {
    return fail(error, "%s.slots: must have %d slots", kind, dude->slots);
  }

  int size = dim * dim;
  for (int z = 0; z < dude->slots; z++) {
    const cJSON *slot = cJSON_GetArrayItem(slots, z);
    char field[64];
    int min, max;

    snprintf(field, sizeof field, "%s.slots[%d]", kind, z);
    if (!cJSON_IsObject(slot)) {
      return fail(error, "%s: must be an object", field);
    }
    if (!get_int(slot, "min", &min) || !get_int(slot, "max", &max) ||
        min != z * dude->range || max != (z + 1) * dude->range) {
      return fail(error,
                  "%s: min and max must be %d and %d, no gaps and no overlaps",
                  field, z * dude->range, (z + 1) * dude->range);
    }

    char matrix_field[80];
    snprintf(matrix_field, sizeof matrix_field, "%s.counts", field);
    if (read_matrix(cJSON_GetObjectItemCaseSensitive(slot, "counts"), NULL,
                    dim, matrix_field, error)) {
      return 1;
    }
    snprintf(matrix_field, sizeof matrix_field, "%s.channel", field);
    if (read_matrix(cJSON_GetObjectItemCaseSensitive(slot, "channel"),
                    channel + z * size, dim, matrix_field, error)) {
      return 1;
    }
    snprintf(matrix_field, sizeof matrix_field, "%s.inverse", field);
    if (read_matrix(cJSON_GetObjectItemCaseSensitive(slot, "inverse"),
                    inverse + z * size, dim, matrix_field, error)) {
      return 1;
    }
    if (validate_matrices(channel + z * size, inverse + z * size, dim, field,
                          error)) {
      return 1;
    }

    set_schur_product(channel + z * size, loss, schur + z * size * dim, dim);
  }
  return 0;
}

int parse_matrices(const cJSON *root, dude_t *dude, char *error) {
  return parse_section(root, "substitutions", SUBS_MATRIX_DIM, SUBS_SYMBOLS,
                       dude, &dude->subs_channel[0][0][0],
                       &dude->subs_channel_i[0][0][0],
                       &dude->subs_schur_product[0][0][0][0],
                       &dude->subs_loss[0][0], error) ||
         parse_section(root, "homopolymers", HOM_MATRIX_DIM, NULL, dude,
                       &dude->hom_channel[0][0][0],
                       &dude->hom_channel_i[0][0][0],
                       &dude->hom_schur_product[0][0][0][0],
                       &dude->hom_loss[0][0], error);
}

matrix_meta_t read_matrix_meta(const char *path) {
  matrix_meta_t meta;
  char error[MATRIX_ERROR_SIZE];

  cJSON *root = parse_matrix_file(path);
  int failed = parse_matrix_meta(root, &meta, error);
  cJSON_Delete(root);
  if (failed) {
    matrix_error(path, error);
  }
  return meta;
}

void load_matrices(const char *path, dude_t *dude) {
  char error[MATRIX_ERROR_SIZE];

  printf("Loading and validating the matrices...\n");

  cJSON *root = parse_matrix_file(path);
  int failed = parse_matrices(root, dude, error);
  cJSON_Delete(root);
  if (failed) {
    matrix_error(path, error);
  }

  printf("Validation completed\n");
}

// create a new context helper
context_helper_t *create_helper(void) {
  context_helper_t *helper =
      (context_helper_t *)malloc(sizeof(context_helper_t));
  if (helper == NULL) {
    error_exit();
  }
  return helper;
}

// Create a new DUDE params struct, with the slots of the matrices file
dude_t *create_dude_params(int k, int h, matrix_meta_t meta) {

  dude_t *dude = (dude_t *)malloc(sizeof(dude_t));

  if (dude == NULL) {
    error_exit();
  }

  // substitutions context length
  dude->k = k;

  // indels context length
  dude->h = h;

  dude->type = meta.type;

  // the last slot can go past max_length, ie. 150 and range 20 -> 8 slots
  dude->range = meta.range;
  dude->slots = (meta.max_length + meta.range - 1) / meta.range;

  // ======================================
  //   Substitutions Hamming loss function
  // ======================================
  for (int i = 0; i < SUBS_MATRIX_DIM; i++) {
    for (int j = 0; j < SUBS_MATRIX_DIM; j++) {
      dude->subs_loss[i][j] = i != j;
    }
  }

  // ======================================
  //   Homopolymers Hamming loss function
  // ======================================
  for (int i = 0; i < HOM_MATRIX_DIM; i++) {
    for (int j = 0; j < HOM_MATRIX_DIM; j++) {
      dude->hom_loss[i][j] = i != j;
    }
  }

  int size = dude->slots;

  dude->subs_channel = calloc(size, sizeof(subs_matrix_type));
  dude->subs_channel_i = calloc(size, sizeof(subs_matrix_type));
  dude->subs_schur_product = calloc(size, sizeof(subs_schur_product));

  dude->hom_channel = calloc(size, sizeof(hom_matrix_type));
  dude->hom_channel_i = calloc(size, sizeof(hom_matrix_type));
  dude->hom_schur_product = calloc(size, sizeof(hom_schur_product));

  if (!dude->subs_channel || !dude->subs_channel_i ||
      !dude->subs_schur_product || !dude->hom_channel || !dude->hom_channel_i ||
      !dude->hom_schur_product) {
    error_exit();
  }

  return dude;
}

void destroy_dude_params(dude_t *dude) {
  if (dude == NULL)
    return;
  free(dude->subs_channel);
  free(dude->subs_channel_i);
  free(dude->subs_schur_product);
  free(dude->hom_channel);
  free(dude->hom_channel_i);
  free(dude->hom_schur_product);
  free(dude);
}

int get_index_for_matrix(int position_in_read, dude_t *dude) {
  int p = position_in_read / dude->range;
  if (p >= dude->slots) {
    // positions past the max read length use the last range
    return dude->slots - 1;
  } else {
    return p;
  }
}

void destroy_helper(context_helper_t *helper) {
  if (helper == NULL)
    return;
  // free(helper->counter);
  // free(helper->counters_inv);
  free(helper);
}

void error_exit() {
  puts("Error on malloc");
  exit(-1);
}
