#include <gmp.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <stdatomic.h>

#define NUM_THREADS 1

typedef struct {
  mpz_t x;
  mpz_t y;
} ec_point_t;

typedef struct {
  int thread_id;
  ec_point_t public_key;
  mpz_t p;
  mpz_t start_k;
  mpz_t end_k;
} thread_arg_t;

atomic_uint_least64_t current_step = 0;
atomic_uint_least64_t private_key = 0;
atomic_int found_collision = 0;

pthread_mutex_t collision_mutex = PTHREAD_MUTEX_INITIALIZER;

const char *P_HEX = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F";

void ec_point_init(ec_point_t *point) { mpz_inits(point->x, point->y, NULL); }

void ec_point_clear(ec_point_t *point) { mpz_clears(point->x, point->y, NULL); }

void ec_point_add(ec_point_t *result, ec_point_t *p1, ec_point_t *p2, mpz_t p) {
    mpz_t lambda, temp1, temp2;
    mpz_inits(lambda, temp1, temp2, NULL);

    if (mpz_cmp(p1->x, p2->x) == 0 && mpz_cmp(p1->y, p2->y) == 0) {
        mpz_t three, two;
        mpz_inits(three, two, NULL);
        mpz_set_ui(three, 3);
        mpz_set_ui(two, 2);

        mpz_mul(temp1, p1->x, p1->x);
        mpz_mul(temp1, temp1, three);
        mpz_mod(temp1, temp1, p);

        mpz_mul(temp2, p1->y, two);
        mpz_invert(temp2, temp2, p);
        mpz_mul(lambda, temp1, temp2);
        mpz_mod(lambda, lambda, p);

        mpz_clears(three, two, NULL);
    } else if (mpz_cmp(p1->x, p2->x) == 0) {
        mpz_clears(lambda, temp1, temp2, NULL);
        return;
    } else {
        mpz_sub(temp1, p2->y, p1->y);
        mpz_sub(temp2, p2->x, p1->x);
        mpz_invert(temp2, temp2, p);
        mpz_mul(lambda, temp1, temp2);
        mpz_mod(lambda, lambda, p);
    }

    mpz_mul(temp1, lambda, lambda);
    mpz_sub(temp1, temp1, p1->x);
    mpz_sub(temp1, temp1, p2->x);
    mpz_mod(result->x, temp1, p);
    mpz_sub(temp1, p1->x, result->x);
    mpz_mul(temp1, lambda, temp1);
    mpz_sub(temp1, temp1, p1->y);
    mpz_mod(result->y, temp1, p);

    mpz_clears(lambda, temp1, temp2, NULL);
}

int ec_point_equal(ec_point_t *p1, ec_point_t *p2) {
  return (mpz_cmp(p1->x, p2->x) == 0 && mpz_cmp(p1->y, p2->y) == 0);
}

int y_coordinate_from_compressed(mpz_t y, mpz_t x, const char *prefix, mpz_t p) {
  mpz_t y_squared, beta, temp;
  mpz_inits(y_squared, beta, temp, NULL);

  mpz_powm_ui(y_squared, x, 3, p);
  mpz_add_ui(y_squared, y_squared, 7);
  mpz_mod(y_squared, y_squared, p);

  mpz_add_ui(temp, p, 1);
  mpz_fdiv_q_ui(temp, temp, 4);
  mpz_powm(beta, y_squared, temp, p);

  if ((mpz_tstbit(beta, 0) == 0 && prefix[1] == '2') || (mpz_tstbit(beta, 0) == 1 && prefix[1] == '3')) {
    mpz_set(y, beta);
  } else {
    mpz_sub(y, p, beta);
  }

  mpz_clears(y_squared, beta, temp, NULL);
  return 1;
}

void show_progress(uint64_t current_step, mpz_t max_k, struct timespec *last_update_time) {
  struct timespec current_time;
  clock_gettime(CLOCK_MONOTONIC, &current_time);

  static uint64_t last_step_count = 0;
  double elapsed_sec = (double)(current_time.tv_sec - last_update_time->tv_sec) + (double)(current_time.tv_nsec - last_update_time->tv_nsec) / 1e9;

  if (elapsed_sec >= 1.0) {
    double current_steps_per_sec = (double)(current_step - last_step_count) / elapsed_sec;
    const char *steps_per_second_unit;

    if (current_steps_per_sec >= 1e18) {
      current_steps_per_sec /= 1e18;
      steps_per_second_unit = "Es/s";
    } else if (current_steps_per_sec >= 1e15) {
      current_steps_per_sec /= 1e15;
      steps_per_second_unit = "Ps/s";
    } else if (current_steps_per_sec >= 1e12) {
      current_steps_per_sec /= 1e12;
      steps_per_second_unit = "Ts/s";
    } else if (current_steps_per_sec >= 1e9) {
      current_steps_per_sec /= 1e9;
      steps_per_second_unit = "Gs/s";
    } else if (current_steps_per_sec >= 1e6) {
      current_steps_per_sec /= 1e6;
      steps_per_second_unit = "Ms/s";
    } else if (current_steps_per_sec >= 1e3) {
      current_steps_per_sec /= 1e3;
      steps_per_second_unit = "Ks/s";
    } else {
      steps_per_second_unit = "S/s";
    }

    double key_range_double = mpz_get_d(max_k);
    int bits = (int)round(log2(key_range_double));

    printf("\rSteps: %llu %.2f %s Key Range: %d bits", (unsigned long long)current_step, current_steps_per_sec, steps_per_second_unit, bits);
    fflush(stdout);

    *last_update_time = current_time;
    last_step_count = current_step;
  }
}

void ec_point_set(ec_point_t *dest, ec_point_t *src) {
  mpz_set(dest->x, src->x);
  mpz_set(dest->y, src->y);
}

void *thread_function(void *arg) {
  thread_arg_t *thread_args = (thread_arg_t *)arg;
  ec_point_t *public_key = &thread_args->public_key;
  mpz_t p;
  mpz_init_set(p, thread_args->p);
  ec_point_t temp, result;
  ec_point_init(&temp);
  ec_point_init(&result);

  mpz_set(temp.x, public_key->x);
  mpz_set(temp.y, public_key->y);

  mpz_t current_k, end_k;
  mpz_init_set(current_k, thread_args->start_k);
  mpz_init_set(end_k, thread_args->end_k);

  struct timespec last_update_time;
  clock_gettime(CLOCK_MONOTONIC, &last_update_time);

  while (mpz_cmp(current_k, end_k) < 0 && !atomic_load(&found_collision)) {
      ec_point_add(&result, &temp, public_key, p);
      mpz_set(temp.x, result.x);
      mpz_set(temp.y, result.y);

      mpz_add_ui(current_k, current_k, 1);

      atomic_fetch_add(&current_step, 1);

      show_progress(current_step, thread_args->end_k, &last_update_time);

      if (ec_point_equal(&result, public_key)) {
          pthread_mutex_lock(&collision_mutex);
          if (!atomic_load(&found_collision)) {
              atomic_store(&private_key, mpz_get_ui(current_k));
              atomic_store(&found_collision, 1);
              printf("\rCollision found! Private key: %lu\n", atomic_load(&private_key));
              fflush(stdout);
          }
          pthread_mutex_unlock(&collision_mutex);
          break;
      }
  }

  mpz_clear(current_k);
  mpz_clear(end_k);
  ec_point_clear(&temp);
  ec_point_clear(&result);
  mpz_clear(p);

  return NULL;
}

int pollardsrho(int argc, char *argv[]) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s <public key> <key range>\n", argv[0]);
    return 1;
  }

  mpz_t p;
  mpz_init_set_str(p, P_HEX, 16);

  ec_point_t public_key_a;
  ec_point_init(&public_key_a);

  char *public_key_hex = argv[1];
  if (strlen(public_key_hex) != 66) {
    fprintf(stderr, "Invalid public key length. Expected 66 characters (33 bytes in hex).\n");
    return 1;
  }

  char prefix[3];
  char x_hex[65];
  strncpy(prefix, public_key_hex, 2);
  prefix[2] = '\0';
  strncpy(x_hex, public_key_hex + 2, 64);
  x_hex[64] = '\0';

  mpz_set_str(public_key_a.x, x_hex, 16);
  if (!y_coordinate_from_compressed(public_key_a.y, public_key_a.x, prefix, p)) {
    fprintf(stderr, "Failed to calculate y-coordinate from compressed public key.\n");
    return 1;
  }

  mpz_t key_range;
  mpz_init_set_str(key_range, argv[2], 10);
  mpz_t max_k;
  mpz_init(max_k);
  mpz_ui_pow_ui(max_k, 2, mpz_get_ui(key_range));

  pthread_t threads[NUM_THREADS];
  thread_arg_t thread_args[NUM_THREADS];

  mpz_t range_size, thread_range, start_k;
  mpz_inits(range_size, thread_range, start_k, NULL);
  mpz_set(range_size, max_k);
  mpz_div_ui(thread_range, range_size, NUM_THREADS);

  for (int i = 0; i < NUM_THREADS; i++) {
    thread_args[i].thread_id = i;
    thread_args[i].public_key = public_key_a;
    mpz_init_set(thread_args[i].p, p);

    mpz_init_set(thread_args[i].start_k, start_k);
    if (i == NUM_THREADS - 1) {
      mpz_init_set(thread_args[i].end_k, max_k);
    } else {
      mpz_init(thread_args[i].end_k);
      mpz_add(thread_args[i].end_k, start_k, thread_range);
    }

    if (pthread_create(&threads[i], NULL, thread_function, (void *)&thread_args[i])) {
      perror("pthread_create");
      exit(1);
    }

    mpz_add(start_k, start_k, thread_range);
  }

  for (int i = 0; i < NUM_THREADS; i++) {
    if (pthread_join(threads[i], NULL)) {
      perror("pthread_join");
      exit(1);
    }
    mpz_clear(thread_args[i].p);
    mpz_clear(thread_args[i].start_k);
    mpz_clear(thread_args[i].end_k);
  }

  ec_point_clear(&public_key_a);
  mpz_clears(p, key_range, max_k, range_size, thread_range, start_k, NULL);

  pthread_mutex_lock(&collision_mutex);
  if (!atomic_load(&found_collision)) {
    printf("\rNo collision found within the given range.\n");
    fflush(stdout);
  }
  pthread_mutex_unlock(&collision_mutex);

  return 0;
}

int main(int argc, char *argv[]) {
  return pollardsrho(argc, argv);
}