#include <gmp.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <stdatomic.h>

#define NUM_THREADS 8

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

static bool loading_points = true;

const char *P_HEX = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F";

int key_range;

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

void show_progress(uint64_t current_step, struct timespec *last_update_time) {
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

    printf("\rSteps: %llu %.2f %s Key Range: %d bits", (unsigned long long)current_step, current_steps_per_sec, steps_per_second_unit, key_range);
    fflush(stdout);

    *last_update_time = current_time;
    last_step_count = current_step;
  }
}

void ec_point_set(ec_point_t *dest, ec_point_t *src) {
  mpz_set(dest->x, src->x);
  mpz_set(dest->y, src->y);
}

void points(ec_point_t *derived_points, ec_point_t *A, mpz_t Gx, mpz_t Gy, mpz_t p, int num_points) {
    mpz_t increment;
    mpz_init(increment);
    mpz_tdiv_q_ui(increment, p, num_points);

    ec_point_t temp, G;
    ec_point_init(&temp);
    ec_point_init(&G);

    mpz_set(G.x, Gx);
    mpz_set(G.y, Gy);

    ec_point_t x1, x2;
    ec_point_init(&x1);
    ec_point_init(&x2);

    mpz_set(x1.x, A->x);
    mpz_set(x1.y, A->y);
    mpz_set(x2.x, A->x);
    mpz_set(x2.y, A->y);

    mpz_t a1, a2, b1, b2;
    mpz_inits(a1, a2, b1, b2, NULL);
    mpz_set_ui(a1, 0);
    mpz_set_ui(a2, 1);
    mpz_set_ui(b1, 0);
    mpz_set_ui(b2, 1);

    ec_point_t tortoise, hare;
    ec_point_init(&tortoise);
    ec_point_init(&hare);
    ec_point_set(&tortoise, &x1);
    ec_point_set(&hare, &x1);

    int steps = 2;
    int step_size = 1;

    while (steps < num_points) {
        for (int j = 0; j < step_size && steps < num_points; j++) {
            ec_point_add(&hare, &hare, &G, p);
            //printf("Step %d: Hare - X: ", steps);
            //mpz_out_str(stdout, 16, hare.x);
            //printf(", Y: ");
            //mpz_out_str(stdout, 16, hare.y);
            //printf("\n");
            steps++;
        }

        for (int j = 0; j < step_size && steps < num_points; j++) {
            ec_point_add(&tortoise, &tortoise, &G, p);
            ec_point_add(&hare, &hare, &G, p);
            //printf("Step %d: Tortoise - X: ", steps);
            //mpz_out_str(stdout, 16, tortoise.x);
            //printf(", Y: ");
            //mpz_out_str(stdout, 16, tortoise.y);
            //printf("\n");
            steps++;

            if (ec_point_equal(&tortoise, &hare)) {
                printf("Collision detected at step %d\n", steps);
                fflush(stdout);
                ec_point_set(&derived_points[0], &tortoise);
                mpz_clears(increment, a1, a2, b1, b2, NULL);
                ec_point_clear(&temp);
                ec_point_clear(&G);
                ec_point_clear(&x1);
                ec_point_clear(&x2);
                ec_point_clear(&tortoise);
                ec_point_clear(&hare);
                return;
            }
        }

        step_size *= 2;
    }

    for (int i = 1; i < num_points; i++) {
        ec_point_add(&x1, &x1, &G, p);
        ec_point_add(&x2, &x2, &G, p);
        ec_point_add(&x2, &x2, &G, p);
        ec_point_set(&derived_points[i], &x1);

        //Invert the point to mirror the (x) symmetry >> (x, y) for (x, -y)
        mpz_set(derived_points[num_points + i].x, derived_points[i].x);
        mpz_sub(derived_points[num_points + i].y, p, derived_points[i].y);
    }

    mpz_clear(increment);
    ec_point_clear(&temp);
    ec_point_clear(&G);
    ec_point_clear(&x1);
    ec_point_clear(&x2);
    mpz_clears(a1, a2, b1, b2, NULL);
}

void *thread_function(void *arg) {
  thread_arg_t *thread_args = (thread_arg_t *)arg;
  ec_point_t *public_key = &thread_args->public_key;
  mpz_t p, Gx, Gy, n, n_half;
  mpz_inits(p, Gx, Gy, n, n_half, NULL);
  mpz_set(p, thread_args->p);

  //SECP25651
  mpz_set_str(Gx, "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
  mpz_set_str(Gy, "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16);
  mpz_set_str(n, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141", 16);
  mpz_tdiv_q_ui(n_half, n, 2);

  uint_least64_t num_derived_points;

  if (key_range <= 30) { num_derived_points = (uint_least64_t)pow(2, key_range / 2.0); } 
  else { num_derived_points = 1ULL << key_range; }

  if (num_derived_points > 32000) {
      num_derived_points = 32000;
  }

  //Derive the two symmetric parts of secp256k1:
  ec_point_t derived_points[2 * num_derived_points];
  for (int i = 0; i < 2 * num_derived_points; i++) {
    ec_point_init(&derived_points[i]);
  }

  points(derived_points, public_key, Gx, Gy, n_half, num_derived_points);

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

      if (mpz_cmp_ui(temp.y, 0) < 0) {
          mpz_sub(temp.y, p, temp.y);
      }

      mpz_add_ui(current_k, current_k, 1);
      atomic_fetch_add(&current_step, 1);

      show_progress(current_step, &last_update_time);

      if (ec_point_equal(&result, public_key)) {
          pthread_mutex_lock(&collision_mutex);
          if (!atomic_load(&found_collision)) {
              atomic_store(&private_key, mpz_get_ui(current_k));
              atomic_store(&found_collision, 1);
              mpz_t private_key_mpz;
              mpz_init_set_ui(private_key_mpz, atomic_load(&private_key));

              printf("\rCollision found! Private key: ");
              gmp_printf("%ZX\n", private_key_mpz);
              printf("derived points: %d\n", num_derived_points);
              fflush(stdout);

              FILE *file = fopen("KeysFound.txt", "a");
              if (file == NULL) {
                  perror("Error opening KeysFound.txt, is null!");
              } else {
                  gmp_fprintf(file, "%ZX\n", private_key_mpz);
                  fclose(file);
              }
              mpz_clear(private_key_mpz);
          }
          pthread_mutex_unlock(&collision_mutex);
          break;
      }
      
      for (int i = 0; i < 2 * num_derived_points; i++) {
           if (ec_point_equal(&result, &derived_points[i])) {
                pthread_mutex_lock(&collision_mutex);
                if (!atomic_load(&found_collision)) {
                     atomic_store(&private_key, mpz_get_ui(current_k));
                     atomic_store(&found_collision, 1);
                     mpz_t private_key_mpz;
                     mpz_init_set_ui(private_key_mpz, atomic_load(&private_key));

                     printf("\rCollision found with derived point! Private key: ");
                     gmp_printf("%ZX\n", private_key_mpz);
                     fflush(stdout);

                     FILE *file = fopen("KeysFound.txt", "a");
                     if (file == NULL) { perror("Error opening KeysFound.txt, is null!");}
                     else {
                     gmp_fprintf(file, "%ZX\n", private_key_mpz);
                     fclose(file); }
                     mpz_clear(private_key_mpz);
                }
                pthread_mutex_unlock(&collision_mutex);
                break;
           }
      }
  }

  mpz_clear(current_k);
  mpz_clear(end_k);
  mpz_clear(n_half);
  ec_point_clear(&temp);
  ec_point_clear(&result);
  for (int i = 0; i < 2 * num_derived_points; i++) {
    ec_point_clear(&derived_points[i]);
  }
  mpz_clears(p, Gx, Gy, n, NULL);

  return NULL;
}

int pollardsrho(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <public key> <key range>\n", argv[0]);
        return 1;
    }

    if(loading_points) { 
       printf("\rLoading Points...");
       fflush(stdout);

       loading_points = false; 
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

    key_range = atoi(argv[2]);

    if (key_range <= 0 || key_range > 256) {
        fprintf(stderr, "Invalid key_range. Must be between 1 and 256.\n");
        return 1;
    }

    mpz_t max_k;
    mpz_init(max_k);
    mpz_ui_pow_ui(max_k, 2, key_range);

    pthread_t threads[NUM_THREADS];
    thread_arg_t thread_args[NUM_THREADS];

    mpz_t thread_range, start_k;
    mpz_inits(thread_range, start_k, NULL);
    mpz_tdiv_q_ui(thread_range, max_k, NUM_THREADS);
    mpz_set_ui(start_k, 0);

    for (int i = 0; i < NUM_THREADS; i++) {
        thread_args[i].thread_id = i;
        mpz_init_set(thread_args[i].p, p);
        thread_args[i].public_key = public_key_a;
        mpz_init_set_ui(thread_args[i].start_k, mpz_get_ui(start_k));
        mpz_init(thread_args[i].end_k);
        mpz_add(thread_args[i].end_k, start_k, thread_range);
        if (i == NUM_THREADS - 1) { mpz_set(thread_args[i].end_k, max_k); }   
        if (pthread_create(&threads[i], NULL, thread_function, (void *)&thread_args[i])) { perror("pthread_create");exit(1); }
        mpz_add_ui(start_k, start_k, mpz_get_ui(thread_range));
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
    mpz_clears(p, max_k, thread_range, start_k, NULL);

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
