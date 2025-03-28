#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>

#define PI 3.141592
#define RRR 6378.388
#define time_limit 1.35

double to_geographical_radians(double degree) {
    double deg = floor(degree + 0.5);
    double min = degree - deg;
    return PI * (deg + 5.0 * min / 3.0) / 180.0;
}

double calculate_geographical_distance(double x1, double y1, double x2, double y2) {
    double lat1 = to_geographical_radians(x1);
    double long1 = to_geographical_radians(y1);
    double lat2 = to_geographical_radians(x2);
    double long2 = to_geographical_radians(y2);

    double q1 = cos(long1 - long2);
    double q2 = cos(lat1 - lat2);
    double q3 = cos(lat1 + lat2);

    double tmp = RRR * acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0;
    return floor(tmp);
}

double calculate_euclidean_distance(double x1, double y1, double x2, double y2) {
    return (int)(sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2)));
}

void read_coordinates(const char* filename, double** coordinates, int* num_points) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Failed to open file.\n");
        exit(1);
    }

    fscanf(file, "%d", num_points);
    *coordinates = (double*)malloc(*num_points * 2 * sizeof(double));

    int vertex;
    for (int i = 0; i < *num_points; i++) {
        fscanf(file, "%d %lf %lf", &vertex, &(*coordinates)[2 * i], &(*coordinates)[2 * i + 1]);
    }

    fclose(file);
}

void nearest_insertion_tsp(double* coordinates, int num_points, int* tour, double* total_euclidean_distance, double* total_geographical_distance) {
    int* visited = (int*)calloc(num_points, sizeof(int));
    int num_tour = 0;

    int first = 0;
    int second = 1;
    double min_dist = calculate_euclidean_distance(
        coordinates[2 * first], coordinates[2 * first + 1],
        coordinates[2 * second], coordinates[2 * second + 1]
    );

    for (int i = 2; i < num_points; i++) {
        double dist = calculate_euclidean_distance(
            coordinates[2 * first], coordinates[2 * first + 1],
            coordinates[2 * i], coordinates[2 * i + 1]
        ) + calculate_euclidean_distance(
            coordinates[2 * second], coordinates[2 * second + 1],
            coordinates[2 * i], coordinates[2 * i + 1]
        ) - calculate_euclidean_distance(
            coordinates[2 * first], coordinates[2 * first + 1],
            coordinates[2 * second], coordinates[2 * second + 1]
        );

        if (dist < min_dist) {
            min_dist = dist;
            second = i;
        }
    }

    tour[num_tour++] = first;
    tour[num_tour++] = second;
    visited[first] = 1;
    visited[second] = 1;
    *total_euclidean_distance = calculate_euclidean_distance(
        coordinates[2 * first], coordinates[2 * first + 1],
        coordinates[2 * second], coordinates[2 * second + 1]
    );
    *total_geographical_distance = calculate_geographical_distance(
        coordinates[2 * first], coordinates[2 * first + 1],
        coordinates[2 * second], coordinates[2 * second + 1]
    );

    while (num_tour < num_points) {
        double min_insertion_dist = INT_MAX;
        int insert_index = -1;
        int new_city = -1;

        for (int i = 0; i < num_points; i++) {
            if (!visited[i]) {
                for (int j = 0; j < num_tour; j++) {
                    int next_index = (j + 1) % num_tour;
                    double dist = calculate_euclidean_distance(
                        coordinates[2 * tour[j]], coordinates[2 * tour[j] + 1],
                        coordinates[2 * i], coordinates[2 * i + 1]
                    ) + calculate_euclidean_distance(
                        coordinates[2 * i], coordinates[2 * i + 1],
                        coordinates[2 * tour[next_index]], coordinates[2 * tour[next_index] + 1]
                    ) - calculate_euclidean_distance(
                        coordinates[2 * tour[j]], coordinates[2 * tour[j] + 1],
                        coordinates[2 * tour[next_index]], coordinates[2 * tour[next_index] + 1]
                    );

                    if (dist < min_insertion_dist) {
                        min_insertion_dist = dist;
                        insert_index = j;
                        new_city = i;
                    }
                }
            }
        }

        for (int i = num_tour; i > insert_index + 1; i--) {
            tour[i] = tour[i - 1];
        }
        tour[insert_index + 1] = new_city;
        visited[new_city] = 1;
        num_tour++;
        *total_euclidean_distance += min_insertion_dist;
        *total_geographical_distance += calculate_geographical_distance(
            coordinates[2 * tour[insert_index]], coordinates[2 * tour[insert_index] + 1],
            coordinates[2 * new_city], coordinates[2 * new_city + 1]
        ) + calculate_geographical_distance(
            coordinates[2 * new_city], coordinates[2 * new_city + 1],
            coordinates[2 * tour[(insert_index + 2) % num_tour]], coordinates[2 * tour[(insert_index + 2) % num_tour] + 1]
        ) - calculate_geographical_distance(
            coordinates[2 * tour[insert_index]], coordinates[2 * tour[insert_index] + 1],
            coordinates[2 * tour[(insert_index + 2) % num_tour]], coordinates[2 * tour[(insert_index + 2) % num_tour] + 1]
        );
    }

    *total_euclidean_distance += calculate_euclidean_distance(
        coordinates[2 * tour[num_tour - 1]], coordinates[2 * tour[num_tour - 1] + 1],
        coordinates[2 * tour[0]], coordinates[2 * tour[0] + 1]
    );
    *total_geographical_distance += calculate_geographical_distance(
        coordinates[2 * tour[num_tour - 1]], coordinates[2 * tour[num_tour - 1] + 1],
        coordinates[2 * tour[0]], coordinates[2 * tour[0] + 1]
    );

    free(visited);
}

void two_opt_swap(int* tour, int i, int k, int num_points) {
    while (i < k) {
        int temp = tour[i];
        tour[i] = tour[k];
        tour[k] = temp;
        i++;
        k--;
    }
}

void local_search_2opt(double* coordinates, int num_points, int* tour, double* total_euclidean_distance, double* total_geographical_distance) {
    clock_t start_time = clock();
    int improvement = 1;

    while (improvement) {
        improvement = 0;
        for (int i = 1; i < num_points - 2; i++) {
            for (int k = i + 1; k < num_points - 1; k++) {
                clock_t end_time = clock();
                double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
                if (elapsed_time > time_limit) {
                    return;
                }

                two_opt_swap(tour, i, k, num_points - 1);

                double new_euclidean_distance = 0.0;
                double new_geographical_distance = 0.0;
                for (int j = 0; j < num_points - 1; j++) {
                    new_euclidean_distance += calculate_euclidean_distance(
                        coordinates[2 * tour[j]], coordinates[2 * tour[j] + 1],
                        coordinates[2 * tour[j + 1]], coordinates[2 * tour[j + 1] + 1]
                    );
                    new_geographical_distance += calculate_geographical_distance(
                        coordinates[2 * tour[j]], coordinates[2 * tour[j] + 1],
                        coordinates[2 * tour[j + 1]], coordinates[2 * tour[j + 1] + 1]
                    );
                }
                new_euclidean_distance += calculate_euclidean_distance(
                    coordinates[2 * tour[num_points - 1]], coordinates[2 * tour[num_points - 1] + 1],
                    coordinates[2 * tour[0]], coordinates[2 * tour[0] + 1]
                );
                new_geographical_distance += calculate_geographical_distance(
                    coordinates[2 * tour[num_points - 1]], coordinates[2 * tour[num_points - 1] + 1],
                    coordinates[2 * tour[0]], coordinates[2 * tour[0] + 1]
                );

                if (new_euclidean_distance < *total_euclidean_distance && new_geographical_distance < *total_geographical_distance) {
                    *total_euclidean_distance = new_euclidean_distance;
                    *total_geographical_distance = new_geographical_distance;
                    improvement = 1;
                } else {
                    two_opt_swap(tour, i, k, num_points - 1);
                }
            }
        }
    }
}

int main(void) {
    double* coordinates;
    int num_points;
    char filename[100];

    printf("Enter the filename: ");
    scanf("%s", filename);

    read_coordinates(filename, &coordinates, &num_points);

    int* tour = (int*)malloc(num_points * sizeof(int));
    double total_euclidean_distance, total_geographical_distance;

    clock_t start_t, end_t;
    double utime;

    start_t = clock();

    nearest_insertion_tsp(coordinates, num_points, tour, &total_euclidean_distance, &total_geographical_distance);
    local_search_2opt(coordinates, num_points, tour, &total_euclidean_distance, &total_geographical_distance);

    end_t = clock();
    utime = (double)(end_t - start_t) / CLOCKS_PER_SEC;

    printf("Time = %lf seconds\n", utime);
    printf("Total Euclidean distance: %lf\n", total_euclidean_distance);
    printf("Total Geographical distance: %lf\n", total_geographical_distance);
    printf("Number of vertices: %d\n", num_points);
    printf("Tour: ");
    for (int i = 0; i < num_points; i++) {
        printf("%d ", tour[i] + 1);
    }
    printf("\n");

    free(coordinates);
    free(tour);

    return 0;
}
