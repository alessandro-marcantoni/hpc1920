#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define FULL_ANGLE 360.0

typedef struct {
    double x, y;
} point_t;

typedef struct {
    point_t *p;
    int n;
} points_t;

void print_circle(FILE *f, const points_t *circle) {
    fprintf(f, "%d\n%d\n", 2, circle->n);
    for (int i=0; i<circle->n; i++) {
        fprintf(f, "%f %f\n", circle->p[i].x, circle->p[i].y);
    }
}

int main(int argc, char *argv[]) {

    int size = atoi(argv[1]);

    points_t circle;
    circle.n = 0;
    circle.p = (point_t *)malloc(size * sizeof(point_t));

    double angle = FULL_ANGLE / size;
    for (int i=0; i<size; i++) {
        circle.p[i].x = cos(i * angle);
        circle.p[i].y = sin(i * angle);
        circle.n++;
    }

    print_circle(stdout, &circle);

    return EXIT_SUCCESS;
}
