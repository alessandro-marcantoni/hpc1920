/****************************************************************************
 *
 * convex-hull.c
 *
 * Compute the convex hull of a set of points in 2D
 *
 * Copyright (C) 2019 Moreno Marzolla <moreno.marzolla(at)unibo.it>
 * Last updated on 2019-11-25
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ****************************************************************************
 *
 * Questo programma calcola l'inviluppo convesso (convex hull) di un
 * insieme di punti 2D letti da standard input usando l'algoritmo
 * "gift wrapping". Le coordinate dei vertici dell'inviluppo sono
 * stampate su standard output.  Per una descrizione completa del
 * problema si veda la specifica del progetto sul sito del corso:
 *
 * http://moreno.marzolla.name/teaching/HPC/
 *
 * Per compilare:
 *
 * gcc -D_XOPEN_SOURCE=600 -std=c99 -Wall -Wpedantic -O2 convex-hull.c -o convex-hull -lm
 *
 * (il flag -D_XOPEN_SOURCE=600 e' superfluo perche' viene settato
 * nell'header "hpc.h", ma definirlo tramite la riga di comando fa si'
 * che il programma compili correttamente anche se non si include
 * "hpc.h", o per errore non lo si include come primo file).
 *
 * Per eseguire il programma si puo' usare la riga di comando:
 *
 * ./convex-hull < ace.in > ace.hull
 * 
 * Per visualizzare graficamente i punti e l'inviluppo calcolato è
 * possibile usare lo script di gnuplot (http://www.gnuplot.info/)
 * incluso nella specifica del progetto:
 *
 * gnuplot -c plot-hull.gp ace.in ace.hull ace.png
 *
 ****************************************************************************/
#include "hpc.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* A single point */
typedef struct {
    double x, y;
} point_t;

/* An array of n points */
typedef struct {
    int n;      /* number of points     */
    point_t *p; /* array of points      */
} points_t;

enum {
    LEFT = -1,
    COLLINEAR,
    RIGHT
};

int n_threads;

/**
 * Read input from file f, and store the set of points into the
 * structure pset.
 */
void read_input( FILE *f, points_t *pset )
{
    char buf[1024];
    int i, dim, npoints;
    point_t *points;
    
    if ( 1 != fscanf(f, "%d", &dim) ) {
        fprintf(stderr, "FATAL: can not read dimension\n");
        exit(EXIT_FAILURE);
    }
    if (dim != 2) {
        fprintf(stderr, "FATAL: This program supports dimension 2 only (got dimension %d instead)\n", dim);
        exit(EXIT_FAILURE);
    }
    if (NULL == fgets(buf, sizeof(buf), f)) { /* ignore rest of the line */
        fprintf(stderr, "FATAL: failed to read rest of first line\n");
        exit(EXIT_FAILURE);
    }
    if (1 != fscanf(f, "%d", &npoints)) {
        fprintf(stderr, "FATAL: can not read number of points\n");
        exit(EXIT_FAILURE);
    }
    assert(npoints > 2);
    points = (point_t*)malloc( npoints * sizeof(*points) );
    assert(points);
    for (i=0; i<npoints; i++) {
        if (2 != fscanf(f, "%lf %lf", &(points[i].x), &(points[i].y))) {
            fprintf(stderr, "FATAL: failed to get coordinates of point %d\n", i);
            exit(EXIT_FAILURE);
        }
    }
    pset->n = npoints;
    pset->p = points;
}

/**
 * Free the memory allocated by structure pset.
 */
void free_pointset( points_t *pset )
{
    pset->n = 0;
    free(pset->p);
    pset->p = NULL;
}

/**
 * Dump the convex hull to file f. The first line is the number of
 * dimensione (always 2); the second line is the number of vertices of
 * the hull PLUS ONE; the next (n+1) lines are the vertices of the
 * hull, in clockwise order. The first point is repeated twice, in
 * order to be able to plot the result using gnuplot as a closed
 * polygon
 */
void write_hull( FILE *f, const points_t *hull )
{
    int i;
    fprintf(f, "%d\n%d\n", 2, hull->n + 1);
    for (i=0; i<hull->n; i++) {
        fprintf(f, "%f %f\n", hull->p[i].x, hull->p[i].y);
    }
    /* write again the coordinates of the first point */
    fprintf(f, "%f %f\n", hull->p[0].x, hull->p[0].y);    
}

/**
 * Return LEFT, RIGHT or COLLINEAR depending on the shape
 * of the vectors p0p1 and p1p2
 *
 * LEFT            RIGHT           COLLINEAR
 * 
 *  p2              p1----p2            p2
 *    \            /                   /
 *     \          /                   /
 *      p1       p0                  p1
 *     /                            /
 *    /                            /
 *  p0                            p0
 *
 * See Cormen, Leiserson, Rivest and Stein, "Introduction to Algorithms",
 * 3rd ed., MIT Press, 2009, Section 33.1 "Line-Segment properties"
 */
int turn(const point_t p0, const point_t p1, const point_t p2)
{
    /*
      This function returns the correct result (COLLINEAR) also in the
      following cases:
      
      - p0==p1==p2
      - p0==p1
      - p1==p2
    */
    const double cross = (p1.x-p0.x)*(p2.y-p0.y) - (p2.x-p0.x)*(p1.y-p0.y);
    if (cross > 0.0) {
        return LEFT;
    } else {
        if (cross < 0.0) {
            return RIGHT;
        } else {
            return COLLINEAR;
        }
    }
}

/**
 * Get the clockwise angle between the line p0p1 and the vector p1p2 
 *
 *         .
 *        . 
 *       .--+ (this angle) 
 *      .   |    
 *     .    V
 *    p1--------------p2
 *    /
 *   /
 *  /
 * p0
 *
 * The function is not used in this program, but it might be useful.
 */
double cw_angle(const point_t p0, const point_t p1, const point_t p2)
{
    const double x1 = p2.x - p1.x;
    const double y1 = p2.y - p1.y;    
    const double x2 = p1.x - p0.x;
    const double y2 = p1.y - p0.y;
    const double dot = x1*x2 + y1*y2;
    const double det = x1*y2 - y1*x2;
    const double result = atan2(det, dot);
    return (result >= 0 ? result : 2*M_PI + result);
}

/**
 * Compute the convex hull of all points in pset using the "Gift
 * Wrapping" algorithm. The vertices are stored in the hull data
 * structure, that does not need to be initialized by the caller.
 */
void convex_hull(const points_t *pset, points_t *hull)
{
    const int n = pset->n;
    const point_t *p = pset->p;
    points_t global_set;
    int i, j;
    int global_cur, global_next, global_leftmost, global_rightmost;
    double size;
    
    global_set.n = 0;
    global_set.p = (point_t *)malloc(n * sizeof(point_t));

    hull->n = 0;
    /* There can be at most n points in the convex hull. At the end of
       this function we trim the excess space. */
    hull->p = (point_t*)malloc(n * sizeof(*(hull->p))); assert(hull->p);
    
    /* Identify the leftmost point p[leftmost] */
    global_leftmost = 0;
    global_rightmost = 0;
    for (i = 1; i<n; i++) {
        if (p[i].x < p[global_leftmost].x) {
            global_leftmost = i;
        }
        if (p[i].x > p[global_rightmost].x) {
			global_rightmost = i;
		}
    }
    size = p[global_rightmost].x - p[global_leftmost].x;
    
#pragma omp parallel default(none) shared(n_threads, p, global_leftmost, size, global_set) private(i, j)
{
	n_threads = omp_get_num_threads();
	points_t local_set, local_hull;
	int local_leftmost, local_cur, local_next;
	double local_size, local_start, local_end;
	
	/* Define the limits of the local set. */
	local_size = size / omp_get_num_threads();
	local_start = local_size * omp_get_thread_num() + p[global_leftmost].x;
	local_end = local_size * (omp_get_thread_num() + 1) + p[global_leftmost].x;

	/* Create and fill the local set of points on which compute the hull. */
	local_set.n = 0;
	local_set.p = (point_t *)malloc(n * sizeof(point_t));

	for (i=0; i<n; i++) {
		if (p[i].x > local_start && p[i].x <= local_end) {
			local_set.p[local_set.n].x = p[i].x;
			local_set.p[local_set.n].y = p[i].y;
			local_set.n++;
		}
		if (0 == omp_get_thread_num()) {
			if (p[i].x == local_start) {
				local_set.p[local_set.n].x = p[i].x;
				local_set.p[local_set.n].y = p[i].y;
				local_set.n++;
			}
		}
	}
    
    local_leftmost = 0;
    for (i=1; i<local_set.n; i++) {
		if (local_set.p[i].x < local_set.p[local_leftmost].x) {
			local_leftmost = i;
		}
	}
	local_cur = local_leftmost;

	local_hull.n = 0;
	local_hull.p = (point_t *)malloc(local_set.n * sizeof(point_t));

	printf("local(%d)\n", local_set.n);

    /* Main loop of the Gift Wrapping algorithm. This is where most of
       the time is spent; therefore, this is the block of code that
       must be parallelized. */
    if (local_set.n > 0) {
		do {
			/* Add the current vertex to the hull */
			assert(local_hull.n < local_set.n);
			local_hull.p[local_hull.n].x = local_set.p[local_cur].x;
			local_hull.p[local_hull.n].y = local_set.p[local_cur].y;
			local_hull.n++;
        
			/* Search for the next vertex */
			local_next = (local_cur + 1) % local_set.n;
			for (j=0; j<local_set.n; j++) {
				if (LEFT == turn(local_set.p[local_cur], local_set.p[local_next], local_set.p[j])) {
					local_next = j;
				}
			}
			//assert(local_cur != local_next);
			local_cur = local_next;
		} while (local_cur != local_leftmost);
    }

    /* Trim the excess space in the convex hull array */
    local_hull.p = (point_t *)realloc(local_hull.p, local_hull.n * sizeof(point_t));
    //assert(local_hull.p);

#pragma omp critical
{
	for (i=global_set.n; i<global_set.n+local_hull.n; i++) {
		global_set.p[i].x = local_hull.p[i - global_set.n].x;
		global_set.p[i].y = local_hull.p[i - global_set.n].y;
	}
	global_set.n += local_hull.n;
}
}

	global_leftmost = 0;
	for (i=0; i<global_set.n; i++) {
		if (global_set.p[i].x < global_set.p[global_leftmost].x) {
			global_leftmost = i;
		}
	}
	global_cur = global_leftmost;

	do {
		/* Add the current vertex to the hull */
        assert(hull->n < global_set.n);
        hull->p[hull->n] = global_set.p[global_cur];
        hull->n++;
        
        /* Search for the next vertex */
        global_next = (global_cur + 1) % global_set.n;
        for (j=0; j<global_set.n; j++) {
            if (LEFT == turn(global_set.p[global_cur], global_set.p[global_next], global_set.p[j])) {
                global_next = j;
            }
        }
        assert(global_cur != global_next);
        global_cur = global_next;
	} while (global_cur != global_leftmost);

	hull->p = (point_t *)realloc(hull->p, (hull->n) * sizeof(*(hull->p)));
	assert(hull->p);

}

/**
 * Compute the area ("volume", in qconvex terminoloty) of a convex
 * polygon whose vertices are stored in pset using Gauss' area formula
 * (also known as the "shoelace formula"). See:
 *
 * https://en.wikipedia.org/wiki/Shoelace_formula
 *
 * This function does not need to be parallelized.
 */
double hull_volume( const points_t *hull )
{
    const int n = hull->n;
    const point_t *p = hull->p;
    double sum = 0.0;
    int i;
    for (i=0; i<n-1; i++) {
        sum += ( p[i].x * p[i+1].y - p[i+1].x * p[i].y );
    }
    sum += p[n-1].x*p[0].y - p[0].x*p[n-1].y;
    return 0.5*fabs(sum);
}

/**
 * Compute the length of the perimeter ("facet area", in qconvex
 * terminoloty) of a convex polygon whose vertices are stored in pset.
 * This function does not need to be parallelized.
 */
double hull_facet_area( const points_t *hull )
{
    const int n = hull->n;
    const point_t *p = hull->p;
    double length = 0.0;
    int i;
    for (i=0; i<n-1; i++) {
        length += hypot( p[i].x - p[i+1].x, p[i].y - p[i+1].y );
    }
    /* Add the n-th side connecting point n-1 to point 0 */
    length += hypot( p[n-1].x - p[0].x, p[n-1].y - p[0].y );
    return length;
}

int main( void )
{
    points_t pset, hull;
    double tstart, elapsed;
    
    read_input(stdin, &pset);
    tstart = hpc_gettime();
    convex_hull(&pset, &hull);
    elapsed = hpc_gettime() - tstart;
    fprintf(stderr, "\nConvex hull of %d points in 2-d:\n\n", pset.n);
    fprintf(stderr, "  Number of vertices: %d\n", hull.n);
    fprintf(stderr, "  Total facet area: %f\n", hull_facet_area(&hull));
    fprintf(stderr, "  Total volume: %f\n\n", hull_volume(&hull));
    fprintf(stderr, "Elapsed time: %f with %d threads\n\n", elapsed, n_threads);
    write_hull(stdout, &hull);
    free_pointset(&pset);
    free_pointset(&hull);
    return EXIT_SUCCESS;    
}
