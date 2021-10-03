#ifndef _BARY_CENTRIC_H
#define _BARY_CENTRIC_H

/*
 *	Return the bary-centric coordinates of a uniformly sampled triangle.
 *  @param n: the level of sampling. n^2 points are sampled.
 */
void uniTriSample(int n, float (*bary)[3]);
#endif

