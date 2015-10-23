/*
 * hog.h
 *
 *  Created on: Oct 18, 2015
 *      Author: ra146446
 */

#ifndef HOG_H_
#define HOG_H_

#include "ift.h"
#include "iftExtractFeatures.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

/**
 * HoG cell abstraction.
 */
typedef struct {
	iftImage *img;
	int cx;	// center x coordinate.
	int cy; // center y coordinate.

} Cell;

/**
 * Extract image features using HoG image descriptor.
 *
 * @param[in] img 	gray scale image.
 * @return			features structure containing characteristic vector.
 */
iftFeatures *extractHog(iftImage *img);

/**
 * Normalize image
 * @param[in]	img			gray scale image.
 * @return					normalized image.
 */
iftImage *normalize(iftImage *img);

/**
 * Calculates image gradient.
 *
 * @param[in]	img			gray scale image.
 * @param[out]	magnitude	gradient magnitude.
 * @param[out]	direction	gradient direction.
 */
void gradient(iftImage *img, iftImage **magnitude, iftImage **direction);

/**
 * Calculates histrograms of a block.
 *
 * @param[in]	img			gray scale image.
 * @param[out]	magnitude	gradient magnitude.
 * @param[out]	direction	gradient direction.
 * @return					array of int representing histograms
 */
iftHistogram *create_histogram(iftImage *magnitude, iftImage *direction);

/**
 *
 */
Cell create_cell(iftImage *img, int xIni, int xEnd, int yIni, int yEnd);

/**************************************************************
 **************************************************************/

iftFeatures *extractHog(iftImage *img) {
	int nOfCells = 4;
	Cell cells[nOfCells * nOfCells];
	int cellXSize = img->xsize / nOfCells;
	int cellYSize = img->ysize / nOfCells;

	for (int i = 0; i < nOfCells; ++i) {
		;
	}

	return NULL;
}

void gradient(iftImage *img, iftImage **magnitude, iftImage **direction) {
	iftImage* normImg = iftCreateImage(img->xsize, img->ysize, img->zsize);
	iftVoxel v;
	iftVoxel u;
	float r = 5.0;
	float sigma = r / 3;
	iftAdjRel* adj = iftCircular(r);

	int p;
	int q;
	int i;

	float g;
	float gx;
	float gy;
	float factor;
	float distance;
	float xPart;
	float yPart;
	float _const = 180 / PI;

	*magnitude = iftCreateImage(img->xsize, img->ysize, img->zsize);
	*direction = iftCreateImage(img->xsize, img->ysize, img->zsize);

	for (p = 0; p < img->n; ++p) { //foreach pixel in the image
		v = iftGetVoxelCoord(img, p); //gets the multidimensional coordinate using the unidimensional index

		gx = 0;
		gy = 0;
		for (i = 1; i < adj->n; ++i) { //foreach neighbor voxel in the adjacency

			u = iftGetAdjacentVoxel(adj, v, i);

			if (iftValidVoxel(img, u)) { //check if it is a valid voxel (is it inside the image?)
				q = iftGetVoxelIndex(img, u);

				distance = sqrt(pow(v.x - u.x, 2.0) + pow(v.y - u.y, 2.0));
				factor = img->val[p] * img->val[q]
						* exp(-(pow(distance, 2.0) / (2 * pow(sigma, 2.0))));

				xPart = (v.x - u.x) / distance;
				yPart = (v.y - u.y) / distance;

				gx += factor * xPart;
				gy += factor * yPart;

			}
		}

		g = sqrt(pow(gx, 2.0) + pow(gy, 2.0));
		(*magnitude)->val[p] = g;

		if (gy / g >= 0) {
			(*direction)->val[p] = _const * acos(gx / g);
		} else {
			(*direction)->val[p] = 360 - _const * acos(gx / g);
		}

	}

}

iftHistogram *create_histogram(iftImage *magnitude, iftImage *direction) {
	int numOfBins = 9;
	iftHistogram *hist = iftCreateHistogram(numOfBins);

	for (int p = 0; p < magnitude->n; ++p) {
		;
	}
}

iftImage *normalize(iftImage *img) {
	iftImage* normImg = iftCreateImage(img->xsize, img->ysize, img->zsize);
	iftVoxel v;
	iftVoxel u;
	iftAdjRel* adj = iftCircular(5.0);

	int p;
	int q;
	int i;
	double sum;

	for (p = 0; p < img->n; ++p) { //foreach pixel in the image
		v = iftGetVoxelCoord(img, p); //gets the multidimensional coordinate using the unidimensional index

		sum = 0;
		for (i = 0; i < adj->n; ++i) { //foreach neighbor voxel in the adjacency

			u = iftGetAdjacentVoxel(adj, v, i);

			if (iftValidVoxel(img, u)) { //check if it is a valid voxel (is it inside the image?)
				q = iftGetVoxelIndex(img, u);
				sum += (img->val[q] * img->val[q]);
			}
		}

		normImg->val[p] = ((img->val[p] / sqrt(sum))) * 255;
	}

	return normImg;
}

#endif /* HOG_H_ */
