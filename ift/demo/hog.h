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
	iftHistogram * histogram;
	int cx;	// center x coordinate.
	int cy; // center y coordinate.
} Cell;

/**
 * Extract image features using HoG image descriptor.
 *
 * @param[in] img 	gray scale image.
 * @return			features structure containing characteristic vector.
 */
iftFeatures *extractHog(iftImage *window);

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
void calc_histograms(Cell **cells, int row, int col, iftImage *gradMag,
		iftImage *gradDir);

/**
 *
 */
Cell create_cell(iftImage *window, int xIni, int xEnd, int yIni, int yEnd);

/**************************************************************
 **************************************************************/

iftFeatures *extractHog(iftImage *window) {
	int nOfCells = 4;
	Cell **cells = (Cell**) malloc(nOfCells * sizeof(Cell*));
	int cellXSize = window->xsize / nOfCells;
	int cellYSize = window->ysize / nOfCells;
//	iftImage *cellImg = iftCreateImage(cellXSize, cellYSize, window->zsize);
	iftImage *gradMag = iftCreateImage(window->xsize, window->ysize,
			window->zsize);
	iftImage *gradDir = iftCreateImage(window->xsize, window->ysize,
			window->zsize);

	int minY;
	int maxY;
	int minX;
	int maxX;

	int x;
	int y;
	int origp;
	int p;

	/* Create cells */
	for (int i = 0; i < nOfCells; ++i) {
		minY = i * cellYSize;
		maxY = minY + cellYSize - 1;

		cells[i] = (Cell *) malloc(nOfCells * sizeof(Cell));

		for (int j = 0; j < nOfCells; ++j) {
			minX = j * cellXSize;
			maxX = minX + cellXSize - 1;

//			/* Copy pixel values from original window to a image representing
//			 * cell. */
//			for (x = minY; x < (maxY + 1); x++) {
//				for (y = minX; y < (maxX + 1); y++) {
//					origp = x * window->xsize + y;
//					p = (x - minY) * (maxX - minX + 1) + (y - minX);
//					cellImg->val[p] = window->val[origp];
//					// Copy Cb and Cr bands for color images
//					if (window->Cb != NULL) {
//						cellImg->Cb[p] = window->Cb[origp];
//						cellImg->Cr[p] = window->Cr[origp];
//					}
//				}
//			}
//
//			gradient(cellImg, &(cells[i][j].gradMag), &(cells[i][j].gradDir));

			cells[i][j].cx = minX + cellXSize / 2;
			cells[i][j].cy = minY + cellYSize / 2;
		};

	}

	gradient(window, &gradMag, &gradDir);
	calc_histograms(cells, nOfCells, nOfCells, gradMag, gradDir);

	//TODO dealloc cells
	iftDestroyImage(&gradMag);
	iftDestroyImage(&gradDir);

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

typedef struct {
	float val;
	int i;
	int j;
} Distance;

int cmp(const void *a, const void *b) {
	Distance *ad = (Distance*) a;
	Distance *bd = (Distance*) b;
	return ad->val > bd->val;
}

void nearestCells(Cell **cells, int row, int col, int px, int py, int *nearestI,
		int *nearestJ) {

	Distance dist[row * col];
	int aux = 0;

	for (int i = 0; i < row; ++i) {
		for (int j = 0; j < col; ++j) {
			dist[aux].val = sqrt(
					pow(cells[i][j].cx - px, 2.0)
							+ pow(cells[i][j].cx - px, 2.0));
			dist[aux].i = i;
			dist[aux++].j = j;
		}
	}
	qsort(dist, row * col, sizeof(Distance), cmp);

	nearestI[0] = dist[0].i;
	nearestI[1] = dist[1].i;
	nearestI[2] = dist[2].i;
	nearestI[3] = dist[3].i;

	nearestJ[0] = dist[0].j;
	nearestJ[1] = dist[1].j;
	nearestJ[2] = dist[2].j;
	nearestJ[3] = dist[3].j;

}

void calc_histograms(Cell **cells, int row, int col, iftImage *gradMag,
		iftImage *gradDir) {
	int nearestI[4];
	int nearestJ[4];
	int w[8];
	int x[14];
	int y[14];
	int z[14];
	int q;
	int ww;

	int p;
	iftVoxel v;
	int xp;
	int yp;

	for (int i = 0; i < row; ++i) {
		for (int j = 0; j < col; ++j) {

			for (p = 0; p < gradMag->n; ++p) {
				v = iftGetVoxelCoord(gradMag, p);
				nearestCells(cells, row, col, v.x, v.y, nearestI, nearestJ);

				/** Calculate weight parameters **/

				xp = v.x;
				yp = v.y;

				y[0] = xp;

				y[1] = yp;

				x[6] = cells[nearestI[0]][nearestJ[0]].cx;
				y[6] = cells[nearestI[0]][nearestJ[0]].cy;
				v.x = x[6];
				v.y = y[6];
				q = iftGetVoxelIndex(gradDir, v);
				z[6] = gradDir->val[q];

				x[0] = x[6];
				v.x = x[0];
				v.y = y[0];
				q = iftGetVoxelIndex(gradDir, v);
				z[0] = gradDir->val[q];

				x[2] = x[6];
				y[2] = y[6];
				z[2] = z[6];

				x[10] = x[6];
				y[10] = y[6];
				z[10] = z[6];

				x[7] = cells[nearestI[1]][nearestJ[1]].cx;
				y[7] = cells[nearestI[1]][nearestJ[1]].cy;
				v.x = x[7];
				v.y = y[7];
				q = iftGetVoxelIndex(gradDir, v);
				z[7] = gradDir->val[q];

				x[3] = x[7];
				y[3] = y[7];
				z[3] = z[7];

				x[11] = x[7];
				y[11] = y[7];
				z[11] = z[7];

				x[8] = cells[nearestI[2]][nearestJ[2]].cx;
				y[8] = cells[nearestI[2]][nearestJ[2]].cy;
				v.x = x[8];
				v.y = y[8];
				q = iftGetVoxelIndex(gradDir, v);
				z[8] = gradDir->val[q];

				x[1] = x[8];
				v.x = x[1];
				v.y = y[1];
				q = iftGetVoxelIndex(gradDir, v);
				z[1] = gradDir->val[q];

				x[5] = x[8];
				y[5] = y[8];
				z[5] = z[8];

				x[12] = x[8];
				y[12] = y[8];
				z[12] = z[8];

				x[9] = cells[nearestI[3]][nearestJ[3]].cx;
				y[9] = cells[nearestI[3]][nearestJ[3]].cy;
				v.x = x[9];
				v.y = y[9];
				q = iftGetVoxelIndex(gradDir, v);
				z[9] = gradDir->val[q];

				x[4] = x[9];
				y[4] = y[9];
				z[4] = z[9];

				x[13] = x[9];
				y[13] = y[9];
				z[13] = z[9];

				ww = gradMag->val[p];

				/** Now calculate weights **/
			w[0] = ww   * (x[1] - xp) / (x[1] - x[0]);
			w[1] = ww   * (xp - x[0]) / (x[1] - x[0]);
			w[2] = w[1] * (xp - x[0]) / (x[1] - x[0]);

		}

	}

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
