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

#define EPSLON 1e-6

/**
 * HoG cell abstraction.
 */
typedef struct {
	iftHistogram * histogram;
	int cx;	// center x coordinate.
	int cy; // center y coordinate.
} Cell;

/**
 * Block abstraction
 */
typedef struct {
	float *val;
	int sz;
} Block;

Block createBlock(int sz) {
	Block block;
	block.sz = sz;
	block.val = (float *) malloc(sz * sizeof(float));
	return block;
}

void destroyBlock(Block *bl) {
	free(bl->val);
}

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
void calc_histograms(Cell **cells, int row, int col, int cellSzX, int cellSzY,
		iftImage *gradMag, iftImage *gradDir);

/**************************************************************
 **************************************************************/

void normalizeBlockVector(Block *block, int normFactor) {

	for (int i = 0; i < block->sz; ++i)
		block->val[i] = block->val[i] / (sqrt(normFactor) + EPSLON);
}

iftFeatures *extractHog(iftImage *window) {
	int nOfCells = 4;
	Cell **cells = (Cell**) malloc(nOfCells * sizeof(Cell*));
	int cellXSize = window->xsize / nOfCells;
	int cellYSize = window->ysize / nOfCells;

	int nOfBlocks = 2;
	int blockSz = nOfCells / nOfBlocks;
	Block **blocks = (Block **) malloc(nOfBlocks * sizeof(Block*));

	iftImage *gradMag = iftCreateImage(window->xsize, window->ysize,
			window->zsize);
	iftImage *gradDir = iftCreateImage(window->xsize, window->ysize,
			window->zsize);

	int minY;
	int maxY;
	int minX;
	int maxX;

	/* Create cells */
	for (int i = 0; i < nOfCells; ++i) {
		minY = i * cellYSize;
		maxY = minY + cellYSize - 1;

		cells[i] = (Cell *) malloc(nOfCells * sizeof(Cell));

		for (int j = 0; j < nOfCells; ++j) {
			minX = j * cellXSize;
			maxX = minX + cellXSize - 1;
			cells[i][j].cx = minX + cellXSize / 2;
			cells[i][j].cy = minY + cellYSize / 2;
			cells[i][j].histogram = iftCreateHistogram(9);
		};
	}

	/**Calculates window gradient magnitude and direction**/
	gradient(window, &gradMag, &gradDir);
	/**Calculates histograms for each cell**/
	calc_histograms(cells, nOfCells, nOfCells, cellXSize, cellYSize, gradMag, gradDir);

	int valCounter;
	int bi = 0;
	int bj = 0;
	float normFactor;
	float vj;
	/** Calculates characteristic vector for each block **/
	for (int i = 0; i < nOfCells; i += blockSz) {

		blocks[bi] = (Block *) malloc(nOfBlocks * sizeof(Block));

		for (int j = 0; j < nOfCells; j += blockSz) {
			blocks[bi][bj] = createBlock(blockSz * blockSz * 9);

			valCounter = 0;
			normFactor = 0;
			/**Get each cell of a block**/
			for (int ii = i; ii < blockSz; ++ii) {
				for (int jj = j; jj < blockSz; ++jj) {
					/** Get each bin value of a cell.**/
					for (int b = 0; b < 9; ++b) {
						vj = cells[ii][jj].histogram->val[b];
//						printf("%.2f\n", vj);
						blocks[bi][bj].val[valCounter++] = vj;
						normFactor += (vj * vj);
					}
				}
			}
			normalizeBlockVector(&blocks[bi][bj], normFactor);
			bj++;
		}
		bi++;
	}

	iftFeatures *hogFeatures = iftCreateFeatures(nOfCells * nOfCells * 9);

	int featIndex = 0;
	for (int i = 0; i < nOfBlocks; ++i) {
		for (int j = 0; j < nOfBlocks; ++j) {
			for (int c = 0; c < blocks[i][j].sz; ++c) {
				hogFeatures->val[featIndex++] = blocks[i][j].val[c];
				printf("%f\n", hogFeatures->val[featIndex++]);
			}
		}
	}

	//TODO dealloc cells and blocks
	iftDestroyImage(&gradMag);
	iftDestroyImage(&gradDir);

	return hogFeatures;
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
							+ pow(cells[i][j].cy - py, 2.0));
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

void calc_histograms(Cell **cells, int row, int col, int cellSzX, int cellSzY,
		iftImage *gradMag, iftImage *gradDir) {
	int nearestI[4], nearestJ[4];
	int w[14];
	int ww;
	int center[] = { 22, 67, 112, 157, 202, 247, 292, 337 };
	int bin1, bin2;
	int diff1, diff2;

	int p, q;
	iftVoxel v;
	int xp, yp;
	int x1, x2;
	int y1, y2;
	int auxX, auxY;

	for (p = 0; p < gradMag->n; ++p) {
		v = iftGetVoxelCoord(gradMag, p);
		nearestCells(cells, row, col, v.x, v.y, nearestI, nearestJ);

		/** Find corresponding bins **/
		// find first bin
		bin1 = gradDir->val[p] / 45;

		// find second bin
		if (bin1 + 1 < 8)
			diff1 = center[bin1 + 1] - gradDir->val[p];
		else
			diff1 = 361;

		if (bin1 - 1 >= 0)
			diff2 = center[bin1 - 1] - gradDir->val[p];
		else
			diff2 = 361;

		bin1++;
		if (diff1 < diff2)
			bin2 = bin1 + 1;
		else
			bin2 = bin1 - 1;

		/** Calculate weight parameters **/
		xp = v.x;
		yp = v.y;

		x1 = cells[nearestI[0]][nearestJ[0]].cx;
		y1 = cells[nearestI[0]][nearestJ[0]].cy;

		x2 = cells[nearestI[1]][nearestJ[1]].cx;
		y2 = cells[nearestI[1]][nearestJ[1]].cy;

		auxX = cells[nearestI[2]][nearestJ[2]].cx;
		auxY = cells[nearestI[2]][nearestJ[2]].cy;

		if (x2 == x1)
			x2 = auxX;
		if (y2 == y1)
			y2 = auxY;

		auxX = cells[nearestI[3]][nearestJ[3]].cx;
		auxY = cells[nearestI[3]][nearestJ[3]].cy;

		if (x2 == x1)
			x2 = auxX;
		if (y2 == y1)
			y2 = auxY;

		ww = gradMag->val[p];

		if (ww != 0) {

			/** Now calculate weights **/
			w[0] = ww * abs(x1 - xp) / cellSzX;
			w[1] = ww * abs(xp - x1) / cellSzX;
			w[2] = w[0] * (cellSzY / 2) / cellSzY;
			w[3] = w[0] * (cellSzY / 2) / cellSzY;
			w[4] = w[1] * (cellSzY / 2) / cellSzY;
			w[5] = w[1] * (cellSzY / 2) / cellSzY;
			w[6] = w[2]; //* 45 / 45;
			w[10] = w[2] * 45 / (y1 - center[bin1 - 1]);
			w[7] = w[3];
			w[11] = w[3];
			w[9] = w[4];
			w[13] = w[4];
			w[8] = w[5];
			w[12] = w[5];

			cells[nearestI[0]][nearestJ[0]].histogram->val[bin1] += w[6];
			cells[nearestI[1]][nearestJ[1]].histogram->val[bin1] += w[7];
			cells[nearestI[2]][nearestJ[2]].histogram->val[bin1] += w[9];
			cells[nearestI[3]][nearestJ[3]].histogram->val[bin1] += w[8];

			cells[nearestI[0]][nearestJ[0]].histogram->val[bin2] += w[10];
			cells[nearestI[1]][nearestJ[1]].histogram->val[bin2] += w[11];
			cells[nearestI[2]][nearestJ[2]].histogram->val[bin2] += w[13];
			cells[nearestI[3]][nearestJ[3]].histogram->val[bin2] += w[12];

		} else {
			cells[nearestI[0]][nearestJ[0]].histogram->val[0] += 1;
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
