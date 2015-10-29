/*
 * selectUtils.h
 *
 *  Created on: Oct 28, 2015
 *      Author: ra146446
 */

#ifndef SELECTUTILS_H_
#define SELECTUTILS_H_

#include "ift.h"

float findFormFactor(iftImage *img, int label) {
	iftVoxel v;
	int p;
	int minY, maxY;
	int minX, maxX;
	float formFactor;

	v.z = 0;

	/**Get component limits in Y**/
	/****Get minY value****/
	for (v.y = 0; v.y < img->ysize; ++v.y) {
		for (v.x = 0; v.x < img->xsize; ++v.x) {
			p = iftGetVoxelIndex(img, v);
			if (img->val[p] == label) {
				minY = v.y;
				break;
			}
		}
	}
	/****Get maxY value****/
	for (v.y = img->ysize - 1; v.y >= 0; --v.y) {
		for (v.x = 0; v.x < img->xsize; ++v.x) {
			p = iftGetVoxelIndex(img, v);
			if (img->val[p] == label) {
				maxY = v.y;
				break;
			}
		}
	}

	/**Get component limits in X**/
	/****Get minX value****/
	for (v.x = 0; v.x < img->xsize; ++v.x) {
		for (v.y = 0; v.y < img->ysize; ++v.y) {
			p = iftGetVoxelIndex(img, v);
			if (img->val[p] == label) {
				minX = v.x;
				break;
			}
		}
	}

	/****Get maxX value****/
	for (v.x = img->xsize - 1; v.x >= 0; --v.x) {
		for (v.y = 0; v.y < img->ysize; ++v.y) {
			p = iftGetVoxelIndex(img, v);
			if (img->val[p] == label) {
				maxX = v.x;
				break;
			}
		}
	}

	formFactor = (maxX - minX) / (maxY - minY);
	return formFactor;
}
/**
 * @brief Remove components with less than minVolume pixels.
 */
void removeNonRectangular(iftImage *img) {
	int p;
	int nOfElements = iftMaximumValue(img);
	float *formFactor = iftAllocFloatArray(nOfElements + 1);
	int *labels = iftAllocIntArray(nOfElements + 1);

	for (int label = 1; label <= nOfElements; ++label)
		formFactor[label] = findFormFactor(img, label);

	int nlabels = 1;
	for (int label = 1; label <= nOfElements; ++label) {
		if (formFactor[label] > 2)
			labels[label] = nlabels++;
		else
			if (label != nOfElements) //check for biggest probability
				labels[label] = 0;
			else
				labels[label] = nlabels++;
	}

	for (p = 0; p < img->n; ++p) {
		img->val[p] = labels[img->val[p]];
	}

	free(formFactor);
	free(labels);
}

#endif /* SELECTUTILS_H_ */
