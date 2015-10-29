/*
 * testHog.c
 *
 *  Created on: Oct 18, 2015
 *      Author: ra146446
 */

#include <stdio.h>
#include "hog.h"

int main(int argc, char *argv[]) {
	if (argc < 2) {
		printf("usage:\n\t %s <image_path>\n", argv[0]);
		exit(1);
	}

	iftImage *img = iftReadImageByExt(argv[1]);

	iftImage *norm = normalize(img);
	iftWriteImageP2(norm, "norm.pgm");
	iftDestroyImage(&norm);

	iftImage *gradMag;
	iftImage *gradDir;
	gradient(img, &gradMag, &gradDir);
	iftWriteImageP2(gradMag, "grad-mag.pgm");
	iftWriteImageP2(gradDir, "grad-dir.pgm");
	iftDestroyImage(&gradMag);
	iftDestroyImage(&gradDir);

	iftImage *label = iftReadImageByExt(
			"/tmp/placas/LicensePlates/cand/0001.pgm");

	iftImage *box = iftCreateBoundingBox2D(img, label, 1);

	iftFeatures *feat = extractHog(box);

	iftWriteImageP2(box, "box.pgm");

	iftDestroyImage(&img);
	iftDestroyImage(&box);
	iftDestroyFeatures(&feat);

	return 0;

}

