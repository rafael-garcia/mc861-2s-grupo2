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

	FILE *f = fopen("gradDir-escala-correta.pgm", "w");
	FILE *g = fopen("gradMag-escala-correta.pgm", "w");
	if (f == NULL || g == NULL) {
	    printf("Error opening file!\n");
	    exit(1);
	}

	int max = iftMaximumValue(gradDir);
	fprintf(f, "P2\n%d %d\n%d\n", gradDir->xsize, gradDir->ysize, max);
	int i = 0;
	for (int j = 0; j < gradDir->xsize && i < gradDir->n; j++) {
		for (int k = 0; k < gradDir->ysize && i < gradDir->n; k++) {
			fprintf(f, "%d ", gradDir->val[i]);
			i++;
		}
		fprintf(f, "\n");
	}

	max = iftMaximumValue(gradMag);
	fprintf(g, "P2\n%d %d\n%d\n", gradMag->xsize, gradMag->ysize, max);
	i = 0;
	for (int j = 0; j < gradMag->xsize && i < gradMag->n; j++) {
		for (int k = 0; k < gradMag->ysize && i < gradMag->n; k++) {
			fprintf(g, "%d ", gradMag->val[i]);
			i++;
		}
		fprintf(g, "\n");
	}
//		fprintf(f, "gradDir[n] = %d, Cb = %d, Cr = %d onde n = %d \n", gradDir->val[i], gradDir->Cb[i], gradDir->Cr[i], i);
//		fprintf(f, "gradDir[n] = %d, onde n = %d \n", gradDir->val[i], i);

	fclose(f);

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

