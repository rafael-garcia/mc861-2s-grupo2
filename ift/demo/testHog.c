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

//	iftFeatures *feat = extractHog(img);

	iftImage *norm = normalize(img);

	iftWriteImageP2(norm, "normalized.pgm");

	iftDestroyImage(&img);
	iftDestroyImage(&norm);

	return 0;

}

