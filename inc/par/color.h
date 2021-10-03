#pragma once

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

/**
 * RGB representation of color
 */
typedef struct {
	float r;
	float g;
	float b;
} RGB;

/**
 * Hue, Saturation, and Value representation of color
 */
typedef struct {
	float h;
	float s;
	float v;
} HSV;

RGB HSV2RGB(HSV c1);
HSV RGB2HSV(RGB c1);