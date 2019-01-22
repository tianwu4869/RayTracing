#include "stdafx.h"
#include "image.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <math.h>
using namespace std;

/**
 * Image
 **/
Image::Image (int width_, int height_){

    assert(width_ > 0);
    assert(height_ > 0);

    width           = width_;
    height          = height_;
    num_pixels      = width * height;
    sampling_method = IMAGE_SAMPLING_POINT;
    
    data.raw = new uint8_t[num_pixels*4];
	int b = 0; //which byte to write to
	for (int j = 0; j < height; j++){
		for (int i = 0; i < width; i++){
			data.raw[b++] = 0;
			data.raw[b++] = 0;
			data.raw[b++] = 0;
			data.raw[b++] = 0;
		}
	}

    assert(data.raw != NULL);
}

Image::Image (const Image& src){
	
	width           = src.width;
    height          = src.height;
    num_pixels      = width * height;
    sampling_method = IMAGE_SAMPLING_POINT;
    
    data.raw = new uint8_t[num_pixels*4];
    
    //memcpy(data.raw, src.data.raw, num_pixels);
    *data.raw = *src.data.raw;
}

Image::Image (char* fname){

	int numComponents; //(e.g., Y, YA, RGB, or RGBA)
	data.raw = stbi_load(fname, &width, &height, &numComponents, 4);
	
	if (data.raw == NULL){
		printf("Error loading image: %s", fname);
		exit(-1);
	}
	

	num_pixels = width * height;
	sampling_method = IMAGE_SAMPLING_POINT;
	
}

Image::~Image (){
    delete data.raw;
    data.raw = NULL;
}

void Image::Write(char* fname){
	
	int lastc = strlen(fname);

	switch (fname[lastc-1]){
	   case 'g': //jpeg (or jpg) or png
	     if (fname[lastc-2] == 'p' || fname[lastc-2] == 'e') //jpeg or jpg
	        stbi_write_jpg(fname, width, height, 4, data.raw, 95);  //95% jpeg quality
	     else //png
	        stbi_write_png(fname, width, height, 4, data.raw, width*4);
	     break;
	   case 'a': //tga (targa)
	     stbi_write_tga(fname, width, height, 4, data.raw);
	     break;
	   case 'p': //bmp
	   default:
	     stbi_write_bmp(fname, width, height, 4, data.raw);
	}
}
/**
 * Image Sample
 **/
void Image::SetSamplingMethod(int method)
{
    assert((method >= 0) && (method < IMAGE_N_SAMPLING_METHODS));
    sampling_method = method;
}


Pixel Image::Sample (double u, double v){
	Pixel temp;
	switch (sampling_method)
	{
	case 0: {//point sampling
		int ori_x, ori_y;
		ori_x = round(u);
		ori_y = round(v);
		if (ori_x > (Width() - 1)) {
			ori_x = (Width() - 1);
		}
		if (ori_y > (Height() - 1)) {
			ori_y = (Height() - 1);
		}
		temp = GetPixel(ori_x, ori_y);
		break;
	}
	case 1: {//linear sampling
		int upLeft_x, upLeft_y, upRight_x, upRight_y, downLeft_x, downLeft_y, downRight_x, downRight_y;
		upLeft_x = floor(u);
		upLeft_y = floor(v);
		upRight_x = ceil(u);
		upRight_y = floor(v);
		downLeft_x = floor(u);
		downLeft_y = ceil(v);
		downRight_x = ceil(u);
		downRight_y = ceil(v);
		if (upRight_x >= Width()) {
			upRight_x = Width() - 1;
		}
		else if (upRight_x < 0) {
			upRight_x = 0;
		}
		if (upRight_y >= Height()) {
			upRight_y = Height() - 1;
		}
		else if (upRight_y < 0) {
			upRight_y = 0;
		}

		if (upLeft_x >= Width()) {
			upLeft_x = Width() - 1;
		}
		else if (upLeft_x < 0) {
			upLeft_x = 0;
		}
		if (upLeft_y >= Height()) {
			upLeft_y = Height() - 1;
		}
		else if (upLeft_y < 0) {
			upLeft_y = 0;
		}

		if (downRight_x >= Width()) {
			downRight_x = Width() - 1;
		}
		else if (downRight_x < 0) {
			downRight_x = 0;
		}
		if (downRight_y >= Height()) {
			downRight_y = Height() - 1;
		}
		else if (downRight_y < 0) {
			downRight_y = 0;
		}

		if (downLeft_x >= Width()) {
			downLeft_x = Width() - 1;
		}
		else if (downLeft_x < 0) {
			downLeft_x = 0;
		}
		if (downLeft_y >= Height()) {
			downLeft_y = Height() - 1;
		}
		else if (downLeft_y < 0) {
			downLeft_y = 0;
		}
		Pixel upLeft, upRight, downLeft, downRight;
		upLeft = GetPixel(upLeft_x, upLeft_y);
		upRight = GetPixel(upRight_x, upRight_y);
		downLeft = GetPixel(downLeft_x, downLeft_y);
		downRight = GetPixel(downRight_x, downRight_y);
		Pixel up, down;
		float alpha, beta;
		alpha = downRight_x - u;
		beta = downRight_y - v;
		up.r = upRight.r - (upRight.r - upLeft.r) * alpha;
		up.g = upRight.g - (upRight.g - upLeft.g) * alpha;
		up.b = upRight.b - (upRight.b - upLeft.b) * alpha;
		down.r = downRight.r - (downRight.r - downLeft.r) * alpha;
		down.g = downRight.g - (downRight.g - downLeft.g) * alpha;
		down.b = downRight.b - (downRight.b - downLeft.b) * alpha;
		temp.r = down.r - (down.r - up.r) * beta;
		temp.g = down.g - (down.g - up.g) * beta;
		temp.b = down.b - (down.b - up.b) * beta;
		break;
	}
	case 2: {//Gaussian Sampling
		int upLeft_x, upLeft_y, upRight_x, upRight_y, downLeft_x, downLeft_y, downRight_x, downRight_y;
		upLeft_x = floor(u);
		upLeft_y = floor(v);
		upRight_x = ceil(u);
		upRight_y = floor(v);
		downLeft_x = floor(u);
		downLeft_y = ceil(v);
		downRight_x = ceil(u);
		downRight_y = ceil(v);
		if (upRight_x >= Width()) {
			upRight_x = Width() - 1;
		}
		if (downRight_x >= Width()) {
			downRight_x = Width() - 1;
		}
		if (downRight_y >= Height()) {
			downRight_y = Height() - 1;
		}
		if (downLeft_y >= Height()) {
			downLeft_y = Height() - 1;
		}
		Pixel upLeft, upRight, downLeft, downRight;
		upLeft = GetPixel(upLeft_x, upLeft_y);
		upRight = GetPixel(upRight_x, upRight_y);
		downLeft = GetPixel(downLeft_x, downLeft_y);
		downRight = GetPixel(downRight_x, downRight_y);
		float up, down, left, right;
		right = downRight_x - u;
		down = downRight_y - v;
		up = 1 - down;
		left = 1 - right;
		float tr, tg, tb, upLeftWeight, upRightWeight, downRightWeight, downLeftWeight;
		upLeftWeight = 0.3989 * exp((pow(up, 2) + pow(left, 2)) / (-2));
		upRightWeight = 0.3989 * exp((pow(up, 2) + pow(right, 2)) / (-2));
		downLeftWeight = 0.3989 * exp((pow(down, 2) + pow(left, 2)) / (-2));
		downRightWeight = 0.3989 * exp((pow(down, 2) + pow(right, 2)) / (-2));
		tr = upLeft.r * upLeftWeight + upRight.r * upRightWeight + downLeft.r * downLeftWeight + downRight.r * downRightWeight;
		tg = upLeft.g * upLeftWeight + upRight.g * upRightWeight + downLeft.g * downLeftWeight + downRight.g * downRightWeight;
		tb = upLeft.b * upLeftWeight + upRight.b * upRightWeight + downLeft.b * downLeftWeight + downRight.b * downRightWeight;

		float totalWeight = upLeftWeight + upRightWeight + downLeftWeight + downRightWeight;
		
		float tempDistance, tempWeight;
		if ((upLeft_y - 1) >= 0) {
			tempDistance = pow(left, 2) + pow((1 + up), 2);
			tempWeight = 0.3989 * exp(tempDistance / (-2));
			totalWeight = totalWeight + tempWeight;
			tr = tr + GetPixel(upLeft_x, upLeft_y - 1).r * tempWeight;
			tg = tg + GetPixel(upLeft_x, upLeft_y - 1).g * tempWeight;
			tb = tb + GetPixel(upLeft_x, upLeft_y - 1).b * tempWeight;
			tempDistance = pow(right, 2) + pow((1 + up), 2);
			tempWeight = 0.3989 * exp(tempDistance / (-2));
			totalWeight = totalWeight + tempWeight;
			tr = tr + GetPixel(upRight_x, upLeft_y - 1).r * tempWeight;
			tg = tg + GetPixel(upRight_x, upLeft_y - 1).g * tempWeight;
			tb = tb + GetPixel(upRight_x, upLeft_y - 1).b * tempWeight;
		}

		if ((downLeft_y + 1) < Height()) {
			tempDistance = pow(left, 2) + pow((1 + down), 2);
			tempWeight = 0.3989 * exp(tempDistance / (-2));
			totalWeight = totalWeight + tempWeight;
			tr = tr + GetPixel(downLeft_x, downLeft_y + 1).r * tempWeight;
			tg = tg + GetPixel(downLeft_x, downLeft_y + 1).g * tempWeight;
			tb = tb + GetPixel(downLeft_x, downLeft_y + 1).b * tempWeight;
			tempDistance = pow(right, 2) + pow((1 + down), 2);
			tempWeight = 0.3989 * exp(tempDistance / (-2));
			totalWeight = totalWeight + tempWeight;
			tr = tr + GetPixel(downRight_x, downRight_y + 1).r * tempWeight;
			tg = tg + GetPixel(downRight_x, downRight_y + 1).g * tempWeight;
			tb = tb + GetPixel(downRight_x, downRight_y + 1).b * tempWeight;
		}
		
		if ((downLeft_x - 1) >= 0) {
			tempDistance = pow(up, 2) + pow((1 + left), 2);
			tempWeight = 0.3989 * exp(tempDistance / (-2));
			totalWeight = totalWeight + tempWeight;
			tr = tr + GetPixel(upLeft_x - 1, upLeft_y).r * tempWeight;
			tg = tg + GetPixel(upLeft_x - 1, upLeft_y).g * tempWeight;
			tb = tb + GetPixel(upLeft_x - 1, upLeft_y).b * tempWeight;
			tempDistance = pow(down, 2) + pow((1 + left), 2);
			tempWeight = 0.3989 * exp(tempDistance / (-2));
			totalWeight = totalWeight + tempWeight;
			tr = tr + GetPixel(downLeft_x - 1, downLeft_y).r * tempWeight;
			tg = tg + GetPixel(downLeft_x - 1, downLeft_y).g * tempWeight;
			tb = tb + GetPixel(downLeft_x - 1, downLeft_y).b * tempWeight;
		}

		if ((downRight_x + 1) < Width()) {
			tempDistance = pow(up, 2) + pow((1 + right), 2);
			tempWeight = 0.3989 * exp(tempDistance / (-2));
			totalWeight = totalWeight + tempWeight;
			tr = tr + GetPixel(upRight_x + 1, upRight_y).r * tempWeight;
			tg = tg + GetPixel(upRight_x + 1, upRight_y).g * tempWeight;
			tb = tb + GetPixel(upRight_x + 1, upRight_y).b * tempWeight;
			tempDistance = pow(down, 2) + pow((1 + right), 2);
			tempWeight = 0.3989 * exp(tempDistance / (-2));
			totalWeight = totalWeight + tempWeight;
			tr = tr + GetPixel(downRight_x + 1, downRight_y).r * tempWeight;
			tg = tg + GetPixel(downRight_x + 1, downRight_y).g * tempWeight;
			tb = tb + GetPixel(downRight_x + 1, downRight_y).b * tempWeight;
		}
		
		temp.r = tr / totalWeight;
		temp.g = tg / totalWeight;
		temp.b = tb / totalWeight;

		break;
	}
	default:
		break;
	}

	return temp;
}