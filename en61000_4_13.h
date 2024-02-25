/*
-------------------------------------------------------------------------------
EN6100-4-13 Implementation
Copyright (c) 2024 Cesar Ruano Alvarez

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
-------------------------------------------------------------------------------
*/
#ifndef EN61000_4_13_H
#define EN61000_4_13_H

#include <stdbool.h>

typedef enum{
	EN61000_4_13_FLAT_TOP,
	EN61000_4_13_OVER_SWING,
	EN61000_4_13_FREQUENCY_SWEEP,
	EN61000_4_13_INDIVIDUAL_HARMONICS,
	EN61000_4_13_INTERHARMONICS,
	EN61000_4_13_MEISTER_CURVE,
} EN61000_4_13_Test;

typedef enum{
	EN61000_4_13_CLASS_A,
	EN61000_4_13_CLASS_B,
	EN61000_4_13_CLASS_C,
	LAST_CLASS = EN61000_4_13_CLASS_C
} EN61000_4_13_Class;

typedef struct{
	double rms;
	double base_freq;
	EN61000_4_13_Class class;
	double dwell;
} TestParams;

void* xwave_create();

void xwave_start_test(void* obj, EN61000_4_13_Test test, TestParams* params);

void xwave_step(void* obj, double dt);

void xwave_step_t(void* obj, double t, double* V);

bool xwave_finished(void* obj);

#endif
