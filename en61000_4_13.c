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

#include "en61000_4_13.h"
#include "stdbool.h"
#include <math.h>
#include <stdio.h>


#define GET_OBJ TestInstance obj = (TestInstance)_obj
#define PI 3.14159265358979323846
#define RMS_TO_AMP 1.414

#define _180_DEG (PI)

#define SAT_HI(x, max) ((x>max)?(max):(x))
#define SAT_LO(x, min) SAT_HI(min, x)
#define SAT_HI_LO(x, max, min) (SAT_HI(SAT_LO(x, min), max))

#define NUM_PH 3

typedef void (*StepCallback)(void* _obj, double t, double* V, double fade);


typedef struct{
	double amplitude;
	double period;
} DerivedTestParams;

typedef struct{
	double V[3];
	double t;
	double prev_t;
	bool finished;
	bool fade_out_finished;
	bool fade_in_finished;
	TestParams params;
	DerivedTestParams derived_params;
	StepCallback step_callback;
	double fade_acc;
	int istep;
	double dstep;
	double dwell_acc;
} TestInstance;

typedef struct {
	double clip_prop;
} FlatTopEntry;

FlatTopEntry FLAT_TOP_PARAMS[LAST_CLASS+1] = {
		[EN61000_4_13_CLASS_A] = {.clip_prop = 0.95},
		[EN61000_4_13_CLASS_B] = {.clip_prop = 0.90},
		[EN61000_4_13_CLASS_C] = {.clip_prop = 0.80},
};

typedef struct {
	double third_prop;
	double fifth_prop;
} OverSwingEntry;

OverSwingEntry OVER_SWING_PARAMS[LAST_CLASS+1] = {
		[EN61000_4_13_CLASS_A] = {.third_prop = 0.08, .fifth_prop=0.05},
		[EN61000_4_13_CLASS_B] = {.third_prop = 0.08, .fifth_prop=0.05},
		[EN61000_4_13_CLASS_C] = {.third_prop = 0.08, .fifth_prop=0.05},
};

const double INDIVIDUAL_HARMONIC_PARAMS[LAST_CLASS+1][41] = {
	[EN61000_4_13_CLASS_A] =
	{
		[0]  = 0.0  , [1]  = 0.0,   [2]  = 0.05,  [3]  = 0.09,  [4]  = 0.015,
		[5]  = 0.12,  [6]  = 0.015, [7]  = 0.10,  [8]  = 0.015, [9]  = 0.04,
		[10] = 0.015, [11] = 0.07,  [12] = 0.015, [13] = 0.07,  [14] = 0.015,
		[15] = 0.03,  [16] = 0.015, [17] = 0.06,  [18] = 0.015, [19] = 0.06,
		[20] = 0.015, [21] = 0.02,  [22] = 0.015, [23] = 0.06,  [24] = 0.015,
		[25] = 0.06,  [26] = 0.015, [27] = 0.02,  [28] = 0.015, [29] = 0.05,
		[30] = 0.015, [31] = 0.03,  [32] = 0.015, [33] = 0.015, [34] = 0.015,
		[35] = 0.03,  [36] = 0.015, [37] = 0.03,  [38] = 0.015, [39] = 0.02, [40] = 0.015
	},
	[EN61000_4_13_CLASS_B] =
	{
		[0]  = 0.0,   [1]  = 0.0,   [2]  = 0.03,  [3]  = 0.08,  [4]  = 0.015,
		[5]  = 0.09,  [6]  = 0.015, [7]  = 0.075, [8]  = 0.015, [9]  = 0.025,
		[10] = 0.015, [11] = 0.05,  [12] = 0.015, [13] = 0.045, [14] = 0.015,
		[15] = 0.03,  [16] = 0.015, [17] = 0.03,  [18] = 0.015, [19] = 0.02,
		[20] = 0.015, [21] = 0.015, [22] = 0.015, [23] = 0.02,  [24] = 0.015,
		[25] = 0.02,  [26] = 0.015, [27] = 0.02,  [28] = 0.015, [29] = 0.015,
		[30] = 0.015, [31] = 0.015, [32] = 0.015, [33] = 0.02,  [34] = 0.015,
		[35] = 0.015, [36] = 0.015, [37] = 0.015, [38] = 0.015, [39] = 0.02, [40] = 0.015
	},
	[EN61000_4_13_CLASS_C] =
	{
		[0]  = 0.0,   [1]  = 0.0,   [2]  = 0.03,  [3]  = 0.045, [4]  = 0.015,
		[5]  = 0.045, [6]  = 0.015, [7]  = 0.045, [8]  = 0.015, [9]  = 0.02,
		[10] = 0.015, [11] = 0.045, [12] = 0.015, [13] = 0.04,  [14] = 0.015,
		[15] = 0.03,  [16] = 0.015, [17] = 0.03,  [18] = 0.015, [19] = 0.02,
		[20] = 0.015, [21] = 0.015, [22] = 0.015, [23] = 0.02,  [24] = 0.015,
		[25] = 0.02,  [26] = 0.015, [27] = 0.02,  [28] = 0.015, [29] = 0.015,
		[30] = 0.015, [31] = 0.015, [32] = 0.015, [33] = 0.02,  [34] = 0.015,
		[35] = 0.015, [36] = 0.015, [37] = 0.015, [38] = 0.015, [39] = 0.02, [40] = 0.015
	},
};

typedef struct {
	double freq_prop;
	double dfreq;
	double prop;
} SweepEntry;

const SweepEntry FREQUENCY_SWEEP_PARAMS[LAST_CLASS+1][6] = {
		[EN61000_4_13_CLASS_A] =
		{
			[0]  = {.freq_prop = 2.0,  .dfreq = 0.10, .prop = 0.02},
			[1]  = {.freq_prop = 10.0, .dfreq = 0.20, .prop = 0.05},
			[2]  = {.freq_prop = 20.0, .dfreq = 0.20, .prop = 0.04},
			[4]  = {.freq_prop = 30.0, .dfreq = 0.50, .prop = 0.02},
			[5]  = {.freq_prop = 40.0, .dfreq = 0.50, .prop = 0.02},
			[6]  = {.freq_prop = 40.0, .dfreq = 0.0,  .prop = 0.0},
		},
		[EN61000_4_13_CLASS_B] =
		{
			[0]  = {.freq_prop = 2.0,  .dfreq = 0.10, .prop = 0.03},
			[1]  = {.freq_prop = 10.0, .dfreq = 0.20, .prop = 0.09},
			[2]  = {.freq_prop = 20.0, .dfreq = 0.20, .prop = 0.045},
			[4]  = {.freq_prop = 30.0, .dfreq = 0.50, .prop = 0.02},
			[5]  = {.freq_prop = 40.0, .dfreq = 0.50, .prop = 0.02},
			[6]  = {.freq_prop = 40.0, .dfreq = 0.0,  .prop = 0.0},
		},
		[EN61000_4_13_CLASS_C] =
		{
			[0]  = {.freq_prop = 2.0,  .dfreq = 0.10, .prop = 0.045},
			[1]  = {.freq_prop = 10.0, .dfreq = 0.20, .prop = 0.14},
			[2]  = {.freq_prop = 20.0, .dfreq = 0.20, .prop = 0.09},
			[4]  = {.freq_prop = 30.0, .dfreq = 0.50, .prop = 0.06},
			[5]  = {.freq_prop = 40.0, .dfreq = 0.50, .prop = 0.04},
			[6]  = {.freq_prop = 40.0, .dfreq = 0.0,  .prop = 0.0},
		},
};


const SweepEntry INTER_HARMONIC_PARAMS[LAST_CLASS+1][5] = {
		[EN61000_4_13_CLASS_A] =
		{
			[0]  = {.freq_prop = 2.0,  .dfreq = 0.10, .prop = 0.02},
			[1]  = {.freq_prop = 10.0, .dfreq = 0.20, .prop = 0.05},
			[2]  = {.freq_prop = 15.0, .dfreq = 0.20, .prop = 0.04},
			[3]  = {.freq_prop = 20.0, .dfreq = 0.20, .prop = 0.02},
			[4]  = {.freq_prop = 40.0, .dfreq = 0.50, .prop = 0.02},
			[5]  = {.freq_prop = 40.0, .dfreq = 0.0,  .prop = 0.0},
		},
		[EN61000_4_13_CLASS_B] =
		{
			[0]  = {.freq_prop = 2.0,  .dfreq = 0.10, .prop = 0.02},
			[1]  = {.freq_prop = 10.0, .dfreq = 0.20, .prop = 0.05},
			[2]  = {.freq_prop = 15.0, .dfreq = 0.20, .prop = 0.04},
			[3]  = {.freq_prop = 20.0, .dfreq = 0.20, .prop = 0.02},
			[4]  = {.freq_prop = 40.0, .dfreq = 0.50, .prop = 0.02},
			[5]  = {.freq_prop = 40.0, .dfreq = 0.0,  .prop = 0.0},
		},
		[EN61000_4_13_CLASS_C] =
		{
			[0]  = {.freq_prop = 2.0,  .dfreq = 0.10, .prop = 0.04},
			[1]  = {.freq_prop = 10.0, .dfreq = 0.20, .prop = 0.09},
			[2]  = {.freq_prop = 15.0, .dfreq = 0.20, .prop = 0.045},
			[3]  = {.freq_prop = 20.0, .dfreq = 0.20, .prop = 0.035},
			[4]  = {.freq_prop = 40.0, .dfreq = 0.50, .prop = 0.02},
			[5]  = {.freq_prop = 40.0, .dfreq = 0.0,  .prop = 0.0},
		},
};

static void get_wave(double amplitude, double frequency, double t, double phase, double* V){
	double inter_phase_shift = 2 * PI / 3; //120

	V[0] = amplitude * sin(2 * PI * frequency * t + phase);
	V[1] = amplitude * sin(2 * PI * frequency * t - inter_phase_shift + phase);
	V[2] = amplitude * sin(2 * PI * frequency * t + inter_phase_shift + phase);
}

static void get_base_wave(void* _obj, double t, double* V){
	GET_OBJ;
	get_wave(obj->derived_params.amplitude, obj->params.base_freq, t, V, 0.0);
}

static void add_waves(double* V0, double* V1, double* Vout){
	for(int i=0; i<NUM_PH; i++){
		Vout[i] = V0[i] + V1[i];
	}
}

static void flat_curve_step(void* _obj, double t, double* V, double fade){
	GET_OBJ;
	double clip_prop = FLAT_TOP_PARAMS[obj->params.class].clip_prop;
	double amplitude = obj->derived_params.amplitude;

	get_base_wave(_obj, t, V);
	for(int i=0; i<NUM_PH; i++)
		V[i] = SAT_HI_LO(V[i], clip_prop*amplitude, -clip_prop*amplitude);

	if(obj->dwell_acc > obj->params->dwell){
		obj->finished = true;
		return;
	}

}

static void over_swing_step(void* _obj, double t, double* V, double fade){
	GET_OBJ;
	double third_prop = OVER_SWING_PARAMS[obj->params.class].third_prop;
	double fifth_prop = OVER_SWING_PARAMS[obj->params.class].fifth_prop;
	double amplitude = obj->derived_params.amplitude;

	get_base_wave(_obj, t, V);
	double V_harm[NUM_PH];
	get_wave(amplitude*third_prop*fade, obj->params.base_freq*3.0, t, V, _180_DEG);
	add_waves(V, V_harm, V);
	get_wave(amplitude*fifth_prop*fade, obj->params.base_freq*5.0, t, V, 0);
	add_waves(V, V_harm, V);

	if(obj->dwell_acc > obj->params->dwell){
		obj->finished = true;
	}
}

static void individual_harmonic_step(void* _obj, double t, double* V, double fade){
	GET_OBJ;
	static bool anti_phase = false;

	double harm_prop = INDIVIDUAL_HARMONIC_PARAMS[obj->params.class][obj->istep];
	double amplitude = obj->derived_params.amplitude;

	get_base_wave(_obj, t, V);
	double V_harm[NUM_PH];
	get_wave(amplitude*harm_prop*fade, obj->params.base_freq*obj->istep, t, V, (anti_phase)?(_180_DEG):(0.0));
	add_waves(V, V_harm, V);

	if(obj->dwell_acc > obj->params->dwell){
		obj->dwell_acc = 0.0;
		if(obj->istep >= 40){
			obj->finished = true;
			return;
		} else if((obj->istep > 9) || (anti_phase)){
			obj->istep++;
			anti_phase = false;
		} else {
			anti_phase = true;
		}
	}
}

static void sweep_step(void* _obj, double t, double* V, const SweepEntry* sweep_table, double fade){
	GET_OBJ;

	double harm_prop = sweep_table[obj->istep].prop;
	double amplitude = obj->derived_params.amplitude;

	get_base_wave(_obj, t, V);
	double V_harm[NUM_PH];
	get_wave(amplitude*harm_prop*fade, obj->dstep, t, V, 0.0);
	add_waves(V, V_harm, V);

	if(obj->dwell_acc > obj->params->dwell){
		obj->dwell_acc = 0.0;
		obj->dstep += sweep_table[obj->istep].dfreq*obj->params.base_freq;
		if(sweep_table[obj->istep].dfreq == 0.0){
			obj->finished = true;
			return;
		} else if(obj->dstep > sweep_table[obj->istep].freq_prop*obj->params.base_freq)
			obj->istep++;
	}
}

static void frequency_sweep_step(void* _obj, double t, double* V, double fade){
	GET_OBJ;
	sweep_step(_obj, t, V, FREQUENCY_SWEEP_PARAMS[obj->params.class], fade);
}

static void interharmonic_step(void* _obj, double t, double* V, double fade){
	GET_OBJ;
	sweep_step(_obj, t, V, INTER_HARMONIC_PARAMS[obj->params.class], fade);
}

static void assign_step_callback(void* _obj, EN61000_4_13_Test test){
	GET_OBJ;
	switch(test){
	case EN61000_4_13_FLAT_TOP:
		obj->step_callback = flat_curve_step;
		break;
	case EN61000_4_13_OVER_SWING:
		obj->step_callback = over_swing_step;
		break;
	case EN61000_4_13_FREQUENCY_SWEEP:
		obj->step_callback = frequency_sweep_step;
		break;
	case EN61000_4_13_INDIVIDUAL_HARMONICS:
		obj->step_callback = individual_harmonic_step;
		break;
	case EN61000_4_13_INTERHARMONICS:
		obj->step_callback = interharmonic_step;
		break;
	case EN61000_4_13_MEISTER_CURVE:
		obj->step_callback = NULL;
		break;
	}
}

void* EN61000_4_13_create(){
	return malloc(sizeof(TestInstance));
}

void EN61000_4_13_start_test(void* _obj, EN61000_4_13_Test test, TestParams* params){
	GET_OBJ;
	obj->params = *params;
	obj->derived_params.amplitude = params->rms*RMS_TO_AMP;
	obj->derived_params.period = 1/params->base_freq;
	obj->fade_in_finished = false;
	obj->fade_out_finished = false;
	obj->fade_acc = 0.0;
	obj->finished = false;
	obj->istep = 0;
	obj->dstep = obj->params.base_freq;
	assign_step_callback(_obj, test);
}

void EN61000_4_13_step(void* _obj, double dt, double* V){
	GET_OBJ;
	EN61000_4_13_step_t(_obj, obj->t+dt, V);

}

void fade_in(void* _obj, double t, double* fade){
	GET_OBJ;
	obj->fade_acc += t-obj->prev_t;
	if(obj->fade_acc >= obj->derived_params.period){
		obj->fade_in_finished = true;
		obj->fade_acc = 0.0;
	} else
		fade = obj->fade_acc/obj->derived_params.period;
}

void fade_out(void* _obj, double t, double* fade){
	GET_OBJ;
	obj->fade_acc += t-obj->prev_t;
	if(obj->fade_acc >= obj->derived_params.period)
		obj->fade_out_finished = true;
	else
		fade = 1- (obj->fade_acc/obj->derived_params.period);
}

void EN61000_4_13_step_t(void* _obj, double t, double* V){
	GET_OBJ;
	obj->t = t;

	double fade = 1.0;

	if(!obj->fade_in_finished){
		fade_in(_obj, t, &fade);
	} else if((obj->finished)&&(!obj->fade_out_finished)){
		fade_out(_obj, t, &fade);
	} else if(obj->fade_out_finished)
		fade = 0.0;

	obj->dwell_acc += t-obj->prev_t;
	obj->step_callback(_obj, t, V, fade);
	obj->prev_t = t;
}

bool EN61000_4_13_finished(void* _obj){
	GET_OBJ;
	return obj->fade_out_finished;
}
