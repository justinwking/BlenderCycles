/*
 * Adapted from Open Shading Language with this license:
 *
 * Copyright (c) 2009-2010 Sony Pictures Imageworks Inc., et al.
 * All Rights Reserved.
 *
 * Modifications Copyright 2012, Blender Foundation.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 * * Neither the name of Sony Pictures Imageworks nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include "stdlib.h"
#include "math.h"
#include <stdio.h>

#ifndef __BSDF_TEST_CLOSURE_H__
#define __BSDF_TEST_CLOSURE_H__

#define BRDF_SAMPLING_RES_THETA_H       90
#define BRDF_SAMPLING_RES_THETA_D       90
#define BRDF_SAMPLING_RES_PHI_D         360

#define RED_SCALE (1.0/1500.0)
#define GREEN_SCALE (1.15/1500.0)
#define BLUE_SCALE (1.66/1500.0)
//#define M_PI	3.1415926535897932384626433832795
CCL_NAMESPACE_BEGIN

bool readagain = true;
const char *filename = "/Users/niverik2k/Downloads/polyurethane-foam.binary";
double* brdflookup;

// cross product of two vectors
void cross_product (double* v1, double* v2, double* out)
{
	out[0] = v1[1]*v2[2] - v1[2]*v2[1];
	out[1] = v1[2]*v2[0] - v1[0]*v2[2];
	out[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

// normalize vector
void normalize(double* v)
{
	// normalize
	double len = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	v[0] = v[0] / len;
	v[1] = v[1] / len;
	v[2] = v[2] / len;
}

// rotate vector along one axis
void rotate_vector(double* vector, double* axis, double angle, double* out)
{
	double temp;
	double cross[3];
	double cos_ang = cos(angle);
	double sin_ang = sin(angle);
    
	out[0] = vector[0] * cos_ang;
	out[1] = vector[1] * cos_ang;
	out[2] = vector[2] * cos_ang;
    
	temp = axis[0]*vector[0]+axis[1]*vector[1]+axis[2]*vector[2];
	temp = temp*(1.0-cos_ang);
    
	out[0] += axis[0] * temp;
	out[1] += axis[1] * temp;
	out[2] += axis[2] * temp;
    
	cross_product (axis,vector,cross);
	
	out[0] += cross[0] * sin_ang;
	out[1] += cross[1] * sin_ang;
	out[2] += cross[2] * sin_ang;
}


// convert standard coordinates to half vector/difference vector coordinates
void std_coords_to_half_diff_coords(double theta_in, double fi_in, double theta_out, double fi_out,
                                    double& theta_half,double& fi_half,double& theta_diff,double& fi_diff )
{
    
	// compute in vector
	double in_vec_z = cos(theta_in);
	double proj_in_vec = sin(theta_in);
	double in_vec_x = proj_in_vec*cos(fi_in);
	double in_vec_y = proj_in_vec*sin(fi_in);
	double in[3]= {in_vec_x,in_vec_y,in_vec_z};
	normalize(in);
    
    
	// compute out vector
	double out_vec_z = cos(theta_out);
	double proj_out_vec = sin(theta_out);
	double out_vec_x = proj_out_vec*cos(fi_out);
	double out_vec_y = proj_out_vec*sin(fi_out);
	double out[3]= {out_vec_x,out_vec_y,out_vec_z};
	normalize(out);
    
    
	// compute halfway vector
	double half_x = (in_vec_x + out_vec_x)/2.0f;
	double half_y = (in_vec_y + out_vec_y)/2.0f;
	double half_z = (in_vec_z + out_vec_z)/2.0f;
	double half[3] = {half_x,half_y,half_z};
	normalize(half);
    
	// compute  theta_half, fi_half
	theta_half = acos(half[2]);
	fi_half = atan2(half[1], half[0]);
    
    
	double bi_normal[3] = {0.0, 1.0, 0.0};
	double normal[3] = { 0.0, 0.0, 1.0 };
	double temp[3];
	double diff[3];
    
	// compute diff vector
	rotate_vector(in, normal , -fi_half, temp);
	rotate_vector(temp, bi_normal, -theta_half, diff);
	
	// compute  theta_diff, fi_diff
	theta_diff = acos(diff[2]);
	fi_diff = atan2(diff[1], diff[0]);
    
}

// Lookup theta_half index
// This is a non-linear mapping!
// In:  [0 .. pi/2]
// Out: [0 .. 89]
inline int theta_half_index(double theta_half)
{
	if (theta_half <= 0.0)
		return 0;
	double theta_half_deg = ((theta_half / (M_PI/2.0))*BRDF_SAMPLING_RES_THETA_H);
	double temp = theta_half_deg*BRDF_SAMPLING_RES_THETA_H;
	temp = sqrt(temp);
	int ret_val = (int)temp;
	if (ret_val < 0) ret_val = 0;
	if (ret_val >= BRDF_SAMPLING_RES_THETA_H)
		ret_val = BRDF_SAMPLING_RES_THETA_H-1;
	return ret_val;
}

// Lookup theta_diff index
// In:  [0 .. pi/2]
// Out: [0 .. 89]
inline int theta_diff_index(double theta_diff)
{
	int tmp = int(theta_diff / (M_PI * 0.5) * BRDF_SAMPLING_RES_THETA_D);
	if (tmp < 0)
		return 0;
	else if (tmp < BRDF_SAMPLING_RES_THETA_D - 1)
		return tmp;
	else
		return BRDF_SAMPLING_RES_THETA_D - 1;
}


// Lookup phi_diff index
inline int phi_diff_index(double phi_diff)
{
	// Because of reciprocity, the BRDF is unchanged under
	// phi_diff -> phi_diff + M_PI
	if (phi_diff < 0.0)
		phi_diff += M_PI;
    
	// In: phi_diff in [0 .. pi]
	// Out: tmp in [0 .. 179]
	int tmp = int(phi_diff / M_PI * BRDF_SAMPLING_RES_PHI_D / 2);
	if (tmp < 0)
		return 0;
	else if (tmp < BRDF_SAMPLING_RES_PHI_D / 2 - 1)
		return tmp;
	else
		return BRDF_SAMPLING_RES_PHI_D / 2 - 1;
}

// Given a pair of incoming/outgoing angles, look up the BRDF.
void lookup_brdf_val(double* brdflookup, double theta_in, double fi_in,
                     double theta_out, double fi_out,
                     double& red_val,double& green_val,double& blue_val)
{
	// Convert to halfangle / difference angle coordinates
	double theta_half, fi_half, theta_diff, fi_diff;
	
	std_coords_to_half_diff_coords(theta_in, fi_in, theta_out, fi_out,
                                   theta_half, fi_half, theta_diff, fi_diff);
    
    
	// Find index.
	// Note that phi_half is ignored, since isotropic BRDFs are assumed
	int ind = phi_diff_index(fi_diff) +
    theta_diff_index(theta_diff) * BRDF_SAMPLING_RES_PHI_D / 2 +
    theta_half_index(theta_half) * BRDF_SAMPLING_RES_PHI_D / 2 *
    BRDF_SAMPLING_RES_THETA_D;
    
	red_val = brdflookup[ind] * RED_SCALE;
	green_val = brdflookup[ind + BRDF_SAMPLING_RES_THETA_H*BRDF_SAMPLING_RES_THETA_D*BRDF_SAMPLING_RES_PHI_D/2] * GREEN_SCALE;
	blue_val = brdflookup[ind + BRDF_SAMPLING_RES_THETA_H*BRDF_SAMPLING_RES_THETA_D*BRDF_SAMPLING_RES_PHI_D] * BLUE_SCALE;
    
	
	if (red_val < 0.0 || green_val < 0.0 || blue_val < 0.0)
		fprintf(stderr, "Below horizon.\n");
    
}



// Read BRDF data
bool read_brdf(const char *filename, double* &brdflookup)
{
	FILE *f = fopen(filename, "rb");
	if (!f)
        return false;
    
	int dims[3];
	fread(dims, sizeof(int), 3, f);
	int n = dims[0] * dims[1] * dims[2];
	if (n != BRDF_SAMPLING_RES_THETA_H *
        BRDF_SAMPLING_RES_THETA_D *
        BRDF_SAMPLING_RES_PHI_D / 2)
	{
		fprintf(stderr, "Dimensions don't match\n");
		fclose(f);
		return false;
	}
    
	brdflookup = (double*) malloc (sizeof(double)*3*n);
	fread(brdflookup, sizeof(double), 3*n, f);
    
	fclose(f);
	return true;

}

float2 vectortospherical(float3 vector)
{
    float x = vector[0];
    float y = vector[1];
    float z = vector[2];
    float row = sqrt((x * x) + (y * y) + (z * z));
    float S = sqrt((x * x) + (y * y));
    float theta = acosf(z/row);
    float phi;
    if (x <= 0)
    {
        phi = asinf(y/S);
    }
    else
    {
        phi = M_PI - asinf(y/S);
    }
    return make_float2(theta,phi);
}

__device float3 dolookup(float3 omega_in,float3 I)
{
    float2 thetaphi_in = vectortospherical(omega_in);
    float2 thetaphi_out = vectortospherical(I);
    double theta_in = double(thetaphi_in[0]);
    double phi_in = double(thetaphi_in[1]);
    double theta_out = double(thetaphi_out[0]);
    double phi_out = double(thetaphi_out[1]);
    printf("%f,%f,%f,%f\n",(float)thetaphi_in[0],(float)thetaphi_out[0],(float)thetaphi_in[1],(float)thetaphi_out[1]);
    printf("%d*\n",(double*)brdflookup);
    
	// return color value based on lookup
    double red,green,blue;
    lookup_brdf_val(brdflookup, theta_in, phi_in, theta_out, phi_out, red, green, blue);
    /*green = 1.0;
    red = 1.0;
    blue = 1.0;
    free(&theta_in);
    free(&phi_in);
    free(&theta_out);
    free(&phi_out);*/
    printf("%f,%f,%f\n",(float)red,(float)green,(float)blue);
    return make_float3((float)red, (float)green, (float)blue);
}

__device int bsdf_test_closure_setup(ShaderClosure *sc)
{

    if(!read_brdf(filename, brdflookup))
    {
         fprintf(stderr, "Error reading %s\n", filename);
        exit(1);
     }
	sc->type = CLOSURE_BSDF_TEST_CLOSURE_ID;
	return SD_BSDF | SD_BSDF_HAS_EVAL;
}

__device void bsdf_test_closure_blur(ShaderClosure *sc, float roughness)
{
}

__device float3 bsdf_test_closure_eval_reflect(const ShaderClosure *sc, const float3 I, const float3 omega_in, float *pdf)
{
	float3 N = sc->N;
    if(dot(N, omega_in) > 0.0f)
    {
        float cos_pi = fmaxf(dot(N, omega_in), 0.0f);
        *pdf = cos_pi * M_1_PI_F;
        return dolookup(omega_in,I);
    }
	else
		return make_float3(0.0f,0.0f,0.0f);
}

__device float3 bsdf_test_closure_eval_transmit(const ShaderClosure *sc, const float3 I, const float3 omega_in, float *pdf)
{
	return make_float3(0.0f, 0.0f, 0.0f);
}

__device int bsdf_test_closure_sample(const ShaderClosure *sc, float3 Ng, float3 I, float3 dIdx, float3 dIdy, float randu, float randv, float3 *eval, float3 *omega_in, float3 *domega_in_dx, float3 *domega_in_dy, float *pdf)
{
	float3 N = sc->N;

	// distribution over the hemisphere
	sample_cos_hemisphere(N, randu, randv, omega_in, pdf);

	if(dot(Ng, *omega_in) > 0.0f) {
		//*eval = bsdf_test_closure_get_color(sc, colors, *pdf * M_PI_F) * M_1_PI_F;
        *eval = dolookup(*omega_in,I);
#ifdef __RAY_DIFFERENTIALS__
		*domega_in_dx = (2 * dot(N, dIdx)) * N - dIdx;
		*domega_in_dy = (2 * dot(N, dIdy)) * N - dIdy;
#endif
	}
	else
		*pdf = 0.0f;
	
	return LABEL_REFLECT|LABEL_DIFFUSE;
}

CCL_NAMESPACE_END

#endif /* __BSDF_TEST_CLOSURE_H__ */
