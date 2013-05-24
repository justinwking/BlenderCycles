/*
 * Based on an implement for RSL of the BRDF representation of Lafortune et al.
 *
 * DESCRIPTION:
 *	Implement the BRDF representation described in "Non-Linear
 *	Approximation of Reflectance Functions", by Eric P. F. Lafortune,
 *	Sing-Choong Foo, Kenneth E. Torrance, and Donald P. Greenberg, in
 *	SIGGRAPH 97 Proceedings.
 *
 * PARAMETERS:
 *   coeff:	Coefficients for each lobe. There will be
 *		3 coefficients per lobe per wavelength.
 *		for each lobe of Nlobes:
 *		  for each channel of N_WAVES:
 *		    cxy, cz, n
 *		where cxy, cz are the directional and scale components
 *		n is the exponent for the cosine
 *
 *   Cs:	Diffuse color of surface. This provides the diffuse
 *		term that is an integral part of the model. This
 *		will be the diffuse component from the Lafortune
 *		model multiplied by pi, as the RenderMan definition
 *		of reflectance is different.
 * colormatrix:	3x3 color transformation matrix, if needed.
 *		This is here because measured data may be in the
 *		RGB space of some sensor (e.g. a digital camera)
 *		and may need correction for display. Since the three
 *		channels are represented independently (i.e. each
 *		channel will have different directionality), this
 *		correction must be done on a per-pixel basis.
 *		Since the correction may be different for each
 *		shader, we must do it here. Default is the identity
 *		matrix.
 *
 * ANTIALIASING: should antialias itself fairly well
 *
 * AUTHOR: written by Stephen H. Westin, Cornell University Program of
 *	Computer Graphics
 *
 * HISTORY:
 *
 * initial version 22 March 1999 S. H. Westin
 *
 * 22 March 1999 S. H. Westin
 *		Flipped sign on Cx to agree with description in SIGGRAPH paper
 * 25 October 2000 S. H. Westin
 *		Multiply lobes by incident cosine.
 * 27 October 2000 S. H. Westin
 *		Multiply lobes by incident cosine. Before, I was trying to
 *		use diffuse(), which doesn't work.
 * 24 May 2013 Justin W. King
        Ported to Open Shading Language
 */

/* Number of coefficients per lobe: 3 for an isotropic surface */
#define LOBESIZE 3
/* Number of wavelengths (wired in) */
#define N_WAVES 3
/* Number of lobes (wired in) */
#define Nlobes 3

#define COEFFLEN (LOBESIZE*N_WAVES*Nlobes)


#ifndef __BSDF_LAFORTUNE_H__
#define __BSDF_LAFORTUNE_H__

CCL_NAMESPACE_BEGIN
 
__device float3 bsdf_lafortune_colormatrix(const ShaderClosure *sc, const float colormatrix[9], const float3 ColorIn)
{
	float ColorInR = ColorIn[0];
	float ColorInG = ColorIn[1];
    float ColorInB = ColorIn[2];
    
    float fr = colormatrix[0] * ColorInR + colormatrix[1] * ColorInG  + colormatrix[2] * ColorInB;
    float fg = colormatrix[3] * ColorInR + colormatrix[4] * ColorInG  + colormatrix[5] * ColorInB;
    float fb = colormatrix[6] * ColorInR + colormatrix[7] * ColorInG  + colormatrix[8] * ColorInB;
    
    return make_float3(fr,fg,fb);

}

__device float3 bsdf_lafortune_get_color(const ShaderClosure *sc, const float coeff[27], const float3 I, const float3 omega_in)
{
    float3 local_z;
    float3 local_x, local_y;  /*Unit vector in "u" and "v" directions */
    float3 V;                   /* Normalized eye vector */
    float3 Ln;          /* Normalized vector to light */
    float x, z, f;      /*subterms */
    float fr, fg, fb;   /* RGB components of the non-Lambertian term */
    int basepointer, j; /* loop counters */
    
    /* Get unit vector in "u" parameter direction */
    local_x = normalize(sc->T);
    local_z = normalize(sc->N);
    V = normalize(I);
    
    /* Get y component of the local coordinate system, with the cross product of local_z and local_x. */
    local_y = cross(local_z,local_x);
    
    Ln = normalize(omega_in);
    
    // We will compute the following terms
    // x = x_in * x_view + y_in * y_view
    // z = z_in * z_view
    
    x = dot(local_x,V) * dot(local_x,Ln) + dot(local_y,V) * dot(local_y,Ln);
    z = -(dot(local_z,V) * dot(local_z,Ln));
    
    /*Coefficient structure:
     for each lobe of Nlobes:
     for each channel of N_WAVES:
     cxy, cz, n (where cxy, cz are the directional and scale components
     n is the exponent for the cosine
     */
    
    for ( basepointer=0; basepointer<COEFFLEN; basepointer += LOBESIZE * N_WAVES ) {
        float rexponent = coeff[basepointer + 2];
        float gexponent = coeff[basepointer + LOBESIZE + 2];
        float bexponent = coeff[basepointer + (2 * LOBESIZE) + 2];
        
        fr = 0.0f;
        fg = 0.0f;
        fb = 0.0f;
        
        f = -x * coeff[basepointer] + z * coeff[basepointer + 1];
        
        if ( f > 0.001 * rexponent )
        {
            fr = pow ( f, rexponent);
        }
        
        f = -x * coeff[basepointer+LOBESIZE] + z * coeff[basepointer+LOBESIZE+1];
        
        if ( f > 0.001 * gexponent )
        {
            fg = pow ( f, gexponent);
        }
        
        f = -x * coeff[basepointer+2*LOBESIZE] + z * coeff[basepointer+2*LOBESIZE+1];
        
        if ( f > 0.001 * bexponent )
        {
            fb = pow ( f, bexponent);
        }
    }
        return make_float3(fr,fg,fb)  *  dot(local_z,Ln);
    
}

__device int bsdf_lafortune_setup(ShaderClosure *sc)
{
	sc->type = CLOSURE_BSDF_LAFORTUNE_ID;
	return SD_BSDF | SD_BSDF_HAS_EVAL;
}

__device void bsdf_lafortune_blur(ShaderClosure *sc, float roughness)
{
}

__device float3 bsdf_lafortune_eval_reflect(const ShaderClosure *sc, const float coeff[27], const float colormatrix[8], const float3 I, const float3 omega_in, float *pdf)
{
	float3 N = sc->N;
    float3 coeffcolor = bsdf_lafortune_get_color(sc, coeff, I, omega_in);
    
    float cos_pi = fmaxf(dot(N, omega_in), 0.0f);
	*pdf = cos_pi * M_1_PI_F;
    
    return bsdf_lafortune_colormatrix(sc, colormatrix, coeffcolor);
    
    /*
	float cos_pi = fmaxf(dot(N, omega_in), 0.0f);
	*pdf = cos_pi * M_1_PI_F;
	return bsdf_diffuse_ramp_get_color(sc, colors, cos_pi) * M_1_PI_F;
     */
}

__device float3 bsdf_lafortune_eval_transmit(const ShaderClosure *sc, const float3 I, const float3 omega_in, float *pdf)
{
	return make_float3(0.0f, 0.0f, 0.0f);
}

__device int bsdf_lafortune_sample(const ShaderClosure *sc, const float coeff[27], const float colormatrix[8], float3 Ng, float3 I, float3 dIdx, float3 dIdy, float randu, float randv, float3 *eval, float3 *omega_in, float3 *domega_in_dx, float3 *domega_in_dy, float *pdf)
{
	float3 N = sc->N;

	// distribution over the hemisphere
	sample_cos_hemisphere(N, randu, randv, omega_in, pdf);

	if(dot(Ng, *omega_in) > 0.0f) {
		*eval = bsdf_lafortune_get_color(sc, coeff, I, *omega_in);
#ifdef __RAY_DIFFERENTIALS__
		*domega_in_dx = (2 * dot(N, dIdx)) * N - dIdx;
		*domega_in_dy = (2 * dot(N, dIdy)) * N - dIdy;
#endif
	}
	else
		*pdf = 0.0f;
	
	return LABEL_REFLECT|LABEL_GLOSSY;
}

CCL_NAMESPACE_END

#endif /* __BSDF_LAFORTUNE_H__ */
