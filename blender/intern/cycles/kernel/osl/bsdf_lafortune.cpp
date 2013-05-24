/*
 * Adapted from Open Shading Language with this license:
 *
 * Copyright (c) 2009-2010 Sony Pictures Imageworks Inc., et al.
 * All Rights Reserved.
 *
 * Modifications Copyright 2011, Blender Foundation.
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
  * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <OpenImageIO/fmath.h>

#include <OSL/genclosure.h>

#include "osl_closures.h"

#include "kernel_types.h"
#include "kernel_montecarlo.h"
#include "closure/bsdf_lafortune.h"

CCL_NAMESPACE_BEGIN

using namespace OSL;

class LafortuneClosure : public CBSDFClosure {
public:
	LafortuneClosure() : CBSDFClosure(LABEL_GLOSSY) {}
	float colormatrix[9];
	float coeff[27];

	size_t memsize() const { return sizeof(*this); }
	const char *name() const { return "lafortune"; }

	void setup()
	{
		sc.prim = this;
		m_shaderdata_flag = bsdf_lafortune_setup(&sc);
	}

	bool mergeable(const ClosurePrimitive *other) const
	{
		return false;
	}

	void blur(float roughness)
	{
		bsdf_lafortune_blur(&sc, roughness);
	}

	void print_on(std::ostream &out) const
	{
		out << name() << " ((" << sc.N[0] << ", " << sc.N[1] << ", " << sc.N[2] << "))";
	}

	float3 eval_reflect(const float3 &omega_out, const float3 &omega_in, float& pdf) const
	{
		return bsdf_lafortune_eval_reflect(&sc, colormatrix, coeff, omega_out, omega_in, &pdf);
	}

	float3 eval_transmit(const float3 &omega_out, const float3 &omega_in, float& pdf) const
	{
		return bsdf_lafortune_eval_transmit(&sc, omega_out, omega_in, &pdf);
	}

	int sample(const float3 &Ng,
	           const float3 &omega_out, const float3 &domega_out_dx, const float3 &domega_out_dy,
	           float randu, float randv,
	           float3 &omega_in, float3 &domega_in_dx, float3 &domega_in_dy,
	           float &pdf, float3 &eval) const
	{
		return bsdf_lafortune_sample(&sc, colormatrix, coeff, Ng, omega_out, domega_out_dx, domega_out_dy,
			randu, randv, &eval, &omega_in, &domega_in_dx, &domega_in_dy, &pdf);
	}
};

ClosureParam *closure_bsdf_lafortune_params()
{
	static ClosureParam params[] = {
		CLOSURE_FLOAT3_PARAM(LafortuneClosure, sc.N),
		CLOSURE_COLOR_ARRAY_PARAM(LafortuneClosure, colormatrix, 9),
        CLOSURE_FLOAT_ARRAY_PARAM(LafortuneClosure, coeff, 27),
		CLOSURE_STRING_KEYPARAM("label"),
	    CLOSURE_FINISH_PARAM(LafortuneClosure)
	};
	return params;
}

CLOSURE_PREPARE(closure_bsdf_lafortune_prepare, LafortuneClosure)

CCL_NAMESPACE_END

