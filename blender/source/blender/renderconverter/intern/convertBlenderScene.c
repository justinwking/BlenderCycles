/**
 * $Id$
 *
 * ***** BEGIN GPL/BL DUAL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version. The Blender
 * Foundation also sells licenses for use in proprietary software under
 * the Blender License.  See http://www.blender.org/BL/ for information
 * about this.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 * The Original Code is Copyright (C) 2001-2002 by NaN Holding BV.
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): none yet.
 *
 * ***** END GPL/BL DUAL LICENSE BLOCK *****
 * Interface to transform the Blender scene into renderable data.
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>          /* for INT_MAX                 */

#include "blendef.h"
#include "MTC_matrixops.h"

#include "MEM_guardedalloc.h"

#include "BLI_arithb.h"
#include "BLI_blenlib.h"
#include "BLI_rand.h"
#include "BLI_memarena.h"
#include "BLI_ghash.h"

#include "DNA_armature_types.h"
#include "DNA_camera_types.h"
#include "DNA_material_types.h"
#include "DNA_curve_types.h"
#include "DNA_effect_types.h"
#include "DNA_group_types.h"
#include "DNA_lamp_types.h"
#include "DNA_lattice_types.h"
#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_meta_types.h"
#include "DNA_object_types.h"
#include "DNA_object_force.h"
#include "DNA_scene_types.h"
#include "DNA_texture_types.h"
#include "DNA_view3d_types.h"

#include "BKE_anim.h"
#include "BKE_armature.h"
#include "BKE_action.h"
#include "BKE_curve.h"
#include "BKE_constraint.h"
#include "BKE_displist.h"
#include "BKE_deform.h"
#include "BKE_DerivedMesh.h"
#include "BKE_effect.h"
#include "BKE_global.h"
#include "BKE_group.h"
#include "BKE_key.h"
#include "BKE_ipo.h"
#include "BKE_lattice.h"
#include "BKE_material.h"
#include "BKE_main.h"
#include "BKE_mball.h"
#include "BKE_mesh.h"
#include "BKE_object.h"
#include "BKE_scene.h"
#include "BKE_subsurf.h"
#include "BKE_texture.h"
#include "BKE_utildefines.h"
#include "BKE_world.h"

#include "render_types.h"
#include "rendercore.h"
#include "RE_renderconverter.h"

#include "BIF_space.h"
#include "BIF_screen.h"
#include "BIF_editkey.h"

#include "BSE_sequence.h"

#include "nla.h"

#include "BPY_extern.h"

#include "butspace.h"

#include "radio.h"
#include "YafRay_Api.h"

extern Render R;
extern void error (char *fmt, ...);  /* defined in BIF_toolbox.h, but we dont need to include the rest */ 

/* yafray: Identity transform 'hack' removed, exporter now transforms vertices back to world.
 * Same is true for lamp coords & vec.
 * Duplicated data objects & dupliframe/duplivert objects are only stored once,
 * only the matrix is stored for all others, in yafray these objects are instances of the original.
 * The main changes are in RE_rotateBlenderScene().
 */

/* ------------------------------------------------------------------------- */
/* Local functions                                                           */
/* ------------------------------------------------------------------------- */
static Material *give_render_material(Object *ob, int nr);
static short test_for_displace(Object *ob);
static void do_displacement(Object *ob, int startface, int numface, int startvert, int numvert );

/* ------------------------------------------------------------------------- */
/* tool functions/defines for ad hoc simplification and possible future 
   cleanup      */
/* ------------------------------------------------------------------------- */

#define UVTOINDEX(u,v) (startvlak + (u) * sizev + (v))
/*

NOTE THAT U/V COORDINATES ARE SOMETIMES SWAPPED !!
	
^	()----p4----p3----()
|	|     |     |     |
u	|     |  F1 |  F2 |
	|     |     |     |
	()----p1----p2----()
	       v ->
*/

/* ------------------------------------------------------------------------- */

/* Stuff for stars. This sits here because it uses gl-things. Part of
   this code may move down to the converter.  */
/* ------------------------------------------------------------------------- */
/* this is a bad beast, since it is misused by the 3d view drawing as well. */

/* more star stuff, here used to be a cliptest, removed for envmap render or panorama... */
static HaloRen *initstar(float *vec, float hasize)
{
	HaloRen *har;
	float hoco[4];
	
	RE_projectverto(vec, hoco);
	
	har= RE_findOrAddHalo(R.tothalo++);
	
	/* projectvert is done in function zbufvlaggen again, because of parts */
	VECCOPY(har->co, vec);
	har->hasize= hasize;
	
	har->zd= 0.0;
	
	return har;
}

extern unsigned char hash[512];

/* there must be a 'fixed' amount of stars generated between
 *         near and far
 * all stars must by preference lie on the far and solely
 *        differ in clarity/color
 */

void RE_make_stars(void (*initfunc)(void),
				   void (*vertexfunc)(float*),
				   void (*termfunc)(void))
{
	HaloRen *har;
	double dblrand, hlfrand;
	float vec[4], fx, fy, fz;
	float fac, starmindist, clipend;
	float mat[4][4], stargrid, maxrand, maxjit, force, alpha;
/* 	float loc_far_var, loc_near_var; */
	int x, y, z, sx, sy, sz, ex, ey, ez, done = 0;
	Camera * camera;

	if(initfunc) R.wrld= *(G.scene->world);

	stargrid = R.wrld.stardist;		/* distance between stars */
	maxrand = 2.0;						/* amount a star can be shifted (in grid units) */
	maxjit = (R.wrld.starcolnoise);			/* amount a color is being shifted */

/* 	loc_far_var = R.far; */
/* 	loc_near_var = R.near; */


	/* size of stars */
	force = ( R.wrld.starsize );

	/* minimal free space (starting at camera) */
	starmindist= R.wrld.starmindist;

	if (stargrid <= 0.10) return;

	if (!initfunc) R.flag |= R_HALO;
	else stargrid *= 1.0;				/* then it draws fewer */


	MTC_Mat4Invert(mat, R.viewmat);

	/* BOUNDING BOX CALCULATION
	 * bbox goes from z = loc_near_var | loc_far_var,
	 * x = -z | +z,
	 * y = -z | +z
	 */

	camera = G.scene->camera->data;
	clipend = camera->clipend;

	/* convert to grid coordinates */

	sx = ((mat[3][0] - clipend) / stargrid) - maxrand;
	sy = ((mat[3][1] - clipend) / stargrid) - maxrand;
	sz = ((mat[3][2] - clipend) / stargrid) - maxrand;

	ex = ((mat[3][0] + clipend) / stargrid) + maxrand;
	ey = ((mat[3][1] + clipend) / stargrid) + maxrand;
	ez = ((mat[3][2] + clipend) / stargrid) + maxrand;

	dblrand = maxrand * stargrid;
	hlfrand = 2.0 * dblrand;

	if (initfunc) {
		initfunc();	
	}

	for (x = sx, fx = sx * stargrid; x <= ex; x++, fx += stargrid) {
		for (y = sy, fy = sy * stargrid; y <= ey ; y++, fy += stargrid) {
			for (z = sz, fz = sz * stargrid; z <= ez; z++, fz += stargrid) {

				BLI_srand((hash[z & 0xff] << 24) + (hash[y & 0xff] << 16) + (hash[x & 0xff] << 8));
				vec[0] = fx + (hlfrand * BLI_drand()) - dblrand;
				vec[1] = fy + (hlfrand * BLI_drand()) - dblrand;
				vec[2] = fz + (hlfrand * BLI_drand()) - dblrand;
				vec[3] = 1.0;

				if (vertexfunc) {
					if(done & 1) vertexfunc(vec);
					done++;
				}
				else {
					MTC_Mat4MulVecfl(R.viewmat, vec);

					/* in vec are global coordinates
					 * calculate distance to camera
					 * and using that, define the alpha
					 */

					{
						float tx, ty, tz;

						tx = vec[0];
						ty = vec[1];
						tz = vec[2];

						alpha = sqrt(tx * tx + ty * ty + tz * tz);

						if (alpha >= clipend) alpha = 0.0;
						else if (alpha <= starmindist) alpha = 0.0;
						else if (alpha <= 2.0 * starmindist) {
							alpha = (alpha - starmindist) / starmindist;
						} else {
							alpha -= 2.0 * starmindist;
							alpha /= (clipend - 2.0 * starmindist);
							alpha = 1.0 - alpha;
						}
					}


					if (alpha != 0.0) {
						fac = force * BLI_drand();

						har = initstar(vec, fac);

						if (har) {
							har->alfa = sqrt(sqrt(alpha));
							har->add= 255;
							har->r = har->g = har->b = 1.0;
							if (maxjit) {
								har->r += ((maxjit * BLI_drand()) ) - maxjit;
								har->g += ((maxjit * BLI_drand()) ) - maxjit;
								har->b += ((maxjit * BLI_drand()) ) - maxjit;
							}
							har->hard = 32;

							har->type |= HA_ONLYSKY;
							done++;
						}
					}
				}
			}
			/* do not call blender_test_break() here, since it is used in UI as well, confusing the callback system */
			/* main cause is G.afbreek of course, a global again... (ton) */
		}
	}
	if (termfunc) termfunc();
}

/* ------------------------------------------------------------------------- */

static VertRen *RE_duplicate_vertren(VertRen *ver)
{
	VertRen *v1= RE_findOrAddVert(R.totvert++);
	int index= v1->index;
	*v1= *ver;
	v1->index= index;
	return v1;
}

static void split_v_renderfaces(int startvlak, int startvert, int usize, int vsize, int uIndex, int cyclu, int cyclv)
{
	int vLen = vsize-1+(!!cyclv);
	int uLen = usize-1+(!!cyclu);
	int v;

	for (v=0; v<vLen; v++) {
		VlakRen *vlr = RE_findOrAddVlak(startvlak + vLen*uIndex + v);
		VertRen *vert = RE_duplicate_vertren(vlr->v2);

		if (cyclv) {
			vlr->v2 = vert;

			if (v==vLen-1) {
				VlakRen *vlr = RE_findOrAddVlak(startvlak + vLen*uIndex + 0);
				vlr->v1 = vert;
			} else {
				VlakRen *vlr = RE_findOrAddVlak(startvlak + vLen*uIndex + v+1);
				vlr->v1 = vert;
			}
		} else {
			vlr->v2 = vert;

			if (v<vLen-1) {
				VlakRen *vlr = RE_findOrAddVlak(startvlak + vLen*uIndex + v+1);
				vlr->v1 = vert;
			}

			if (v==0) {
				vlr->v1 = RE_duplicate_vertren(vlr->v1);
			} 
		}
	}
}

/* ------------------------------------------------------------------------- */

static int contrpuntnormr(float *n, float *puno)
{
	float inp;

	inp=n[0]*puno[0]+n[1]*puno[1]+n[2]*puno[2];
	if(inp<0.0) return 1;
	return 0;
}

/* ------------------------------------------------------------------------- */

static void calc_edge_stress_add(float *accum, VertRen *v1, VertRen *v2)
{
	float len= VecLenf(v1->co, v2->co)/VecLenf(v1->orco, v2->orco);
	float *acc;
	
	acc= accum + 2*v1->index;
	acc[0]+= len;
	acc[1]+= 1.0f;
	
	acc= accum + 2*v2->index;
	acc[0]+= len;
	acc[1]+= 1.0f;
}

static void calc_edge_stress(Mesh *me, int startvert, int startvlak)
{
	float loc[3], size[3], *accum, *acc, *accumoffs, *stress;
	int a;
	
	if(startvert==R.totvert) return;
	
	mesh_get_texspace(me, loc, NULL, size);
	
	accum= MEM_callocN(2*sizeof(float)*(R.totvert-startvert), "temp accum for stress");
	
	/* de-normalize orco */
	for(a=startvert; a<R.totvert; a++, acc+=2) {
		VertRen *ver= RE_findOrAddVert(a);
		if(ver->orco) {
			ver->orco[0]= ver->orco[0]*size[0] +loc[0];
			ver->orco[1]= ver->orco[1]*size[1] +loc[1];
			ver->orco[2]= ver->orco[2]*size[2] +loc[2];
		}
	}
	
	/* add stress values */
	accumoffs= accum - 2*startvert;	/* so we can use vertex index */
	for(a=startvlak; a<R.totvlak; a++) {
		VlakRen *vlr= RE_findOrAddVlak(a);

		if(vlr->v1->orco && vlr->v4) {
			calc_edge_stress_add(accumoffs, vlr->v1, vlr->v2);
			calc_edge_stress_add(accumoffs, vlr->v2, vlr->v3);
			calc_edge_stress_add(accumoffs, vlr->v3, vlr->v1);
			if(vlr->v4) {
				calc_edge_stress_add(accumoffs, vlr->v3, vlr->v4);
				calc_edge_stress_add(accumoffs, vlr->v4, vlr->v1);
				calc_edge_stress_add(accumoffs, vlr->v2, vlr->v4);
			}
		}
	}
	
	for(a=startvert; a<R.totvert; a++) {
		VertRen *ver= RE_findOrAddVert(a);
		if(ver->orco) {
			/* find stress value */
			acc= accumoffs + 2*ver->index;
			if(acc[1]!=0.0f)
				acc[0]/= acc[1];
			stress= RE_vertren_get_stress(ver, 1);
			*stress= *acc;
			
			/* restore orcos */
			ver->orco[0] = (ver->orco[0]-loc[0])/size[0];
			ver->orco[1] = (ver->orco[1]-loc[1])/size[1];
			ver->orco[2] = (ver->orco[2]-loc[2])/size[2];
		}
	}
	
	MEM_freeN(accum);
}

static void calc_tangent_vector(VlakRen *vlr, float fac1, float fac2, float fac3, float fac4)
{
	TFace *tface= vlr->tface;
	
	if(tface) {
		VertRen *v1=vlr->v1, *v2=vlr->v2, *v3=vlr->v3, *v4=vlr->v4;
		float *uv1= tface->uv[0], *uv2= tface->uv[1], *uv3= tface->uv[2], *uv4= tface->uv[3];
		float tang[3], *tav;
		float s1, s2, t1, t2, det;
		
		/* we calculate quads as two triangles, so weight for diagonal gets halved */
		if(v4) {
			fac1*= 0.5f;
			fac3*= 0.5f;
		}
		
		/* first tria, we use the V now */
		s1= uv2[0] - uv1[0];
		s2= uv3[0] - uv1[0];
		t1= uv2[1] - uv1[1];
		t2= uv3[1] - uv1[1];
		det= 1.0f / (s1 * t2 - s2 * t1);
		
		/* normals in render are inversed... */
		tang[0]= (t2 * (v1->co[0]-v2->co[0]) - t1 * (v1->co[0]-v3->co[0])); 
		tang[1]= (t2 * (v1->co[1]-v2->co[1]) - t1 * (v1->co[1]-v3->co[1]));
		tang[2]= (t2 * (v1->co[2]-v2->co[2]) - t1 * (v1->co[2]-v3->co[2]));
		
		tav= RE_vertren_get_tangent(v1, 1);
		VECADDFAC(tav, tav, tang, fac1);
		tav= RE_vertren_get_tangent(v2, 1);
		VECADDFAC(tav, tav, tang, fac2);
		tav= RE_vertren_get_tangent(v3, 1);
		VECADDFAC(tav, tav, tang, fac3);
		
		if(v4) {
			/* 2nd tria, we use the V now */
			s1= uv3[0] - uv1[0];
			s2= uv4[0] - uv1[0];
			t1= uv3[1] - uv1[1];
			t2= uv4[1] - uv1[1];
			det= 1.0f / (s1 * t2 - s2 * t1);
			
			/* normals in render are inversed... */
			tang[0]= (t2 * (v1->co[0]-v3->co[0]) - t1 * (v1->co[0]-v4->co[0])); 
			tang[1]= (t2 * (v1->co[1]-v3->co[1]) - t1 * (v1->co[1]-v4->co[1]));
			tang[2]= (t2 * (v1->co[2]-v3->co[2]) - t1 * (v1->co[2]-v4->co[2]));
			
			Normalise(tang);
			
			tav= RE_vertren_get_tangent(v1, 1);
			VECADDFAC(tav, tav, tang, fac1);
			tav= RE_vertren_get_tangent(v3, 1);
			VECADDFAC(tav, tav, tang, fac3);
			tav= RE_vertren_get_tangent(v4, 1);
			VECADDFAC(tav, tav, tang, fac4);
		}
	}	
}

static void calc_vertexnormals(int startvert, int startvlak, int do_tangent)
{
	int a;

		/* clear all vertex normals */
	for(a=startvert; a<R.totvert; a++) {
		VertRen *ver= RE_findOrAddVert(a);
		ver->n[0]=ver->n[1]=ver->n[2]= 0.0;
	}

		/* calculate cos of angles and point-masses, use as weight factor to
		   add face normal to vertex */
	for(a=startvlak; a<R.totvlak; a++) {
		VlakRen *vlr= RE_findOrAddVlak(a);
		if(vlr->flag & ME_SMOOTH) {
			VertRen *adrve1= vlr->v1;
			VertRen *adrve2= vlr->v2;
			VertRen *adrve3= vlr->v3;
			VertRen *adrve4= vlr->v4;
			float n1[3], n2[3], n3[3], n4[3];
			float fac1, fac2, fac3, fac4=0.0f;

			VecSubf(n1, adrve2->co, adrve1->co);
			Normalise(n1);
			VecSubf(n2, adrve3->co, adrve2->co);
			Normalise(n2);
			if(adrve4==NULL) {
				VecSubf(n3, adrve1->co, adrve3->co);
				Normalise(n3);

				fac1= saacos(-n1[0]*n3[0]-n1[1]*n3[1]-n1[2]*n3[2]);
				fac2= saacos(-n1[0]*n2[0]-n1[1]*n2[1]-n1[2]*n2[2]);
				fac3= saacos(-n2[0]*n3[0]-n2[1]*n3[1]-n2[2]*n3[2]);
			}
			else {
				VecSubf(n3, adrve4->co, adrve3->co);
				Normalise(n3);
				VecSubf(n4, adrve1->co, adrve4->co);
				Normalise(n4);

				fac1= saacos(-n4[0]*n1[0]-n4[1]*n1[1]-n4[2]*n1[2]);
				fac2= saacos(-n1[0]*n2[0]-n1[1]*n2[1]-n1[2]*n2[2]);
				fac3= saacos(-n2[0]*n3[0]-n2[1]*n3[1]-n2[2]*n3[2]);
				fac4= saacos(-n3[0]*n4[0]-n3[1]*n4[1]-n3[2]*n4[2]);

				if(!(vlr->flag & R_NOPUNOFLIP)) {
					if( contrpuntnormr(vlr->n, adrve4->n) ) fac4= -fac4;
				}

				adrve4->n[0] +=fac4*vlr->n[0];
				adrve4->n[1] +=fac4*vlr->n[1];
				adrve4->n[2] +=fac4*vlr->n[2];
			}

			if(!(vlr->flag & R_NOPUNOFLIP)) {
				if( contrpuntnormr(vlr->n, adrve1->n) ) fac1= -fac1;
				if( contrpuntnormr(vlr->n, adrve2->n) ) fac2= -fac2;
				if( contrpuntnormr(vlr->n, adrve3->n) ) fac3= -fac3;
			}

			adrve1->n[0] +=fac1*vlr->n[0];
			adrve1->n[1] +=fac1*vlr->n[1];
			adrve1->n[2] +=fac1*vlr->n[2];

			adrve2->n[0] +=fac2*vlr->n[0];
			adrve2->n[1] +=fac2*vlr->n[1];
			adrve2->n[2] +=fac2*vlr->n[2];

			adrve3->n[0] +=fac3*vlr->n[0];
			adrve3->n[1] +=fac3*vlr->n[1];
			adrve3->n[2] +=fac3*vlr->n[2];
			
			if(do_tangent)
				calc_tangent_vector(vlr, fac1, fac2, fac3, fac4);
		}
	}

		/* do solid faces */
	for(a=startvlak; a<R.totvlak; a++) {
		VlakRen *vlr= RE_findOrAddVlak(a);
		if((vlr->flag & ME_SMOOTH)==0) {
			float *f1= vlr->v1->n;
			if(f1[0]==0.0 && f1[1]==0.0 && f1[2]==0.0) VECCOPY(f1, vlr->n);
			f1= vlr->v2->n;
			if(f1[0]==0.0 && f1[1]==0.0 && f1[2]==0.0) VECCOPY(f1, vlr->n);
			f1= vlr->v3->n;
			if(f1[0]==0.0 && f1[1]==0.0 && f1[2]==0.0) VECCOPY(f1, vlr->n);
			if(vlr->v4) {
				f1= vlr->v4->n;
				if(f1[0]==0.0 && f1[1]==0.0 && f1[2]==0.0) VECCOPY(f1, vlr->n);
			}			
		}
	}
	
		/* normalise vertex normals */
	for(a=startvert; a<R.totvert; a++) {
		VertRen *ver= RE_findOrAddVert(a);
		Normalise(ver->n);
		if(do_tangent) {
			float *tav= RE_vertren_get_tangent(ver, 0);
			if(tav) Normalise(tav);
		}
	}

		/* vertex normal (puno) switch flags for during render */
	for(a=startvlak; a<R.totvlak; a++) {
		VlakRen *vlr= RE_findOrAddVlak(a);

		if((vlr->flag & R_NOPUNOFLIP)==0) {
			VertRen *adrve1= vlr->v1;
			VertRen *adrve2= vlr->v2;
			VertRen *adrve3= vlr->v3;
			VertRen *adrve4= vlr->v4;
			vlr->puno &= ~15;
			if ((vlr->n[0]*adrve1->n[0]+vlr->n[1]*adrve1->n[1]+vlr->n[2]*adrve1->n[2])<0.0) vlr->puno= 1;
			if ((vlr->n[0]*adrve2->n[0]+vlr->n[1]*adrve2->n[1]+vlr->n[2]*adrve2->n[2])<0.0) vlr->puno+= 2;
			if ((vlr->n[0]*adrve3->n[0]+vlr->n[1]*adrve3->n[1]+vlr->n[2]*adrve3->n[2])<0.0) vlr->puno+= 4;
			if(adrve4) {
				if((vlr->n[0]*adrve4->n[0]+vlr->n[1]*adrve4->n[1]+vlr->n[2]*adrve4->n[2])<0.0) vlr->puno+= 8;
			}
		}
	}
}

/* ------------------------------------------------------------------------- */
/* Autosmoothing:                                                            */
/* ------------------------------------------------------------------------- */

typedef struct ASvert {
	int totface;
	ListBase faces;
} ASvert;

typedef struct ASface {
	struct ASface *next, *prev;
	VlakRen *vlr[4];
	VertRen *nver[4];
} ASface;

static void as_addvert(ASvert *asv, VertRen *v1, VlakRen *vlr)
{
	ASface *asf;
	int a;
	
	if(v1 == NULL) return;
	
	if(asv->faces.first==NULL) {
		asf= MEM_callocN(sizeof(ASface), "asface");
		BLI_addtail(&asv->faces, asf);
	}
	
	asf= asv->faces.last;
	for(a=0; a<4; a++) {
		if(asf->vlr[a]==NULL) {
			asf->vlr[a]= vlr;
			asv->totface++;
			break;
		}
	}
	
	/* new face struct */
	if(a==4) {
		asf= MEM_callocN(sizeof(ASface), "asface");
		BLI_addtail(&asv->faces, asf);
		asf->vlr[0]= vlr;
		asv->totface++;
	}
}

static int as_testvertex(VlakRen *vlr, VertRen *ver, ASvert *asv, float thresh) 
{
	/* return 1: vertex needs a copy */
	ASface *asf;
	float inp;
	int a;
	
	if(vlr==0) return 0;
	
	asf= asv->faces.first;
	while(asf) {
		for(a=0; a<4; a++) {
			if(asf->vlr[a] && asf->vlr[a]!=vlr) {
				inp= fabs( vlr->n[0]*asf->vlr[a]->n[0] + vlr->n[1]*asf->vlr[a]->n[1] + vlr->n[2]*asf->vlr[a]->n[2] );
				if(inp < thresh) return 1;
			}
		}
		asf= asf->next;
	}
	
	return 0;
}

static VertRen *as_findvertex(VlakRen *vlr, VertRen *ver, ASvert *asv, float thresh) 
{
	/* return when new vertex already was made */
	ASface *asf;
	float inp;
	int a;
	
	asf= asv->faces.first;
	while(asf) {
		for(a=0; a<4; a++) {
			if(asf->vlr[a] && asf->vlr[a]!=vlr) {
				/* this face already made a copy for this vertex! */
				if(asf->nver[a]) {
					inp= fabs( vlr->n[0]*asf->vlr[a]->n[0] + vlr->n[1]*asf->vlr[a]->n[1] + vlr->n[2]*asf->vlr[a]->n[2] );
					if(inp >= thresh) {
						return asf->nver[a];
					}
				}
			}
		}
		asf= asf->next;
	}
	
	return NULL;
}

static void autosmooth(int startvert, int startvlak, int degr)
{
	ASvert *asv, *asverts, *asvertoffs;
	ASface *asf;
	VertRen *ver, *v1;
	VlakRen *vlr;
	float thresh;
	int a, b, totvert;
	
	if(startvert==R.totvert) return;
	asverts= MEM_callocN(sizeof(ASvert)*(R.totvert-startvert), "all smooth verts");
	asvertoffs= asverts-startvert;	 /* se we can use indices */
	
	thresh= cos( M_PI*((float)degr)/180.0 );
	
	/* step one: construct listbase of all vertices and pointers to faces */
	for(a=startvlak; a<R.totvlak; a++) {
		vlr= RE_findOrAddVlak(a);
		
		as_addvert(asvertoffs+vlr->v1->index, vlr->v1, vlr);
		as_addvert(asvertoffs+vlr->v2->index, vlr->v2, vlr);
		as_addvert(asvertoffs+vlr->v3->index, vlr->v3, vlr);
		if(vlr->v4) 
			as_addvert(asvertoffs+vlr->v4->index, vlr->v4, vlr);
	}
	
	/* we now test all vertices, when faces have a normal too much different: they get a new vertex */
	totvert= R.totvert;
	for(a=startvert, asv=asverts; a<totvert; a++, asv++) {
		if(asv && asv->totface>1) {
			ver= RE_findOrAddVert(a);
			
			asf= asv->faces.first;
			while(asf) {
				for(b=0; b<4; b++) {
				
					/* is there a reason to make a new vertex? */
					vlr= asf->vlr[b];
					if( as_testvertex(vlr, ver, asv, thresh) ) {
						
						/* already made a new vertex within threshold? */
						v1= as_findvertex(vlr, ver, asv, thresh);
						if(v1==NULL) {
							/* make a new vertex */
							v1= RE_duplicate_vertren(ver);
						}
						asf->nver[b]= v1;
						if(vlr->v1==ver) vlr->v1= v1;
						if(vlr->v2==ver) vlr->v2= v1;
						if(vlr->v3==ver) vlr->v3= v1;
						if(vlr->v4==ver) vlr->v4= v1;
					}
				}
				asf= asf->next;
			}
		}
	}
	
	/* free */
	for(a=0; a<totvert-startvert; a++) {
		BLI_freelistN(&asverts[a].faces);
	}
	MEM_freeN(asverts);
}

/* ------------------------------------------------------------------------- */
/* End of autosmoothing:                                                     */
/* ------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------- */
/* Orco hash																 */
/* ------------------------------------------------------------------------- */


static GHash *g_orco_hash = NULL;

static float *get_object_orco(Object *ob)
{
	float *orco;
	
	if (!g_orco_hash)
		g_orco_hash = BLI_ghash_new(BLI_ghashutil_ptrhash, BLI_ghashutil_ptrcmp);
	
	orco = BLI_ghash_lookup(g_orco_hash, ob);
	
	if (!orco) {
		if (ob->type==OB_MESH) {
			orco = mesh_create_orco_render(ob);
		} else if (ELEM(ob->type, OB_CURVE, OB_FONT)) {
			orco = make_orco_curve(ob);
		} else if (ob->type==OB_SURF) {
			orco = make_orco_surf(ob);
		}
		
		if (orco)
			BLI_ghash_insert(g_orco_hash, ob, orco);
	}
	
	return orco;
}
static void free_mesh_orco_hash(void) 
{
	if (g_orco_hash) {
		BLI_ghash_free(g_orco_hash, NULL, (GHashValFreeFP)MEM_freeN);
		g_orco_hash = NULL;
	}
}

/* ******************** END ORCO HASH ***************** */


static void make_render_halos(Object *ob, Mesh *me, int totvert, MVert *mvert, Material *ma, float *orco)
{
	HaloRen *har;
	float xn, yn, zn, nor[3], view[3];
	float vec[3], hasize, mat[4][4], imat[3][3];
	int a, ok, seed= ma->seed1;

	MTC_Mat4MulMat4(mat, ob->obmat, R.viewmat);
	MTC_Mat3CpyMat4(imat, ob->imat);

	R.flag |= R_HALO;

	for(a=0; a<totvert; a++, mvert++) {
		ok= 1;

		if(ok) {
			hasize= ma->hasize;

			VECCOPY(vec, mvert->co);
			MTC_Mat4MulVecfl(mat, vec);

			if(ma->mode & MA_HALOPUNO) {
				xn= mvert->no[0];
				yn= mvert->no[1];
				zn= mvert->no[2];

				/* transpose ! */
				nor[0]= imat[0][0]*xn+imat[0][1]*yn+imat[0][2]*zn;
				nor[1]= imat[1][0]*xn+imat[1][1]*yn+imat[1][2]*zn;
				nor[2]= imat[2][0]*xn+imat[2][1]*yn+imat[2][2]*zn;
				Normalise(nor);

				VECCOPY(view, vec);
				Normalise(view);

				zn= nor[0]*view[0]+nor[1]*view[1]+nor[2]*view[2];
				if(zn>=0.0) hasize= 0.0;
				else hasize*= zn*zn*zn*zn;
			}

			if(orco) har= RE_inithalo(ma, vec, NULL, orco, hasize, 0.0, seed);
			else har= RE_inithalo(ma, vec, NULL, mvert->co, hasize, 0.0, seed);
			if(har) har->lay= ob->lay;
		}
		if(orco) orco+= 3;
		seed++;
	}
}

/* ------------------------------------------------------------------------- */


static void render_particle_system(Object *ob, PartEff *paf)
{
	Particle *pa=0;
	HaloRen *har=0;
	Material *ma=0;
	float xn, yn, zn, imat[3][3], tmat[4][4], mat[4][4], hasize, stime, ptime, ctime, vec[3], vec1[3], view[3], nor[3];
	int a, mat_nr=1, seed;

	pa= paf->keys;
	if(pa==NULL || paf->disp!=100) {
		build_particle_system(ob);
		pa= paf->keys;
		if(pa==NULL) return;
	}

	ma= give_render_material(ob, paf->omat);
	
	MTC_Mat4MulMat4(mat, ob->obmat, R.viewmat);
	MTC_Mat4Invert(ob->imat, mat);	/* this is correct, for imat texture */

	/* enable duplicators to work */
	Mat4MulMat4(tmat, paf->imat, ob->obmat);
	MTC_Mat4MulMat4(mat, tmat, R.viewmat);
	
	MTC_Mat4Invert(tmat, mat);
	MTC_Mat3CpyMat4(imat, tmat);

	R.flag |= R_HALO;

	if(ob->ipoflag & OB_OFFS_PARTICLE) ptime= ob->sf;
	else ptime= 0.0;
	ctime= bsystem_time(ob, 0, (float)G.scene->r.cfra, ptime);
	seed= ma->seed1;

	for(a=0; a<paf->totpart; a++, pa+=paf->totkey, seed++) {

		/* offset time for calculating normal */
		stime= ctime;
		ptime= ctime+1.0f;
		if(ctime < pa->time) {
			if(paf->flag & PAF_UNBORN)
				ptime= pa->time+1.0f;
			else
				continue;
		}
		if(ctime > pa->time+pa->lifetime) {
			if(paf->flag & PAF_DIED)
				stime= pa->time+pa->lifetime-1.0f;
			else
				continue;
		}
		
		/* watch it: also calculate the normal of a particle */
		if(paf->stype==PAF_VECT || ma->mode & MA_HALO_SHADE) {
			where_is_particle(paf, pa, stime, vec);
			MTC_Mat4MulVecfl(mat, vec);
			where_is_particle(paf, pa, ptime, vec1);
			MTC_Mat4MulVecfl(mat, vec1);
		}
		else {
			where_is_particle(paf, pa, ctime, vec);
			MTC_Mat4MulVecfl(mat, vec);
		}

		if(pa->mat_nr != mat_nr) {
			mat_nr= pa->mat_nr;
			ma= give_render_material(ob, mat_nr);
		}

		if(ma->ipo) {
			/* correction for lifetime */
			ptime= 100.0*(ctime-pa->time)/pa->lifetime;
			calc_ipo(ma->ipo, ptime);
			execute_ipo((ID *)ma, ma->ipo);
		}

		hasize= ma->hasize;

		if(ma->mode & MA_HALOPUNO) {
			xn= pa->no[0];
			yn= pa->no[1];
			zn= pa->no[2];

			/* transpose ! */
			nor[0]= imat[0][0]*xn+imat[0][1]*yn+imat[0][2]*zn;
			nor[1]= imat[1][0]*xn+imat[1][1]*yn+imat[1][2]*zn;
			nor[2]= imat[2][0]*xn+imat[2][1]*yn+imat[2][2]*zn;
			Normalise(nor);

			VECCOPY(view, vec);
			Normalise(view);

			zn= nor[0]*view[0]+nor[1]*view[1]+nor[2]*view[2];
			if(zn>=0.0) hasize= 0.0;
			else hasize*= zn*zn*zn*zn;
		}

		if(paf->stype==PAF_VECT) har= RE_inithalo(ma, vec, vec1, pa->co, hasize, paf->vectsize, seed);
		else {
			har= RE_inithalo(ma, vec, NULL, pa->co, hasize, 0.0, seed);
			if(har && ma->mode & MA_HALO_SHADE) {
				VecSubf(har->no, vec, vec1);
				Normalise(har->no);
			}
		}
		if(har) har->lay= ob->lay;
	}

	/* restore material */
	for(a=1; a<=ob->totcol; a++) {
		ma= give_render_material(ob, a);
		if(ma) do_mat_ipo(ma);
	}
	
	if(paf->disp!=100) {
		MEM_freeN(paf->keys);
		paf->keys= NULL;
	}
}


/* ------------------------------------------------------------------------- */

/* old call... for when duplicators were using temporal objects */
static Object *vlr_set_ob(Object *ob)
{
	return ob;
}

/* ------------------------------------------------------------------------- */

static void static_particle_strand(Object *ob, Material *ma, float *orco, float *vec, float *vec1, float ctime, int first)
{
	static VertRen *v1= NULL, *v2= NULL;
	VlakRen *vlr;
	float nor[3], cross[3], w, dx, dy;
	int flag;
	
	VecSubf(nor, vec, vec1);
	Normalise(nor);		// nor needed as tangent 
	Crossf(cross, vec, nor);
	
	/* turn cross in pixelsize */
	w= vec[2]*R.winmat[2][3] + R.winmat[3][3];
	dx= R.rectx*cross[0]*R.winmat[0][0]/w;
	dy= R.recty*cross[1]*R.winmat[1][1]/w;
	w= sqrt(dx*dx + dy*dy);
	if(w!=0.0f) {
		float fac;
		if(ma->strand_ease!=0.0f) {
			if(ma->strand_ease<0.0f)
				fac= pow(ctime, 1.0+ma->strand_ease);
			else
				fac= pow(ctime, 1.0/(1.0f-ma->strand_ease));
		}
		else fac= ctime;
		
		VecMulf(cross, ((1.0f-fac)*ma->strand_sta + (fac)*ma->strand_end)/w);
	}
	
	if(ma->mode & MA_TANGENT_STR)
		flag= R_SMOOTH|R_NOPUNOFLIP|R_STRAND|R_TANGENT;
	else
		flag= R_SMOOTH|R_STRAND;
	
	/* first two vertices */
	if(first) {
		v1= RE_findOrAddVert(R.totvert++);
		v2= RE_findOrAddVert(R.totvert++);
		
		VECCOPY(v1->co, vec);
		VecAddf(v1->co, v1->co, cross);
		VECCOPY(v1->n, nor);
		v1->orco= orco;
		v1->accum= -1.0f;	// accum abuse for strand texco
		
		VECCOPY(v2->co, vec);
		VecSubf(v2->co, v2->co, cross);
		VECCOPY(v2->n, nor);
		v2->orco= orco;
		v2->accum= v1->accum;
	}
	else {
		
		vlr= RE_findOrAddVlak(R.totvlak++);
		vlr->flag= flag;
		vlr->ob= vlr_set_ob(ob);
		vlr->v1= v1;
		vlr->v2= v2;
		vlr->v3= RE_findOrAddVert(R.totvert++);
		vlr->v4= RE_findOrAddVert(R.totvert++);
		
		v1= vlr->v4; // cycle
		v2= vlr->v3; // cycle
		
		VECCOPY(vlr->v4->co, vec);
		VecAddf(vlr->v4->co, vlr->v4->co, cross);
		VECCOPY(vlr->v4->n, nor);
		vlr->v4->orco= orco;
		vlr->v4->accum= -1.0f + 2.0f*ctime;	// accum abuse for strand texco
		
		VECCOPY(vlr->v3->co, vec);
		VecSubf(vlr->v3->co, vlr->v3->co, cross);
		VECCOPY(vlr->v3->n, nor);
		vlr->v3->orco= orco;
		vlr->v3->accum= vlr->v4->accum;
		
		CalcNormFloat4(vlr->v4->co, vlr->v3->co, vlr->v2->co, vlr->v1->co, vlr->n);
		
		vlr->mat= ma;
		vlr->ec= ME_V2V3;
		vlr->lay= ob->lay;
	}
}

static void render_static_particle_system(Object *ob, PartEff *paf)
{
	Particle *pa=0;
	HaloRen *har=0;
	Material *ma=0;
	VertRen *v1= NULL, *v2= NULL;
	VlakRen *vlr;
	float xn, yn, zn, imat[3][3], mat[4][4], hasize;
	float mtime, ptime, ctime, vec[3], vec1[3], view[3], nor[3];
	float *orco= NULL, loc_tex[3], size_tex[3];
	int a, mat_nr=1, seed, totvlako, totverto, first;

	pa= paf->keys;
	if(pa==NULL || (paf->flag & PAF_ANIMATED) || paf->disp!=100) {
		build_particle_system(ob);
		pa= paf->keys;
		if(pa==NULL) return;
	}

	totvlako= R.totvlak;
	totverto= R.totvert;
	
	ma= give_render_material(ob, paf->omat);
	if(ma->mode & MA_HALO)
		R.flag |= R_HALO;

	MTC_Mat4MulMat4(mat, ob->obmat, R.viewmat);
	MTC_Mat4Invert(ob->imat, mat);	/* need to be that way, for imat texture */

	MTC_Mat3CpyMat4(imat, ob->imat);
	
	/* orcos */
	if(!(ma->mode & (MA_HALO|MA_WIRE))) {
		orco= MEM_mallocN(3*sizeof(float)*paf->totpart, "static particle orcos");
		if (!g_orco_hash)
			g_orco_hash = BLI_ghash_new(BLI_ghashutil_ptrhash, BLI_ghashutil_ptrcmp);
		BLI_ghash_insert(g_orco_hash, paf, orco);	/* pointer is particles, otherwise object uses it */
	}
	
	mesh_get_texspace(ob->data, loc_tex, NULL, size_tex);
	
	if(ob->ipoflag & OB_OFFS_PARTICLE) ptime= ob->sf;
	else ptime= 0.0;
	ctime= bsystem_time(ob, 0, (float)G.scene->r.cfra, ptime);
	seed= ma->seed1;

	for(a=0; a<paf->totpart; a++, pa+=paf->totkey) {
		
		where_is_particle(paf, pa, pa->time, vec1);
		if(orco) {
			orco[0] = (vec1[0]-loc_tex[0])/size_tex[0];
			orco[1] = (vec1[1]-loc_tex[1])/size_tex[1];
			orco[2] = (vec1[2]-loc_tex[2])/size_tex[2];
		}
		MTC_Mat4MulVecfl(mat, vec1);
		mtime= pa->time+pa->lifetime+paf->staticstep-1;
		
		first= 1;
		for(ctime= pa->time; ctime<mtime; ctime+=paf->staticstep) {
			
			/* make sure hair grows until the end.. */
			if(ctime>pa->time+pa->lifetime) ctime= pa->time+pa->lifetime;
			
			/* watch it: also calc the normal of a particle */
			if(paf->stype==PAF_VECT || ma->mode & MA_HALO_SHADE) {
				where_is_particle(paf, pa, ctime+1.0, vec);
				MTC_Mat4MulVecfl(mat, vec);
			}
			else {
				where_is_particle(paf, pa, ctime, vec);
				MTC_Mat4MulVecfl(mat, vec);
			}

			if(pa->mat_nr != mat_nr) {
				mat_nr= pa->mat_nr;
				ma= give_render_material(ob, mat_nr);
			}
			
			/* wires */
			if(ma->mode & MA_WIRE) {
				if(ctime == pa->time) {
					v1= RE_findOrAddVert(R.totvert++);
					VECCOPY(v1->co, vec);
				}
				else {
					vlr= RE_findOrAddVlak(R.totvlak++);
					vlr->ob= vlr_set_ob(ob);
					vlr->v1= v1;
					vlr->v2= RE_findOrAddVert(R.totvert++);
					vlr->v3= vlr->v2;
					vlr->v4= NULL;
					
					v1= vlr->v2; // cycle
					VECCOPY(v1->co, vec);
					
					VecSubf(vlr->n, vec, vec1);
					Normalise(vlr->n);
					VECCOPY(v1->n, vlr->n);
					
					vlr->mat= ma;
					vlr->ec= ME_V1V2;
					vlr->lay= ob->lay;
				}
			}
			else {
				if(ma->ipo) {
					/* correction for lifetime */
					ptime= 100.0*(ctime-pa->time)/pa->lifetime;
					calc_ipo(ma->ipo, ptime);
					execute_ipo((ID *)ma, ma->ipo);
				}
				
				if(ma->mode & MA_HALO) {
					hasize= ma->hasize;

					if(ma->mode & MA_HALOPUNO) {
						xn= pa->no[0];
						yn= pa->no[1];
						zn= pa->no[2];

						/* transpose ! */
						nor[0]= imat[0][0]*xn+imat[0][1]*yn+imat[0][2]*zn;
						nor[1]= imat[1][0]*xn+imat[1][1]*yn+imat[1][2]*zn;
						nor[2]= imat[2][0]*xn+imat[2][1]*yn+imat[2][2]*zn;
						Normalise(nor);

						VECCOPY(view, vec);
						Normalise(view);

						zn= nor[0]*view[0]+nor[1]*view[1]+nor[2]*view[2];
						if(zn>=0.0) hasize= 0.0;
						else hasize*= zn*zn*zn*zn;
					}

					if(paf->stype==PAF_VECT) har= RE_inithalo(ma, vec, vec1, pa->co, hasize, paf->vectsize, seed);
					else {
						har= RE_inithalo(ma, vec, NULL, pa->co, hasize, 0.0, seed);
						if(har && (ma->mode & MA_HALO_SHADE)) {
							VecSubf(har->no, vec, vec1);
							Normalise(har->no);
							har->lay= ob->lay;
						}
					}
					if(har) har->lay= ob->lay;
				}
				else {	/* generate pixel sized hair strand */
					static_particle_strand(ob, ma, orco, vec, vec1, (ctime-pa->time)/(mtime-pa->time), first);
				}
			}
			
			VECCOPY(vec1, vec);
			first= 0;
		}
		
		seed++;
		if(orco) orco+=3;
	}

	if(paf->disp!=100) {
		MEM_freeN(paf->keys);
		paf->keys= NULL;
	}

	if((ma->mode & MA_TANGENT_STR)==0)
		calc_vertexnormals(totverto, totvlako, 0);
}


/* ------------------------------------------------------------------------- */

static int verghalo(const void *a1, const void *a2)
{
	const struct halosort *x1=a1, *x2=a2;
	
	if( x1->z < x2->z ) return 1;
	else if( x1->z > x2->z) return -1;
	return 0;
}

/* ------------------------------------------------------------------------- */
extern int rblohalen;
static void sort_halos(void)
{
	struct halosort *hablock, *haso;
	HaloRen *har = NULL, **bloha;
	int a;

	if(R.tothalo==0) return;

	/* make datablock with halo pointers, sort */
	haso= hablock= MEM_mallocN(sizeof(struct halosort)*R.tothalo, "hablock");

	for(a=0; a<R.tothalo; a++) {
		if((a & 255)==0) har= R.bloha[a>>8];
		else har++;
		haso->har= har;
		haso->z= har->zs;
		haso++;
	}

	qsort(hablock, R.tothalo, sizeof(struct halosort), verghalo);

	/* re-assamble R.bloha */

	bloha= R.bloha;
	R.bloha= (HaloRen **)MEM_callocN(sizeof(void *)*(rblohalen),"Bloha");

	haso= hablock;
	for(a=0; a<R.tothalo; a++) {
		har= RE_findOrAddHalo(a);
		*har= *(haso->har);

		haso++;
	}

	/* free */
	a= 0;
	while(bloha[a]) {
		MEM_freeN(bloha[a]);
		a++;
	}
	MEM_freeN(bloha);
	MEM_freeN(hablock);

}



static Material *give_render_material(Object *ob, int nr)
{
	extern Material defmaterial;	// initrender.c
	Object *temp;
	Material *ma;

	ma= give_current_material(ob, nr);
	if(ma==NULL) ma= &defmaterial;

	return ma;
}

/* ------------------------------------------------------------------------- */
static void init_render_mball(Object *ob)
{
	DispList *dl, *dlo;
	VertRen *ver;
	VlakRen *vlr, *vlr1;
	Material *ma;
	float *data, *nors, mat[4][4], imat[3][3], xn, yn, zn;
	int a, need_orco, startvert, *index;

	if (ob!=find_basis_mball(ob))
		return;

	MTC_Mat4MulMat4(mat, ob->obmat, R.viewmat);
	MTC_Mat4Invert(ob->imat, mat);
	MTC_Mat3CpyMat4(imat, ob->imat);

	ma= give_render_material(ob, 1);

	need_orco= 0;
	if(ma->texco & TEXCO_ORCO) {
		need_orco= 1;
	}

	dlo= ob->disp.first;
	if(dlo) BLI_remlink(&ob->disp, dlo);

	makeDispListMBall(ob);
	dl= ob->disp.first;
	if(dl==0) return;

	startvert= R.totvert;
	data= dl->verts;
	nors= dl->nors;

	for(a=0; a<dl->nr; a++, data+=3, nors+=3) {

		ver= RE_findOrAddVert(R.totvert++);
		VECCOPY(ver->co, data);
		MTC_Mat4MulVecfl(mat, ver->co);

		xn= nors[0];
		yn= nors[1];
		zn= nors[2];

		/* transpose ! */
		ver->n[0]= imat[0][0]*xn+imat[0][1]*yn+imat[0][2]*zn;
		ver->n[1]= imat[1][0]*xn+imat[1][1]*yn+imat[1][2]*zn;
		ver->n[2]= imat[2][0]*xn+imat[2][1]*yn+imat[2][2]*zn;
		Normalise(ver->n);
		//if(ob->transflag & OB_NEG_SCALE) VecMulf(ver->n. -1.0);
		
		if(need_orco) ver->orco= data;
	}

	index= dl->index;
	for(a=0; a<dl->parts; a++, index+=4) {

		vlr= RE_findOrAddVlak(R.totvlak++);
		vlr->ob= vlr_set_ob(ob);
		vlr->v1= RE_findOrAddVert(startvert+index[0]);
		vlr->v2= RE_findOrAddVert(startvert+index[1]);
		vlr->v3= RE_findOrAddVert(startvert+index[2]);
		vlr->v4= 0;

		if(ob->transflag & OB_NEG_SCALE) 
			CalcNormFloat(vlr->v1->co, vlr->v2->co, vlr->v3->co, vlr->n);
		else
			CalcNormFloat(vlr->v3->co, vlr->v2->co, vlr->v1->co, vlr->n);

		vlr->mat= ma;
		vlr->flag= ME_SMOOTH+R_NOPUNOFLIP;
		vlr->ec= 0;
		vlr->lay= ob->lay;

		/* mball -too bad- always has triangles, because quads can be non-planar */
		if(index[3]) {
			vlr1= RE_findOrAddVlak(R.totvlak++);
			*vlr1= *vlr;
			vlr1->v2= vlr1->v3;
			vlr1->v3= RE_findOrAddVert(startvert+index[3]);
			if(ob->transflag & OB_NEG_SCALE) 
				CalcNormFloat(vlr1->v1->co, vlr1->v2->co, vlr1->v3->co, vlr1->n);
			else
				CalcNormFloat(vlr1->v3->co, vlr1->v2->co, vlr1->v1->co, vlr1->n);
		}
	}

	if(need_orco) {
		/* store displist and scale */
		make_orco_mball(ob);
		if(dlo) BLI_addhead(&ob->disp, dlo);

	}
	else {
		freedisplist(&ob->disp);
		if(dlo) BLI_addtail(&ob->disp, dlo);
	}
}
/* ------------------------------------------------------------------------- */
/* convert */

struct edgesort {
	int v1, v2;
	int has_mcol;
	TFace *tface;
	float uv1[2], uv2[2];
	unsigned int mcol1, mcol2;
};

/* edges have to be added with lowest index first for sorting */
static void to_edgesort(struct edgesort *ed, int i1, int i2, int v1, int v2, unsigned int *mcol, TFace *tface)
{
	if(v1<v2) {
		ed->v1= v1; ed->v2= v2;
	}
	else {
		ed->v1= v2; ed->v2= v1;
		SWAP(int, i1, i2);
	}
	/* copy color and tface, edges use different ordering */
	ed->tface= tface;
	if(tface) {
		ed->uv1[0]= tface->uv[i1][0];
		ed->uv1[1]= tface->uv[i1][1];
		ed->uv2[0]= tface->uv[i2][0];
		ed->uv2[1]= tface->uv[i2][1];
		
		ed->mcol1= tface->col[i1];
		ed->mcol2= tface->col[i2];
	}	
	ed->has_mcol= mcol!=NULL;
	if(mcol) {
		ed->mcol1= mcol[i1];
		ed->mcol2= mcol[i2];
	}
}

static int vergedgesort(const void *v1, const void *v2)
{
	const struct edgesort *x1=v1, *x2=v2;
	
	if( x1->v1 > x2->v1) return 1;
	else if( x1->v1 < x2->v1) return -1;
	else if( x1->v2 > x2->v2) return 1;
	else if( x1->v2 < x2->v2) return -1;
	
	return 0;
}

static struct edgesort *make_mesh_edge_lookup(Mesh *me, DispListMesh *dlm, int *totedgesort)
{
	MFace *mf, *mface;
	TFace *tface=NULL;
	struct edgesort *edsort, *ed;
	unsigned int *mcol=NULL;
	int a, totedge=0, totface;
	
	if (dlm) {
		mface= dlm->mface;
		totface= dlm->totface;
		if (dlm->tface) 
			tface= dlm->tface;
		else if (dlm->mcol)
			mcol= (unsigned int *)dlm->mcol;
	} else {
		mface= me->mface;
		totface= me->totface;
		if (me->tface) 
			tface= me->tface;
		else if (me->mcol)
			mcol= (unsigned int *)me->mcol;
	}
	
	if(mcol==NULL && tface==NULL) return NULL;
	
	/* make sorted table with edges and and tface/mcol pointers in it */
	for(a= totface, mf= mface; a>0; a--, mf++) {
		if(mf->v4) totedge+=4;
		else if(mf->v3) totedge+=3;
	}
	if(totedge==0) return NULL;
	
	ed= edsort= MEM_mallocN(totedge*sizeof(struct edgesort), "edgesort");
	
	for(a= me->totface, mf= mface; a>0; a--, mf++) {
		if(mface->v4 || mface->v3) {
			to_edgesort(ed++, 0, 1, mf->v1, mf->v2, mcol, tface);
			to_edgesort(ed++, 1, 2, mf->v2, mf->v3, mcol, tface);
			if(mf->v4) {
				to_edgesort(ed++, 2, 3, mf->v3, mf->v4, mcol, tface);
				to_edgesort(ed++, 3, 0, mf->v4, mf->v1, mcol, tface);
			}
			else if(mf->v3) {
				to_edgesort(ed++, 2, 3, mf->v3, mf->v1, mcol, tface);
			}
		}
		if(mcol) mcol+=4;
		if(tface) tface++;
	}
	
	qsort(edsort, totedge, sizeof(struct edgesort), vergedgesort);
	
	*totedgesort= totedge;
	return edsort;
}

static void use_mesh_edge_lookup(Mesh *me, DispListMesh *dlm, MEdge *medge, VlakRen *vlr, struct edgesort *edgetable, int totedge)
{
	struct edgesort ed, *edp;
	
	if(medge->v1 < medge->v2) {
		ed.v1= medge->v1; ed.v2= medge->v2;
	}
	else {
		ed.v1= medge->v2; ed.v2= medge->v1;
	}
	
	edp= bsearch(&ed, edgetable, totedge, sizeof(struct edgesort), vergedgesort);
	if(edp) {
		/* since edges have different index ordering, we have to duplicate mcol and tface */
		if(edp->tface) {
			vlr->tface= BLI_memarena_alloc(R.memArena, sizeof(TFace));
			vlr->vcol= vlr->tface->col;
			memcpy(vlr->tface, edp->tface, sizeof(TFace));
			
			if(edp->v1==medge->v1) {
				vlr->vcol[0]= edp->mcol1;
				vlr->vcol[1]= edp->mcol2;
			}
			else {
				vlr->vcol[0]= edp->mcol2;
				vlr->vcol[1]= edp->mcol1;
			}
			vlr->vcol[2]= vlr->vcol[1];
			vlr->vcol[3]= vlr->vcol[1];
			
			if(edp->v1==medge->v1) {
				memcpy(vlr->tface->uv[0], edp->uv1, 2*sizeof(float));
				memcpy(vlr->tface->uv[1], edp->uv2, 2*sizeof(float));
			}
			else {
				memcpy(vlr->tface->uv[0], edp->uv2, 2*sizeof(float));
				memcpy(vlr->tface->uv[1], edp->uv1, 2*sizeof(float));
			}
			memcpy(vlr->tface->uv[2], vlr->tface->uv[1], 2*sizeof(float));
			memcpy(vlr->tface->uv[3], vlr->tface->uv[1], 2*sizeof(float));
		} 
		else if(edp->has_mcol) {
			vlr->vcol= BLI_memarena_alloc(R.memArena, sizeof(MCol)*4);
			vlr->vcol[0]= edp->mcol1;
			vlr->vcol[1]= edp->mcol2;
			vlr->vcol[2]= vlr->vcol[1];
			vlr->vcol[3]= vlr->vcol[1];
		}
	}
}

static void init_render_mesh(Object *ob)
{
	Mesh *me;
	MVert *mvert = NULL;
	MFace *mface;
	VlakRen *vlr; //, *vlr1;
	VertRen *ver;
	Material *ma;
	MSticky *ms = NULL;
	PartEff *paf;
	DispListMesh *dlm = NULL;
	DerivedMesh *dm;
	unsigned int *vertcol;
	float xn, yn, zn,  imat[3][3], mat[4][4];  //nor[3],
	float *orco=0;
	int a, a1, ok, need_orco=0, need_stress=0, need_tangent=0, totvlako, totverto, vertofs;
	int end, do_autosmooth=0, totvert = 0, dm_needsfree;
	
	me= ob->data;

	paf = give_parteff(ob);
	if(paf) {
		/* warning; build_particle_system does modifier calls itself */
		if(paf->flag & PAF_STATIC) render_static_particle_system(ob, paf);
		else render_particle_system(ob, paf);
		if((paf->flag & PAF_SHOWE)==0) return;
	}

	MTC_Mat4MulMat4(mat, ob->obmat, R.viewmat);
	MTC_Mat4Invert(ob->imat, mat);
	MTC_Mat3CpyMat4(imat, ob->imat);

	if(me->totvert==0) {
		return;
	}
	
	totvlako= R.totvlak;
	totverto= R.totvert;

	need_orco= 0;
	for(a=1; a<=ob->totcol; a++) {
		ma= give_render_material(ob, a);
		if(ma) {
			if(ma->texco & (TEXCO_ORCO|TEXCO_STRESS))
				need_orco= 1;
			if(ma->texco & TEXCO_STRESS)
				need_stress= 1;
			if(ma->mode & MA_TANGENT_V)
				need_tangent= 1;
		}
	}
	
	if(need_orco) orco = get_object_orco(ob);
	
	dm = mesh_create_derived_render(ob);
	dm_needsfree= 1;
	
	if(dm==NULL) return;	/* in case duplicated object fails? */
	
	dlm = dm->convertToDispListMesh(dm, 1);

	mvert= dlm->mvert;
	totvert= dlm->totvert;

	ms = (totvert==me->totvert)?me->msticky:NULL;
	
	ma= give_render_material(ob, 1);

	if(ma->mode & MA_HALO) {
		make_render_halos(ob, me, totvert, mvert, ma, orco);
	}
	else {

		for(a=0; a<totvert; a++, mvert++) {
			ver= RE_findOrAddVert(R.totvert++);
			VECCOPY(ver->co, mvert->co);
			MTC_Mat4MulVecfl(mat, ver->co);

			if(orco) {
				ver->orco= orco;
				orco+=3;
			}
			if(ms) {
				float *sticky= RE_vertren_get_sticky(ver, 1);
				sticky[0]= ms->co[0];
				sticky[1]= ms->co[1];
				ms++;
			}
		}
		/* still to do for keys: the correct local texture coordinate */

		/* faces in order of color blocks */
		vertofs= R.totvert - totvert;
		for(a1=0; (a1<ob->totcol || (a1==0 && ob->totcol==0)); a1++) {

			ma= give_render_material(ob, a1+1);
			
			/* test for 100% transparant */
			ok= 1;
			if(ma->alpha==0.0 && ma->spectra==0.0) {
				ok= 0;
				/* texture on transparency? */
				for(a=0; a<MAX_MTEX; a++) {
					if(ma->mtex[a] && ma->mtex[a]->tex) {
						if(ma->mtex[a]->mapto & MAP_ALPHA) ok= 1;
					}
				}
			}
			
			/* if wire material, and we got edges, don't do the faces */
			if(ma->mode & MA_WIRE) {
				end= dlm?dlm->totedge:me->totedge;
				if(end) ok= 0;
			}

			if(ok) {
				TFace *tface= NULL;

				/* radio faces need autosmooth, to separate shared vertices in corners */
				if(R.r.mode & R_RADIO)
					if(ma->mode & MA_RADIO) 
						do_autosmooth= 1;
				
				end= dlm?dlm->totface:me->totface;
				if (dlm) {
					mface= dlm->mface;
					if (dlm->tface) {
						tface= dlm->tface;
						vertcol= NULL;
					} else if (dlm->mcol) {
						vertcol= (unsigned int *)dlm->mcol;
					} else {
						vertcol= NULL;
					}
				} else {
					mface= me->mface;
					if (me->tface) {
						tface= me->tface;
						vertcol= NULL;
					} else if (me->mcol) {
						vertcol= (unsigned int *)me->mcol;
					} else {
						vertcol= NULL;
					}
				}

				for(a=0; a<end; a++) {
					int v1, v2, v3, v4, flag;
					
					if( mface->mat_nr==a1 ) {
						float len;
							
						v1= mface->v1;
						v2= mface->v2;
						v3= mface->v3;
						v4= mface->v4;
						flag= mface->flag & ME_SMOOTH;
						
						vlr= RE_findOrAddVlak(R.totvlak++);
						vlr->ob= vlr_set_ob(ob);
						vlr->v1= RE_findOrAddVert(vertofs+v1);
						vlr->v2= RE_findOrAddVert(vertofs+v2);
						vlr->v3= RE_findOrAddVert(vertofs+v3);
						if(v4) vlr->v4= RE_findOrAddVert(vertofs+v4);
						else vlr->v4= 0;

						/* render normals are inverted in render */
						if(vlr->v4) len= CalcNormFloat4(vlr->v4->co, vlr->v3->co, vlr->v2->co,
							vlr->v1->co, vlr->n);
						else len= CalcNormFloat(vlr->v3->co, vlr->v2->co, vlr->v1->co,
							vlr->n);

						vlr->mat= ma;
						vlr->flag= flag;
						if((me->flag & ME_NOPUNOFLIP) ) {
							vlr->flag |= R_NOPUNOFLIP;
						}
						vlr->ec= 0; /* mesh edges rendered separately */
						vlr->lay= ob->lay;

						if(len==0) R.totvlak--;
						else {
							if(dlm) {
								if(tface) {
									vlr->tface= BLI_memarena_alloc(R.memArena, sizeof(TFace));
									vlr->vcol= vlr->tface->col;
									memcpy(vlr->tface, tface, sizeof(TFace));
								} 
								else if (vertcol) {
									vlr->vcol= BLI_memarena_alloc(R.memArena, sizeof(int)*4);
									memcpy(vlr->vcol, vertcol+4*a, sizeof(int)*4);
								}
							} else {
								if(tface) {
									vlr->vcol= tface->col;
									vlr->tface= tface;
								} 
								else if (vertcol) {
									vlr->vcol= vertcol+4*a;
								}
							}
						}
					}

					mface++;
					if(tface) tface++;
				}
			}
		}
		
		/* exception... we do edges for wire mode. potential conflict when faces exist... */
		end= dlm?dlm->totedge:me->totedge;
		mvert= dlm?dlm->mvert:me->mvert;
		ma= give_render_material(ob, 1);
		if(end && (ma->mode & MA_WIRE)) {
			MEdge *medge;
			struct edgesort *edgetable;
			int totedge;
			
			medge= dlm?dlm->medge:me->medge;
			
			/* we want edges to have UV and vcol too... */
			edgetable= make_mesh_edge_lookup(me, dlm, &totedge);
			
			for(a1=0; a1<end; a1++, medge++) {
				if (medge->flag&ME_EDGERENDER) {
					MVert *v0 = &mvert[medge->v1];
					MVert *v1 = &mvert[medge->v2];

					vlr= RE_findOrAddVlak(R.totvlak++);
					vlr->ob= vlr_set_ob(ob);
					vlr->v1= RE_findOrAddVert(vertofs+medge->v1);
					vlr->v2= RE_findOrAddVert(vertofs+medge->v2);
					vlr->v3= vlr->v2;
					vlr->v4= NULL;
					
					if(edgetable) {
						use_mesh_edge_lookup(me, dlm, medge, vlr, edgetable, totedge);
					}
					
					xn= (v0->no[0]+v1->no[0]);
					yn= (v0->no[1]+v1->no[1]);
					zn= (v0->no[2]+v1->no[2]);
					/* transpose ! */
					vlr->n[0]= imat[0][0]*xn+imat[0][1]*yn+imat[0][2]*zn;
					vlr->n[1]= imat[1][0]*xn+imat[1][1]*yn+imat[1][2]*zn;
					vlr->n[2]= imat[2][0]*xn+imat[2][1]*yn+imat[2][2]*zn;
					Normalise(vlr->n);
					
					vlr->mat= ma;
					vlr->flag= 0;
					vlr->ec= ME_V1V2;
					vlr->lay= ob->lay;
				}
			}
			if(edgetable)
				MEM_freeN(edgetable);
		}
	}
	
	if (test_for_displace( ob ) ) {
		calc_vertexnormals(totverto, totvlako, 0);
		do_displacement(ob, totvlako, R.totvlak-totvlako, totverto, R.totvert-totverto);
	}

	if(do_autosmooth || (me->flag & ME_AUTOSMOOTH)) {
		autosmooth(totverto, totvlako, me->smoothresh);
	}

	calc_vertexnormals(totverto, totvlako, need_tangent);

	if(need_stress)
		calc_edge_stress(me, totverto, totvlako);
	
	if(dlm) displistmesh_free(dlm);
	if(dm_needsfree) dm->release(dm);
}

/* ------------------------------------------------------------------------- */

static void area_lamp_vectors(LampRen *lar)
{
	float xsize= 0.5*lar->area_size, ysize= 0.5*lar->area_sizey;

	/* corner vectors */
	lar->area[0][0]= lar->co[0] - xsize*lar->mat[0][0] - ysize*lar->mat[1][0];
	lar->area[0][1]= lar->co[1] - xsize*lar->mat[0][1] - ysize*lar->mat[1][1];
	lar->area[0][2]= lar->co[2] - xsize*lar->mat[0][2] - ysize*lar->mat[1][2];	

	/* corner vectors */
	lar->area[1][0]= lar->co[0] - xsize*lar->mat[0][0] + ysize*lar->mat[1][0];
	lar->area[1][1]= lar->co[1] - xsize*lar->mat[0][1] + ysize*lar->mat[1][1];
	lar->area[1][2]= lar->co[2] - xsize*lar->mat[0][2] + ysize*lar->mat[1][2];	

	/* corner vectors */
	lar->area[2][0]= lar->co[0] + xsize*lar->mat[0][0] + ysize*lar->mat[1][0];
	lar->area[2][1]= lar->co[1] + xsize*lar->mat[0][1] + ysize*lar->mat[1][1];
	lar->area[2][2]= lar->co[2] + xsize*lar->mat[0][2] + ysize*lar->mat[1][2];	

	/* corner vectors */
	lar->area[3][0]= lar->co[0] + xsize*lar->mat[0][0] - ysize*lar->mat[1][0];
	lar->area[3][1]= lar->co[1] + xsize*lar->mat[0][1] - ysize*lar->mat[1][1];
	lar->area[3][2]= lar->co[2] + xsize*lar->mat[0][2] - ysize*lar->mat[1][2];	
	/* only for correction button size, matrix size works on energy */
	lar->areasize= lar->dist*lar->dist/(4.0*xsize*ysize);
}

/* If lar takes more lamp data, the decoupling will be better. */
void RE_add_render_lamp(Object *ob, int actual_render)
{
	Lamp *la= ob->data;
	LampRen *lar;
	GroupObject *go;
	float mat[4][4], hoek, xn, yn;
	int c;

	/* prevent only shadow from rendering light, but only return on render, not preview */
	if(actual_render) {
		if(la->mode & LA_ONLYSHADOW)
			if((R.r.mode & R_SHADOW)==0)
				return;
	}
	
	go= MEM_callocN(sizeof(GroupObject), "groupobject");
	BLI_addtail(&R.lights, go);
	R.totlamp++;
	lar= (LampRen *)MEM_callocN(sizeof(LampRen),"lampren");
	go->lampren= lar;
	go->ob= ob;

	MTC_Mat4MulMat4(mat, ob->obmat, R.viewmat);
	MTC_Mat4Invert(ob->imat, mat);

	MTC_Mat3CpyMat4(lar->mat, mat);
	MTC_Mat3CpyMat4(lar->imat, ob->imat);

	lar->bufsize = la->bufsize;
	lar->samp = la->samp;
	lar->soft = la->soft;
	lar->shadhalostep = la->shadhalostep;
	lar->clipsta = la->clipsta;
	lar->clipend = la->clipend;
	lar->bias = la->bias;

	lar->type= la->type;
	lar->mode= la->mode;

	lar->energy= la->energy;
	lar->energy= la->energy;
	if(la->mode & LA_NEG) lar->energy= -lar->energy;

	lar->vec[0]= -mat[2][0];
	lar->vec[1]= -mat[2][1];
	lar->vec[2]= -mat[2][2];
	Normalise(lar->vec);
	lar->co[0]= mat[3][0];
	lar->co[1]= mat[3][1];
	lar->co[2]= mat[3][2];
	lar->dist= la->dist;
	lar->haint= la->haint;
	lar->distkw= lar->dist*lar->dist;
	lar->r= lar->energy*la->r;
	lar->g= lar->energy*la->g;
	lar->b= lar->energy*la->b;
	lar->k= la->k;

	// area
	lar->ray_samp= la->ray_samp;
	lar->ray_sampy= la->ray_sampy;
	lar->ray_sampz= la->ray_sampz;

	lar->area_size= la->area_size;
	lar->area_sizey= la->area_sizey;
	lar->area_sizez= la->area_sizez;

	lar->area_shape= la->area_shape;
	lar->ray_samp_type= la->ray_samp_type;

	if(lar->type==LA_AREA) {
		switch(lar->area_shape) {
		case LA_AREA_SQUARE:
			lar->ray_totsamp= lar->ray_samp*lar->ray_samp;
			lar->ray_sampy= lar->ray_samp;
			lar->area_sizey= lar->area_size;
			break;
		case LA_AREA_RECT:
			lar->ray_totsamp= lar->ray_samp*lar->ray_sampy;
			break;
		case LA_AREA_CUBE:
			lar->ray_totsamp= lar->ray_samp*lar->ray_samp*lar->ray_samp;
			lar->ray_sampy= lar->ray_samp;
			lar->ray_sampz= lar->ray_samp;
			lar->area_sizey= lar->area_size;
			lar->area_sizez= lar->area_size;
			break;
		case LA_AREA_BOX:
			lar->ray_totsamp= lar->ray_samp*lar->ray_sampy*lar->ray_sampz;
			break;
		}

		area_lamp_vectors(lar);
	}
	else lar->ray_totsamp= 0;
	
	/* yafray: photonlight and other params */
	if (R.r.renderer==R_YAFRAY) {
		lar->YF_numphotons = la->YF_numphotons;
		lar->YF_numsearch = la->YF_numsearch;
		lar->YF_phdepth = la->YF_phdepth;
		lar->YF_useqmc = la->YF_useqmc;
		lar->YF_causticblur = la->YF_causticblur;
		lar->YF_ltradius = la->YF_ltradius;
		lar->YF_bufsize = la->YF_bufsize;
		lar->YF_glowint = la->YF_glowint;
		lar->YF_glowofs = la->YF_glowofs;
		lar->YF_glowtype = la->YF_glowtype;
	}

	lar->spotsi= la->spotsize;
	if(lar->mode & LA_HALO) {
		if(lar->spotsi>170.0) lar->spotsi= 170.0;
	}
	lar->spotsi= cos( M_PI*lar->spotsi/360.0 );
	lar->spotbl= (1.0-lar->spotsi)*la->spotblend;

	memcpy(lar->mtex, la->mtex, MAX_MTEX*sizeof(void *));

	lar->lay= ob->lay & 0xFFFFFF;	// higher 8 bits are localview layers

	lar->ld1= la->att1;
	lar->ld2= la->att2;

	if(lar->type==LA_SPOT) {

		Normalise(lar->imat[0]);
		Normalise(lar->imat[1]);
		Normalise(lar->imat[2]);

		xn= saacos(lar->spotsi);
		xn= sin(xn)/cos(xn);
		lar->spottexfac= 1.0/(xn);

		if(lar->mode & LA_ONLYSHADOW) {
			if((lar->mode & (LA_SHAD|LA_SHAD_RAY))==0) lar->mode -= LA_ONLYSHADOW;
		}

	}

	/* set flag for spothalo en initvars */
	if(la->type==LA_SPOT && (la->mode & LA_HALO)) {
		if(la->haint>0.0) {
			R.flag |= R_LAMPHALO;

			/* camera position (0,0,0) rotate around lamp */
			lar->sh_invcampos[0]= -lar->co[0];
			lar->sh_invcampos[1]= -lar->co[1];
			lar->sh_invcampos[2]= -lar->co[2];
			MTC_Mat3MulVecfl(lar->imat, lar->sh_invcampos);

			/* z factor, for a normalized volume */
			hoek= saacos(lar->spotsi);
			xn= lar->spotsi;
			yn= sin(hoek);
			lar->sh_zfac= yn/xn;
			/* pre-scale */
			lar->sh_invcampos[2]*= lar->sh_zfac;

		}
	}

	for(c=0; c<MAX_MTEX; c++) {
		if(la->mtex[c] && la->mtex[c]->tex) {
			lar->mode |= LA_TEXTURE;

			if(R.flag & R_RENDERING) {
				if(R.osa) {
					if(la->mtex[c]->tex->type==TEX_IMAGE) lar->mode |= LA_OSATEX;
				}
			}
		}
	}

	/* yafray: shadowbuffers and jitter only needed for internal render */
	if (actual_render && R.r.renderer==R_INTERN) {
		if(R.r.mode & R_SHADOW) {
			if (la->type==LA_SPOT && (lar->mode & LA_SHAD) ) {
				/* Per lamp, one shadow buffer is made. */
				Mat4CpyMat4(mat, ob->obmat);
				RE_initshadowbuf(lar, mat);	// mat is altered
			}
			else if(la->type==LA_AREA && (lar->mode & LA_SHAD_RAY) ) {
				init_jitter_plane(lar);
			}
		}
	}
	
	/* yafray: shadow flag should not be cleared, only used with internal renderer */
	if (R.r.renderer==R_INTERN) {
		/* to make sure we can check ray shadow easily in the render code */
		if(lar->mode & LA_SHAD_RAY) {
			if( (R.r.mode & R_RAYTRACE)==0)
				lar->mode &= ~LA_SHAD_RAY;
		}
	}
}

/* ------------------------------------------------------------------------- */
static void init_render_surf(Object *ob)
{
	extern Material defmaterial;	// initrender.c
	Nurb *nu=0;
	Curve *cu;
	ListBase displist;
	DispList *dl;
	VertRen *ver, *v1, *v2, *v3, *v4;
	VlakRen *vlr;
	Material *matar[32];
	float *data, *orco=NULL, *orcobase=NULL, n1[3], flen, mat[4][4];
	int a, need_orco=0, startvlak, startvert, p1, p2, p3, p4;
	int u, v;
	int sizeu, sizev;
	VlakRen *vlr1, *vlr2, *vlr3;
	float  vn[3]; // n2[3],

	cu= ob->data;
	nu= cu->nurb.first;
	if(nu==0) return;

	MTC_Mat4MulMat4(mat, ob->obmat, R.viewmat);
	MTC_Mat4Invert(ob->imat, mat);

	/* material array */
	memset(matar, 0, 4*32);
	matar[0]= &defmaterial;
	for(a=0; a<ob->totcol; a++) {
		matar[a]= give_render_material(ob, a+1);
		if(matar[a] && matar[a]->texco & TEXCO_ORCO) {
			need_orco= 1;
		}
	}

	if(ob->parent && (ob->parent->type==OB_LATTICE)) need_orco= 1;

	if(need_orco) orcobase= orco= get_object_orco(ob);

	displist.first= displist.last= 0;
	makeDispListSurf(ob, &displist, 1);

	dl= displist.first;
	/* walk along displaylist and create rendervertices/-faces */
	while(dl) {
			/* watch out: u ^= y, v ^= x !! */
		if(dl->type==DL_SURF) {
			int nsizeu, nsizev;

			startvert= R.totvert;
			nsizeu = sizeu = dl->parts; nsizev = sizev = dl->nr; 

			data= dl->verts;
			for (u = 0; u < sizeu; u++) {
				v1 = RE_findOrAddVert(R.totvert++); /* save this for possible V wrapping */
				VECCOPY(v1->co, data); data += 3;
				if(orco) {
					v1->orco= orco; orco+= 3;
				}	
				MTC_Mat4MulVecfl(mat, v1->co);

				for (v = 1; v < sizev; v++) {
					ver= RE_findOrAddVert(R.totvert++);
					VECCOPY(ver->co, data); data += 3;
					if(orco) {
						ver->orco= orco; orco+= 3;
					}	
					MTC_Mat4MulVecfl(mat, ver->co);
				}
				/* if V-cyclic, add extra vertices at end of the row */
				if (dl->flag & DL_CYCL_U) {
					ver= RE_findOrAddVert(R.totvert++);
					VECCOPY(ver->co, v1->co);
					if(orco) {
						ver->orco= orcobase + 3*(u*sizev + 0);
					}
				}	
			}	

				/* Done before next loop to get corner vert */
			if (dl->flag & DL_CYCL_U) nsizev++;
			if (dl->flag & DL_CYCL_V) nsizeu++;

			/* if U cyclic, add extra row at end of column */
			if (dl->flag & DL_CYCL_V) {
				for (v = 0; v < nsizev; v++) {
					v1= RE_findOrAddVert(startvert + v);
					ver= RE_findOrAddVert(R.totvert++);
					VECCOPY(ver->co, v1->co);
					if(orco) {
						ver->orco= orcobase + 3*(0*sizev + v);
					}
				}
			}
			
			sizeu = nsizeu;
			sizev = nsizev;

			startvlak= R.totvlak;

			for(u = 0; u < sizeu - 1; u++) {
				p1 = startvert + u * sizev; /* walk through face list */
				p2 = p1 + 1;
				p3 = p2 + sizev;
				p4 = p3 - 1;

				for(v = 0; v < sizev - 1; v++) {
					v1= RE_findOrAddVert(p1);
					v2= RE_findOrAddVert(p2);
					v3= RE_findOrAddVert(p3);
					v4= RE_findOrAddVert(p4);

					vlr= RE_findOrAddVlak(R.totvlak++);
					vlr->ob= vlr_set_ob(ob);
					vlr->v1= v1; vlr->v2= v2; vlr->v3= v3; vlr->v4= v4;
					
					flen= CalcNormFloat4(vlr->v4->co, vlr->v3->co, vlr->v2->co, vlr->v1->co, n1);
					VECCOPY(vlr->n, n1);
					
					vlr->lay= ob->lay;
					vlr->mat= matar[ dl->col];
					vlr->ec= ME_V1V2+ME_V2V3;
					vlr->flag= dl->rt;
					if( (cu->flag & CU_NOPUNOFLIP) ) {
						vlr->flag |= R_NOPUNOFLIP;
					}

					VecAddf(v1->n, v1->n, n1);
					VecAddf(v2->n, v2->n, n1);
					VecAddf(v3->n, v3->n, n1);
					VecAddf(v4->n, v4->n, n1);

					p1++; p2++; p3++; p4++;
				}
			}	
			/* fix normals for U resp. V cyclic faces */
			sizeu--; sizev--;  /* dec size for face array */
			if (dl->flag & DL_CYCL_V) {

				for (v = 0; v < sizev; v++)
				{
					/* optimize! :*/
					vlr= RE_findOrAddVlak(UVTOINDEX(sizeu - 1, v));
					vlr1= RE_findOrAddVlak(UVTOINDEX(0, v));
					VecAddf(vlr1->v1->n, vlr1->v1->n, vlr->n);
					VecAddf(vlr1->v2->n, vlr1->v2->n, vlr->n);
					VecAddf(vlr->v3->n, vlr->v3->n, vlr1->n);
					VecAddf(vlr->v4->n, vlr->v4->n, vlr1->n);
				}
			}
			if (dl->flag & DL_CYCL_U) {

				for (u = 0; u < sizeu; u++)
				{
					/* optimize! :*/
					vlr= RE_findOrAddVlak(UVTOINDEX(u, 0));
					vlr1= RE_findOrAddVlak(UVTOINDEX(u, sizev-1));
					VecAddf(vlr1->v2->n, vlr1->v2->n, vlr->n);
					VecAddf(vlr1->v3->n, vlr1->v3->n, vlr->n);
					VecAddf(vlr->v1->n, vlr->v1->n, vlr1->n);
					VecAddf(vlr->v4->n, vlr->v4->n, vlr1->n);
				}
			}
			/* last vertex is an extra case: 

			^	()----()----()----()
			|	|     |     ||     |
			u	|     |(0,n)||(0,0)|
				|     |     ||     |
			 	()====()====[]====()
			 	|     |     ||     |
			 	|     |(m,n)||(m,0)|
				|     |     ||     |
				()----()----()----()
				       v ->

			vertex [] is no longer shared, therefore distribute
			normals of the surrounding faces to all of the duplicates of []
			*/

			if ((dl->flag & DL_CYCL_V) && (dl->flag & DL_CYCL_U))
			{
				vlr= RE_findOrAddVlak(UVTOINDEX(sizeu - 1, sizev - 1)); /* (m,n) */
				vlr1= RE_findOrAddVlak(UVTOINDEX(0,0));  /* (0,0) */
				VecAddf(vn, vlr->n, vlr1->n);
				vlr2= RE_findOrAddVlak(UVTOINDEX(0, sizev-1)); /* (0,n) */
				VecAddf(vn, vn, vlr2->n);
				vlr3= RE_findOrAddVlak(UVTOINDEX(sizeu-1, 0)); /* (m,0) */
				VecAddf(vn, vn, vlr3->n);
				VECCOPY(vlr->v3->n, vn);
				VECCOPY(vlr1->v1->n, vn);
				VECCOPY(vlr2->v2->n, vn);
				VECCOPY(vlr3->v4->n, vn);
			}
			for(a = startvert; a < R.totvert; a++) {
				ver= RE_findOrAddVert(a);
				Normalise(ver->n);
			}


		}

		dl= dl->next;
	}
	freedisplist(&displist);
}

static void init_render_curve(Object *ob)
{
	extern Material defmaterial;	// initrender.c
	Curve *cu;
	VertRen *ver;
	VlakRen *vlr;
	DispList *dl;
	Material *matar[32];
	float len, *data, *fp, *orco=NULL;
	float n[3], mat[4][4];
	int nr, startvert, startvlak, a, b;
	int frontside, need_orco=0;

	cu= ob->data;
	if(cu->nurb.first==NULL) return;

	/* no modifier call here, is in makedisp */

	/* test displist */
	if(cu->disp.first==0) makeDispListCurveTypes(ob, 0);
	dl= cu->disp.first;
	if(cu->disp.first==0) return;

	MTC_Mat4MulMat4(mat, ob->obmat, R.viewmat);
	MTC_Mat4Invert(ob->imat, mat);

	/* material array */
	memset(matar, 0, 4*32);
	matar[0]= &defmaterial;
	for(a=0; a<ob->totcol; a++) {
		matar[a]= give_render_material(ob, a+1);
		if(matar[a]->texco & TEXCO_ORCO) {
			need_orco= 1;
		}
	}

	if(need_orco) orco= get_object_orco(ob);

	dl= cu->disp.first;
	while(dl) {
		if(dl->type==DL_INDEX3) {
			int *index;

			startvert= R.totvert;
			data= dl->verts;

			n[0]= ob->imat[0][2];
			n[1]= ob->imat[1][2];
			n[2]= ob->imat[2][2];
			Normalise(n);

			/* copy first, rotate later for comparision trick */
			for(a=0; a<dl->nr; a++, data+=3) {
				ver= RE_findOrAddVert(R.totvert++);
				VECCOPY(ver->co, data);
				MTC_Mat4MulVecfl(mat, ver->co);

				if(ver->co[2] < 0.0) {
					VECCOPY(ver->n, n);
					ver->flag = 1;
				}
				else {
					ver->n[0]= -n[0]; ver->n[1]= -n[1]; ver->n[2]= -n[2];
					ver->flag = 0;
				}

				if (orco) {
					ver->orco = orco;
					orco += 3;
				}
			}

			startvlak= R.totvlak;
			index= dl->index;
			for(a=0; a<dl->parts; a++, index+=3) {

				vlr= RE_findOrAddVlak(R.totvlak++);
				vlr->ob = vlr_set_ob(ob);	/* yafray: correction for curve rendering, obptr was not set */
				vlr->v1= RE_findOrAddVert(startvert+index[0]);
				vlr->v2= RE_findOrAddVert(startvert+index[1]);
				vlr->v3= RE_findOrAddVert(startvert+index[2]);
				vlr->v4= NULL;
				
				if(vlr->v1->flag) {
					VECCOPY(vlr->n, n);
				}
				else {
					vlr->n[0]= -n[0]; vlr->n[1]= -n[1]; vlr->n[2]= -n[2];
				}
				
				vlr->mat= matar[ dl->col ];
				vlr->flag= 0;
				if( (cu->flag & CU_NOPUNOFLIP) ) {
					vlr->flag |= R_NOPUNOFLIP;
				}
				vlr->ec= 0;
				vlr->lay= ob->lay;
			}
		}
		else if (dl->type==DL_SURF) {
			int p1,p2,p3,p4;
			float *surforco = orco;

			fp= dl->verts;
			startvert= R.totvert;
			nr= dl->nr*dl->parts;

			while(nr--) {
				ver= RE_findOrAddVert(R.totvert++);
					
				VECCOPY(ver->co, fp);
				MTC_Mat4MulVecfl(mat, ver->co);
				fp+= 3;

				if (orco) {
					ver->orco = orco;
					orco += 3;
				}
			}

			startvlak= R.totvlak;

			for(a=0; a<dl->parts; a++) {

				frontside= (a >= dl->nr/2);

				DL_SURFINDEX(dl->flag & DL_CYCL_U, dl->flag & DL_CYCL_V, dl->nr, dl->parts);
				p1+= startvert;
				p2+= startvert;
				p3+= startvert;
				p4+= startvert;

				for(; b<dl->nr; b++) {
					vlr= RE_findOrAddVlak(R.totvlak++);
					vlr->ob= vlr_set_ob(ob);
					vlr->v1= RE_findOrAddVert(p2);
					vlr->v2= RE_findOrAddVert(p1);
					vlr->v3= RE_findOrAddVert(p3);
					vlr->v4= RE_findOrAddVert(p4);
					vlr->ec= ME_V2V3+ME_V3V4;
					if(a==0) vlr->ec+= ME_V1V2;

					vlr->flag= dl->rt;
					vlr->lay= ob->lay;

					/* this is not really scientific: the vertices
						* 2, 3 en 4 seem to give better vertexnormals than 1 2 3:
						* front and backside treated different!!
						*/

					if(frontside)
						CalcNormFloat(vlr->v2->co, vlr->v3->co, vlr->v4->co, vlr->n);
					else 
						CalcNormFloat(vlr->v1->co, vlr->v2->co, vlr->v3->co, vlr->n);

					vlr->mat= matar[ dl->col ];

					p4= p3;
					p3++;
					p2= p1;
					p1++;
				}
			}

			if (dl->bevelSplitFlag) {
				for(a=0; a<dl->parts-1+!!(dl->flag&DL_CYCL_V); a++)
					if(dl->bevelSplitFlag[a>>5]&(1<<(a&0x1F)))
						split_v_renderfaces(startvlak, startvert, dl->parts, dl->nr, a, dl->flag&DL_CYCL_V, dl->flag&DL_CYCL_U);
			}

			/* vertex normals */
			for(a= startvlak; a<R.totvlak; a++) {
				vlr= RE_findOrAddVlak(a);

				VecAddf(vlr->v1->n, vlr->v1->n, vlr->n);
				VecAddf(vlr->v3->n, vlr->v3->n, vlr->n);
				VecAddf(vlr->v2->n, vlr->v2->n, vlr->n);
				VecAddf(vlr->v4->n, vlr->v4->n, vlr->n);
			}
			for(a=startvert; a<R.totvert; a++) {
				ver= RE_findOrAddVert(a);
				len= Normalise(ver->n);
				if(len==0.0) ver->flag= 1;	/* flag use, its only used in zbuf now  */
				else ver->flag= 0;
			}
			for(a= startvlak; a<R.totvlak; a++) {
				vlr= RE_findOrAddVlak(a);
				if(vlr->v1->flag) VECCOPY(vlr->v1->n, vlr->n);
				if(vlr->v2->flag) VECCOPY(vlr->v2->n, vlr->n);
				if(vlr->v3->flag) VECCOPY(vlr->v3->n, vlr->n);
				if(vlr->v4->flag) VECCOPY(vlr->v4->n, vlr->n);
			}
		}

		dl= dl->next;
	}
}

/* prevent phong interpolation for giving ray shadow errors (terminator problem) */
static void set_phong_threshold(Object *ob, int startface, int numface, int startvert, int numvert )
{
//	VertRen *ver;
	VlakRen *vlr;
	float thresh= 0.0, dot;
	int tot=0, i;
	
	/* Added check for 'pointy' situations, only dotproducts of 0.9 and larger 
	   are taken into account. This threshold is meant to work on smooth geometry, not
	   for extreme cases (ton) */
	
	for(i=startface; i<startface+numface; i++) {
		vlr= RE_findOrAddVlak(i);
		if(vlr->flag & R_SMOOTH) {
			dot= INPR(vlr->n, vlr->v1->n);
			dot= ABS(dot);
			if(dot>0.9) {
				thresh+= dot; tot++;
			}
			dot= INPR(vlr->n, vlr->v2->n);
			dot= ABS(dot);
			if(dot>0.9) {
				thresh+= dot; tot++;
			}

			dot= INPR(vlr->n, vlr->v3->n);
			dot= ABS(dot);
			if(dot>0.9) {
				thresh+= dot; tot++;
			}

			if(vlr->v4) {
				dot= INPR(vlr->n, vlr->v4->n);
				dot= ABS(dot);
				if(dot>0.9) {
					thresh+= dot; tot++;
				}
			}
		}
	}
	
	if(tot) {
		thresh/= (float)tot;
		ob->smoothresh= cos(0.5*M_PI-acos(thresh));
	}
}

static void init_render_object(Object *ob)
{
	float mat[4][4];
	int startface, startvert;
	
	startface=R.totvlak;
	startvert=R.totvert;

	ob->flag |= OB_DONE;

	if(ob->type==OB_LAMP)
		RE_add_render_lamp(ob, 1);
	else if ELEM(ob->type, OB_FONT, OB_CURVE)
		init_render_curve(ob);
	else if(ob->type==OB_SURF)
		init_render_surf(ob);
	else if(ob->type==OB_MESH)
		init_render_mesh(ob);
	else if(ob->type==OB_MBALL)
		init_render_mball(ob);
	else {
		MTC_Mat4MulMat4(mat, ob->obmat, R.viewmat);
		MTC_Mat4Invert(ob->imat, mat);
	}
	
	/* generic post process here */
	if(startvert!=R.totvert) {
	
		/* the exception below is because displace code now is in init_render_mesh call, 
		I will look at means to have autosmooth enabled for all object types 
		and have it as general postprocess, like displace */
		if (ob->type!=OB_MESH && test_for_displace( ob ) ) 
			do_displacement(ob, startface, R.totvlak-startface, startvert, R.totvert-startvert);
	
		/* phong normal interpolation can cause error in tracing (terminator prob) */
		ob->smoothresh= 0.0;
		if( (R.r.mode & R_RAYTRACE) && (R.r.mode & R_SHADOW) ) 
			set_phong_threshold(ob, startface, R.totvlak-startface, startvert, R.totvert-startvert);
	}
}

void RE_freeRotateBlenderScene(void)
{
	ShadBuf *shb;
	Object *ob = NULL;
	GroupObject *go;
	unsigned long *ztile;
	int a, b, v;
	char *ctile;

	/* FREE */
	
	BLI_memarena_free(R.memArena);
	R.memArena = NULL;
	
	for(go= R.lights.first; go; go= go->next) {
		struct LampRen *lar= go->lampren;
		if(lar->shb) {
			shb= lar->shb;
			v= (shb->size*shb->size)/256;
			ztile= shb->zbuf;
			ctile= shb->cbuf;
			for(b=0; b<v; b++, ztile++, ctile++) {
				if(*ctile) MEM_freeN((void *) *ztile);
			}
			
			MEM_freeN(shb->zbuf);
			MEM_freeN(shb->cbuf);
			MEM_freeN(lar->shb);
		}
		if(lar->jitter) MEM_freeN(lar->jitter);
		MEM_freeN(lar);
	}
	
	BLI_freelistN(&R.lights);

	/* note; these pointer arrays were allocated, with last element NULL to stop loop */
	RE_free_vertex_tables();
	
	a=0;
	while(R.blovl[a]) {
		MEM_freeN(R.blovl[a]);
		R.blovl[a]= NULL;
		a++;
	}
	a=0;
	while(R.bloha[a]) {
		MEM_freeN(R.bloha[a]);
		R.bloha[a]= NULL;
		a++;
	}

	/* free orco. check all objects because of duplis and sets */
	ob= G.main->object.first;
	while(ob) {
		if(ob->type==OB_MBALL) {
			if(ob->disp.first && ob->disp.first!=ob->disp.last) {
				DispList *dl= ob->disp.first;
				BLI_remlink(&ob->disp, dl);
				freedisplist(&ob->disp);
				BLI_addtail(&ob->disp, dl);
			}
		}
		ob= ob->id.next;
	}

	free_mesh_orco_hash();

	end_radio_render();
	end_render_materials();
	
	if(R.wrld.aosphere) {
		MEM_freeN(R.wrld.aosphere);
		R.wrld.aosphere= NULL;
		G.scene->world->aosphere= NULL;
	}
	
	R.totvlak=R.totvert=R.totlamp=R.tothalo= 0;
}

/* per face check if all samples should be taken.
   if raytrace, do always for raytraced material, or when material full_osa set */
static void set_fullsample_flag(void)
{
	VlakRen *vlr;
	int a, trace;

	trace= R.r.mode & R_RAYTRACE;
	
	for(a=R.totvlak-1; a>=0; a--) {
		vlr= RE_findOrAddVlak(a);
		
		if(vlr->mat->mode & MA_FULL_OSA) vlr->flag |= R_FULL_OSA;
		else if(trace) {
			if(vlr->mat->mode & MA_SHLESS);
			else if(vlr->mat->mode & (MA_RAYTRANSP|MA_RAYMIRROR|MA_SHADOW))
				vlr->flag |= R_FULL_OSA;
		}
	}
}

/* 10 times larger than normal epsilon, test it on default nurbs sphere with ray_transp */
#ifdef FLT_EPSILON
#undef FLT_EPSILON
#endif
#define FLT_EPSILON 1.19209290e-06F


static void check_non_flat_quads(void)
{
	VlakRen *vlr, *vlr1;
	VertRen *v1, *v2, *v3, *v4;
	float nor[3], xn, flen;
	int a;

	for(a=R.totvlak-1; a>=0; a--) {
		vlr= RE_findOrAddVlak(a);
		
		/* test if rendering as a quad or triangle, skip wire */
		if(vlr->v4 && (vlr->flag & R_STRAND)==0 && (vlr->mat->mode & MA_WIRE)==0) {
			
			/* check if quad is actually triangle */
			v1= vlr->v1;
			v2= vlr->v2;
			v3= vlr->v3;
			v4= vlr->v4;
			VECSUB(nor, v1->co, v2->co);
			if( ABS(nor[0])<FLT_EPSILON &&  ABS(nor[1])<FLT_EPSILON && ABS(nor[2])<FLT_EPSILON ) {
				vlr->v1= v2;
				vlr->v2= v3;
				vlr->v3= v4;
				vlr->v4= NULL;
			}
			else {
				VECSUB(nor, v2->co, v3->co);
				if( ABS(nor[0])<FLT_EPSILON &&  ABS(nor[1])<FLT_EPSILON && ABS(nor[2])<FLT_EPSILON ) {
					vlr->v2= v3;
					vlr->v3= v4;
					vlr->v4= NULL;
				}
				else {
					VECSUB(nor, v3->co, v4->co);
					if( ABS(nor[0])<FLT_EPSILON &&  ABS(nor[1])<FLT_EPSILON && ABS(nor[2])<FLT_EPSILON ) {
						vlr->v4= NULL;
					}
					else {
						VECSUB(nor, v4->co, v1->co);
						if( ABS(nor[0])<FLT_EPSILON &&  ABS(nor[1])<FLT_EPSILON && ABS(nor[2])<FLT_EPSILON ) {
							vlr->v4= NULL;
						}
					}
				}
			}
			
			if(vlr->v4) {
				
				/* Face is divided along edge with the least gradient 		*/
				/* Flagged with R_DIVIDE_24 if divide is from vert 2 to 4 	*/
				/* 		4---3		4---3 */
				/*		|\ 1|	or  |1 /| */
				/*		|0\ |		|/ 0| */
				/*		1---2		1---2 	0 = orig face, 1 = new face */
				
				/* render normals are inverted in render! we calculate normal of single tria here */
				flen= CalcNormFloat(vlr->v4->co, vlr->v3->co, vlr->v1->co, nor);
				if(flen==0.0) CalcNormFloat(vlr->v4->co, vlr->v2->co, vlr->v1->co, nor);
				
				xn= nor[0]*vlr->n[0] + nor[1]*vlr->n[1] + nor[2]*vlr->n[2];
				if(ABS(xn) < 0.99995 ) {	// checked on noisy fractal grid
					float d1, d2;
					
					vlr1= RE_findOrAddVlak(R.totvlak++);
					*vlr1= *vlr;
					vlr1->flag |= R_FACE_SPLIT;
					
					/* split direction based on vnorms */
					CalcNormFloat(vlr->v1->co, vlr->v2->co, vlr->v3->co, nor);
					d1= nor[0]*vlr->v1->n[0] + nor[1]*vlr->v1->n[1] + nor[2]*vlr->v1->n[2];

					CalcNormFloat(vlr->v2->co, vlr->v3->co, vlr->v4->co, nor);
					d2= nor[0]*vlr->v2->n[0] + nor[1]*vlr->v2->n[1] + nor[2]*vlr->v2->n[2];
					
					if( fabs(d1) < fabs(d2) ) vlr->flag |= R_DIVIDE_24;
					else vlr->flag &= ~R_DIVIDE_24;
					
					/* new vertex pointers */
					if (vlr->flag & R_DIVIDE_24) {
						vlr1->v1= vlr->v2;
						vlr1->v2= vlr->v3;
						vlr1->v3= vlr->v4;

						vlr->v3 = vlr->v4;
						
						vlr1->flag |= R_DIVIDE_24;
					}
					else {
						vlr1->v1= vlr->v1;
						vlr1->v2= vlr->v3;
						vlr1->v3= vlr->v4;
						
						vlr1->flag &= ~R_DIVIDE_24;
					}
					vlr->v4 = vlr1->v4 = NULL;
					
					/* new normals */
					CalcNormFloat(vlr->v3->co, vlr->v2->co, vlr->v1->co, vlr->n);
					CalcNormFloat(vlr1->v3->co, vlr1->v2->co, vlr1->v1->co, vlr1->n);
					
					/* so later UV can be pulled from original tface, look for R_DIVIDE_24 for direction */
					vlr1->tface=vlr->tface; 

				}
				/* clear the flag when not divided */
				else vlr->flag &= ~R_DIVIDE_24;
			}
		}
	}
}

static void set_material_lightgroups(void)
{
	GroupObject *go, *gol;
	Material *ma;
	
	/* it's a bit too many loops in loops... but will survive */
	for(ma= G.main->mat.first; ma; ma=ma->id.next) {
		if(ma->group) {
			for(go= ma->group->gobject.first; go; go= go->next) {
				for(gol= R.lights.first; gol; gol= gol->next) {
					if(gol->ob==go->ob) {
						go->lampren= gol->lampren;
						break;
					}
				}
			}
		}
	}
}

extern int slurph_opt;	/* key.c */
extern ListBase duplilist;
void RE_rotateBlenderScene(void)
{
	Base *base;
	Object *ob, *obd;
	Scene *sce, *setscene;
	unsigned int lay;
	float mat[4][4];

	if(G.scene->camera==NULL) return;

	R.memArena = BLI_memarena_new(BLI_MEMARENA_STD_BUFSIZE);
	R.totvlak=R.totvert=R.totlamp=R.tothalo= 0;
	R.lights.first= R.lights.last= NULL;
	
	slurph_opt= 0;

	/* in localview, lamps are using normal layers, objects only local bits */
	if(G.scene->lay & 0xFF000000) lay= G.scene->lay & 0xFF000000;
	else lay= G.scene->lay;
	
	/* applies changes fully */
	scene_update_for_newframe(G.scene, lay);

	MTC_Mat4CpyMat4(R.viewinv, G.scene->camera->obmat);
	MTC_Mat4Ortho(R.viewinv);
	MTC_Mat4Invert(R.viewmat, R.viewinv);

	RE_setwindowclip(1,-1); /*  no jit:(-1) */

	/* clear imat flags */
	ob= G.main->object.first;
	while(ob) {
		ob->flag &= ~OB_DO_IMAT;
		ob= ob->id.next;
	}

	init_render_world();	/* do first, because of ambient. also requires R.osa set correct */
	if( (R.wrld.mode & WO_AMB_OCC) && (R.r.mode & R_RAYTRACE) ) {
		R.wrld.aosphere= MEM_mallocN(2*3*R.wrld.aosamp*R.wrld.aosamp*sizeof(float), "AO sphere");
		/* we make twice the amount of samples, because only a hemisphere is used */
		init_ao_sphere(R.wrld.aosphere, 2*R.wrld.aosamp*R.wrld.aosamp, 16);
		
		/* bah... init_render_world writes this over, and that is called/needed in envmap. */
		G.scene->world->aosphere= R.wrld.aosphere;
	}
	init_render_textures();
	init_render_materials();
	set_node_shader_lamp_loop(shade_material_loop);

	/* imat objects, OB_DO_IMAT can be set in init_render_materials
	   has to be done here, since displace can have texture using Object map-input */
	ob= G.main->object.first;
	while(ob) {
		if(ob->flag & OB_DO_IMAT) {
			ob->flag &= ~OB_DO_IMAT;
			MTC_Mat4MulMat4(mat, ob->obmat, R.viewmat);
			MTC_Mat4Invert(ob->imat, mat);
		}
		ob= ob->id.next;
	}

	/* MAKE RENDER DATA */

	/* each object should only be rendered once */
	ob= G.main->object.first;
	while(ob) {
		ob->flag &= ~OB_DONE;
		ob= ob->id.next;
	}
	
	sce= G.scene;
	setscene= G.scene->set;
	base= G.scene->base.first;
	while(base) {

		ob= base->object;

		if(ob->flag & OB_DONE) {
			/* yafray: this object needs to be included in renderlist for duplivert instancing.
				 This only works for dupliverts, dupliframes handled below.
				 This is based on the assumption that OB_DONE is only set for duplivert objects,
				 before scene conversion, there are no other flags set to indicate it's use as far as I know...
				 NOT done for lamps, these are included as separate objects, see below.
				 correction: also ignore lattices, armatures and camera's (.....) */
			if ((ob->type!=OB_LATTICE) && (ob->type!=OB_ARMATURE) &&
					(ob->type!=OB_LAMP) && (ob->type!=OB_CAMERA) && (R.r.renderer==R_YAFRAY))
			{
				printf("Adding %s to renderlist\n", ob->id.name);
				ob->flag &= ~OB_DONE;
				init_render_object(ob);
				ob->flag |= OB_DONE;
			}
		}
		else {
			if( (base->lay & lay) || (ob->type==OB_LAMP && (base->lay & G.scene->lay)) ) {

				if(ob->transflag & OB_DUPLI) {
					
					/* exception: mballs! */
					/* yafray: Include at least one copy of a dupliframe object for yafray in the renderlist.
					   mballs comment above true as well for yafray, they are not included, only all other object types */
					if (R.r.renderer==R_YAFRAY) {
						if ((ob->type!=OB_MBALL) && ((ob->transflag & OB_DUPLIFRAMES)!=0)) {
							printf("Object %s has OB_DUPLIFRAMES set, adding to renderlist\n", ob->id.name);
							init_render_object(ob);
						}
					}
					/* before make duplis, update particle for current frame */
					if(ob->transflag & OB_DUPLIVERTS) {
						PartEff *paf= give_parteff(ob);
						if(paf) {
							if(paf->flag & PAF_ANIMATED) build_particle_system(ob);
						}
					}
					
					if(ob->type==OB_MBALL) {
						init_render_object(ob);
					}
					else {
						DupliObject *dob;
						ListBase *lb= object_duplilist(sce, ob);
						
						for(dob= lb->first; dob; dob= dob->next) {
							Object *obd= dob->ob;
							Mat4CpyMat4(obd->obmat, dob->mat);
							
							if(obd->type!=OB_MBALL) {
								/* yafray: special handling of duplivert objects for yafray:
								   only the matrix is stored, together with the source object name.
									 Since the original object is needed as well, it is included in the renderlist (see above)
									 NOT done for lamps, these need to be included as normal lamps separately
									 correction: also ignore lattices, armatures and cameras (....) */
								if ((obd->type!=OB_LATTICE) && (obd->type!=OB_ARMATURE) &&
										(obd->type!=OB_LAMP) && (obd->type!=OB_CAMERA) && (R.r.renderer==R_YAFRAY))
								{
									printf("Adding dupli matrix for object %s\n", obd->id.name);
									YAF_addDupliMtx(obd);
								}
								else init_render_object(obd);
							}
							Mat4CpyMat4(obd->obmat, dob->omat);
						}
						BLI_freelistN(lb);
					}
				}
				else {
					/* yafray: if there are linked data objects (except lamps, empties or armatures),
					   yafray only needs to know about one, the rest can be instanciated.
					   The dupliMtx list is used for this purpose.
					   Exception: objects which have object linked materials, these cannot be instanciated. */
					if ((R.r.renderer==R_YAFRAY) && (ob->colbits==0))
					{
						/* Special case, parent object dupli's: ignore if object itself is lamp or parent is lattice or empty */
						if (ob->parent) {
							if ((ob->type!=OB_LAMP) && (ob->parent->type!=OB_EMPTY) &&
									(ob->parent->type!=OB_LATTICE) && YAF_objectKnownData(ob))
								printf("From parent: Added dupli matrix for linked data object %s\n", ob->id.name);
							else
								init_render_object(ob);
						}
						else if ((ob->type!=OB_EMPTY) && (ob->type!=OB_LAMP) && (ob->type!=OB_ARMATURE) && YAF_objectKnownData(ob))
							printf("Added dupli matrix for linked data object %s\n", ob->id.name);
						else
							init_render_object(ob);
					}
					else init_render_object(ob);
				}

			}
			else {
				MTC_Mat4MulMat4(mat, ob->obmat, R.viewmat);
				MTC_Mat4Invert(ob->imat, mat);
			}

		}
		if(blender_test_break()) break;

		base= base->next;
		if(base==0 && setscene) {
			sce= setscene;
			base= setscene->base.first;
			setscene= setscene->set;
		}
		/*if(base->next==0 && G.scene->set && base==G.scene->base.last) {
			base= G.scene->set->base.first;
			sce= G.scene->set;
		}
		else base= base->next;*/

	}

	sort_halos();
	
	set_material_lightgroups();
	
	if(R.wrld.mode & WO_STARS) RE_make_stars(NULL, NULL, NULL);

	slurph_opt= 1;

	if(blender_test_break()) return;

	set_fullsample_flag();
	check_non_flat_quads();
	set_normalflags();
}

/* **************************************************************** */
/*                sticky texture coords                             */
/* **************************************************************** */

void RE_make_sticky(void)
{
	Object *ob;
	Base *base;
	MVert *mvert;
	Mesh *me;
	MSticky *ms;
	float ho[4], mat[4][4];
	int a;
	
	if(G.scene->camera==0) return;
	
	if(G.obedit) {
		error("Unable to make sticky in Edit Mode");
		return;
	}
	base= FIRSTBASE;
	while(base) {
		if TESTBASELIB(base) {
			if(base->object->type==OB_MESH) {
				ob= base->object;
				
				me= ob->data;
				mvert= me->mvert;
				if(me->msticky) MEM_freeN(me->msticky);
				me->msticky= MEM_mallocN(me->totvert*sizeof(MSticky), "sticky");
				
				/* like convert to render data */		
				R.r= G.scene->r;
				R.r.xsch= (R.r.size*R.r.xsch)/100;
				R.r.ysch= (R.r.size*R.r.ysch)/100;
				
				R.afmx= R.r.xsch/2;
				R.afmy= R.r.ysch/2;
				
				R.ycor= ( (float)R.r.yasp)/( (float)R.r.xasp);
				
				R.rectx= R.r.xsch; 
				R.recty= R.r.ysch;
				R.xstart= -R.afmx; 
				R.ystart= -R.afmy;
				R.xend= R.xstart+R.rectx-1;
				R.yend= R.ystart+R.recty-1;
				
				where_is_object(G.scene->camera);
				Mat4CpyMat4(R.viewinv, G.scene->camera->obmat);
				Mat4Ortho(R.viewinv);
				Mat4Invert(R.viewmat, R.viewinv);
				
				RE_setwindowclip(1, -1);
				
				where_is_object(ob);
				Mat4MulMat4(mat, ob->obmat, R.viewmat);
				
				ms= me->msticky;
				for(a=0; a<me->totvert; a++, ms++, mvert++) {
					VECCOPY(ho, mvert->co);
					Mat4MulVecfl(mat, ho);
					RE_projectverto(ho, ho);
					ms->co[0]= ho[0]/ho[3];
					ms->co[1]= ho[1]/ho[3];
				}
			}
		}
		base= base->next;
	}
}


/* **************************************************************** */
/*                Displacement mapping                              */
/* **************************************************************** */
static short test_for_displace(Object *ob)
{
	/* return 1 when this object uses displacement textures. */
	Material *ma;
	int i;
	
	for (i=1; i<=ob->totcol; i++) {
		ma=give_render_material(ob, i);
		/* ma->mapto is ORed total of all mapto channels */
		if(ma && (ma->mapto & MAP_DISPLACE)) return 1;
	}
	return 0;
}

static void displace_render_vert(ShadeInput *shi, VertRen *vr, float *scale)
{
	short texco= shi->mat->texco;
	float sample=0;
	/* shi->co is current render coord, just make sure at least some vector is here */
	VECCOPY(shi->co, vr->co);
	/* vertex normal is used for textures type 'col' and 'var' */
	VECCOPY(shi->vn, vr->n);
	
	/* set all rendercoords, 'texco' is an ORed value for all textures needed */
	if ((texco & TEXCO_ORCO) && (vr->orco)) {
		VECCOPY(shi->lo, vr->orco);
	}
	if (texco & TEXCO_STICKY) {
		float *sticky= RE_vertren_get_sticky(vr, 0);
		if(sticky) {
			shi->sticky[0]= sticky[0];
			shi->sticky[1]= sticky[1];
			shi->sticky[2]= 0.0f;
		}
	}
	if (texco & TEXCO_GLOB) {
		VECCOPY(shi->gl, shi->co);
		MTC_Mat4MulVecfl(R.viewinv, shi->gl);
	}
	if (texco & TEXCO_NORM) {
		VECCOPY(shi->orn, shi->vn);
	}
	if(texco & TEXCO_REFL) {
		/* not (yet?) */
	}
	
	shi->displace[0]= shi->displace[1]= shi->displace[2]= 0.0;
	
	do_material_tex(shi);
	
	//printf("no=%f, %f, %f\nbefore co=%f, %f, %f\n", vr->n[0], vr->n[1], vr->n[2], 
	//vr->co[0], vr->co[1], vr->co[2]);
	
	/* 0.5 could become button once?  */
	vr->co[0] +=  shi->displace[0] * scale[0] ; 
	vr->co[1] +=  shi->displace[1] * scale[1] ; 
	vr->co[2] +=  shi->displace[2] * scale[2] ; 
	
	//printf("after co=%f, %f, %f\n", vr->co[0], vr->co[1], vr->co[2]); 
	
	/* we just don't do this vertex again, bad luck for other face using same vertex with
		different material... */
	vr->flag |= 1;
	
	/* Pass sample back so displace_face can decide which way to split the quad */
	sample  = shi->displace[0]*shi->displace[0];
	sample += shi->displace[1]*shi->displace[1];
	sample += shi->displace[2]*shi->displace[2];
	
	vr->accum=sample; 
	/* Should be sqrt(sample), but I'm only looking for "bigger".  Save the cycles. */
	return;
}

static void displace_render_face(VlakRen *vlr, float *scale)
{
	ShadeInput shi;
	//	VertRen vr;
	//	float samp1,samp2, samp3, samp4, xn;
	short hasuv=0;
	/* set up shadeinput struct for multitex() */
	
	shi.osatex= 0;		/* signal not to use dx[] and dy[] texture AA vectors */
	shi.vlr= vlr;		/* current render face */
	shi.mat= vlr->mat;		/* current input material */
	
	
	/* UV coords must come from face */
	hasuv = vlr->tface && (shi.mat->texco & TEXCO_UV);
	if (hasuv) shi.uv[2]=0.0f; 
	/* I don't think this is used, but seting it just in case */
	
	/* Displace the verts, flag is set when done */
	if (! (vlr->v1->flag)){		
		if (hasuv)	{
			shi.uv[0] = 2*vlr->tface->uv[0][0]-1.0f; /* shi.uv and tface->uv are */
			shi.uv[1]=  2*vlr->tface->uv[0][1]-1.0f; /* scalled differently 	 */
		}
		displace_render_vert(&shi, vlr->v1, scale);
	}
	
	if (! (vlr->v2->flag)) {
		if (hasuv)	{
			shi.uv[0] = 2*vlr->tface->uv[1][0]-1.0f; 
			shi.uv[1]=  2*vlr->tface->uv[1][1]-1.0f;
		}
		displace_render_vert(&shi, vlr->v2, scale);
	}
	
 	if (! (vlr->v3->flag)) {
		if (hasuv)	{
			shi.uv[0] = 2*vlr->tface->uv[2][0]-1.0f; 
			shi.uv[1]=  2*vlr->tface->uv[2][1]-1.0f;
		}	
		displace_render_vert(&shi, vlr->v3, scale);
	}
	
	if (vlr->v4) {
		if (! (vlr->v4->flag)) {
		 	if (hasuv)	{
				shi.uv[0] = 2*vlr->tface->uv[3][0]-1.0f; 
				shi.uv[1]=  2*vlr->tface->uv[3][1]-1.0f;
			}	
			displace_render_vert(&shi, vlr->v4, scale);
		}
		/* We want to split the quad along the opposite verts that are */
		/*	closest in displace value.  This will help smooth edges.   */ 
		if ( fabs(vlr->v1->accum - vlr->v3->accum) > fabs(vlr->v2->accum - vlr->v4->accum)) 
			vlr->flag |= R_DIVIDE_24;
		else vlr->flag &= ~R_DIVIDE_24;	// E: typo?, was missing '='
	}
	
	/* Recalculate the face normal  - if flipped before, flip now */
	if(vlr->v4) {
		CalcNormFloat4(vlr->v4->co, vlr->v3->co, vlr->v2->co, vlr->v1->co, vlr->n);
	}	
	else {
		CalcNormFloat(vlr->v3->co, vlr->v2->co, vlr->v1->co, vlr->n);
	}
	
}


static void do_displacement(Object *ob, int startface, int numface, int startvert, int numvert )
{
	VertRen *vr;
	VlakRen *vlr;
//	float min[3]={1e30, 1e30, 1e30}, max[3]={-1e30, -1e30, -1e30};
	float scale[3]={1.0f, 1.0f, 1.0f}, temp[3];//, xn
	int i; //, texflag=0;
	Object *obt;
		
	/* Object Size with parenting */
	obt=ob;
	while(obt){
		VecAddf(temp, obt->size, obt->dsize);
		scale[0]*=temp[0]; scale[1]*=temp[1]; scale[2]*=temp[2];
		obt=obt->parent;
	}
	
	/* Clear all flags */
	for(i=startvert; i<startvert+numvert; i++){ 
		vr= RE_findOrAddVert(i);
		vr->flag= 0;
	}
	
	for(i=startface; i<startface+numface; i++){
		vlr=RE_findOrAddVlak(i);
		displace_render_face(vlr, scale);
	}
	
	/* Recalc vertex normals */
	calc_vertexnormals(startvert, startface, 0);
}

