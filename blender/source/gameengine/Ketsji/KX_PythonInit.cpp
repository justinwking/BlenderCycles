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
 * Initialize Python thingies.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef WIN32
#pragma warning (disable : 4786)
#endif //WIN32

#include "KX_PythonInit.h"

#include "SCA_IInputDevice.h"
#include "SCA_PropertySensor.h"
#include "SCA_RandomActuator.h"
#include "KX_ConstraintActuator.h"
#include "KX_IpoActuator.h"
#include "RAS_IRasterizer.h"
#include "RAS_ICanvas.h"
#include "MT_Vector3.h"
#include "MT_Point3.h"
#include "ListValue.h"
#include "KX_Scene.h"
#include "SND_DeviceManager.h"

#include "KX_PyMath.h"

// FIXME: Enable for access to blender python modules.  This is disabled because
// python has dependencies on a lot of other modules and is a pain to link.
//#define USE_BLENDER_PYTHON
#ifdef USE_BLENDER_PYTHON
//#include "BPY_extern.h"
#endif 

static void setSandbox(TPythonSecurityLevel level);


// 'local' copy of canvas ptr, for window height/width python scripts
static RAS_ICanvas* gp_Canvas = NULL;
static KX_Scene*	gp_KetsjiScene = NULL;
static RAS_IRasterizer* gp_Rasterizer = NULL;

/* Macro for building the keyboard translation */
//#define KX_MACRO_addToDict(dict, name) PyDict_SetItemString(dict, #name, PyInt_FromLong(SCA_IInputDevice::KX_##name))
#define KX_MACRO_addToDict(dict, name) PyDict_SetItemString(dict, #name, PyInt_FromLong(name))
/* For the defines for types from logic bricks, we do stuff explicitly... */
#define KX_MACRO_addTypesToDict(dict, name, name2) PyDict_SetItemString(dict, #name, PyInt_FromLong(name2))


// temporarily python stuff, will be put in another place later !
#include "KX_Python.h"
#include "SCA_PythonController.h"
// List of methods defined in the module

static PyObject* ErrorObject;
STR_String gPyGetRandomFloat_doc="getRandomFloat returns a random floating point value in the range [0..1)";

static PyObject* gPyGetRandomFloat(PyObject* self,
										 PyObject* args, 
										 PyObject* kwds)
{
	return PyFloat_FromDouble(MT_random());
}





void GlobalConvertPythonVectorArg(PyObject* args, MT_Vector3 &pos)
{
	PyObject* pylist;
	PyArg_ParseTuple(args,"O",&pylist);

	pos = MT_Vector3FromPyList(pylist);
}

void GlobalConvertPythonVectorArg(PyObject* args, MT_Vector4 &vec)
{
	PyObject* pylist;
	PyArg_ParseTuple(args,"O",&pylist);

	vec = MT_Vector4FromPyList(pylist);
}


static PyObject* gPySetGravity(PyObject* self,
										 PyObject* args, 
										 PyObject* kwds)
{
	MT_Vector3 vec = MT_Vector3(0., 0., 0.);
	GlobalConvertPythonVectorArg(args, vec);

	if (gp_KetsjiScene)
		gp_KetsjiScene->SetGravity(vec);
	
	Py_Return;
}


static bool usedsp = false;

// this gets a pointer to an array filled with floats
static PyObject* gPyGetSpectrum(PyObject* self,
								PyObject* args, 
								PyObject* kwds)
{
	SND_IAudioDevice* audiodevice = SND_DeviceManager::Instance();

	PyObject* resultlist = PyList_New(512);

	if (audiodevice)
	{
		if (!usedsp)
		{
			audiodevice->StartUsingDSP();
			usedsp = true;
		}
			
		float* spectrum = audiodevice->GetSpectrum();

		for (int index = 0; index < 512; index++)
		{
			PyList_SetItem(resultlist, index, PyFloat_FromDouble(spectrum[index]));
		}
	}

	return resultlist;
}



static PyObject* gPyStartDSP(PyObject* self,
						PyObject* args, 
						PyObject* kwds)
{
	SND_IAudioDevice* audiodevice = SND_DeviceManager::Instance();

	if (audiodevice)
	{
		if (!usedsp)
		{
			audiodevice->StartUsingDSP();
			usedsp = true;
			Py_Return;
		}
	}
	return NULL;
}



static PyObject* gPyStopDSP(PyObject* self,
					   PyObject* args, 
					   PyObject* kwds)
{
	SND_IAudioDevice* audiodevice = SND_DeviceManager::Instance();

	if (audiodevice)
	{
		if (usedsp)
		{
			audiodevice->StopUsingDSP();
			usedsp = false;
			Py_Return;
		}
	}
	return NULL;
}


static STR_String gPyGetCurrentScene_doc =  
"getCurrentScene()\n"
"Gets a reference to the current scene.\n";
static PyObject* gPyGetCurrentScene(PyObject* self,
					   PyObject* args, 
					   PyObject* kwds)
{
	Py_INCREF(gp_KetsjiScene);
	return (PyObject*) gp_KetsjiScene;
}
					   


static struct PyMethodDef game_methods[] = {
	{"getCurrentController",
	(PyCFunction) SCA_PythonController::sPyGetCurrentController,
	METH_VARARGS, SCA_PythonController::sPyGetCurrentController__doc__},
	{"getCurrentScene", (PyCFunction) gPyGetCurrentScene,
	METH_VARARGS, gPyGetCurrentScene_doc.Ptr()},
	{"addActiveActuator",(PyCFunction) SCA_PythonController::sPyAddActiveActuator,
	METH_VARARGS, SCA_PythonController::sPyAddActiveActuator__doc__},
	{"getRandomFloat",(PyCFunction) gPyGetRandomFloat,
	METH_VARARGS,gPyGetRandomFloat_doc.Ptr()},
	{"setGravity",(PyCFunction) gPySetGravity, METH_VARARGS,"set Gravitation"},
	{"getSpectrum",(PyCFunction) gPyGetSpectrum, METH_VARARGS,"get audio spectrum"},
	{"stopDSP",(PyCFunction) gPyStopDSP, METH_VARARGS,"stop using the audio dsp (for performance reasons)"},
	{NULL, (PyCFunction) NULL, 0, NULL }
};


static PyObject* gPyGetWindowHeight(PyObject* self, 
										 PyObject* args, 
										 PyObject* kwds)
{
	int height = (gp_Canvas ? gp_Canvas->GetHeight() : 0);

		PyObject* heightval = PyInt_FromLong(height);
		return heightval;
}



static PyObject* gPyGetWindowWidth(PyObject* self, 
										 PyObject* args, 
										 PyObject* kwds)
{
		

	int width = (gp_Canvas ? gp_Canvas->GetWidth() : 0);
	
		PyObject* widthval = PyInt_FromLong(width);
		return widthval;
}



// temporarility visibility thing, will be moved to rasterizer/renderer later
bool gUseVisibilityTemp = false;

static PyObject* gPyEnableVisibility(PyObject* self, 
										 PyObject* args, 
										 PyObject* kwds)
{
	int visible;
	if (PyArg_ParseTuple(args,"i",&visible))
	{
	    gUseVisibilityTemp = (visible != 0);
	}
	else
	{
	  Py_Return;	     
	}
   Py_Return;
}



static PyObject* gPyShowMouse(PyObject* self, 
										 PyObject* args, 
										 PyObject* kwds)
{
	int visible;
	if (PyArg_ParseTuple(args,"i",&visible))
	{
	    if (visible)
		{
			if (gp_Canvas)
				gp_Canvas->SetMouseState(RAS_ICanvas::MOUSE_NORMAL);
		} else
		{
			if (gp_Canvas)
				gp_Canvas->SetMouseState(RAS_ICanvas::MOUSE_INVISIBLE);
		}
	}
	
   Py_Return;
}



static PyObject* gPySetMousePosition(PyObject* self, 
										 PyObject* args, 
										 PyObject* kwds)
{
	int x,y;
	if (PyArg_ParseTuple(args,"ii",&x,&y))
	{
	    if (gp_Canvas)
			gp_Canvas->SetMousePosition(x,y);
	}
	
   Py_Return;
}



static PyObject* gPySetBackgroundColor(PyObject* self, 
										 PyObject* args, 
										 PyObject* kwds)
{
	
	MT_Vector4 vec = MT_Vector4(0., 0., 0.3, 0.);
	GlobalConvertPythonVectorArg(args, vec);

	if (gp_Canvas)
	{
		gp_Rasterizer->SetBackColor(vec[0], vec[1], vec[2], vec[3]);
	}
   Py_Return;
}



static PyObject* gPySetMistColor(PyObject* self, 
										 PyObject* args, 
										 PyObject* kwds)
{
	
	MT_Vector3 vec = MT_Vector3(0., 0., 0.);
	GlobalConvertPythonVectorArg(args, vec);

	if (gp_Rasterizer)
	{
		gp_Rasterizer->SetFogColor(vec[0], vec[1], vec[2]);
	}
   Py_Return;
}



static PyObject* gPySetMistStart(PyObject* self, 
										 PyObject* args, 
										 PyObject* kwds)
{

	float miststart;
	if (PyArg_ParseTuple(args,"f",&miststart))
	{
		if (gp_Rasterizer)
		{
			gp_Rasterizer->SetFogStart(miststart);
		}
	}
   Py_Return;
}



static PyObject* gPySetMistEnd(PyObject* self, 
										 PyObject* args, 
										 PyObject* kwds)
{

	float mistend;
	if (PyArg_ParseTuple(args,"f",&mistend))
	{
		if (gp_Rasterizer)
		{
			gp_Rasterizer->SetFogEnd(mistend);
		}
	}
   Py_Return;
}



static PyObject* gPyMakeScreenshot(PyObject* self,
									PyObject* args,
									PyObject* kwds)
{
	char* filename;
	if (PyArg_ParseTuple(args,"s",&filename))
	{
		if (gp_Canvas)
		{
			gp_Canvas->MakeScreenShot(filename);
		}
	}
	Py_Return;
}



STR_String	gPyGetWindowHeight__doc__="getWindowHeight doc";
STR_String	gPyGetWindowWidth__doc__="getWindowWidth doc";
STR_String	gPyEnableVisibility__doc__="enableVisibility doc";
STR_String	gPyMakeScreenshot__doc__="make Screenshot doc";
STR_String	gPyShowMouse__doc__="showMouse(bool visible)";
STR_String	gPySetMousePosition__doc__="setMousePosition(int x,int y)";

static struct PyMethodDef rasterizer_methods[] = {
  {"getWindowWidth",(PyCFunction) gPyGetWindowWidth,
   METH_VARARGS, gPyGetWindowWidth__doc__.Ptr()},
   {"getWindowHeight",(PyCFunction) gPyGetWindowHeight,
   METH_VARARGS, gPyGetWindowHeight__doc__.Ptr()},
  {"makeScreenshot",(PyCFunction)gPyMakeScreenshot,
	METH_VARARGS, gPyMakeScreenshot__doc__.Ptr()},
   {"enableVisibility",(PyCFunction) gPyEnableVisibility,
   METH_VARARGS, gPyEnableVisibility__doc__.Ptr()},
	{"showMouse",(PyCFunction) gPyShowMouse,
   METH_VARARGS, gPyShowMouse__doc__.Ptr()},
   {"setMousePosition",(PyCFunction) gPySetMousePosition,
   METH_VARARGS, gPySetMousePosition__doc__.Ptr()},
  {"setBackgroundColor",(PyCFunction)gPySetBackgroundColor,METH_VARARGS,"set Background Color (rgb)"},
  {"setMistColor",(PyCFunction)gPySetMistColor,METH_VARARGS,"set Mist Color (rgb)"},
  {"setMistStart",(PyCFunction)gPySetMistStart,METH_VARARGS,"set Mist Start(rgb)"},
  {"setMistEnd",(PyCFunction)gPySetMistEnd,METH_VARARGS,"set Mist End(rgb)"},
  
  { NULL, (PyCFunction) NULL, 0, NULL }
};



// Initialization function for the module (*must* be called initGameLogic)

static char GameLogic_module_documentation[] =
"This is the Python API for the game engine of GameLogic"
;

static char Rasterizer_module_documentation[] =
"This is the Python API for the game engine of Rasterizer"
;



PyObject* initGameLogic(KX_Scene* scene) // quick hack to get gravity hook
{
	PyObject* m;
	PyObject* d;

	gp_KetsjiScene = scene;

	gUseVisibilityTemp=false;

	// Create the module and add the functions
	m = Py_InitModule4("GameLogic", game_methods,
					   GameLogic_module_documentation,
					   (PyObject*)NULL,PYTHON_API_VERSION);

	// Add some symbolic constants to the module
	d = PyModule_GetDict(m);

	ErrorObject = PyString_FromString("GameLogic.error");
	PyDict_SetItemString(d, "error", ErrorObject);

	// XXXX Add constants here
	/* To use logic bricks, we need some sort of constants. Here, we associate */
	/* constants and sumbolic names. Add them to dictionary d.                 */

	/* 1. true and false: needed for everyone                                  */
	KX_MACRO_addTypesToDict(d, KX_TRUE,  SCA_ILogicBrick::KX_TRUE);
	KX_MACRO_addTypesToDict(d, KX_FALSE, SCA_ILogicBrick::KX_FALSE);

	/* 2. Property sensor                                                      */
	KX_MACRO_addTypesToDict(d, KX_PROPSENSOR_EQUAL,      SCA_PropertySensor::KX_PROPSENSOR_EQUAL);
	KX_MACRO_addTypesToDict(d, KX_PROPSENSOR_NOTEQUAL,   SCA_PropertySensor::KX_PROPSENSOR_NOTEQUAL);
	KX_MACRO_addTypesToDict(d, KX_PROPSENSOR_INTERVAL,   SCA_PropertySensor::KX_PROPSENSOR_INTERVAL);
	KX_MACRO_addTypesToDict(d, KX_PROPSENSOR_CHANGED,    SCA_PropertySensor::KX_PROPSENSOR_CHANGED);
	KX_MACRO_addTypesToDict(d, KX_PROPSENSOR_EXPRESSION, SCA_PropertySensor::KX_PROPSENSOR_EXPRESSION);

	/* 3. Constraint actuator                                                  */
	KX_MACRO_addTypesToDict(d, KX_CONSTRAINTACT_LOCX, KX_ConstraintActuator::KX_ACT_CONSTRAINT_LOCX);
	KX_MACRO_addTypesToDict(d, KX_CONSTRAINTACT_LOCY, KX_ConstraintActuator::KX_ACT_CONSTRAINT_LOCY);
	KX_MACRO_addTypesToDict(d, KX_CONSTRAINTACT_LOCZ, KX_ConstraintActuator::KX_ACT_CONSTRAINT_LOCZ);
	KX_MACRO_addTypesToDict(d, KX_CONSTRAINTACT_ROTX, KX_ConstraintActuator::KX_ACT_CONSTRAINT_ROTX);
	KX_MACRO_addTypesToDict(d, KX_CONSTRAINTACT_ROTY, KX_ConstraintActuator::KX_ACT_CONSTRAINT_ROTY);
	KX_MACRO_addTypesToDict(d, KX_CONSTRAINTACT_ROTZ, KX_ConstraintActuator::KX_ACT_CONSTRAINT_ROTZ);

	/* 4. Ipo actuator, simple part                                            */
	KX_MACRO_addTypesToDict(d, KX_IPOACT_PLAY,     KX_IpoActuator::KX_ACT_IPO_PLAY);
	KX_MACRO_addTypesToDict(d, KX_IPOACT_PINGPONG, KX_IpoActuator::KX_ACT_IPO_PINGPONG);
	KX_MACRO_addTypesToDict(d, KX_IPOACT_FLIPPER,  KX_IpoActuator::KX_ACT_IPO_FLIPPER);
	KX_MACRO_addTypesToDict(d, KX_IPOACT_LOOPSTOP, KX_IpoActuator::KX_ACT_IPO_LOOPSTOP);
	KX_MACRO_addTypesToDict(d, KX_IPOACT_LOOPEND,  KX_IpoActuator::KX_ACT_IPO_LOOPEND);

	/* 5. Random distribution types                                            */
	KX_MACRO_addTypesToDict(d, KX_RANDOMACT_BOOL_CONST,      SCA_RandomActuator::KX_RANDOMACT_BOOL_CONST);
	KX_MACRO_addTypesToDict(d, KX_RANDOMACT_BOOL_UNIFORM,    SCA_RandomActuator::KX_RANDOMACT_BOOL_UNIFORM);
	KX_MACRO_addTypesToDict(d, KX_RANDOMACT_BOOL_BERNOUILLI, SCA_RandomActuator::KX_RANDOMACT_BOOL_BERNOUILLI);
	KX_MACRO_addTypesToDict(d, KX_RANDOMACT_INT_CONST,       SCA_RandomActuator::KX_RANDOMACT_INT_CONST);
	KX_MACRO_addTypesToDict(d, KX_RANDOMACT_INT_UNIFORM,     SCA_RandomActuator::KX_RANDOMACT_INT_UNIFORM);
	KX_MACRO_addTypesToDict(d, KX_RANDOMACT_INT_POISSON,     SCA_RandomActuator::KX_RANDOMACT_INT_POISSON);
	KX_MACRO_addTypesToDict(d, KX_RANDOMACT_FLOAT_CONST,     SCA_RandomActuator::KX_RANDOMACT_FLOAT_CONST);
	KX_MACRO_addTypesToDict(d, KX_RANDOMACT_FLOAT_UNIFORM,   SCA_RandomActuator::KX_RANDOMACT_FLOAT_UNIFORM);
	KX_MACRO_addTypesToDict(d, KX_RANDOMACT_FLOAT_NORMAL,    SCA_RandomActuator::KX_RANDOMACT_FLOAT_NORMAL);
	KX_MACRO_addTypesToDict(d, KX_RANDOMACT_FLOAT_NEGATIVE_EXPONENTIAL, SCA_RandomActuator::KX_RANDOMACT_FLOAT_NEGATIVE_EXPONENTIAL);

	// Check for errors
	if (PyErr_Occurred())
    {
		Py_FatalError("can't initialize module GameLogic");
    }

	return d;
}



// Python Sandbox code
// override builtin functions import() and open()


PyObject *KXpy_open(PyObject *self, PyObject *args)
{
	PyErr_SetString(PyExc_RuntimeError, "Sandbox: open() function disabled!\nGame Scripts should not use this function.");
	return NULL;
}



PyObject *KXpy_import(PyObject *self, PyObject *args)
{
	char *name;
	PyObject *globals = NULL;
	PyObject *locals = NULL;
	PyObject *fromlist = NULL;
	PyObject *l, *m, *n;

	if (!PyArg_ParseTuple(args, "s|OOO:m_import",
	        &name, &globals, &locals, &fromlist))
	    return NULL;

	/* check for builtin modules */
	m = PyImport_AddModule("sys");
	l = PyObject_GetAttrString(m, "builtin_module_names");
	n = PyString_FromString(name);
	
	if (PySequence_Contains(l, n)) {
		return PyImport_ImportModuleEx(name, globals, locals, fromlist);
	}

	/* quick hack for GamePython modules 
		TODO: register builtin modules properly by ExtendInittab */
	if (!strcmp(name, "GameLogic") || !strcmp(name, "GameKeys") ||
		!strcmp(name, "Rasterizer")) {
		return PyImport_ImportModuleEx(name, globals, locals, fromlist);
	}
		
	PyErr_Format(PyExc_ImportError,
		 "Import of external Module %.20s not allowed.", name);
	return NULL;

}



static PyMethodDef meth_open[] = {
	{ "open", KXpy_open, METH_VARARGS,
		"(disabled)"}
};


static PyMethodDef meth_import[] = {
	{ "import", KXpy_import, METH_VARARGS,
		"our own import"}
};



//static PyObject *g_oldopen = 0;
//static PyObject *g_oldimport = 0;
//static int g_security = 0;


void setSandbox(TPythonSecurityLevel level)
{
    PyObject *m = PyImport_AddModule("__builtin__");
    PyObject *d = PyModule_GetDict(m);
	PyObject *meth = PyCFunction_New(meth_open, NULL);

	switch (level) {
	case psl_Highest:
		//if (!g_security) {
			//g_oldopen = PyDict_GetItemString(d, "open");
			PyDict_SetItemString(d, "open", meth);
			meth = PyCFunction_New(meth_import, NULL);
			PyDict_SetItemString(d, "__import__", meth);
			//g_security = level;
		//}
		break;
	/*
	case psl_Lowest:
		if (g_security) {
			PyDict_SetItemString(d, "open", g_oldopen);
			PyDict_SetItemString(d, "__import__", g_oldimport);
			g_security = level;
		}
	*/
	default:
		break;
	}
}

/**
 * Python is not initialised.
 */
PyObject* initGamePlayerPythonScripting(const STR_String& progname, TPythonSecurityLevel level)
{
	STR_String pname = progname;
	Py_SetProgramName(pname.Ptr());
	Py_NoSiteFlag=1;
	Py_FrozenFlag=1;
	Py_Initialize();

	//importBlenderModules()
	
	setSandbox(level);

	PyObject* moduleobj = PyImport_AddModule("__main__");
	return PyModule_GetDict(moduleobj);
}

void exitGamePlayerPythonScripting()
{
	Py_Finalize();
}

/**
 * Python is already initialized.
 */
PyObject* initGamePythonScripting(const STR_String& progname, TPythonSecurityLevel level)
{
	STR_String pname = progname;
	Py_SetProgramName(pname.Ptr());
	Py_NoSiteFlag=1;
	Py_FrozenFlag=1;

	setSandbox(level);

	PyObject* moduleobj = PyImport_AddModule("__main__");
	return PyModule_GetDict(moduleobj);
}



void exitGamePythonScripting()
{
}



PyObject* initRasterizer(RAS_IRasterizer* rasty,RAS_ICanvas* canvas)
{
	gp_Canvas = canvas;
	gp_Rasterizer = rasty;


  PyObject* m;
  PyObject* d;

  // Create the module and add the functions
  m = Py_InitModule4("Rasterizer", rasterizer_methods,
		     Rasterizer_module_documentation,
		     (PyObject*)NULL,PYTHON_API_VERSION);

  // Add some symbolic constants to the module
  d = PyModule_GetDict(m);
  ErrorObject = PyString_FromString("Rasterizer.error");
  PyDict_SetItemString(d, "error", ErrorObject);

  // XXXX Add constants here

  // Check for errors
  if (PyErr_Occurred())
    {
      Py_FatalError("can't initialize module Rasterizer");
    }

  return d;
}



/* ------------------------------------------------------------------------- */
/* GameKeys: symbolic constants for key mapping                              */
/* ------------------------------------------------------------------------- */

static char GameKeys_module_documentation[] =
"This modules provides defines for key-codes"
;



static struct PyMethodDef gamekeys_methods[] = {
	{ NULL, (PyCFunction) NULL, 0, NULL }
};



PyObject* initGameKeys()
{
	PyObject* m;
	PyObject* d;

	// Create the module and add the functions
	m = Py_InitModule4("GameKeys", gamekeys_methods,
					   GameKeys_module_documentation,
					   (PyObject*)NULL,PYTHON_API_VERSION);

	// Add some symbolic constants to the module
	d = PyModule_GetDict(m);

	// XXXX Add constants here

	KX_MACRO_addTypesToDict(d, AKEY, SCA_IInputDevice::KX_AKEY);
	KX_MACRO_addTypesToDict(d, BKEY, SCA_IInputDevice::KX_BKEY);
	KX_MACRO_addTypesToDict(d, CKEY, SCA_IInputDevice::KX_CKEY);
	KX_MACRO_addTypesToDict(d, DKEY, SCA_IInputDevice::KX_DKEY);
	KX_MACRO_addTypesToDict(d, EKEY, SCA_IInputDevice::KX_EKEY);
	KX_MACRO_addTypesToDict(d, FKEY, SCA_IInputDevice::KX_FKEY);
	KX_MACRO_addTypesToDict(d, GKEY, SCA_IInputDevice::KX_GKEY);
	KX_MACRO_addTypesToDict(d, HKEY, SCA_IInputDevice::KX_HKEY);
	KX_MACRO_addTypesToDict(d, IKEY, SCA_IInputDevice::KX_IKEY);
	KX_MACRO_addTypesToDict(d, JKEY, SCA_IInputDevice::KX_JKEY);
	KX_MACRO_addTypesToDict(d, KKEY, SCA_IInputDevice::KX_KKEY);
	KX_MACRO_addTypesToDict(d, LKEY, SCA_IInputDevice::KX_LKEY);
	KX_MACRO_addTypesToDict(d, MKEY, SCA_IInputDevice::KX_MKEY);
	KX_MACRO_addTypesToDict(d, NKEY, SCA_IInputDevice::KX_NKEY);
	KX_MACRO_addTypesToDict(d, OKEY, SCA_IInputDevice::KX_OKEY);
	KX_MACRO_addTypesToDict(d, PKEY, SCA_IInputDevice::KX_PKEY);
	KX_MACRO_addTypesToDict(d, QKEY, SCA_IInputDevice::KX_QKEY);
	KX_MACRO_addTypesToDict(d, RKEY, SCA_IInputDevice::KX_RKEY);
	KX_MACRO_addTypesToDict(d, SKEY, SCA_IInputDevice::KX_SKEY);
	KX_MACRO_addTypesToDict(d, TKEY, SCA_IInputDevice::KX_TKEY);
	KX_MACRO_addTypesToDict(d, UKEY, SCA_IInputDevice::KX_UKEY);
	KX_MACRO_addTypesToDict(d, VKEY, SCA_IInputDevice::KX_VKEY);
	KX_MACRO_addTypesToDict(d, WKEY, SCA_IInputDevice::KX_WKEY);
	KX_MACRO_addTypesToDict(d, XKEY, SCA_IInputDevice::KX_XKEY);
	KX_MACRO_addTypesToDict(d, YKEY, SCA_IInputDevice::KX_YKEY);
	KX_MACRO_addTypesToDict(d, ZKEY, SCA_IInputDevice::KX_ZKEY);
	
	KX_MACRO_addTypesToDict(d, ZEROKEY, SCA_IInputDevice::KX_ZEROKEY);		
	KX_MACRO_addTypesToDict(d, ONEKEY, SCA_IInputDevice::KX_ONEKEY);		
	KX_MACRO_addTypesToDict(d, TWOKEY, SCA_IInputDevice::KX_TWOKEY);		
	KX_MACRO_addTypesToDict(d, THREEKEY, SCA_IInputDevice::KX_THREEKEY);
	KX_MACRO_addTypesToDict(d, FOURKEY, SCA_IInputDevice::KX_FOURKEY);		
	KX_MACRO_addTypesToDict(d, FIVEKEY, SCA_IInputDevice::KX_FIVEKEY);		
	KX_MACRO_addTypesToDict(d, SIXKEY, SCA_IInputDevice::KX_SIXKEY);		
	KX_MACRO_addTypesToDict(d, SEVENKEY, SCA_IInputDevice::KX_SEVENKEY);
	KX_MACRO_addTypesToDict(d, EIGHTKEY, SCA_IInputDevice::KX_EIGHTKEY);
	KX_MACRO_addTypesToDict(d, NINEKEY, SCA_IInputDevice::KX_NINEKEY);		
		
	KX_MACRO_addTypesToDict(d, CAPSLOCKKEY, SCA_IInputDevice::KX_CAPSLOCKKEY);
		
	KX_MACRO_addTypesToDict(d, LEFTCTRLKEY, SCA_IInputDevice::KX_LEFTCTRLKEY);	
	KX_MACRO_addTypesToDict(d, LEFTALTKEY, SCA_IInputDevice::KX_LEFTALTKEY); 		
	KX_MACRO_addTypesToDict(d, RIGHTALTKEY, SCA_IInputDevice::KX_RIGHTALTKEY); 	
	KX_MACRO_addTypesToDict(d, RIGHTCTRLKEY, SCA_IInputDevice::KX_RIGHTCTRLKEY); 	
	KX_MACRO_addTypesToDict(d, RIGHTSHIFTKEY, SCA_IInputDevice::KX_RIGHTSHIFTKEY);	
	KX_MACRO_addTypesToDict(d, LEFTSHIFTKEY, SCA_IInputDevice::KX_LEFTSHIFTKEY);
		
	KX_MACRO_addTypesToDict(d, ESCKEY, SCA_IInputDevice::KX_ESCKEY);
	KX_MACRO_addTypesToDict(d, TABKEY, SCA_IInputDevice::KX_TABKEY);
	KX_MACRO_addTypesToDict(d, RETKEY, SCA_IInputDevice::KX_RETKEY);
	KX_MACRO_addTypesToDict(d, SPACEKEY, SCA_IInputDevice::KX_SPACEKEY);
	KX_MACRO_addTypesToDict(d, LINEFEEDKEY, SCA_IInputDevice::KX_LINEFEEDKEY);		
	KX_MACRO_addTypesToDict(d, BACKSPACEKEY, SCA_IInputDevice::KX_BACKSPACEKEY);
	KX_MACRO_addTypesToDict(d, DELKEY, SCA_IInputDevice::KX_DELKEY);
	KX_MACRO_addTypesToDict(d, SEMICOLONKEY, SCA_IInputDevice::KX_SEMICOLONKEY);
	KX_MACRO_addTypesToDict(d, PERIODKEY, SCA_IInputDevice::KX_PERIODKEY);		
	KX_MACRO_addTypesToDict(d, COMMAKEY, SCA_IInputDevice::KX_COMMAKEY);		
	KX_MACRO_addTypesToDict(d, QUOTEKEY, SCA_IInputDevice::KX_QUOTEKEY);		
	KX_MACRO_addTypesToDict(d, ACCENTGRAVEKEY, SCA_IInputDevice::KX_ACCENTGRAVEKEY);	
	KX_MACRO_addTypesToDict(d, MINUSKEY, SCA_IInputDevice::KX_MINUSKEY);		
	KX_MACRO_addTypesToDict(d, SLASHKEY, SCA_IInputDevice::KX_SLASHKEY);		
	KX_MACRO_addTypesToDict(d, BACKSLASHKEY, SCA_IInputDevice::KX_BACKSLASHKEY);
	KX_MACRO_addTypesToDict(d, EQUALKEY, SCA_IInputDevice::KX_EQUALKEY);		
	KX_MACRO_addTypesToDict(d, LEFTBRACKETKEY, SCA_IInputDevice::KX_LEFTBRACKETKEY);	
	KX_MACRO_addTypesToDict(d, RIGHTBRACKETKEY, SCA_IInputDevice::KX_RIGHTBRACKETKEY);	
		
	KX_MACRO_addTypesToDict(d, LEFTARROWKEY, SCA_IInputDevice::KX_LEFTARROWKEY);
	KX_MACRO_addTypesToDict(d, DOWNARROWKEY, SCA_IInputDevice::KX_DOWNARROWKEY);
	KX_MACRO_addTypesToDict(d, RIGHTARROWKEY, SCA_IInputDevice::KX_RIGHTARROWKEY);	
	KX_MACRO_addTypesToDict(d, UPARROWKEY, SCA_IInputDevice::KX_UPARROWKEY);		
	
	KX_MACRO_addTypesToDict(d, PAD2	, SCA_IInputDevice::KX_PAD2);
	KX_MACRO_addTypesToDict(d, PAD4	, SCA_IInputDevice::KX_PAD4);
	KX_MACRO_addTypesToDict(d, PAD6	, SCA_IInputDevice::KX_PAD6);
	KX_MACRO_addTypesToDict(d, PAD8	, SCA_IInputDevice::KX_PAD8);
		
	KX_MACRO_addTypesToDict(d, PAD1	, SCA_IInputDevice::KX_PAD1);
	KX_MACRO_addTypesToDict(d, PAD3	, SCA_IInputDevice::KX_PAD3);
	KX_MACRO_addTypesToDict(d, PAD5	, SCA_IInputDevice::KX_PAD5);
	KX_MACRO_addTypesToDict(d, PAD7	, SCA_IInputDevice::KX_PAD7);
	KX_MACRO_addTypesToDict(d, PAD9	, SCA_IInputDevice::KX_PAD9);
		
	KX_MACRO_addTypesToDict(d, PADPERIOD, SCA_IInputDevice::KX_PADPERIOD);
	KX_MACRO_addTypesToDict(d, PADSLASHKEY, SCA_IInputDevice::KX_PADSLASHKEY);
	KX_MACRO_addTypesToDict(d, PADASTERKEY, SCA_IInputDevice::KX_PADASTERKEY);
		
		
	KX_MACRO_addTypesToDict(d, PAD0, SCA_IInputDevice::KX_PAD0);
	KX_MACRO_addTypesToDict(d, PADMINUS, SCA_IInputDevice::KX_PADMINUS);
	KX_MACRO_addTypesToDict(d, PADENTER, SCA_IInputDevice::KX_PADENTER);
	KX_MACRO_addTypesToDict(d, PADPLUSKEY, SCA_IInputDevice::KX_PADPLUSKEY);
		
		
	KX_MACRO_addTypesToDict(d, F1KEY , SCA_IInputDevice::KX_F1KEY);
	KX_MACRO_addTypesToDict(d, F2KEY , SCA_IInputDevice::KX_F2KEY);
	KX_MACRO_addTypesToDict(d, F3KEY , SCA_IInputDevice::KX_F3KEY);
	KX_MACRO_addTypesToDict(d, F4KEY , SCA_IInputDevice::KX_F4KEY);
	KX_MACRO_addTypesToDict(d, F5KEY , SCA_IInputDevice::KX_F5KEY);
	KX_MACRO_addTypesToDict(d, F6KEY , SCA_IInputDevice::KX_F6KEY);
	KX_MACRO_addTypesToDict(d, F7KEY , SCA_IInputDevice::KX_F7KEY);
	KX_MACRO_addTypesToDict(d, F8KEY , SCA_IInputDevice::KX_F8KEY);
	KX_MACRO_addTypesToDict(d, F9KEY , SCA_IInputDevice::KX_F9KEY);
	KX_MACRO_addTypesToDict(d, F10KEY, SCA_IInputDevice::KX_F10KEY);
	KX_MACRO_addTypesToDict(d, F11KEY, SCA_IInputDevice::KX_F11KEY);
	KX_MACRO_addTypesToDict(d, F12KEY, SCA_IInputDevice::KX_F12KEY);
		
	KX_MACRO_addTypesToDict(d, PAUSEKEY, SCA_IInputDevice::KX_PAUSEKEY);
	KX_MACRO_addTypesToDict(d, INSERTKEY, SCA_IInputDevice::KX_INSERTKEY);
	KX_MACRO_addTypesToDict(d, HOMEKEY , SCA_IInputDevice::KX_HOMEKEY);
	KX_MACRO_addTypesToDict(d, PAGEUPKEY, SCA_IInputDevice::KX_PAGEUPKEY);
	KX_MACRO_addTypesToDict(d, PAGEDOWNKEY, SCA_IInputDevice::KX_PAGEDOWNKEY);
	KX_MACRO_addTypesToDict(d, ENDKEY, SCA_IInputDevice::KX_ENDKEY);


	// Check for errors
	if (PyErr_Occurred())
    {
		Py_FatalError("can't initialize module GameKeys");
    }

	return d;
}

void PHY_SetActiveScene(class KX_Scene* scene)
{
	gp_KetsjiScene = scene;
}
