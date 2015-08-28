/*
// Copyright (c) 2015 Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef __C_PD__
#define __C_PD__

#include "m_pd.h"

/*!
 * \fn          void cpd_initialize(t_method print, char verbose, char debug, char loadbang, char gui, char usestdpath, char usestderr)
 * \brief       Intializes the Pure Data's environment.
 * \details     The function set up all the static variables and methods.
 * \param print         The print method to use (t_printhook).
 * \param verbose       If the system should post the debug informations.
 * \param debuglvl      The level of the debug informations. 0 for fatal, 1 for error, 2 for\n
                        normal, 3 for debug, 4 for all.
 * \param loadbang      If the system should send a loadbang message when a patch is open.
 * \param gui           If the system uses its internal Tcl/Tk GUI system (sys_gui, etc.).
 * \param usestdpath    If the system uses the standard path ("extra" folder, etc.).
 * \param usestderr     If the system use err sprintf in case the post method is no set up.
 * \see sys_vgui, sys_gui
 */
void cpd_initialize(t_method print, char verbose, char debug, char loadbang, char gui, char usestdpath, char usestderr);

/*!
 * \fn          void epd_add_folder(const char* name, const char* folder)
 * \brief       Adds a subfolder to library folder.
 * \details     The function initializes the environment of Pure Data.
 * \param print         The print method to use (t_printhook).
 * \param verbose       If the system should post the debug informations.
 * \param debuglvl      The level of the debug informations. 0 for fatal, 1 for error, 2 for\n
 normal, 3 for debug, 4 for all.
 * \param loadbang      If the system should send a loadbang message when a patch is open.
 * \param gui           If the system uses its internal Tcl/Tk GUI system (sys_gui, etc.).
 * \param usestdpath    If the system uses the standard path ("extra" folder, etc.).
 * \param usestderr     If the system use err sprintf in case the post method is no set up.
 * \see sys_vgui, sys_gui
 */
void cpd_initialize_dsp();

t_canvas *pdinstance_newcanvas(t_symbol *name, t_symbol *dir);



#endif
