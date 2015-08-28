/*
// Copyright (c) 2015 Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef __CAMOMILE_PD_UGEN__
#define __CAMOMILE_PD_UGEN__

#include "m_pd.h"

void pdinstance_sched_prepare(t_pdinstance* instance, int nins, int nouts, double blocksize, double samplerate);

void pdinstance_sched_tick(t_pdinstance* instance, double blocksize, double samplerate);

void pdinstance_start_dsp(t_pdinstance* instance);
void pdinstance_top_dsp(t_pdinstance* instance);

#endif
