/*
// Copyright (c) 2015 Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include "pdugen.h"
#include "m_imp.h"

typedef void (*t_instanceclockmethod)(void *client);

typedef struct _instanceclock
{
    double                  c_settime;
    void*                   c_owner;
    t_instanceclockmethod   c_fn;
    struct _clock*          c_next;
    t_float                 c_unit;
}t_instanceclock;

static void pdinstance_dsp_tick(t_pdinstance* instance)
{
    if(instance->pd_dspchain)
    {
        t_int *ip;
        for(ip = instance->pd_dspchain; ip; )
        {
            ip = (*(t_perfroutine)(*ip))(ip);
        }
    }
}

void pdinstance_sched_tick(t_pdinstance* instance)
{
    int todo;
    const double next_sys_time = pd_this->pd_systime;// + sys_time_per_dsp_tick;
    while(pd_this->pd_clock_setlist && ((t_instanceclock *)pd_this->pd_clock_setlist)->c_settime < next_sys_time)
    {
        t_instanceclock *c = (t_instanceclock *)pd_this->pd_clock_setlist;
        pd_this->pd_systime = c->c_settime;
        clock_unset(pd_this->pd_clock_setlist);
        outlet_setstacklim();
        (*c->c_fn)(c->c_owner);
    }
    pd_this->pd_systime = next_sys_time;
    pdinstance_dsp_tick(pd_this);
}


