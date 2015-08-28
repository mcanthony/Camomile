/*
// Copyright (c) 2015 Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include "c_pd.h"
#include "s_stuff.h"

EXTERN void pd_init(void);

void cpd_initialize(t_method print, char verbose, char debug, char loadbang, char gui, char usestdpath, char usestderr)
{
    static char initialized = 0;
    if(!initialized)
    {
        // Default signal and midi to NULL
        sys_soundin     = NULL;
        sys_soundout    = NULL;
        sys_nmidiin     = 0;
        sys_nmidiout    = 0;
        sys_searchpath  = NULL;
        
        sys_printhook       = (t_printhook)print;
        sys_verbose         = (verbose) ? 1 : 0;
        sys_debuglevel      = (debug >= 0 && debug < 5) ? debug : 0;
        sys_noloadbang      = (loadbang == 0) ? 1 : 0;
        sys_nogui           = (gui == 0) ? 1 : 0;
        sys_usestdpath      = (usestdpath) ? 1 : 0;
        sys_printtostderr   = (usestderr) ? 1 : 0;
        
        sys_set_audio_api(API_DUMMY);
        sys_init_fdpoll();
        pd_init();
        
        // Don't really understand the function
        sys_schedblocksize      = DEFDACBLKSIZE;
        sys_hipriority          = 0;
        sys_externalschedlib    = 0;
#ifdef HAVE_SCHED_TICK_ARG
        sys_time                = 0;
#endif
        //signal(SIGFPE, SIG_IGN);
        initialized = 1;
    }
}

void cpd_initialize_dsp()
{
    int indev[MAXAUDIOINDEV], inch[MAXAUDIOINDEV],
    outdev[MAXAUDIOOUTDEV], outch[MAXAUDIOOUTDEV];
    indev[0] = outdev[0] = DEFAULTAUDIODEV;
    inch[0] = 16;
    outch[0] = 16;
    sys_set_audio_settings(1, indev, 1, inch,
                           1, outdev, 1, outch, 44100, -1, 1, DEFDACBLKSIZE);
    sched_set_using_audio(SCHED_AUDIO_CALLBACK);
    sys_reopen_audio();
}

void ibinbuf_evalfile(t_symbol *name, t_symbol *dir)
{
    t_binbuf *b = binbuf_new();
    glob_setfilename(0, name, dir);
    if(binbuf_read(b, name->s_name, dir->s_name, 0))
    {
        error("%s: read failed.", name->s_name);
    }
    else
    {
        t_pd *bounda = gensym("#A")->s_thing, *boundn = s__N.s_thing;
        gensym("#A")->s_thing = 0;
        s__N.s_thing = &pd_canvasmaker;
        binbuf_eval(b, 0, 0, 0);
        gensym("#A")->s_thing = bounda;
        s__N.s_thing = boundn;
    }
    glob_setfilename(0, &s_, &s_);
    binbuf_free(b);
}

t_canvas *pdinstance_newcanvas(t_symbol *name, t_symbol *dir)
{
    t_canvas *x = NULL;
    t_pd *boundx = s__X.s_thing;
    s__X.s_thing = NULL;
    ibinbuf_evalfile(name, dir);
    while ((x != (t_canvas *)s__X.s_thing) && s__X.s_thing)
    {
        x = (t_canvas *)s__X.s_thing;
        vmess((t_pd *)x, gensym("pop"), "i", 1);
    }
    pd_vmess((t_pd *)x, gensym("loadbang"), "");
    s__X.s_thing = boundx;
    return x;
}

