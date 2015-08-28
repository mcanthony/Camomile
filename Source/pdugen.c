/*
// Copyright (c) 2015 Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include "pdugen.h"
#include "m_imp.h"
#include "s_stuff.h"
#include "g_canvas.h"

static int idsp_phase;
struct _vinlet;
struct _voutlet;
extern t_class *vinlet_class, *voutlet_class, *canvas_class;
t_float *obj_findsignalscalar(t_object *x, int m);

typedef void (*t_iclockmethod)(void *client);

typedef struct _vinlet
{
    t_object x_obj;
    t_canvas *x_canvas;
    t_inlet *x_inlet;
    int x_bufsize;
    t_float *x_buf;         /* signal buffer; zero if not a signal */
    t_float *x_endbuf;
    t_float *x_fill;
    t_float *x_read;
    int x_hop;
    /* if not reblocking, the next slot communicates the parent's inlet
     signal from the prolog to the DSP routine: */
    t_signal *x_directsignal;
    
    t_resample x_updown;
} t_vinlet;

void vinlet_dspprolog(struct _vinlet *x, t_signal **parentsigs,
                      int myvecsize, int calcsize, int phase, int period, int frequency,
                      int downsample, int upsample,  int reblock, int switched);
void voutlet_dspprolog(struct _voutlet *x, t_signal **parentsigs,
                       int myvecsize, int calcsize, int phase, int period, int frequency,
                       int downsample, int upsample, int reblock, int switched);
void voutlet_dspepilog(struct _voutlet *x, t_signal **parentsigs,
                       int myvecsize, int calcsize, int phase, int period, int frequency,
                       int downsample, int upsample, int reblock, int switched);

typedef struct _block
{
    t_object x_obj;
    int x_vecsize;      /* size of audio signals in this block */
    int x_calcsize;     /* number of samples actually to compute */
    int x_overlap;
    int x_phase;        /* from 0 to period-1; when zero we run the block */
    int x_period;       /* submultiple of containing canvas */
    int x_frequency;    /* supermultiple of comtaining canvas */
    int x_count;        /* number of times parent block has called us */
    int x_chainonset;   /* beginning of code in DSP chain */
    int x_blocklength;  /* length of dspchain for this block */
    int x_epiloglength; /* length of epilog */
    char x_switched;    /* true if we're acting as a a switch */
    char x_switchon;    /* true if we're switched on */
    char x_reblock;     /* true if inlets and outlets are reblocking */
    int x_upsample;     /* upsampling-factor */
    int x_downsample;   /* downsampling-factor */
    int x_return;       /* stop right after this block (for one-shots) */
} t_block;

typedef struct _iclock
{
    double                  c_settime;
    void*                   c_owner;
    t_iclockmethod   c_fn;
    struct _iclock*  c_next;
    t_float                 c_unit;
}t_iclock;

typedef struct _ugenbox
{
    struct _siginlet *u_in;
    int u_nin;
    struct _sigoutlet *u_out;
    int u_nout;
    int u_phase;
    struct _ugenbox *u_next;
    t_object *u_obj;
    int u_done;
} t_ugenbox;

typedef struct _siginlet
{
    int i_nconnect;
    int i_ngot;
    t_signal *i_signal;
} t_siginlet;

typedef struct _sigoutconnect
{
    t_ugenbox *oc_who;
    int oc_inno;
    struct _sigoutconnect *oc_next;
} t_sigoutconnect;

typedef struct _sigoutlet
{
    int o_nconnect;
    int o_nsent;
    t_signal *o_signal;
    t_sigoutconnect *o_connections;
} t_sigoutlet;

typedef struct _dspcontext
{
    struct _ugenbox *dc_ugenlist;
    struct _dspcontext *dc_parentcontext;
    int dc_ninlets;
    int dc_noutlets;
    t_signal **dc_iosigs;
    t_float dc_srate;
    int dc_vecsize;         /* vector size, power of two */
    int dc_calcsize;        /* number of elements to calculate */
    char dc_toplevel;       /* true if "iosigs" is invalid. */
    char dc_reblock;        /* true if we have to reblock inlets/outlets */
    char dc_switched;       /* true if we're switched */
    
} t_dspcontext;

static t_signal *signal_freelist[MAXLOGSIG+1];
static t_signal *signal_freeborrowed;
static int ugen_sortno = 0;
static t_dspcontext *ugen_currentcontext;

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
    idsp_phase++;
}

static void iclock_unset(t_pdinstance* instance, t_iclock *x)
{
    if (x->c_settime >= 0)
    {
        if (x == (t_iclock *)instance->pd_clock_setlist)
            instance->pd_clock_setlist = (t_clock *)x->c_next;
        else
        {
            t_iclock *x2 = (t_iclock *)instance->pd_clock_setlist;
            while (x2->c_next != x) x2 = x2->c_next;
            x2->c_next = x->c_next;
        }
        x->c_settime = -1;
    }
}

void pdinstance_sched_prepare(t_pdinstance* instance, int nins, int nouts, double blocksize, double samplerate)
{
    /*
    t_canvas *x;
    if (instance->pd_dspstate) ugen_stop();
    else sys_gui("pdtk_pd_dsp ON\n");
    ugen_start();
    
    for(x = instance->pd_canvaslist; x; x = x->gl_next)
    {
        canvas_dodsp(x, 1, 0);
    }
    
    canvas_dspstate = instance->pd_dspstate = 1;
     */
}

void pdinstance_sched_tick(t_pdinstance* instance, double blocksize, double samplerate)
{
    const double next_sys_time = instance->pd_systime + 320. * blocksize;
    while(instance->pd_clock_setlist && ((t_iclock *)instance->pd_clock_setlist)->c_settime < next_sys_time)
    {
        t_iclock *c = (t_iclock *)instance->pd_clock_setlist;
        instance->pd_systime = c->c_settime;
        iclock_unset(instance, (t_iclock *)instance->pd_clock_setlist);
        outlet_setstacklim();
        (*c->c_fn)(c->c_owner);
    }
    instance->pd_systime = next_sys_time;
    pdinstance_dsp_tick(instance);
}

///// CANVAS ////

#define PROLOGCALL 2
#define EPILOGCALL 2

static int signal_compatible(t_signal *s1, t_signal *s2)
{
    return (s1->s_n == s2->s_n && s1->s_sr == s2->s_sr);
}

static t_int *block_epilog(t_int *w)
{
    t_block *x = (t_block *)w[1];
    int count = x->x_count - 1;
    if (x->x_return)
        return (0);
    if (!x->x_reblock)
        return (w + x->x_epiloglength + EPILOGCALL);
    if (count)
    {
        x->x_count = count;
        return (w - (x->x_blocklength -
                     (PROLOGCALL + EPILOGCALL)));   /* go to ugen after prolog */
    }
    else return (w + EPILOGCALL);
}
static t_int *iblock_prolog(t_int *w)
{
    t_block *x = (t_block *)w[1];
    int phase = x->x_phase;
    /* if we're switched off, jump past the epilog code */
    if (!x->x_switchon)
        return (w + x->x_blocklength);
    if (phase)
    {
        phase++;
        if (phase == x->x_period) phase = 0;
        x->x_phase = phase;
        return (w + x->x_blocklength);  /* skip block; jump past epilog */
    }
    else
    {
        x->x_count = x->x_frequency;
        x->x_phase = (x->x_period > 1 ? 1 : 0);
        return (w + PROLOGCALL);        /* beginning of block is next ugen */
    }
}

static void isignal_setborrowed(t_signal *sig, t_signal *sig2)
{
    if (!sig->s_isborrowed || sig->s_borrowedfrom)
        bug("signal_setborrowed");
    if (sig == sig2)
        bug("signal_setborrowed 2");
    sig->s_borrowedfrom = sig2;
    sig->s_vec = sig2->s_vec;
    sig->s_n = sig2->s_n;
    sig->s_vecsize = sig2->s_vecsize;
}

void isignal_makereusable(t_signal *sig)
{
    int logn = ilog2(sig->s_vecsize);
#if 1
    t_signal *s5;
    for (s5 = signal_freeborrowed; s5; s5 = s5->s_nextfree)
    {
        if (s5 == sig)
        {
            bug("signal_free 3");
            return;
        }
    }
    for (s5 = signal_freelist[logn]; s5; s5 = s5->s_nextfree)
    {
        if (s5 == sig)
        {
            bug("signal_free 4");
            return;
        }
    }
#endif
    if (sig->s_isborrowed)
    {
        /* if the signal is borrowed, decrement the borrowed-from signal's
         reference count, possibly marking it reusable too */
        t_signal *s2 = sig->s_borrowedfrom;
        if ((s2 == sig) || !s2)
            bug("signal_free");
        s2->s_refcount--;
        if (!s2->s_refcount)
            isignal_makereusable(s2);
        sig->s_nextfree = signal_freeborrowed;
        signal_freeborrowed = sig;
    }
    else
    {
        /* if it's a real signal (not borrowed), put it on the free list
         so we can reuse it. */
        if (signal_freelist[logn] == sig) bug("signal_free 2");
        sig->s_nextfree = signal_freelist[logn];
        signal_freelist[logn] = sig;
    }
}

t_signal *isignal_new(t_pdinstance* instance, int n, t_float sr)
{
    int logn, n2, vecsize = 0;
    t_signal *ret, **whichlist;
    t_sample *fp;
    logn = ilog2(n);
    if (n)
    {
        if ((vecsize = (1<<logn)) != n)
            vecsize *= 2;
        if (logn > MAXLOGSIG)
            bug("signal buffer too large");
        whichlist = signal_freelist + logn;
    }
    else
        whichlist = &signal_freeborrowed;
    
    /* first try to reclaim one from the free list */
    if ((ret = *whichlist))
        *whichlist = ret->s_nextfree;
    else
    {
        /* LATER figure out what to do for out-of-space here! */
        ret = (t_signal *)t_getbytes(sizeof *ret);
        if (n)
        {
            ret->s_vec = (t_sample *)getbytes(vecsize * sizeof (*ret->s_vec));
            ret->s_isborrowed = 0;
        }
        else
        {
            ret->s_vec = 0;
            ret->s_isborrowed = 1;
        }
        ret->s_nextused = instance->pd_signals;
        instance->pd_signals = ret;
    }
    ret->s_n = n;
    ret->s_vecsize = vecsize;
    ret->s_sr = sr;
    ret->s_refcount = 0;
    ret->s_borrowedfrom = 0;
    
    return (ret);
}

static t_signal *isignal_newlike(t_pdinstance* instance, const t_signal *sig)
{
    return (isignal_new(instance,sig->s_n, sig->s_sr));
}

static t_dspcontext *pdinstance_ugen_start_graph(t_pdinstance* instance, int toplevel, t_signal **sp, int ninlets, int noutlets)
{
    t_dspcontext *dc = (t_dspcontext *)getbytes(sizeof(*dc));
    t_float parent_srate, srate;
    int parent_vecsize, vecsize;
    
    if(toplevel)
        ninlets=noutlets=0;
    
    dc->dc_ugenlist = 0;
    dc->dc_toplevel = toplevel;
    dc->dc_iosigs = sp;
    dc->dc_ninlets = ninlets;
    dc->dc_noutlets = noutlets;
    dc->dc_parentcontext = ugen_currentcontext;
    ugen_currentcontext = dc;
    return (dc);
}

static void pdinstance_ugen_add(t_pdinstance* instance, t_dspcontext *dc, t_object *obj)
{
    t_ugenbox *x = (t_ugenbox *)getbytes(sizeof *x);
    int i;
    t_sigoutlet *uout;
    t_siginlet *uin;
    
    x->u_next = dc->dc_ugenlist;
    dc->dc_ugenlist = x;
    x->u_obj = obj;
    x->u_nin = obj_nsiginlets(obj);
    x->u_in = getbytes(x->u_nin * sizeof (*x->u_in));
    for (uin = x->u_in, i = x->u_nin; i--; uin++)
        uin->i_nconnect = 0;
    x->u_nout = obj_nsigoutlets(obj);
    x->u_out = getbytes(x->u_nout * sizeof (*x->u_out));
    for (uout = x->u_out, i = x->u_nout; i--; uout++)
        uout->o_connections = 0, uout->o_nconnect = 0;
}

static void pdinstance_ugen_connect(t_pdinstance* instance, t_dspcontext *dc, t_object *x1, int outno, t_object *x2,
                  int inno)
{
    t_ugenbox *u1, *u2;
    t_sigoutlet *uout;
    t_siginlet *uin;
    t_sigoutconnect *oc;
    int sigoutno = obj_sigoutletindex(x1, outno);
    int siginno = obj_siginletindex(x2, inno);

    for (u1 = dc->dc_ugenlist; u1 && u1->u_obj != x1; u1 = u1->u_next);
    for (u2 = dc->dc_ugenlist; u2 && u2->u_obj != x2; u2 = u2->u_next);
    if (!u1 || !u2 || siginno < 0)
    {
        pd_error(u1->u_obj,
                 "signal outlet connect to nonsignal inlet (ignored)");
        return;
    }
    if (sigoutno < 0 || sigoutno >= u1->u_nout || siginno >= u2->u_nin)
    {
        bug("ugen_connect %s %s %d %d (%d %d)",
            class_getname(x1->ob_pd),
            class_getname(x2->ob_pd), sigoutno, siginno, u1->u_nout,
            u2->u_nin);
    }
    uout = u1->u_out + sigoutno;
    uin = u2->u_in + siginno;
    
    /* add a new connection to the outlet's list */
    oc = (t_sigoutconnect *)getbytes(sizeof *oc);
    oc->oc_next = uout->o_connections;
    uout->o_connections = oc;
    oc->oc_who = u2;
    oc->oc_inno = siginno;
    /* update inlet and outlet counts  */
    uout->o_nconnect++;
    uin->i_nconnect++;
}

static void pdinstance_ugen_doit(t_pdinstance* instance, t_dspcontext *dc, t_ugenbox *u)
{
    t_sigoutlet *uout;
    t_siginlet *uin;
    t_sigoutconnect *oc, *oc2;
    t_class *class = pd_class(&u->u_obj->ob_pd);
    int i, n;
    /* suppress creating new signals for the outputs of signal
     inlets and subpatchs; except in the case we're an inlet and "blocking"
     is set.  We don't yet know if a subcanvas will be "blocking" so there
     we delay new signal creation, which will be handled by calling
     signal_setborrowed in the ugen_done_graph routine below. */
    int nonewsigs = (class == canvas_class ||
                     ((class == vinlet_class) && !(dc->dc_reblock)));
    /* when we encounter a subcanvas or a signal outlet, suppress freeing
     the input signals as they may be "borrowed" for the super or sub
     patch; same exception as above, but also if we're "switched" we
     have to do a copy rather than a borrow.  */
    int nofreesigs = (class == canvas_class ||
                      ((class == voutlet_class) &&  !(dc->dc_reblock || dc->dc_switched)));
    t_signal **insig, **outsig, **sig, *s1, *s2, *s3;
    t_ugenbox *u2;
    
    for (i = 0, uin = u->u_in; i < u->u_nin; i++, uin++)
    {
        if (!uin->i_nconnect)
        {
            t_float *scalar;
            s3 = isignal_new(instance, dc->dc_calcsize, dc->dc_srate);
            /* post("%s: unconnected signal inlet set to zero",
             class_getname(u->u_obj->ob_pd)); */
            if ((scalar = obj_findsignalscalar(u->u_obj, i)))
                dsp_add_scalarcopy(scalar, s3->s_vec, s3->s_n);
            else
                dsp_add_zero(s3->s_vec, s3->s_n);
            uin->i_signal = s3;
            s3->s_refcount = 1;
        }
    }
    insig = (t_signal **)getbytes((u->u_nin + u->u_nout) * sizeof(t_signal *));
    outsig = insig + u->u_nin;
    for (sig = insig, uin = u->u_in, i = u->u_nin; i--; sig++, uin++)
    {
        int newrefcount;
        *sig = uin->i_signal;
        newrefcount = --(*sig)->s_refcount;
        /* if the reference count went to zero, we free the signal now,
         unless it's a subcanvas or outlet; these might keep the
         signal around to send to objects connected to them.  In this
         case we increment the reference count; the corresponding decrement
         is in sig_makereusable(). */
        if (nofreesigs)
            (*sig)->s_refcount++;
        else if (!newrefcount)
            isignal_makereusable(*sig);
    }
    for (sig = outsig, uout = u->u_out, i = u->u_nout; i--; sig++, uout++)
    {
        /* similarly, for outlets of subcanvases we delay creating
         them; instead we create "borrowed" ones so that the refcount
         is known.  The subcanvas replaces the fake signal with one showing
         where the output data actually is, to avoid having to copy it.
         For any other object, we just allocate a new output vector;
         since we've already freed the inputs the objects might get called
         "in place." */
        if (nonewsigs)
        {
            *sig = uout->o_signal =
            isignal_new(instance, 0, dc->dc_srate);
        }
        else
            *sig = uout->o_signal = isignal_new(instance, dc->dc_calcsize, dc->dc_srate);
        (*sig)->s_refcount = uout->o_nconnect;
    }
    /* now call the DSP scheduling routine for the ugen.  This
     routine must fill in "borrowed" signal outputs in case it's either
     a subcanvas or a signal inlet. */
    mess1(&u->u_obj->ob_pd, gensym("dsp"), insig);
    
    /* if any output signals aren't connected to anyone, free them
     now; otherwise they'll either get freed when the reference count
     goes back to zero, or even later as explained above. */
    
    for (sig = outsig, uout = u->u_out, i = u->u_nout; i--; sig++, uout++)
    {
        if (!(*sig)->s_refcount)
            isignal_makereusable(*sig);
    }
    
    /* pass it on and trip anyone whose last inlet was filled */
    for (uout = u->u_out, i = u->u_nout; i--; uout++)
    {
        s1 = uout->o_signal;
        for (oc = uout->o_connections; oc; oc = oc->oc_next)
        {
            u2 = oc->oc_who;
            uin = &u2->u_in[oc->oc_inno];
            /* if there's already someone here, sum the two */
            if ((s2 = uin->i_signal))
            {
                s1->s_refcount--;
                s2->s_refcount--;
                if (!signal_compatible(s1, s2))
                {
                    pd_error(u->u_obj, "%s: incompatible signal inputs",
                             class_getname(u->u_obj->ob_pd));
                    return;
                }
                s3 = isignal_newlike(instance, s1);
                dsp_add_plus(s1->s_vec, s2->s_vec, s3->s_vec, s1->s_n);
                uin->i_signal = s3;
                s3->s_refcount = 1;
                if (!s1->s_refcount) isignal_makereusable(s1);
                if (!s2->s_refcount) isignal_makereusable(s2);
            }
            else uin->i_signal = s1;
            uin->i_ngot++;
            /* if we didn't fill this inlet don't bother yet */
            if (uin->i_ngot < uin->i_nconnect)
                goto notyet;
            /* if there's more than one, check them all */
            if (u2->u_nin > 1)
            {
                for (uin = u2->u_in, n = u2->u_nin; n--; uin++)
                    if (uin->i_ngot < uin->i_nconnect) goto notyet;
            }
            /* so now we can schedule the ugen.  */
            pdinstance_ugen_doit(instance, dc, u2);
        notyet: ;
        }
    }
    t_freebytes(insig,(u->u_nin + u->u_nout) * sizeof(t_signal *));
    u->u_done = 1;
}

void pdinstance_ugen_done_graph(t_pdinstance* instance, t_dspcontext *dc)
{
    t_ugenbox *u, *u2;
    t_sigoutlet *uout;
    t_siginlet *uin;
    t_sigoutconnect *oc, *oc2;
    int i, n;
    t_block *blk;
    t_dspcontext *parent_context = dc->dc_parentcontext;
    t_float parent_srate;
    int parent_vecsize;
    int period, frequency, phase, vecsize, calcsize;
    t_float srate;
    int chainblockbegin;    /* DSP chain onset before block prolog code */
    int chainblockend;      /* and after block epilog code */
    int chainafterall;      /* and after signal outlet epilog */
    int reblock = 0, switched;
    int downsample = 1, upsample = 1;

    /* search for an object of class "block~" */
    for (u = dc->dc_ugenlist, blk = 0; u; u = u->u_next)
    {
        t_pd *zz = &u->u_obj->ob_pd;
        if (pd_class(zz)->c_name == gensym("block~"))
        {
            if(blk)
                pd_error(blk, "conflicting block~ objects in same page");
            else blk = (t_block *)zz;
        }
    }
    
    /* figure out block size, calling frequency, sample rate */
    if (parent_context)
    {
        parent_srate = parent_context->dc_srate;
        parent_vecsize = parent_context->dc_vecsize;
    }
    else
    {
        parent_srate = sys_getsr();
        parent_vecsize = sys_getblksize();
    }
    if (blk)
    {
        int realoverlap;
        vecsize = blk->x_vecsize;
        if (vecsize == 0)
            vecsize = parent_vecsize;
        calcsize = blk->x_calcsize;
        if (calcsize == 0)
            calcsize = vecsize;
        realoverlap = blk->x_overlap;
        if (realoverlap > vecsize) realoverlap = vecsize;
        downsample = blk->x_downsample;
        upsample   = blk->x_upsample;
        if (downsample > parent_vecsize)
            downsample = parent_vecsize;
        period = (vecsize * downsample)/
        (parent_vecsize * realoverlap * upsample);
        frequency = (parent_vecsize * realoverlap * upsample)/
        (vecsize * downsample);
        phase = blk->x_phase;
        srate = parent_srate * realoverlap * upsample / downsample;
        if (period < 1) period = 1;
        if (frequency < 1) frequency = 1;
        blk->x_frequency = frequency;
        blk->x_period = period;
        blk->x_phase = idsp_phase & (period - 1);
        if (! parent_context || (realoverlap != 1) ||
            (vecsize != parent_vecsize) ||
            (downsample != 1) || (upsample != 1))
            reblock = 1;
        switched = blk->x_switched;
    }
    else
    {
        srate = parent_srate;
        vecsize = parent_vecsize;
        calcsize = (parent_context ? parent_context->dc_calcsize : vecsize);
        downsample = upsample = 1;
        period = frequency = 1;
        phase = 0;
        if (!parent_context) reblock = 1;
        switched = 0;
    }
    dc->dc_reblock = reblock;
    dc->dc_switched = switched;
    dc->dc_srate = srate;
    dc->dc_vecsize = vecsize;
    dc->dc_calcsize = calcsize;
    
    /* if we're reblocking or switched, we now have to create output
     signals to fill in for the "borrowed" ones we have now.  This
     is also possibly true even if we're not blocked/switched, in
     the case that there was a signal loop.  But we don't know this
     yet.  */
    
    if (dc->dc_iosigs && (switched || reblock))
    {
        t_signal **sigp;
        for (i = 0, sigp = dc->dc_iosigs + dc->dc_ninlets; i < dc->dc_noutlets;
             i++, sigp++)
        {
            if ((*sigp)->s_isborrowed && !(*sigp)->s_borrowedfrom)
            {
                isignal_setborrowed(*sigp, isignal_new(instance, parent_vecsize, parent_srate));
                (*sigp)->s_refcount++;
            }
        }
    }
    
    /* schedule prologs for inlets and outlets.  If the "reblock" flag
     is set, an inlet will put code on the DSP chain to copy its input
     into an internal buffer here, before any unit generators' DSP code
     gets scheduled.  If we don't "reblock", inlets will need to get
     pointers to their corresponding inlets/outlets on the box we're inside,
     if any.  Outlets will also need pointers, unless we're switched, in
     which case outlet epilog code will kick in. */
    
    for (u = dc->dc_ugenlist; u; u = u->u_next)
    {
        t_pd *zz = &u->u_obj->ob_pd;
        t_signal **insigs = dc->dc_iosigs, **outsigs = dc->dc_iosigs;
        if (outsigs) outsigs += dc->dc_ninlets;
        
        if (pd_class(zz) == vinlet_class)
            vinlet_dspprolog((struct _vinlet *)zz,
                             dc->dc_iosigs, vecsize, calcsize, idsp_phase, period, frequency,
                             downsample, upsample, reblock, switched);
        else if (pd_class(zz) == voutlet_class)
            voutlet_dspprolog((struct _voutlet *)zz,
                              outsigs, vecsize, calcsize, idsp_phase, period, frequency,
                              downsample, upsample, reblock, switched);
    }
    chainblockbegin = instance->pd_dspchainsize;
    
    if (blk && (reblock || switched))   /* add the block DSP prolog */
    {
        dsp_add(iblock_prolog, 1, blk);
        blk->x_chainonset = instance->pd_dspchainsize - 1;
    }
    /* Initialize for sorting */
    for (u = dc->dc_ugenlist; u; u = u->u_next)
    {
        u->u_done = 0;
        for (uout = u->u_out, i = u->u_nout; i--; uout++)
            uout->o_nsent = 0;
        for (uin = u->u_in, i = u->u_nin; i--; uin++)
            uin->i_ngot = 0, uin->i_signal = 0;
    }
    
    /* Do the sort */
    
    for (u = dc->dc_ugenlist; u; u = u->u_next)
    {
        /* check that we have no connected signal inlets */
        if (u->u_done) continue;
        for (uin = u->u_in, i = u->u_nin; i--; uin++)
            if (uin->i_nconnect) goto next;
        
        pdinstance_ugen_doit(instance, dc, u);
    next: ;
    }
    
    /* check for a DSP loop, which is evidenced here by the presence
     of ugens not yet scheduled. */
    
    for (u = dc->dc_ugenlist; u; u = u->u_next)
        if (!u->u_done)
        {
            t_signal **sigp;
            pd_error(u->u_obj,
                     "DSP loop detected (some tilde objects not scheduled)");
            /* this might imply that we have unfilled "borrowed" outputs
             which we'd better fill in now. */
            for (i = 0, sigp = dc->dc_iosigs + dc->dc_ninlets; i < dc->dc_noutlets;
                 i++, sigp++)
            {
                if ((*sigp)->s_isborrowed && !(*sigp)->s_borrowedfrom)
                {
                    t_signal *s3 = isignal_new(instance, parent_vecsize, parent_srate);
                    isignal_setborrowed(*sigp, s3);
                    (*sigp)->s_refcount++;
                    dsp_add_zero(s3->s_vec, s3->s_n);
                }
            }
            break;   /* don't need to keep looking. */
        }
    
    if (blk && (reblock || switched))    /* add block DSP epilog */
        dsp_add(block_epilog, 1, blk);
    chainblockend = instance->pd_dspchainsize;
    
    /* add epilogs for outlets.  */
    
    for (u = dc->dc_ugenlist; u; u = u->u_next)
    {
        t_pd *zz = &u->u_obj->ob_pd;
        if (pd_class(zz) == voutlet_class)
        {
            t_signal **iosigs = dc->dc_iosigs;
            if (iosigs) iosigs += dc->dc_ninlets;
            voutlet_dspepilog((struct _voutlet *)zz,
                              iosigs, vecsize, calcsize, idsp_phase, period, frequency,
                              downsample, upsample, reblock, switched);
        }
    }
    
    chainafterall = instance->pd_dspchainsize;
    if (blk)
    {
        blk->x_blocklength = chainblockend - chainblockbegin;
        blk->x_epiloglength = chainafterall - chainblockend;
        blk->x_reblock = reblock;
    }
    
    /* now delete everything. */
    while (dc->dc_ugenlist)
    {
        for (uout = dc->dc_ugenlist->u_out, n = dc->dc_ugenlist->u_nout;
             n--; uout++)
        {
            oc = uout->o_connections;
            while (oc)
            {
                oc2 = oc->oc_next;
                freebytes(oc, sizeof *oc);
                oc = oc2;
            }
        }
        freebytes(dc->dc_ugenlist->u_out, dc->dc_ugenlist->u_nout *
                  sizeof (*dc->dc_ugenlist->u_out));
        freebytes(dc->dc_ugenlist->u_in, dc->dc_ugenlist->u_nin *
                  sizeof(*dc->dc_ugenlist->u_in));
        u = dc->dc_ugenlist;
        dc->dc_ugenlist = u->u_next;
        freebytes(u, sizeof *u);
    }
    if (ugen_currentcontext == dc)
        ugen_currentcontext = dc->dc_parentcontext;
    else bug("ugen_currentcontext");
    freebytes(dc, sizeof(*dc));
    
}

static void pdinstance_canvas_dsp(t_pdinstance* instance, t_canvas *x, int toplevel, t_signal **sp)
{
    t_linetraverser t;
    t_outconnect *oc;
    t_gobj *y;
    t_object *ob;
    t_symbol *dspsym = gensym("dsp");
    t_dspcontext *dc;
    
    dc = pdinstance_ugen_start_graph(instance, toplevel, sp, obj_nsiginlets(&x->gl_obj), obj_nsigoutlets(&x->gl_obj));
    for (y = x->gl_list; y; y = y->g_next)
    {
        if ((ob = pd_checkobject(&y->g_pd)) && zgetfn(&y->g_pd, dspsym))
        {
            pdinstance_ugen_add(instance, dc, ob);
        }
    }
    
    linetraverser_start(&t, x);
    while((oc = linetraverser_next(&t)))
    {
        if(obj_issignaloutlet(t.tr_ob, t.tr_outno))
        {
            pdinstance_ugen_connect(instance, dc, t.tr_ob, t.tr_outno, t.tr_ob2, t.tr_inno);
        }
    }
    
    pdinstance_ugen_done_graph(instance, dc);
}

///// UGEN ////

static t_int pdinstance_dsp_done(t_int *w)
{
    return (0);
}

static void pdinstance_signal_cleanup(t_pdinstance* instance)
{
    t_signal **svec, *sig, *sig2;
    int i;
    sig = instance->pd_signals;
    while(sig)
    {
        sig = sig->s_nextused;
        if (!sig->s_isborrowed)
        {
            t_freebytes(sig->s_vec, sig->s_vecsize * sizeof (*sig->s_vec));
        }
        t_freebytes(sig, sizeof *sig);
    }
    for(i = 0; i <= MAXLOGSIG; i++)
    {
        signal_freelist[i] = 0;
    }
    signal_freeborrowed = 0;
}

static void pdinstance_ugen_stop(t_pdinstance* instance)
{
    t_signal *s;
    int i;
    if(instance->pd_dspchain)
    {
        freebytes(instance->pd_dspchain, instance->pd_dspchainsize * sizeof (t_int));
        instance->pd_dspchain = NULL;
    }
    pdinstance_signal_cleanup(instance);
}

static void pdinstance_ugen_start(t_pdinstance* instance)
{
    pdinstance_ugen_stop(instance);
    ugen_sortno++;
    instance->pd_dspchain = (t_int *)getbytes(sizeof(*instance->pd_dspchain));
    instance->pd_dspchain[0] = (t_int)pdinstance_dsp_done;
    instance->pd_dspchainsize = 1;
    if(ugen_currentcontext)
    {
        bug("ugen_start");
    }
}

void pdinstance_start_dsp(t_pdinstance* instance)
{
    t_canvas *x;
    if(instance->pd_dspstate)
    {
        pdinstance_ugen_stop(instance);
    }
    else
    {
        pdinstance_ugen_start(instance);
    }
    
    for(x = instance->pd_canvaslist; x; x = x->gl_next)
    {
        pdinstance_canvas_dsp(instance, x, 1, NULL);
    }
    
    canvas_dspstate = instance->pd_dspstate = 1;
}

 void pdinstance_top_dsp(t_pdinstance* instance)
{
    if(instance->pd_dspstate)
    {
        pdinstance_ugen_stop(instance);
        canvas_dspstate = instance->pd_dspstate = 0;
    }
}


