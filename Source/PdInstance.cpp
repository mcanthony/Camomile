/*
// Copyright (c) 2015 Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include "PdInstance.h"
#include "PdPatch.h"

extern "C"
{
#include "m_pd.h"
#include "m_imp.h"
EXTERN  void pd_init(void);
#include "pdugen.h"
#include "c_pd.h"
}

namespace pd
{
    // ==================================================================================== //
    //                                          INSTANCE                                    //
    // ==================================================================================== //
    
    int Instance::s_sample_rate;
    std::mutex Instance::s_mutex;
    std::string Instance::s_console;
    
    Instance::Internal::Internal(std::string const& _name) :
    instance(pdinstance_new()),
    counter(1),
    name(_name),
    inputs(nullptr),
    outputs(nullptr)
    {
        cpd_initialize((t_method)print, true, 2, true, false, false, false);
        static int initialized = 0;
        if(!initialized)
        {
            cpd_initialize_dsp();
            libpd_loadcream();
            s_sample_rate = sys_getsr();
            initialized = 1;
        }
    }
    
    Instance::Internal::~Internal()
    {
        patcher.clear();
        std::lock_guard<std::mutex> guard(s_mutex);
        if(instance)
        {
            pdinstance_free(instance);
        }
    }
    
    Instance::Instance() noexcept : m_internal(nullptr)
    {
        
    }
    
    Instance::Instance(std::string const& name) noexcept : m_internal(new Internal(name))
    {
        
    }
    
    Instance::Instance(Instance const& other) noexcept : m_internal(other.m_internal)
    {
        if(m_internal)
        {
            ++m_internal->counter;
        }
    }
    
    Instance::Instance(Instance&& other) noexcept : m_internal(other.m_internal)
    {
        other.m_internal = nullptr;
    }
    
    Instance& Instance::operator=(Instance const& other) noexcept
    {
        if(other.m_internal)
        {
            other.m_internal->counter++;
            m_internal = other.m_internal;
        }
        else
        {
            m_internal = nullptr;
        }
        return *this;
    }
    
    Instance& Instance::operator=(Instance&& other) noexcept
    {
        std::swap(m_internal, other.m_internal);
        return *this;
    }
    
    Instance::~Instance() noexcept
    {
        if(m_internal && m_internal->counter)
        {
            if(!(--m_internal->counter))
            {
                releaseDsp();
                delete m_internal;
            }
        }
    }
    
    void Instance::prepareDsp(const int nins, const int nouts, const int samplerate, const int nsamples) noexcept
    {
        releaseDsp();
        std::lock_guard<std::mutex> guard(m_internal->mutex);
        std::lock_guard<std::mutex> guard2(s_mutex);
        pd_setinstance(m_internal->instance);
        t_sample *tempin = nullptr, *tempout = nullptr;
        if(m_internal->inputs)
        {
            free(m_internal->inputs);
        }
        if(m_internal->outputs)
        {
            free(m_internal->outputs);
        }
        m_internal->inputs = (t_sample *)malloc((size_t)(16 * nsamples) * sizeof(t_sample));
        memset(m_internal->inputs, 0, (size_t)(nins * nsamples) * sizeof(t_sample));
        m_internal->outputs = (t_sample *)malloc((size_t)(16 * nsamples) * sizeof(t_sample));
        memset(m_internal->outputs, 0, (size_t)(nins * nsamples) * sizeof(t_sample));
        
        pdinstance_top_dsp(m_internal->instance);

        
        if(s_sample_rate != samplerate)
        {
            int indev[MAXAUDIOINDEV], inch[MAXAUDIOINDEV],
            outdev[MAXAUDIOOUTDEV], outch[MAXAUDIOOUTDEV];
            indev[0] = outdev[0] = DEFAULTAUDIODEV;
            inch[0] = s_max_channels;
            outch[0] = s_max_channels;
            sys_set_audio_settings(1, indev, 1, inch,
                                   1, outdev, 1, outch, samplerate, -1, 1, DEFDACBLKSIZE);
            sched_set_using_audio(SCHED_AUDIO_CALLBACK);
            sys_reopen_audio();
            s_sample_rate = sys_getsr();
        }
        tempin  = sys_soundin;
        tempout = sys_soundout;
        sys_soundin     = m_internal->inputs;
        sys_soundout    = m_internal->outputs;
        
        post("prepareDsp start -------------\n");
        post("prepareDsp instance %ld", (long)m_internal->instance);
        post("prepareDsp inputs %ld", (long)m_internal->inputs);
        post("prepareDsp outputs %ld", (long)m_internal->outputs);
        post("prepareDsp real inputs %ld", (long)tempin);
        post("prepareDsp real outputs %ld", (long)tempout);
        
        pdinstance_start_dsp(m_internal->instance);
        
        sys_soundin     = tempin;
        sys_soundout    = tempout;
        post("prepareDsp end -------------\n");
        
    }
    
    void Instance::performDsp(int nsamples, const int nins, const float** inputs, const int nouts, float** outputs) noexcept
    {
        std::lock_guard<std::mutex> guard2(s_mutex);
        std::lock_guard<std::mutex> guard(m_internal->mutex);
        t_sample* ins = m_internal->inputs;
        t_sample* outs= m_internal->outputs;
        const int blksize = sys_getblksize();
        for(int i = 0; i < nsamples; i += blksize)
        {
            for(int j = 0; j < nins; j++)
            {
                memcpy(ins+j*blksize, inputs[j]+i, blksize * sizeof(t_sample));
            }
            memset(outs, 0, blksize * sizeof(t_sample) * nouts);
            pdinstance_sched_tick(m_internal->instance, (double)blksize, (double)sys_getsr());
            for(int j = 0; j < nouts; j++)
            {
                memcpy(outputs[j]+i, outs+j*blksize, blksize * sizeof(t_sample));
            }
        }
    }
    
    void Instance::releaseDsp() noexcept
    {
        std::lock_guard<std::mutex> guard2(s_mutex);
        std::lock_guard<std::mutex> guard(m_internal->mutex);
        if(m_internal->inputs)
        {
            //free(m_internal->inputs);
        }
        if(m_internal->outputs)
        {
            //free(m_internal->outputs);
        }
        m_internal->inputs = nullptr;
        m_internal->outputs = nullptr;
    }
    
    void Instance::addToSearchPath(std::string const& path) noexcept
    {
        std::lock_guard<std::mutex> guard(s_mutex);
        sys_searchpath = namelist_append(sys_searchpath, path.c_str(), 0);
    }
    
    
    void Instance::clearSearchPath() noexcept
    {
        std::lock_guard<std::mutex> guard(s_mutex);
        namelist_free(sys_searchpath);
        sys_searchpath = NULL;
    }
    
    void Instance::setConsole(std::string const& text) noexcept
    {
        s_console = text;
    }
    
    std::string Instance::getConsole() noexcept
    {
        return s_console;
    }
    
    void Instance::print(const char* s)
    {
        s_console.append(s);
        t_symbol* send = gensym("camo-console");
        if(send->s_thing)
        {
            pd_bang(send->s_thing);
        }
        std::cout << s;
    }
}



