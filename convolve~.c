/*

convolve~

Copyright 2010 William Brent

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

version 0.0.2, Dec 7, 2010

¥Êsee note in perform about gain reduction and stability
¥Êalso, this crashes kind of randomly when DSP is turned on after analysis

¥ 0.0.2, shooting for a small buffered block size (N).  0.0.1 requires that irlen==N in order to work.  Here, we make a longer nonoverlapped history so that the latest final_output can be summed based on that history.

*/

#include "m_pd.h"
#include <math.h>
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

static t_class *convolve_tilde_class;

typedef struct _convolve_tilde 
{
    t_object x_obj;
    t_symbol *arrayname;
    t_word *x_vec;
    int arraysize;
    float sr;
    float n;
    int dsp_tick;
    int buffer_limit;
    int overlap;
    int overlap_chunk;
    int buf_hist_size;
    int window;
    int window_half;
    t_float *ir_spec;
    t_float *hann;
	t_sample *signal_buf;
	t_sample *signal_buf_windowed;
	
	t_sample *nonoverlapped_output;
	t_sample *final_output;
	
    float x_f;
    
} convolve_tilde;


/* ------------------------ convolve~ -------------------------------- */

static void convolve_tilde_analyze(convolve_tilde *x, t_symbol *s)
{
    t_garray *a;
    int i;

    x->arrayname = s;
    if (!(a = (t_garray *)pd_findbyclass(x->arrayname, garray_class)))
    {
        if (*s->s_name) pd_error(x, "convolve~: %s: no such array",
            x->arrayname->s_name);
        x->x_vec = 0;
    }
    else if (!garray_getfloatwords(a, &x->arraysize, &x->x_vec))
    {
        pd_error(x, "%s: bad template for convolve~", x->arrayname->s_name);
        x->x_vec = 0;
    }
    else garray_usedindsp(a);
    
    // NOT CHECKING FOR POWER OF 2 TABLE SIZE - IT'S ASSUMED FOR NOW
    
    x->buf_hist_size = x->arraysize/x->overlap_chunk;
	post("buf hist size: %i", x->buf_hist_size);
    
    x->ir_spec = (t_float *)t_resizebytes(x->ir_spec, 0, x->arraysize*sizeof(t_float));
    x->signal_buf_windowed = (t_sample *)t_resizebytes(x->signal_buf_windowed, 0, x->arraysize*sizeof(t_sample));
    x->nonoverlapped_output = (t_sample *)t_resizebytes(x->nonoverlapped_output, 0, x->buf_hist_size*x->arraysize*sizeof(t_sample));

    // copy array internally for upcoming in-place FFT, and init other buffers
    for(i=0; i<x->arraysize; i++)
    {
    	x->ir_spec[i] = x->x_vec[i].w_float;
    	x->signal_buf_windowed[i] = 0.0;
    }

    for(i=0; i<x->buf_hist_size*x->arraysize; i++)
    	x->nonoverlapped_output[i] = 0.0;
    
    // take FT, not Hann windowed for now
	mayer_realfft(x->arraysize, x->ir_spec);
	
	post("analysis complete: %i samples", x->arraysize);
}


static void *convolve_tilde_new(t_floatarg window, t_floatarg overlap)
{
    convolve_tilde *x = (convolve_tilde *)pd_new(convolve_tilde_class);
	int i, isPow2;
	
	outlet_new(&x->x_obj, &s_signal);

	isPow2 = (int)window && !( ((int)window-1) & (int)window );
	overlap = overlap;
	
	x->window = 2048;
	x->overlap = 4;
	
	x->arraysize = 0;
	x->buf_hist_size = 0;
	x->sr = 44100;
	x->n = 64;
	x->dsp_tick = 0; 
	x->window_half = x->window*0.5;
	x->overlap_chunk = x->window/x->overlap;
	
	x->buffer_limit = x->window/x->n/x->overlap; // divide by overlap since we want to push out buffered audio to a window when the main buffer has been updated by 1/overlap the window size.


	// these will be resized when analysis occurs
	x->ir_spec = (t_float *)getbytes(0);
	x->signal_buf_windowed = (t_sample *)getbytes(0);
	x->nonoverlapped_output = (t_sample *)getbytes(0);
	
	x->hann = (t_float *)getbytes(x->window*sizeof(t_float));
	x->signal_buf = (t_sample *)getbytes(x->window*sizeof(t_sample));
	x->final_output = (t_sample *)getbytes(x->overlap_chunk*sizeof(t_sample));


	// initialize hann window
 	for(i=0; i<x->window; i++)
 		x->hann[i] = 0.5 * (1 - cos(2*M_PI*i/x->window));

	// init signal buffer
 	for(i=0; i<x->window; i++)
 		x->signal_buf[i] = 0.0;

	// init output buffer
	for(i=0; i<x->overlap_chunk; i++)
		x->final_output[i] = 0.0;
	
    post("convolve~: window size: %i, overlap: %i", x->window, x->overlap);
    
    return (x);
}


static t_int *convolve_tilde_perform(t_int *w)
{
    int i, j, n, window, window_half, ir_window, ir_window_half, buf_hist_size, overlap, overlap_chunk;
	float amp_scalar;
	
    convolve_tilde *x = (convolve_tilde *)(w[1]);

    t_sample *in = (t_float *)(w[2]);
    t_sample *out = (t_float *)(w[3]);
    
    n = w[4];
	window = x->window;
	window_half = x->window_half;
	ir_window = x->arraysize;
	ir_window_half = x->arraysize*0.5;
	buf_hist_size = x->buf_hist_size;
	overlap = x->overlap;
	overlap_chunk = x->overlap_chunk;
		
	// shift previous contents back
	for(i=0; i<(window-n); i++)
		x->signal_buf[i] = x->signal_buf[n+i];
		
	// buffer most recent block
	for(i=0; i<n; i++)
		x->signal_buf[window-n+i] = in[i];
		
	
	if(x->dsp_tick==x->buffer_limit)
	{
		x->dsp_tick = 0;
		amp_scalar = (2.0/3.0)/ir_window; // appropriate for an overlap of 4
		
		// window the signal
		for(i=0; i<window; i++)
			x->signal_buf_windowed[i] = x->signal_buf[i] * x->hann[i];

		// pad the latest x->window samples with zeros out to the length of the IR
		for(; i<ir_window; i++)
			x->signal_buf_windowed[i] = 0.0;
		
		// take FT of the padded window
		mayer_realfft(ir_window, x->signal_buf_windowed);

		// multiply against IR spectrum
		// this needs to be a complex multiply
		
		// DC and Nyquist are real only, hence normal multiplies
		x->signal_buf_windowed[0] *= x->ir_spec[0];
		x->signal_buf_windowed[ir_window_half] *= x->ir_spec[ir_window_half];
		
		for(i=1, j=ir_window-1; i<ir_window_half; i++, j--)
		{
			float real, imag;
			
			// MINUS the imag part because i^2 = -1
			real = (x->signal_buf_windowed[i] * x->ir_spec[i]) - (x->signal_buf_windowed[j] * x->ir_spec[j]);
			
			imag = (x->signal_buf_windowed[i] * x->ir_spec[j]) + (x->signal_buf_windowed[j] * x->ir_spec[i]);
			
			x->signal_buf_windowed[i] = real;
			x->signal_buf_windowed[j] = imag;
		}
		
		// resynth
		mayer_realifft(ir_window, x->signal_buf_windowed);
		
		// reduce gain
		for(i=0; i<ir_window; i++)
			x->signal_buf_windowed[i] *= amp_scalar;
			
		// shift nonoverlapped blocks back
		for(i=0; i<(buf_hist_size-1); i++)
			for(j=0; j<ir_window; j++)
				x->nonoverlapped_output[(i*ir_window)+j] = x->nonoverlapped_output[((i+1)*ir_window)+j];
		
		// write the new block
		for(i=0; i<ir_window; i++)
			x->nonoverlapped_output[((buf_hist_size-1)*ir_window)+i] = x->signal_buf_windowed[i];
			
		// init this chunk of the final output so it can be summed in the for() below
		for(i=0; i<overlap_chunk; i++)
			x->final_output[i] = 0.0;
			
		// do the overlap/add
		for(i=0; i<buf_hist_size; i++)
			for(j=0; j<overlap_chunk; j++)
				x->final_output[j] += x->nonoverlapped_output[(i*ir_window)+((buf_hist_size-i-1)*overlap_chunk)+j];
	};
	
	// output
	for(i=0; i<n; i++, out++)
		*out = x->final_output[(x->dsp_tick*n)+i];

	x->dsp_tick++;

    return (w+5);
}


static void convolve_tilde_dsp(convolve_tilde *x, t_signal **sp)
{
	dsp_add(
		convolve_tilde_perform,
		4,
		x,
		sp[0]->s_vec,
		sp[1]->s_vec,
		sp[0]->s_n
	);
};

static void convolve_tilde_free(convolve_tilde *x)
{	
    t_freebytes(x->ir_spec, x->arraysize*sizeof(t_float));
    t_freebytes(x->hann, x->window*sizeof(t_float));
    t_freebytes(x->signal_buf, x->window*sizeof(t_sample));
    t_freebytes(x->signal_buf_windowed, x->arraysize*sizeof(t_sample));
    t_freebytes(x->nonoverlapped_output, x->buf_hist_size*x->arraysize*sizeof(t_sample));
    t_freebytes(x->final_output, x->overlap_chunk*sizeof(t_sample));
};

void convolve_tilde_setup(void)
{
    convolve_tilde_class = 
    class_new(
    	gensym("convolve~"),
    	(t_newmethod)convolve_tilde_new,
    	(t_method)convolve_tilde_free,
        sizeof(convolve_tilde),
        CLASS_DEFAULT,
        A_DEFFLOAT,
        A_DEFFLOAT,
		0
    );

    CLASS_MAINSIGNALIN(convolve_tilde_class, convolve_tilde, x_f);
	
	class_addmethod(
		convolve_tilde_class,
		(t_method)convolve_tilde_analyze,
		gensym("analyze"),
		A_SYMBOL,
		0
	);
	
    class_addmethod(
    	convolve_tilde_class,
    	(t_method)convolve_tilde_dsp,
    	gensym("dsp"),
    	0
    );
}
