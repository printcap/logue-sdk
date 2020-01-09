// Two-pole low pass filter using Biquad Filter from dsp::biquad
//
// Transfer function of biquad filter (two quadratic functions, 
// one in numerator and one denominator):
//
//         Y(z)     f_0 + f_1*z^{-1} + f_2*z^{-2}         (f: forward)
// H(z) = ------ = -------------------------------
//         X(z)       1 + b_1*z^{-1} + b_2*z^{-2}         (b: backward)
//
// Y(z)[1 + b_1*z^{-1} + b_2*z^{-2}] = X(z)[f_0 + f_1*z^{-1} + f_2*z^{-2}]
// 
// y_k = f_0*x_k + f_1*x_{k-1} + f_2*x_{k-2} - b_1*y_{k-1} - b_2*y_{k-2}
// 
// dsp::BiQuad implements filter in Direct Transposed Form 2:
// 
// x_k -+->[ f_0]--->[ + ]-----------+-> y_k
//      |              ^             |        
//      |              |             |         y_k  := f_0*x_k + Z1_{k-1}
//      |           [z^{-1}]         |         Z1_k := f_1*x_k + Z2_{k-1} - b_1*y_k
//      |              |             |         Z2_k := f_2*x_k - b_2*y_k
//      +->[ f_1]--->[ + ]<---[-b_1]-+
//      |              ^             |        ==
//      |              |             |
//      |           [z^{-1}]         |        acc   := f_0*x_k + Z1
//      |              |             |        Z1    := f_1*x_k + Z2 - b_1*acc
//      +->[ f_2]--->[ + ]<---[-b_2]-+        Z2    := f_2*x_k - b_2*acc
//                                            y_k   := acc
// 
//
//


#include <usermodfx.h>

#include <biquad.hpp>

static dsp::BiQuad s_bq_l, s_bq_r;
static dsp::BiQuad s_bqs_l, s_bqs_r;

static int s_param_changed = 0;
static float s_wc;
static float s_q;
static const float s_fs_recip = 1.f / 48000.f;

void MODFX_INIT(uint32_t platform, uint32_t api)
{
  s_param_changed = 0;
  s_wc = 0.49f;
  s_q = 1.4142f;

  s_bq_l.flush();
  s_bq_r.flush();
  s_bq_l.mCoeffs.setSOLP(fx_tanpif(s_wc), s_q);
  s_bq_r.mCoeffs = s_bq_l.mCoeffs;

  s_bqs_l.flush();
  s_bqs_r.flush();
  s_bqs_l.mCoeffs = s_bqs_r.mCoeffs = s_bq_l.mCoeffs;
}

void MODFX_PROCESS(const float *main_xn, float *main_yn,
                   const float *sub_xn,  float *sub_yn,
                   uint32_t frames)
{
  const float * mx = main_xn;
  float * __restrict my = main_yn;
  const float * my_e = my + 2*frames;

  const float *sx = sub_xn;
  float * __restrict sy = sub_yn;
  
  const int param_changed = s_param_changed;
  if (param_changed) {
    
    s_bq_l.mCoeffs.setSOLP(fx_tanpif(s_wc), s_q);
    s_bq_r.mCoeffs = s_bq_l.mCoeffs;
    s_bqs_l.mCoeffs = s_bq_l.mCoeffs;
    s_bqs_r.mCoeffs = s_bq_l.mCoeffs;
    s_param_changed = 0;
  }
  
  for (; my != my_e; ) {
    *(my++) = s_bq_l.process_so(*(mx++));
    *(my++) = s_bq_r.process_so(*(mx++));
    *(sy++) = s_bqs_l.process_so(*(sx++));
    *(sy++) = s_bqs_r.process_so(*(sx++));
  }
}


void MODFX_PARAM(uint8_t index, int32_t value)
{
  const float valf = q31_to_f32(value); // valf in [0, 1)
  switch (index) {
  case k_user_modfx_param_time:
    // convert time knob parameter to cut-off frequency 
    // (max cut-off < Nyquist frequency)
    // valf = 0     -> 0.001 * Fs = 0.001 * 48 kHz = 48 Hz
    // valf = 0.999 -> 0.490 * Fs = 0.490 * 48 KHz = 23.520 KHz
    // scale of valf to frequency is exponential 
    s_wc = 0.001f * fasterpow2f(valf * 8.9456f);    // 2*log2(490)/0.999
    s_param_changed = 1;
    break;
  case k_user_modfx_param_depth:
    // convert depth knob parameter to resonance q
    // exponential behavior
    // valf = 0     -> 0 dB resonance peak
    // valf = 0.999 -> +20 dB resonance peak
    s_q = fasterpow2f(valf * 4.0f) * M_1_SQRT2;   
    s_param_changed = 1;
    break;
  default:
    break;
  }
}
