// generated by Fast Light User Interface Designer (fluid) version 1.0308

#ifndef rxmon_h
#define rxmon_h
#include <FL/Fl.H>
#include "configuration.h"
#include "combo.h"
#include "flinput2.h"
#include "flslider2.h"
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Check_Button.H>
extern Fl_Check_Button *btn_mon_xcvr_audio;
extern Fl_Value_Slider2 *sldrRxFilt_vol;
extern Fl_Check_Button *btn_rxgain_x10;
extern Fl_Check_Button *btn_mon_dsp_audio;
extern Fl_Value_Slider2 *sldrRxFilt_bw;
extern Fl_Value_Slider2 *sldrRxFilt_mid;
extern Fl_Value_Slider2 *sldrRxFilt_low;
extern Fl_Value_Slider2 *sldrRxFilt_high;
extern Fl_Check_Button *btn_RxFilt_at_track;
extern Fl_Check_Button *btn_mon_wf_display;
extern Fl_Check_Button *btn_mon_xmt_audio;
extern Fl_Value_Slider2 *sldr_tx_vol;
Fl_Double_Window* make_rxaudio_dialog();
#endif