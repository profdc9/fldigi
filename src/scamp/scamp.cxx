// ----------------------------------------------------------------------------
// scamp.cxx  --  SCAMP modem
//
// Copyright (C) 2012
//		Dave Freese, W1HKJ
//		Stefan Fendt, DL1SMF
//      Daniel Marks, KW4TI
//
// This file is part of fldigi.
//
// Fldigi is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Fldigi is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with fldigi.  If not, see <http://www.gnu.org/licenses/>.
// ----------------------------------------------------------------------------

#include <config.h>
#include <iostream>
#include <fstream>

#include "view_scamp.h"
#include "fl_digi.h"
#include "digiscope.h"
#include "misc.h"
#include "waterfall.h"
#include "confdialog.h"
#include "configuration.h"
#include "status.h"
#include "digiscope.h"
#include "trx.h"
#include "debug.h"
#include "synop.h"
#include "main.h"
#include "modem.h"

#include "threads.h"

#include "scamp.h"

#define FILTER_DEBUG 0

#define SHAPER_BAUD 150

//=====================================================================
// Baudot support
//=====================================================================

static char msg1[20];

void scamp::tx_init()
{
	transmit_phase = 0;
}

void scamp::rx_init()
{
}

void scamp::init()
{
	bool wfrev = wf->Reverse();
	bool wfsb = wf->USB();
	// Probably not necessary because similar to modem::set_reverse
	reverse = wfrev ^ !wfsb;

    if (progStatus.carrier != 0) {
		set_freq(progStatus.carrier);
#if !BENCHMARK_MODE
		progStatus.carrier = 0;
#endif
	} else
		set_freq(wf->Carrier());

	rx_init();
	put_MODEstatus(mode);
    snprintf(msg1, sizeof(msg1), "SCAMP");
	put_Status1(msg1);
    set_scope_mode(Digiscope::BLANK);
}

scamp::~scamp()
{
	if (scampviewer) delete scampviewer;
}

void scamp::reset_filters()
{
	circbuffer_head_tail = 0;
	sample_count = 0;
	circbuffer1_ac_re = 0;
	circbuffer1_ac_im = 0;
	circbuffer2_ac_re = 0;
	circbuffer2_ac_im = 0;
	for (int i=0;i<SCAMP_CIRCBUFFER_MAX;i++)
	{
		circbuffer1_re[i] = 0;
		circbuffer1_im[i] = 0;
		circbuffer2_re[i] = 0;
		circbuffer2_im[i] = 0;
	}
	circ_phase = 0.0;
	carrier_phase = 0.0;
}

void scamp::restart()
{
	uint8_t protocol;
	
    sample_count_check = SCAMP_SampleRate / SCAMP_SAMPLE_COUNT_CHECK_RATE;
    switch (mode)
    {
	   case MODE_SCAMPFSK:
			protocol = PROTOCOL_SCAMP_FSK;
	        scamp_fsk_mode = true;
	        circbuffer_samples = 60*(samplerate/2000);
	        shift_freq = (0.5*TWOPI)*66.6666666/((double)samplerate);
	        mode_bandwidth = 66.66666666 + 2*33.33333333;
	        shift = 66.66666666;
	        channel_bandwidth = 33.33333333;
			break;
	   case MODE_SCAMPOOK:
			protocol = PROTOCOL_SCAMP_OOK;
	        scamp_fsk_mode = false;
	        circbuffer_samples = 64*(samplerate/2000);
	        shift_freq = 0;
	        mode_bandwidth = 2*31.25;
	        channel_bandwidth = 31.25;
	        shift = 0.0;
			break;
	   case MODE_SCFSKFST:
			protocol = PROTOCOL_SCAMP_FSK_FAST;
	        scamp_fsk_mode = true;
	        circbuffer_samples = 24*(samplerate/2000);
	        shift_freq = (0.5*TWOPI)*166.6666666/((double)samplerate);
	        mode_bandwidth = 166.66666666 + 2*83.33333333;
	        shift = 83.3333333333;
	        channel_bandwidth = 83.333333333;
			break;
	   case MODE_SCFSKSLW:
			protocol = PROTOCOL_SCAMP_FSK_SLOW;
	        scamp_fsk_mode = true;
	        circbuffer_samples = 144*(samplerate/2000);
	        shift_freq = (0.5*TWOPI)*41.66666666/((double)samplerate);
	        mode_bandwidth = 41.66666666 + 2*13.88888888;
	        shift = 41.6666666666;
	        channel_bandwidth = 41.666666666;
	        break;
	   case MODE_SCOOKSLW:
			protocol = PROTOCOL_SCAMP_OOK_SLOW;
	        scamp_fsk_mode = false;
	        circbuffer_samples = 144*(samplerate/2000);
	        shift_freq = 0;
	        mode_bandwidth = 2*13.88888888;
	        channel_bandwidth = 41.666666666;
	        shift = 0.0;
	        break;
	   case MODE_SCFSKVSL:
			protocol = PROTOCOL_SCAMP_FSK_VSLW;
	        scamp_fsk_mode = true;
	        circbuffer_samples = 288*(samplerate/2000);
	        shift_freq = (0.25*TWOPI)*41.66666666/((double)samplerate);
	        mode_bandwidth = 0.5*(41.66666666 + 2*13.88888888);
	        shift = 0.5*41.6666666666;
	        channel_bandwidth = 0.5*41.666666666;
            sample_count_check = sample_count_check * 2;
	        break;
	}

    transmit_scale = 1.0/((double)circbuffer_samples);
	scamp_protocol.init(protocol);
	
	set_bandwidth(mode_bandwidth);
	
	wf->redraw_marker();

	reset_filters();

    snprintf(msg1, sizeof(msg1), "SCAMP");
	put_Status1(msg1);
	put_MODEstatus(mode);

	if (scampviewer) scampviewer->restart();
}

scamp::scamp(trx_mode tty_mode)
{
	scampviewer = 0;
	
	mode = tty_mode;

	samplerate = SCAMP_SampleRate;

	restart();

}

void scamp::searchDown()
{
	double srchfreq = frequency - shift -100;
	double minfreq = shift * 2 + 100;
	double spwrlo, spwrhi, spwrct;
	while (srchfreq > minfreq) {
		spwrlo = wf->powerDensity(srchfreq - shift/2, 2*channel_bandwidth) + 1e-10;
		spwrhi = wf->powerDensity(srchfreq + shift/2, 2*channel_bandwidth) + 1e-10;
		spwrct = wf->powerDensity(srchfreq, 2*channel_bandwidth) + 1e-10;
		if (((scamp_fsk_mode) && ((spwrlo / spwrct > 10.0) && (spwrhi / spwrct > 10.0)))
		    || ((!scamp_fsk_mode) && ((spwrct / spwrlo > 10.0) && (spwrct / spwrhi > 10.0))))
		{
			set_freq(frequency);
			sigsearch = SIGSEARCH;
			break;
		}
		srchfreq -= 5.0;
	}
}

void scamp::searchUp()
{
	double srchfreq = frequency + shift +100;
	double maxfreq = IMAGE_WIDTH - shift * 2 - 100;
	double spwrlo, spwrhi, spwrct;
	while (srchfreq < maxfreq) {
		spwrlo = wf->powerDensity(srchfreq - shift/2, 2*channel_bandwidth) + 1e-10;
		spwrhi = wf->powerDensity(srchfreq + shift/2, 2*channel_bandwidth) + 1e-10;
		spwrct = wf->powerDensity(srchfreq, 2*channel_bandwidth) + 1e-10;
		if (((scamp_fsk_mode) && ((spwrlo / spwrct > 10.0) && (spwrhi / spwrct > 10.0)))
		    || ((!scamp_fsk_mode) && ((spwrct / spwrlo > 10.0) && (spwrct / spwrhi > 10.0))))
		{
			frequency = srchfreq;
			set_freq(frequency);
			sigsearch = SIGSEARCH;
			break;
		}
		srchfreq += 5.0;
	}
}

#if FILTER_DEBUG == 1
int snum = 0;
int mnum = 0;
#define ook(sp) \
{ \
	value = sin(2.0*M_PI*( \
		(((sp / symbollen) % 2 == 0) ? (frequency + shift/2.0) : (frequency - shift/2.0))\
		/samplerate)*sp); \
}

std::fstream ook_signal("ook_signal.csv", std::ios::out );
#endif

char lsnrmsg[80];
void scamp::Metric()
{
	double delta = channel_bandwidth/8.0;
	double np = wf->powerDensity(frequency, delta) * 3000 / delta + 1e-8;
	double sp =
		wf->powerDensity(frequency - shift/2, delta) +
		wf->powerDensity(frequency + shift/2, delta) + 1e-8;
	double snr = 0;

	if (np < 1e-6) np = sp * 100;

	sigpwr = decayavg( sigpwr, sp, sp > sigpwr ? 2 : 8);
	noisepwr = decayavg( noisepwr, np, 16 );

	snr = 10*log10(sigpwr / noisepwr);

	snprintf(lsnrmsg, sizeof(lsnrmsg), "s/n %-3.0f dB", snr);
	put_Status2(lsnrmsg);
	metric = CLAMP((3000 / delta) * (sigpwr/noisepwr), 0.0, 100.0);
	display_metric(metric);
}

int scamp::rx_process(const double *buf, int len)
{
	const double DECAY_FUDGE_FACTOR = 0.999999999999; 
			/* Fudge factor to prevent roundoff error accumulation */
	const double *buffer = buf;
	int length = len;
	double phaseinc = frequency*TWOPI/((double)samplerate);

    for (int samp=0;samp<length;samp++)
    {
		double val = buffer[samp];
		double mag1, mag2;
		if (scamp_fsk_mode)
		{
		   double buf_re = val * cos(carrier_phase + circ_phase);
		   double buf_im = val * sin(carrier_phase + circ_phase);
		   circbuffer1_ac_re += buf_re - circbuffer1_re[circbuffer_head_tail];
		   circbuffer1_ac_re *= DECAY_FUDGE_FACTOR;
		   circbuffer1_ac_im += buf_im - circbuffer1_im[circbuffer_head_tail];
		   circbuffer1_ac_im *= DECAY_FUDGE_FACTOR;
		   circbuffer1_re[circbuffer_head_tail] = buf_re;
		   circbuffer1_im[circbuffer_head_tail] = buf_im;
		   mag1 = sqrt(circbuffer1_ac_re*circbuffer1_ac_re+circbuffer1_ac_im*circbuffer1_ac_im);
		   
		   buf_re = val * cos(carrier_phase - circ_phase);
		   buf_im = val * sin(carrier_phase - circ_phase);
		   circbuffer2_ac_re += buf_re - circbuffer2_re[circbuffer_head_tail];
		   circbuffer2_ac_re *= DECAY_FUDGE_FACTOR;
		   circbuffer2_ac_im += buf_im - circbuffer2_im[circbuffer_head_tail];
		   circbuffer2_ac_im *= DECAY_FUDGE_FACTOR;
		   circbuffer2_re[circbuffer_head_tail] = buf_re;
		   circbuffer2_im[circbuffer_head_tail] = buf_im;
		   mag2 = sqrt(circbuffer2_ac_re*circbuffer2_ac_re+circbuffer2_ac_im*circbuffer2_ac_im);
		   
		   circ_phase += shift_freq;
		   if (circ_phase >= TWOPI)
				circ_phase -= TWOPI;
		} else
		{
		   double buf_re = val * cos(carrier_phase);
		   double buf_im = val * sin(carrier_phase);
		   circbuffer1_ac_re += buf_re - circbuffer1_re[circbuffer_head_tail];
		   circbuffer1_ac_re *= DECAY_FUDGE_FACTOR;
		   circbuffer1_ac_im += buf_im - circbuffer1_im[circbuffer_head_tail];
		   circbuffer1_ac_im *= DECAY_FUDGE_FACTOR;
		   circbuffer1_re[circbuffer_head_tail] = buf_re;
		   circbuffer1_im[circbuffer_head_tail] = buf_im;
		   mag1 = sqrt(circbuffer1_ac_re*circbuffer1_ac_re+circbuffer1_ac_im*circbuffer1_ac_im);
		   mag2 = 0.0;
		}
		if ((++circbuffer_head_tail) >= circbuffer_samples)
			circbuffer_head_tail = 0;
	    if ((++sample_count) >= sample_count_check)
		{
			int recv_chars[2];
			sample_count = 0;
			/* call SCAMP code */
			scamp_protocol.decode_process(mag1*transmit_scale,mag2*transmit_scale,recv_chars);
			if (recv_chars[0] != -1) 
			{
				put_rx_char(recv_chars[0]);
				if (recv_chars[0] == '\r') put_rx_char('\n');
			}
			if (recv_chars[1] != -1) 
			{
				put_rx_char(recv_chars[1]);
				if (recv_chars[1] == '\r') put_rx_char('\n');
			}
		}
		carrier_phase += phaseinc;
		if (carrier_phase >= TWOPI)
			carrier_phase -= TWOPI;
	}

	if ( !progdefaults.report_when_visible ||
		 dlgViewer->visible() || progStatus.show_channels )
		if (!bHistory && scampviewer) scampviewer->rx_process(buf, len);

	{
		reverse = wf->Reverse() ^ !wf->USB();
	}

	Metric();
	return 0;
}

void scamp::send_frame(uint32_t frame)
{
  double phaseinc = frequency*TWOPI/((double)samplerate);
  for (int bitno=0;bitno<30;bitno++)
  {
    int bitv = (frame & 0x20000000) != 0;
    frame <<= 1;
    if (scamp_fsk_mode)
    {
		double freq = phaseinc + (bitv ? shift_freq : -shift_freq);
		for (int i=0;i<circbuffer_samples;i++)
		{
			transbuffer[i] = sin(transmit_phase);
			transmit_phase += freq;
			if (transmit_phase > TWOPI)
				transmit_phase -= TWOPI;
		}
	} else
	{
		double amp = bitv ? 1.0 : 0.0;
		for (int i=0;i<circbuffer_samples;i++)
		{
			transbuffer[i] = amp*sin(transmit_phase);
			transmit_phase += phaseinc;
			if (transmit_phase > TWOPI)
				transmit_phase -= TWOPI;
		}
	}
    ModulateXmtr(transbuffer, circbuffer_samples);
  }
}

double scamp::scamp_now()
{
	static struct timeval t1;
	gettimeofday(&t1, NULL);
	return t1.tv_sec + t1.tv_usec / 1e6;
}

// sub millisecond accurate sleep function
// sleep_time in seconds
int scamp::scamp_sleep (double sleep_time)
{
	struct timespec tv;
	double start_at = scamp_now();
	double end_at = start_at + sleep_time;

	double delay = sleep_time - 0.005;

	tv.tv_sec = (time_t) delay;
	tv.tv_nsec = (long) ((delay - tv.tv_sec) * 1000000000L);
	int rval = 0;

#ifdef __WIN32__
	timeBeginPeriod(1);
#endif
//	while (1) {
		rval = nanosleep (&tv, &tv);
		if (errno == EINTR) {
//			continue
std::cout << "EINTR error in scamp_sleep" << std::endl;
		}
//		break;
//	}
	rval = 0;
	while (scamp_now() < end_at) rval++;
//std::cout << "scamp_sleep( " << sleep_time << ") : " << rval << std::endl;
#ifdef __WIN32__
	timeEndPeriod(1);
#endif

	return 0;
}

#define TX_PROCESS_FRAMES 100

int scamp::tx_process()
{
	int ret = 0;
	uint32_t frame_array[TX_PROCESS_FRAMES];
	uint8_t frames = 0;
	modem::tx_process();

	scamp_protocol.set_resync_repeat_frames(progdefaults.ScampResync,progdefaults.ScampRepeat);

	int c = get_tx_char();
	
	if (c == GET_TX_CHAR_ETX || stopflag) {
		stopflag = false;
		frames = scamp_protocol.send_char(-1, TX_PROCESS_FRAMES, frame_array);
		ret = -1;
	} else if (c == GET_TX_CHAR_NODATA) {
		frames = scamp_protocol.send_char(-2, TX_PROCESS_FRAMES, frame_array);
	} else {
		frames = scamp_protocol.send_char(c, TX_PROCESS_FRAMES, frame_array);
		put_echo_char(c);
	}
	for (int fr=0;fr<frames;fr++)
		send_frame(frame_array[fr]);
	return ret;
}
