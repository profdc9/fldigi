// ----------------------------------------------------------------------------
// scamp.cxx  --  SCAMP modem
//
// Copyright (C) 2012
//		Dave Freese, W1HKJ
//		Stefan Fendt, DL1SMF
//
// This file is part of fldigi.
//
// This code bears some resemblance to code contained in gmscampsk from which
// it originated.  Much has been changed, but credit should still be
// given to Tomi Manninen (oh2bns@sral.fi), who so graciously distributed
// his gmscampsk modem under the GPL.
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

static char letters[32] = {
	'\0',	'E',	'\n',	'A',	' ',	'S',	'I',	'U',
	'\r',	'D',	'R',	'J',	'N',	'F',	'C',	'K',
	'T',	'Z',	'L',	'W',	'H',	'Y',	'P',	'Q',
	'O',	'B',	'G',	' ',	'M',	'X',	'V',	' '
};

/*
 * U.S. version of the figures case.
 */
static char figures[32] = {
	'\0',	'3',	'\n',	'-',	' ',	'\a',	'8',	'7',
	'\r',	'$',	'4',	'\'',	',',	'!',	':',	'(',
	'5',	'"',	')',	'2',	'#',	'6',	'0',	'1',
	'9',	'?',	'&',	' ',	'.',	'/',	';',	' '
};

int dspcnsc = 0;

static char msg1[20];

const double	scamp::SHIFT[] = {23, 85, 160, 170, 182, 200, 240, 350, 425, 850};
// FILTLEN must be same size as BAUD
const double	scamp::BAUD[]  = {45, 45.45, 50, 56, 75, 100, 110, 150, 200, 300};
const int		scamp::FILTLEN[] = { 512, 512, 512, 512, 512, 512, 512, 256, 128, 64};
const int		scamp::BITS[]  = {5, 7, 8};
const int		scamp::numshifts = (int)(sizeof(SHIFT) / sizeof(*SHIFT));
const int		scamp::numbauds = (int)(sizeof(BAUD) / sizeof(*BAUD));

void scamp::tx_init()
{
/*	phaseacc = 0;
	preamble = true;
	videoText();

	symbols = 0;
	acc_symbols = 0;
	ovhd_symbols = 0;
*/
}

void scamp::rx_init()
{
	
/*
	rxstate = SCAMP_RX_STATE_IDLE;
	rxmode = LETTERS;
	phaseacc = 0;
	SCAMPSKphaseacc = 0;

	for (int i = 0; i < MAXBITS; i++ ) bit_buf[i] = 0.0;

	mark_phase = 0;
	space_phase = 0;
	xy_phase = 0.0;

	mark_mag = 0;
	space_mag = 0;
	mark_env = 0;
	space_env = 0;

	inp_ptr = 0;

	lastchar = 0;
	*/
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
	sample_count_check = SCAMP_SampleRate / SCAMP_SAMPLE_COUNT_CHECK_RATE;
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
	        break;
	}

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

cmplx scamp::mixer(double &phase, double f, cmplx in)
{
	cmplx z = cmplx( cos(phase), sin(phase)) * in;

	phase -= TWOPI * f / samplerate;
	if (phase < -TWOPI) phase += TWOPI;

	return z;
}

void scamp::searchDown()
{
	double srchfreq = frequency - shift -100;
	double minfreq = shift * 2 + 100;
	double spwrlo, spwrhi, npwr;
	while (srchfreq > minfreq) {
		spwrlo = wf->powerDensity(srchfreq - shift/2, 2*channel_bandwidth);
		spwrhi = wf->powerDensity(srchfreq + shift/2, 2*channel_bandwidth);
		npwr = wf->powerDensity(srchfreq + shift, 2*scamp_baud) + 1e-10;
		if ((spwrlo / npwr > 10.0) && (spwrhi / npwr > 10.0)) {
			frequency = srchfreq;
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
	double spwrhi, spwrlo, npwr;
	while (srchfreq < maxfreq) {
		spwrlo = wf->powerDensity(srchfreq - shift/2, 2*channel_bandwidth);
		spwrhi = wf->powerDensity(srchfreq + shift/2, 2*channel_bandwidth);
		npwr = wf->powerDensity(srchfreq - shift, 2*scamp_baud) + 1e-10;
		if ((spwrlo / npwr > 10.0) && (spwrhi / npwr > 10.0)) {
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
		   circbuffer1_ac_im += buf_im - circbuffer1_im[circbuffer_head_tail];
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
			scamp_protocol.decode_process(mag1,mag2,recv_chars);
			if (recv_chars[0] != -1) put_rx_char(recv_chars[0]);
			if (recv_chars[1] != -1) put_rx_char(recv_chars[1]);
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

//=====================================================================
// SCAMP transmit
//=====================================================================
//double freq1;
double maxamsc = 0;

double scamp::nco(double freq)
{
	phaseacc += TWOPI * freq / samplerate;

	if (phaseacc > TWOPI) phaseacc -= TWOPI;

	return cos(phaseacc);
}

double scamp::SCAMPSKnco()
{
	SCAMPSKphaseacc += TWOPI * 1000 / samplerate;

	if (SCAMPSKphaseacc > TWOPI) SCAMPSKphaseacc -= TWOPI;

	return sin(SCAMPSKphaseacc);

}

extern Cserial CW_KEYLINE_serial;
extern bool CW_KEYLINE_isopen;

void scamp::send_symbol(int symbol, int len)
{
	if (reverse) symbol = !symbol;

	acc_symbols += len;

		double freq;

		if (symbol)
			freq = get_txfreq_woffset() + shift / 2.0;
		else
			freq = get_txfreq_woffset() - shift / 2.0;

		for (int i = 0; i < len; i++) {
			outbuf[i] = nco(freq);
			if (symbol)
				SCAMPSKbuf[i] = SCAMPSKnco();
			else
				SCAMPSKbuf[i] = 0.0;
		}

    ModulateXmtr(outbuf, symbollen);
}

void scamp::send_stop()
{
		double freq;
		bool invert = reverse;

		if (invert)
			freq = get_txfreq_woffset() - shift / 2.0;
		else
			freq = get_txfreq_woffset() + shift / 2.0;

		for (int i = 0; i < stoplen; i++) {
			outbuf[i] = nco(freq);
			if (invert)
				SCAMPSKbuf[i] = 0.0;
			else
				SCAMPSKbuf[i] = SCAMPSKnco();
		}

		ModulateXmtr(outbuf, stoplen);

}

void scamp::flush_stream()
{
	double const freq1 = get_txfreq_woffset() + shift / 2.0;
	double const freq2 = get_txfreq_woffset() - shift / 2.0;
	double mark = 0, space = 0, signal = 0;

	for( int i = 0; i < symbollen * 6; ++i ) {
		mark  = m_SymShaper1->Update(0)*m_Osc1->Update( freq1 );
		space = m_SymShaper2->Update(0)*m_Osc2->Update( freq2 );
		signal = mark + space;

		if (maxamsc < fabs(signal)) maxamsc = fabs(signal);
		
		outbuf[i] = maxamsc ? (signal / maxamsc) : 0.0;

		SCAMPSKbuf[i] = 0.0;
	}

	sig_stop = true;
	
    ModulateXmtr(outbuf, symbollen * 6);

}

void scamp::send_char(int c)
{
	int i;
	if (nbits == 5) {
		if (c == LETTERS)
			c = 0x1F;
		if (c == FIGURES)
			c = 0x1B;
	}

// start bit
	send_symbol(0, symbollen);
// data bits
	for (i = 0; i < nbits; i++) {
		send_symbol((c >> i) & 1, symbollen);
	}
// parity bit
//	if (scamp_parity != SCAMP_PARITY_NONE)
//		send_symbol(scampparity(c, nbits), symbollen);
// stop bit(s)
	send_stop();

	if (nbits == 5) {
		if (c == 0x1F || c == 0x1B)
			return;
		if (shift_state == LETTERS)
			c = letters[c];
		else
			c = figures[c];
		if (c)
			put_echo_char(progdefaults.rx_lowercase ? tolower(c) : c);
	}
	else
		put_echo_char(c);

}

void scamp::send_idle()
{
	if (nbits == 5) {
		send_char(LETTERS);
		shift_state = LETTERS;
	} else
		send_char(0);
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

static int line_char_count = 0;
// 1 start, 5 data, 1.5/2.0 stopbits
#define wait_one_byte(baud, stopbits) \
scamp_sleep( ((6 + (stopbits))*1.0 / (baud)));

int scamp::tx_process()
{
	modem::tx_process();

	int c = get_tx_char();

	if (preamble) {
		sig_start = true;
		for (int i = 0; i < progdefaults.TTY_LTRS; i++)
			send_char(LETTERS);
		preamble = false;
	}

// TX buffer empty
	if (c == GET_TX_CHAR_ETX || stopflag) {
		stopflag = false;
		line_char_count = 0;
		if (nbits != 5) {
			if (progdefaults.rtty_crcrlf) send_char('\r');
			send_char('\r');
			send_char('\n');
		} else {
			if (progdefaults.rtty_crcrlf) send_char(0x08);
			send_char(0x08);
			send_char(0x02);
		}
		flush_stream();
		return -1;
	}
// send idle character if c == -1
	if (c == GET_TX_CHAR_NODATA) {
		send_idle();
		return 0;
	}

// if NOT Baudot
	if (nbits != 5) {
		acc_symbols = 0;
		send_char(c);
		xmt_samples = char_samples = acc_symbols;
		return 0;
	}

	if (isalpha(c) || isdigit(c) || isblank(c) || ispunct(c)) {
		++line_char_count;
	}

	if (progdefaults.rtty_autocrlf && (c != '\n' && c != '\r') &&
		(line_char_count == progdefaults.rtty_autocount ||
		 (line_char_count > progdefaults.rtty_autocount - 5 && c == ' '))) {
		line_char_count = 0;
		if (progdefaults.rtty_crcrlf)
			send_char(0x08); // CR-CR-LF triplet
		send_char(0x08);
		send_char(0x02);
		if (c == ' ')
			return 0;
	}
	if (c == '\r') {
		line_char_count = 0;
		send_char(0x08);
		return 0;
	}
	if (c == '\n') {
		line_char_count = 0;
		if (progdefaults.rtty_crcrlf)
			send_char(0x08); // CR-CR-LF triplet
		send_char(0x02);
		return 0;
	}


/* unshift-on-space */
	if (c == ' ') {
		if (progdefaults.UOStx) {
			send_char(LETTERS);
			send_char(0x04); // coded value for a space
			shift_state = LETTERS;
		} else
			send_char(0x04);
		return 0;
	}

	if ((c = baudot_enc(c)) < 0)
		return 0;

// switch case if necessary

	if ((c & 0x300) != shift_state) {
		if (shift_state == FIGURES) {
			send_char(LETTERS);
			shift_state = LETTERS;
		} else {
			send_char(FIGURES);
			shift_state = FIGURES;
		}
	}
///
	acc_symbols = 0;
	send_char(c & 0x1F);
	xmt_samples = char_samples = acc_symbols;

	return 0;
}

int scamp::baudot_enc(unsigned char data)
{
	int i, c, mode;

	mode = 0;
	c = -1;

	if (islower(data))
		data = toupper(data);

	for (i = 0; i < 32; i++) {
		if (data == letters[i]) {
			mode |= LETTERS;
			c = i;
		}
		if (data == figures[i]) {
			mode |= FIGURES;
			c = i;
		}
		if (c != -1)
			return (mode | c);
	}

	return -1;
}

char scamp::baudot_dec(unsigned char data)
{
	int out = 0;

	switch (data) {
	case 0x1F:		/* letters */
		rxmode = LETTERS;
		break;
	case 0x1B:		/* figures */
		rxmode = FIGURES;
		break;
	case 0x04:		/* unshift-on-space */
		if (progdefaults.UOSrx)
			rxmode = LETTERS;
		return ' ';
		break;
	default:
		if (rxmode == LETTERS)
			out = letters[data];
		else
			out = figures[data];
		break;
	}

	return out;
}

//======================================================================
// methods for class Oscillator and class SymbolShaper
//======================================================================

SCAMPOscillator::SCAMPOscillator( double samplerate )
{
	m_phase = 0;
	m_samplerate = samplerate;
//	std::cerr << "samplerate for Oscillator:"<<m_samplerate<<"\n";
}

double SCAMPOscillator::Update( double frequency )
{
	m_phase += frequency/m_samplerate * TWOPI;
	if ( m_phase > TWOPI ) m_phase -= TWOPI;

	return ( sin( m_phase ) );
}

SCAMPSymbolShaper::SCAMPSymbolShaper(double baud, double sr)
{
	m_sinc_table = 0;
	Preset( baud, sr );
}

SCAMPSymbolShaper::~SCAMPSymbolShaper()
{
	delete [] m_sinc_table;
}

void SCAMPSymbolShaper::reset()
{
	m_State = false;
	m_Accumulator = 0.0;
	m_Counter0 = 1024;
	m_Counter1 = 1024;
	m_Counter2 = 1024;
	m_Counter3 = 1024;
	m_Counter4 = 1024;
	m_Counter5 = 1024;
	m_Factor0 = 0.0;
	m_Factor1 = 0.0;
	m_Factor2 = 0.0;
	m_Factor3 = 0.0;
	m_Factor4 = 0.0;
	m_Factor5 = 0.0;
}

void SCAMPSymbolShaper::Preset(double baud, double sr)
{
    double baud_rate = baud;
    double sample_rate = sr;

    LOG_INFO("Shaper::reset( %f, %f )",  baud_rate, sample_rate);

// calculate new table-size for six integrators ----------------------

    m_table_size = sample_rate / baud_rate * 5.49;
    LOG_INFO("Shaper::m_table_size = %d", m_table_size);

// kill old sinc-table and get memory for the new one -----------------

	if (m_sinc_table)
		delete [] m_sinc_table;
    m_sinc_table = new double[m_table_size];

// set up the new sinc-table based on the new parameters --------------

    long double sum = 0.0;

    for( int x=0; x<m_table_size; ++x ) {
        int const offset = m_table_size/2;
        double wfactor = 1.0 / 1.568; // optimal
// symbol-length in samples if wmultiple = 1.0
        double const T = wfactor * sample_rate / (baud_rate*2.0);
// symbol-time relative to zero
        double const t = (x-offset);

        m_sinc_table[x] = rcos( t, T, 1.0 );

// calculate integral
        sum += m_sinc_table[x];
    }

// scale the values in the table so that the integral over it is as
// exactly 1.0000000 as we can do this...

    for( int x=0; x<m_table_size; ++x ) {
        m_sinc_table[x] *= 1.0 / sum;
    }

// reset internal states
    reset();
    maxamsc = 0;
}

double SCAMPSymbolShaper::Update( bool state )
{
	if( m_State != state ) {
		m_State = state;
		if( m_Counter0 >= m_table_size ) {
			m_Counter0 = 0;
			m_Factor0 = (state)? +1.0 : -1.0;
		} else if( m_Counter1 >= m_table_size ) {
			m_Counter1 = 0;
			m_Factor1 = (state)? +1.0 : -1.0;
		} else if( m_Counter2 >= m_table_size ) {
			m_Counter2 = 0;
			m_Factor2 = (state)? +1.0 : -1.0;
		} else if( m_Counter3 >= m_table_size ) {
			m_Counter3 = 0;
			m_Factor3 = (state)? +1.0 : -1.0;
		} else if( m_Counter4 >= m_table_size ) {
			m_Counter4 = 0;
			m_Factor4 = (state)? +1.0 : -1.0;
		} else  if( m_Counter5 >= m_table_size ) {
			m_Counter5 = 0;
			m_Factor5 = (state)? +1.0 : -1.0;
		}
	}

	if( m_Counter0 < m_table_size )
		m_Accumulator += m_Factor0 * m_sinc_table[m_Counter0++];

	if( m_Counter1 < m_table_size )
		m_Accumulator += m_Factor1 * m_sinc_table[m_Counter1++];

	if( m_Counter2 < m_table_size )
		m_Accumulator += m_Factor2 * m_sinc_table[m_Counter2++];

	if( m_Counter3 < m_table_size )
		m_Accumulator += m_Factor3 * m_sinc_table[m_Counter3++];

	if( m_Counter4 < m_table_size )
		m_Accumulator += m_Factor4 * m_sinc_table[m_Counter4++];

	if( m_Counter5 < m_table_size )
		m_Accumulator += m_Factor5 * m_sinc_table[m_Counter5++];

	return ( m_Accumulator / sqrt(2) );
}

void SCAMPSymbolShaper::print_sinc_table()
{
	for (int i = 0; i < 1024; i++) printf("%f\n", m_SincTable[i]);
}

