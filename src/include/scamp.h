// ----------------------------------------------------------------------------
// scamp.h  --  SCAMP modem
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

#ifndef _SCAMP_H
#define _SCAMP_H

#include <iostream>

#include "complex.h"
#include "modem.h"
#include "globals.h"
#include "filters.h"
#include "fftfilt.h"
#include "digiscope.h"
#include "view_scamp.h"
#include "serial.h"

#include "scampsk.h"

#define	SCAMP_SampleRate	8000
//#define SCAMP_SampleRate 11025
//#define SCAMP_SampleRate 12000

#define MAXPIPE			1024
#define MAXBITS			(2 * SCAMP_SampleRate / 23 + 1)

#define	LETTERS	0x100
#define	FIGURES	0x200

#define dispwidth 100

// simple oscillator-class
class SCAMPOscillator
{
public:
	SCAMPOscillator( double samplerate );
	~SCAMPOscillator() {}
	double Update( double frequency );

private:
	double m_phase;
	double m_samplerate;
};

class SCAMPSymbolShaper
{
public:
	SCAMPSymbolShaper(double baud = 45.45, double sr = 8000.0);
	~SCAMPSymbolShaper();
	void reset();
	void Preset(double baud, double sr);
	void print_sinc_table();
	double Update( bool state );

private:
	int		 m_table_size;
	double*	 m_sinc_table;

	bool		m_State;
	double		m_Accumulator;
	long		m_Counter0;
	long		m_Counter1;
	long		m_Counter2;
	long		m_Counter3;
	long		m_Counter4;
	long		m_Counter5;
	double		m_Factor0;
	double		m_Factor1;
	double		m_Factor2;
	double		m_Factor3;
	double		m_Factor4;
	double		m_Factor5;
	double		m_SincTable[1024];

	double		baudrate;
	double		samplerate;
};

//enum TTY_MODE { LETTERS, FIGURES };

class scamp : public modem {
public:
enum SCAMP_RX_STATE {
	SCAMP_RX_STATE_IDLE = 0,
	SCAMP_RX_STATE_START,
	SCAMP_RX_STATE_DATA,
	SCAMP_RX_STATE_PARITY,
	SCAMP_RX_STATE_STOP,
	SCAMP_RX_STATE_STOP2
};
enum SCAMP_PARITY {
	SCAMP_PARITY_NONE = 0,
	SCAMP_PARITY_EVEN,
	SCAMP_PARITY_ODD,
	SCAMP_PARITY_ZERO,
	SCAMP_PARITY_ONE
};

	static const double	SHIFT[];
	static const double	BAUD[];
	static const int		BITS[];
	static const int		FILTLEN[];
	static const int		numshifts;
	static const int		numbauds;

private:

	SCAMPOscillator		*m_Osc1;
	SCAMPOscillator		*m_Osc2;
	SCAMPSymbolShaper	*m_SymShaper1;
	SCAMPSymbolShaper	*m_SymShaper2;

	double shift;
	int symbollen;
	int nbits;
	int stoplen;
	int msb;

	double		phaseacc;
	double		scamp_squelch;
	double		scamp_shift;
	double		scamp_BW;
	double		scamp_baud;
	int 		scamp_bits;
	SCAMP_PARITY	scamp_parity;
	int			scamp_stop;
	bool		scamp_msbfirst;

	double		mark_noise;
	double		space_noise;
	Cmovavg		*bits;
	bool		nubit;
	bool		bit;

	bool		bit_buf[MAXBITS];

	double mark_phase;
	double space_phase;
	fftfilt *mark_filt;
	fftfilt *space_filt;
	int filter_length;

	double *pipe;
	double *dsppipe;
	int pipeptr;

	cmplx mark_history[MAXPIPE];
	cmplx space_history[MAXPIPE];

	SCAMP_RX_STATE rxstate;

	int counter;
	int bitcntr;
	int rxdata;
	double cfreq; // center frequency between MARK/SPACE tones
	double shift_offset; // 1/2 scamp_shift

	double prevsymbol;
	cmplx prevsmpl;

	double xy_phase;
	double rotate;

	cmplx QI[MAXPIPE];
	int inp_ptr;

	cmplx xy;

	bool   clear_zdata;
	double sigpwr;
	double noisepwr;
	double avgsig;

	double mark_mag;
	double space_mag;
	double mark_env;
	double space_env;
	double	noise_floor;

	double SCAMPSKbuf[OUTBUFSIZE];		// signal array for qrq drive
	double SCAMPSKphaseacc;
	double SCAMPSKnco();

	unsigned char lastchar;

	int rxmode;
	int shift_state;
	bool preamble;

	void Clear_syncscope();
	void Update_syncscope();

	double IF_freq;
	inline cmplx mixer(double &phase, double f, cmplx in);

	unsigned char Bit_reverse(unsigned char in, int n);
	int decode_char();
	bool rx(bool bit);

	view_scamp *scampviewer;

// transmit
	double nco(double freq);
	void send_symbol(int symbol, int len);
	void send_stop();
	void send_char(int c);
	void send_idle();
	int scampxprocess();
	int baudot_enc(unsigned char data);
	char baudot_dec(unsigned char data);
	void Metric();

	bool is_mark_space(int &);
	bool is_mark();

//----------------------------------------------------------------------
// SCAMPSK via signal control line [serial DTR or RTS pin]
//----------------------------------------------------------------------
	SCAMPSK  *scampsk_tty;
	void send_SCAMPSK(int);

// SCAMPSK via flrig DTR/RTS
	void flrig_scampsk_send(char c);
	double scamp_now();
	int scamp_sleep (double);

public:
	scamp(trx_mode mode);
	~scamp();
	void init();
	void rx_init();
	void tx_init();
	void restart();
	void reset_filters();
	int rx_process(const double *buf, int len);
	int tx_process();
	void flush_stream();

	void clear_viewer() { scampviewer->clear(); }
	void clear_ch(int n) { scampviewer->clearch(n); }
	int  viewer_get_freq(int n) { return scampviewer->get_freq(n); }

	void searchDown();
	void searchUp();

	void resetSCAMPSK();
};

int scampparity(unsigned int, int);

#endif
