// ----------------------------------------------------------------------------
// scamp.h  --  SCAMP modem
//
// Copyright (C) 2006
//		Dave Freese, W1HKJ
//
// This file is part of fldigi.  Adapted from code contained in gmscampsk source code
// distribution.
//  gmscampsk Copyright (C) 2001, 2002, 2003
//  Tomi Manninen (oh2bns@sral.fi)
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

#ifndef VIEW_SCAMP_H
#define VIEW_SCAMP_H

//#include "scamp.h"
#include "complex.h"
#include "modem.h"
#include "globals.h"
#include "filters.h"
#include "fftfilt.h"
#include "digiscope.h"

#define	VIEW_SCAMP_SampleRate	8000

#define VIEW_MAXPIPE 1024
#define	VIEW_SCAMP_MAXBITS	(2 * VIEW_SCAMP_SampleRate / 23 + 1)

#define MAX_CHANNELS 30

class view_scamp : public modem {
public:
enum CHANNEL_STATE {IDLE, SRCHG, RCVNG, WAITING};
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
	int		filter_length;

private:

struct SCAMP_CHANNEL {

	int				state;

	double			phaseacc;

	fftfilt *mark_filt;
	fftfilt *space_filt;

	Cmovavg		*bits;
	bool		nubit;
	bool		bit;

	bool		bit_buf[VIEW_SCAMP_MAXBITS];

	double mark_phase;
	double space_phase;

	double		metric;

	int			rxmode;
	SCAMP_RX_STATE	rxstate;

	double		frequency;
	double		freqerr;
	double		phase;
	double		posfreq;
	double		negfreq;
	double		freqerrhi;
	double		freqerrlo;
	double		poserr;
	double		negerr;
	int			poscnt;
	int			negcnt;
	int			timeout;

	double		mark_mag;
	double		space_mag;
	double		mark_env;
	double		space_env;
	double		noise_floor;
	double		mark_noise;
	double		space_noise;

	double		sigpwr;
	double		noisepwr;
	double		avgsig;

	double		prevsymbol;
	cmplx		prevsmpl;
	int			counter;
	int			bitcntr;
	int			rxdata;
	int			inp_ptr;

	cmplx		mark_history[VIEW_MAXPIPE];
	cmplx		space_history[VIEW_MAXPIPE];

	int			sigsearch;
};

	double shift;
	int symbollen;
	int nbits;
	int stoplen;
	int msb;
	bool useSCAMPSK;

	SCAMP_CHANNEL	channel[MAX_CHANNELS];

	double		scamp_squelch;
	double		scamp_shift;
	double      scamp_BW;
	double		scamp_baud;
	int 		scamp_bits;
	SCAMP_PARITY	scamp_parity;
	int			scamp_stop;
	bool		scamp_msbfirst;

	int bflen;
	double bp_filt_lo;
	double bp_filt_hi;

	int txmode;
	int preamble;

	void clear_syncscope();
	void update_syncscope();
	cmplx mixer(double &phase, double f, cmplx in);

	unsigned char bitreverse(unsigned char in, int n);
	int decode_char(int ch);
	int scampparity(unsigned int);
	bool rx(int ch, bool bit);

	int scampxprocess();
	char baudot_dec(int ch, unsigned char data);
	void Metric(int ch);
public:
	view_scamp(trx_mode mode);
	~view_scamp();
	void init();
	void rx_init();
	void tx_init() {}
	void restart();
	void reset_filters(int ch);
	int rx_process(const double *buf, int len);
	int tx_process();

	void find_signals();
	void clearch(int ch);
	void clear();
	int get_freq(int n) { return (int)channel[n].frequency;}

	bool is_mark_space(int ch, int &);
	bool is_mark(int ch);

};

//extern view_scamp *scampviewer;

#endif
