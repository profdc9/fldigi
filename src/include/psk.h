// ----------------------------------------------------------------------------
// psk.h  --  psk modem
//
// Copyright (C) 2006-2008
//		Dave Freese, W1HKJ
//
// This file is part of fldigi.  Adapted from code contained in gmfsk source code
// distribution.
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

#ifndef _PSK_H
#define _PSK_H

#include "complex.h"
#include "modem.h"
#include "globals.h"
#include "viterbi.h"
#include "filters.h"
#include "pskcoeff.h"
#include "pskvaricode.h"
#include "viewpsk.h"
#include "pskeval.h"
#include "interleave.h"

//MFSK varicode instead of psk for PSKR modes
#include "mfskvaricode.h"

//=====================================================================
#define PipeLen			(64)

#define SNTHRESHOLD 6.0
#define AFCDECAYSLOW 8

#define NUM_FILTERS 4

#define MAX_CARRIERS 32

//=====================================================================

class psk : public modem {
private:
// tx & rx
	int			symbollen;
	int			symbits;
	bool			_qpsk;
	bool			_pskr;
	bool			_16psk;
	bool			_8psk;
	bool			_xpsk;
	bool			_disablefec;
	bool			_puncturing;
	int			flushlength;
	double 			separation;
	double			phaseacc[MAX_CARRIERS];
	cmplx			prevsymbol[MAX_CARRIERS];
	unsigned int		shreg;
	//FEC: 2nd stream
	unsigned int		shreg2;
	int					numinterleavers; //interleaver size (speed dependant)
	//double 			numcarriers; //Number of parallel carriers for M CAR PSK and PSKR and QPSKR
	int 			numcarriers; //Number of parallel carriers for M CAR PSK and PSKR and QPSKR
	double 			inter_carrier; // Frequency gap betweeb carriers

// rx variables & functions
	C_FIR_filter		*fir1[MAX_CARRIERS];
	C_FIR_filter		*fir2[MAX_CARRIERS];
	double			*fir1c;
	double			*fir2c;

	Cmovavg			*e0_filt;
	Cmovavg			*e1_filt;
	Cmovavg			*e2_filt;
	Cmovavg			*e3_filt;
	bool			displaysn;

	goertzel		*re_Gbin[NUM_FILTERS];
	goertzel		*im_Gbin[NUM_FILTERS];

	encoder 		*enc;
	viterbi 		*dec;
	//PSKR modes - 2nd Viterbi decoder and 2 receive de-interleaver for comparison
	viterbi 		*dec2;
	interleave		*Rxinlv;
	interleave		*Rxinlv2;
	interleave		*Txinlv;
	unsigned int 		bitshreg;
	int 			rxbitstate;
	//PSKR modes - Soft decoding
	unsigned char		symbolpair[2];
	double			fecmet;
	double			fecmet2;

	bool			vestigial;
	int				sfft_size;
	sfft			*vestigial_sfft;
	cmplx			sfft_bins[11];

	double			phase;
	double			freqerr;
	int				bits;
	double 			bitclk;
	double 			syncbuf[16];
	double 			scope_pipe[2*PipeLen];//[PipeLen];
	unsigned int 	pipeptr;
	unsigned int	dcdshreg;
	//PSKR modes - 2nd stream
	unsigned int		dcdshreg2;

	int 			dcd;
	int				dcdbits;
	cmplx			quality;
	int				acquire;
	int             idepth;

	viewpsk*		pskviewer;
	pskeval*		evalpsk;

	void			rx_symbol(cmplx symbol, int car);
	void 			rx_bit(int bit);
	void 			rx_bit2(int bit);
	void			rx_qpsk(int bits);
	void			rx_pskr(unsigned char symbol);
	double 			scopedata[16];
// IMD & s/n variables
	double			k0, k1, k2;
	double			I11, I12, I21, I22, I31, I32;
	double			snratio, s2n, imdratio, imd;
	double			afcmetric;

//PSKR modes
	bool			firstbit;
	bool			startpreamble;

//PSKR & 8PSK modes
	bool			PSKviterbi;
	double			vphase;
	double			maxamp;


//MULTI-CARRIER
	double			sc_bw; // single carrier bandwidth


//	cmplx thirdorder;
// tx variables & functions
	int			accumulated_bits; //JD for multiple carriers
	int				txsymbols[MAX_CARRIERS];

	double			*tx_shape;
	double			*imd_shape;
	int 			preamble;

	void			transmit(double *buf, int len);
	void			tx_carriers();
	void 			tx_symbol(int sym);
	void			tx_bit(int bit);
	void			tx_xpsk(int bit);
	void			tx_char(unsigned char c);
	void			tx_flush();
	void			update_syncscope();
	void			signalquality();
	void			findsignal();
	void			phaseafc();
	void			afc();
	void			coreafc();
	void			vestigial_afc();

	void			initSN_IMD();
	void			resetSN_IMD();
	void			calcSN_IMD(cmplx z);
	//PSKR modes - for Tx interleaver priming
	void 			clearbits();

protected:
	void			s2nreport(void);

public:
	psk(trx_mode mode);
	~psk();
	void init();
	void rx_init();
	void tx_init();
	void restart();
	int rx_process(const double *buf, int len);
	int tx_process();
	void searchDown();
	void searchUp();

	bool viewer_mode();

	void clear_viewer() { 
		if (pskviewer) pskviewer->clear();
//		if ((mode >= MODE_PSK31 && mode <= MODE_PSK125) ||
//			(mode >= MODE_QPSK31 && mode <= MODE_QPSK125)) pskviewer->clear(); 
	}
	void clear_ch(int n) { pskviewer->clearch(n); }
	int  viewer_get_freq(int n) {
		if (pskviewer) pskviewer->get_freq(n);
		return 0;
	}

};

#endif
