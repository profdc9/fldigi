// ----------------------------------------------------------------------------
// anal.h  --  frequency analysis modem
//
// Copyright (C) 2006
//		Dave Freese, W1HKJ
//
// This file is part of fldigi.
//
// Modified for data file creation / analysis JC Gibbons N8OBJ  5/14/19
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

#ifndef _anal_H
#define _anal_H

#include <string>
#include <ctime>
#include <stdio.h>

#include "complex.h"
#include "filters.h"
#include "fftfilt.h"
#include "modem.h"


#define ANAL_SAMPLERATE	8000
#define FILT_LEN 	4		//seconds
#define PIPE_LEN	120
#define DSP_CNT		1		//seconds
#define ANAL_BW		4

class anal : public modem {
private:

	double	phaseacc;

	fftfilt *bpfilt;
	Cmovavg *ffilt;
	Cmovavg *afilt;

	double pipe[PIPE_LEN];

	double prevsymbol;
	cmplx prevsmpl;

	double	fout;
	double	amp;
	long int wf_freq;
	double   dspcnt;
	long int passno;
	int		rxcorr;

	struct timespec start_time;
	struct tm File_Start_Date;

	double elapsed;

	void clear_syncscope();
	inline cmplx mixer(cmplx in);
	int rx(bool bit);

	double nco(double freq);
	void writeFile();
	void createfilename();

	char FileDate [20];
	char FileData [20];

public:
	anal();
	~anal();
	void init();
	void rx_init();
	void tx_init();
	void restart();
	void start_csv();
	void stop_csv();

	int rx_process(const double *buf, int len);
	int tx_process();
	std::string	analysisFilename;
	std::string 	OpenAnalalysisFile;

};

#endif
