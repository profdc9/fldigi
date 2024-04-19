// ----------------------------------------------------------------------------
// scampsk.cxx  --  SCAMPSK signal generator
//
// Copyright (C) 2021
//		Dave Freese, W1HKJ
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#include <iostream>
#include <string>

//#include "io_timer.h"
#include "serial.h"

//#include "io_timer.h"

#ifndef SCAMPSK_H
#define SCAMPSK_H

class SCAMPSK
{
// time (msec) for one symbollen @ 45.45 baud
#define BITLEN 0.022 //22

#define SCAMPSK_UNKNOWN	0x000
#define	SCAMPSK_LETTERS	0x100
#define	SCAMPSK_FIGURES	0x200

#define LTRS 0x1F
#define FIGS 0x1B

#define SCAMPSK_MARK  1
#define SCAMPSK_SPACE 0

public:
	SCAMPSK();
	~SCAMPSK();

	bool open();
	bool close();
	void abort();

	bool sending();

	void send(const char ch);
	void send(std::string s);
	void append(const char ch);
	void append(std::string s);

	void shift_on_space(bool b) { _shift_on_space = b; }
	bool shift_on_space() { return _shift_on_space; }

	void dtr(bool b) { _dtr = b; }
	bool dtr() { return _dtr; }

	void rts(bool b) { _dtr = !b; }
	bool rts() { return !_dtr; }

	void reverse(bool b) { _reverse = b; }
	bool reverse() { return _reverse; }

	void   open_port(std::string device_name);

	void device(std::string device_name) {
		serial_device = device_name;
	}

	void scampsk_shares_port(Cserial *shared_device);

//	size_t io_timer_id;

	Cserial *using_port() { return scampsk_port; }

private:

	Cserial		*scampsk_port;

	std::string serial_device;

	bool   shared_port; // default is false
	static char letters[];
	static char figures[];
	static const char *ascii[];

	int  shift;
	bool _shift_on_space;
	bool _dtr;
	bool _reverse;

	int mode;
	int  shift_state;
	int start_bits;
	int stop_bits;
	int chr_bits;

    std::string str_buff;
	int   chr_out;
    int   scampsk_chr;

	int baudot_enc(int);
	void send_baudot(int ch);
	void scampsk_out (bool);

	double scampsk_now();

	void scampskbit (int, double);

	bool _sending;
public:

	int callback_method();

	int init_scampsk_thread();
	void exit_scampsk_thread();

	bool	scampsk_loop_terminate;

	pthread_t scampsk_thread;

friend
	void *scampsk_loop(void *data);

};

#endif
