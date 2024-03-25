// ----------------------------------------------------------------------------
// Copyright (C) 2014
//              David Freese, W1HKJ
//
// This file is part of fldigi
//
// fldigi is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// fldigi is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ----------------------------------------------------------------------------

#ifndef ADIFIO
#define ADIFIO

#include <cstdio>
#include <cstring>
#include <string>

#include "qso_db.h"

#define ADIF_VERS "2.2.7"

class cAdifIO {
private:
	bool write_all;
	cQsoRec *adifqso;
	FILE *adiFile;
	char *fillfield(int, int, char *);
	static int instances;
public:
	cAdifIO ();
	~cAdifIO ();
	int readAdifRec () {return 0;};
	int  writeAdifRec (cQsoRec *rec, const char *fname);
	void readFile (const char *, cQsoDb *);
	void do_readfile(const char *, cQsoDb *);
	void do_writelog();
	std::string adif_record(cQsoRec *rec);

	int writeFile (const char *, cQsoDb *);
	int writeLog (const char *, cQsoDb *, bool b = true);
	bool log_changed(const char *fname);
};

// crc 16 cycle redundancy check sum

class Ccrc16 {
private:
	unsigned int crcval;
	char ss[5];
public:
	Ccrc16() { crcval = 0xFFFF; }
	~Ccrc16() {};
	void reset() { crcval = 0xFFFF;}
	unsigned int val() {return crcval;}
	std::string sval() {
		snprintf(ss, sizeof(ss), "%04X", crcval);
		return ss;
	}
	void update(char c) {
		crcval ^= c & 255;
		for (int i = 0; i < 8; ++i) {
			if (crcval & 1)
				crcval = (crcval >> 1) ^ 0xA001;
			else
				crcval = (crcval >> 1);
		}
	}
	unsigned int crc16(char c) { 
		update(c); 
		return crcval;
	}
	unsigned int crc16(std::string s) {
		reset();
		for (size_t i = 0; i < s.length(); i++)
			update(s[i]);
		return crcval;
	}
	std::string scrc16(std::string s) {
		crc16(s);
		return sval();
	}
};


#endif
