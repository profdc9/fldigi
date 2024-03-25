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

#include "config.h"

#include "charsetlist.h"
#include "tiniconv.h"

const struct charset_info charset_list[] = {
	{ "ASCII", TINICONV_CHARSET_ASCII },
	{ "CP1250", TINICONV_CHARSET_CP1250 },
	{ "CP1251", TINICONV_CHARSET_CP1251 },
	{ "CP1252", TINICONV_CHARSET_CP1252 },
	{ "CP1253", TINICONV_CHARSET_CP1253 },
	{ "CP1254", TINICONV_CHARSET_CP1254 },
	{ "CP1255", TINICONV_CHARSET_CP1255 },
	{ "CP1256", TINICONV_CHARSET_CP1256 },
	{ "CP1257", TINICONV_CHARSET_CP1257 },
	{ "CP1258", TINICONV_CHARSET_CP1258 },
	{ "CP936", TINICONV_CHARSET_CP936 },
	{ "GB2312", TINICONV_CHARSET_GB2312 },
	{ "GBK", TINICONV_CHARSET_GBK },
	{ "ISO-2022-JP", TINICONV_CHARSET_ISO_2022_JP },
	{ "ISO-8859-1", TINICONV_CHARSET_ISO_8859_1 },
	{ "ISO-8859-2", TINICONV_CHARSET_ISO_8859_2 },
	{ "ISO-8859-3", TINICONV_CHARSET_ISO_8859_3 },
	{ "ISO-8859-4", TINICONV_CHARSET_ISO_8859_4 },
	{ "ISO-8859-5", TINICONV_CHARSET_ISO_8859_5 },
	{ "ISO-8859-6", TINICONV_CHARSET_ISO_8859_6 },
	{ "ISO-8859-7", TINICONV_CHARSET_ISO_8859_7 },
	{ "ISO-8859-8", TINICONV_CHARSET_ISO_8859_8 },
	{ "ISO-8859-9", TINICONV_CHARSET_ISO_8859_9 },
	{ "ISO-8859-10", TINICONV_CHARSET_ISO_8859_10 },
	{ "ISO-8859-11", TINICONV_CHARSET_ISO_8859_11 },
	{ "ISO-8859-13", TINICONV_CHARSET_ISO_8859_13 },
	{ "ISO-8859-14", TINICONV_CHARSET_ISO_8859_14 },
	{ "ISO-8859-15", TINICONV_CHARSET_ISO_8859_15 },
	{ "ISO-8859-16", TINICONV_CHARSET_ISO_8859_16 },
	{ "CP866", TINICONV_CHARSET_CP866 },
	{ "KOI8-R", TINICONV_CHARSET_KOI8_R },
	{ "KOI8-RU", TINICONV_CHARSET_KOI8_RU },
	{ "KOI8-U", TINICONV_CHARSET_KOI8_U },
	{ "MACCYRILLIC", TINICONV_CHARSET_MACCYRILLIC },
	{ "UCS-2", TINICONV_CHARSET_UCS_2 },
	{ "UTF-7", TINICONV_CHARSET_UTF_7 },
	{ "UTF-8", TINICONV_CHARSET_UTF_8 },
	{ "CHINESE", TINICONV_CHARSET_CHINESE },
	{ "BIG5", TINICONV_CHARSET_BIG5 },
};

const unsigned int number_of_charsets = sizeof(charset_list)/sizeof(struct charset_info);
