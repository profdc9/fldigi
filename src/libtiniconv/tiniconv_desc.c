/*
 * This file is part of the TINICONV Library.
 *
 * The TINICONV Library is free software; you can redistribute it
 * and/or modify it under the terms of the GNU Library General Public
 * License version 2 as published by the Free Software Foundation.
 */
// ----------------------------------------------------------------------------
// Copyright (C) 2014
//              David Freese, W1HKJ
//
// This file is part of fldigi
//
// fldigi is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ----------------------------------------------------------------------------

#include "tiniconv_int.h"

#include <assert.h>

#define abort() assert(0)
#define NULL 0

#include "encdec/ascii.h"

#include "encdec/cp1250.h"
#include "encdec/cp1251.h"
#include "encdec/cp1252.h"
#include "encdec/cp1253.h"
#include "encdec/cp1254.h"
#include "encdec/cp1255.h"
#include "encdec/cp1256.h"
#include "encdec/cp1257.h"
#include "encdec/cp1258.h"

#include "encdec/iso8859_1.h"
#include "encdec/iso8859_2.h"
#include "encdec/iso8859_3.h"
#include "encdec/iso8859_4.h"
#include "encdec/iso8859_5.h"
#include "encdec/iso8859_6.h"
#include "encdec/iso8859_7.h"
#include "encdec/iso8859_8.h"
#include "encdec/iso8859_9.h"
#include "encdec/iso8859_10.h"
#include "encdec/iso8859_11.h"
#include "encdec/iso8859_13.h"
#include "encdec/iso8859_14.h"
#include "encdec/iso8859_15.h"
#include "encdec/iso8859_16.h"

#include "encdec/cp866.h"
#include "encdec/koi8_r.h"
#include "encdec/koi8_ru.h"
#include "encdec/koi8_u.h"
#include "encdec/mac_cyrillic.h"

#include "encdec/ucs2.h"
#include "encdec/utf7.h"
#include "encdec/utf8.h"

#include "encdec/gb2312.h" /* is needed for euc_cn.h */
#include "encdec/euc_cn.h"
#include "encdec/gbk.h"
#include "encdec/ces_gbk.h"
#include "encdec/big5.h" /* is needed for ces_big5.h */
#include "encdec/ces_big5.h"
#include "encdec/jisx0208.h"
#include "encdec/jisx0201.h"

#include "encdec/cp936.h"
#include "encdec/iso2022_jp.h"

const struct tiniconv_charset_map_entry_s tiniconv_charset_map[] =
{
  {ascii_mbtowc,        NULL,           ascii_wctomb,        NULL            }, /* 0 */
  {cp1250_mbtowc,       NULL,           cp1250_wctomb,       NULL            }, /* 1 */
  {cp1251_mbtowc,       NULL,           cp1251_wctomb,       NULL            }, /* 2 */
  {cp1252_mbtowc,       NULL,           cp1252_wctomb,       NULL            }, /* 3 */
  {cp1253_mbtowc,       NULL,           cp1253_wctomb,       NULL            }, /* 4 */
  {cp1254_mbtowc,       NULL,           cp1254_wctomb,       NULL            }, /* 5 */
  {cp1255_mbtowc,       cp1255_flushwc, cp1255_wctomb,       NULL            }, /* 6 */
  {cp1256_mbtowc,       NULL,           cp1256_wctomb,       NULL            }, /* 7 */
  {cp1257_mbtowc,       NULL,           cp1257_wctomb,       NULL            }, /* 8 */
  {cp1258_mbtowc,       cp1258_flushwc, cp1258_wctomb,       NULL            }, /* 9 */
  {cp936_mbtowc,        NULL,           cp936_wctomb,        NULL            }, /* 10 */
  {euc_cn_mbtowc,       NULL,           euc_cn_wctomb,       NULL            }, /* 11 */
  {gbk_mbtowc,          NULL,           gbk_wctomb,          NULL            }, /* 12 */
  {iso2022_jp_mbtowc,   NULL,           iso2022_jp_wctomb,   iso2022_jp_reset}, /* 13 */
  {iso8859_1_mbtowc,    NULL,           iso8859_1_wctomb,    NULL            }, /* 14 */
  {iso8859_2_mbtowc,    NULL,           iso8859_2_wctomb,    NULL            }, /* 15 */
  {iso8859_3_mbtowc,    NULL,           iso8859_3_wctomb,    NULL            }, /* 16 */
  {iso8859_4_mbtowc,    NULL,           iso8859_4_wctomb,    NULL            }, /* 17 */
  {iso8859_5_mbtowc,    NULL,           iso8859_5_wctomb,    NULL            }, /* 18 */
  {iso8859_6_mbtowc,    NULL,           iso8859_6_wctomb,    NULL            }, /* 19 */
  {iso8859_7_mbtowc,    NULL,           iso8859_7_wctomb,    NULL            }, /* 20 */
  {iso8859_8_mbtowc,    NULL,           iso8859_8_wctomb,    NULL            }, /* 21 */
  {iso8859_9_mbtowc,    NULL,           iso8859_9_wctomb,    NULL            }, /* 22 */
  {iso8859_10_mbtowc,   NULL,           iso8859_10_wctomb,   NULL            }, /* 23 */
  {iso8859_11_mbtowc,   NULL,           iso8859_11_wctomb,   NULL            }, /* 24 */
  {iso8859_13_mbtowc,   NULL,           iso8859_13_wctomb,   NULL            }, /* 25 */
  {iso8859_14_mbtowc,   NULL,           iso8859_14_wctomb,   NULL            }, /* 26 */
  {iso8859_15_mbtowc,   NULL,           iso8859_15_wctomb,   NULL            }, /* 27 */
  {iso8859_16_mbtowc,   NULL,           iso8859_16_wctomb,   NULL            }, /* 28 */
  {cp866_mbtowc,        NULL,           cp866_wctomb,        NULL            }, /* 29 */
  {koi8_r_mbtowc,       NULL,           koi8_r_wctomb,       NULL            }, /* 30 */
  {koi8_ru_mbtowc,      NULL,           koi8_ru_wctomb,      NULL            }, /* 31 */
  {koi8_u_mbtowc,       NULL,           koi8_u_wctomb,       NULL            }, /* 32 */
  {mac_cyrillic_mbtowc, NULL,           mac_cyrillic_wctomb, NULL            }, /* 33 */
  {ucs2_mbtowc,         NULL,           ucs2_wctomb,         NULL            }, /* 34 */
  {utf7_mbtowc,         NULL,           utf7_wctomb,         utf7_reset      }, /* 35 */
  {utf8_mbtowc,         NULL,           utf8_wctomb,         NULL            }, /* 36 */
  {gb2312_mbtowc,       NULL,           gb2312_wctomb,       NULL            }, /* 37, CHINESE */
  {ces_big5_mbtowc,     NULL,           ces_big5_wctomb,     NULL            }, /* 38 */
  {NULL,                NULL,           NULL,                NULL            }
};
