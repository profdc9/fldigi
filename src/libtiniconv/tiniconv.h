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

#ifndef TINICONV_H_
#define TINICONV_H_

/* For exporting functions under WIN32 */
#ifdef WIN32
#ifdef DLL
#define EXPORT(type) __declspec(dllexport) type
#else 
#define EXPORT(type) /*__declspec(dllimport)*/ type
#endif/*DLL*/
#else
#define EXPORT(type) type
#endif /*WIN32*/

typedef unsigned int ucs4_t;
typedef struct tiniconv_ctx_s * conv_t;

/*
 * int xxx_mb2wc (conv_t conv, ucs4_t *pwc, unsigned char const *s, int n)
 * converts the byte sequence starting at s to a wide character. Up to n bytes
 * are available at s. n is >= 1.
 * Result is number of bytes consumed (if a wide character was read),
 * or -1 if invalid, or -2 if n too small, or -2-(number of bytes consumed)
 * if only a shift sequence was read.
 */
typedef int (*xxx_mb2wc_t) (conv_t conv, ucs4_t *pwc, unsigned char const *s, int n);

/*
 * int xxx_flushwc (conv_t conv, ucs4_t *pwc)
 * returns to the initial state and stores the pending wide character, if any.
 * Result is 1 (if a wide character was read) or 0 if none was pending.
 */
typedef int (*xxx_flushwc_t) (conv_t conv, ucs4_t *pwc);

/*
 * int xxx_wc2mb (conv_t conv, unsigned char *r, ucs4_t wc, int n)
 * converts the wide character wc to the character set xxx, and stores the
 * result beginning at r. Up to n bytes may be written at r. n is >= 1.
 * Result is number of bytes written, or -1 if invalid, or -2 if n too small.
 */
typedef int (*xxx_wc2mb_t) (conv_t conv, unsigned char *r, ucs4_t wc, int n);

/*
 * int xxx_reset (conv_t conv, unsigned char *r, int n)
 * stores a shift sequences returning to the initial state beginning at r.
 * Up to n bytes may be written at r. n is >= 0.
 * Result is number of bytes written, or -2 if n too small.
 */
typedef int (*xxx_reset_t) (conv_t conv, unsigned char *r, int n);

typedef unsigned int conv_state_t;

struct tiniconv_ctx_s {
  conv_state_t istate;
  conv_state_t ostate;
  xxx_mb2wc_t mb2wc;
  xxx_flushwc_t flushwc;
  xxx_wc2mb_t wc2mb;
  xxx_reset_t reset;
  int options;
};

/*
 * tiniconv_init
 */

#define TINICONV_CHARSET_ASCII       0
#define TINICONV_CHARSET_CP1250      1
#define TINICONV_CHARSET_CP1251      2
#define TINICONV_CHARSET_CP1252      3
#define TINICONV_CHARSET_CP1253      4
#define TINICONV_CHARSET_CP1254      5
#define TINICONV_CHARSET_CP1255      6
#define TINICONV_CHARSET_CP1256      7
#define TINICONV_CHARSET_CP1257      8
#define TINICONV_CHARSET_CP1258      9
#define TINICONV_CHARSET_CP936       10
#define TINICONV_CHARSET_GB2312      11
#define TINICONV_CHARSET_GBK         12
#define TINICONV_CHARSET_ISO_2022_JP 13
#define TINICONV_CHARSET_ISO_8859_1  14
#define TINICONV_CHARSET_ISO_8859_2  15
#define TINICONV_CHARSET_ISO_8859_3  16
#define TINICONV_CHARSET_ISO_8859_4  17
#define TINICONV_CHARSET_ISO_8859_5  18
#define TINICONV_CHARSET_ISO_8859_6  19
#define TINICONV_CHARSET_ISO_8859_7  20
#define TINICONV_CHARSET_ISO_8859_8  21
#define TINICONV_CHARSET_ISO_8859_9  22
#define TINICONV_CHARSET_ISO_8859_10 23
#define TINICONV_CHARSET_ISO_8859_11 24
#define TINICONV_CHARSET_ISO_8859_13 25
#define TINICONV_CHARSET_ISO_8859_14 26
#define TINICONV_CHARSET_ISO_8859_15 27
#define TINICONV_CHARSET_ISO_8859_16 28
#define TINICONV_CHARSET_CP866       29
#define TINICONV_CHARSET_KOI8_R      30
#define TINICONV_CHARSET_KOI8_RU     31
#define TINICONV_CHARSET_KOI8_U      32
#define TINICONV_CHARSET_MACCYRILLIC 33
#define TINICONV_CHARSET_UCS_2       34
#define TINICONV_CHARSET_UTF_7       35
#define TINICONV_CHARSET_UTF_8       36
#define TINICONV_CHARSET_CHINESE     37
#define TINICONV_CHARSET_BIG5        38
#define TINICONV_CHARSETSIZE         39

#define TINICONV_OPTION_IGNORE_IN_ILSEQ 1 /*< ignore incorrect input sequences */
#define TINICONV_OPTION_IGNORE_OUT_ILSEQ 2 /*< replace sequence which can't be converted to OUT charset with OUTIL_CHAR */
/* #define TINICONV_OPTION_TRANSLIT 4 */
#define TINICONV_OPTION_OUT_ILSEQ_CHAR(ch) (ch << 8)

#define TINICONV_INIT_OK 0
#define TINICONV_INIT_IN_CHARSET_NA -1
#define TINICONV_INIT_OUT_CHARSET_NA -1


#ifdef __cplusplus
extern "C"
{
#endif
EXPORT(int) tiniconv_init(int in_charset_id, int out_charset_id, int options, struct tiniconv_ctx_s *ctx);


/*
 * tiniconv_convert
 */
#define TINICONV_CONVERT_OK 0
#define TINICONV_CONVERT_IN_TOO_SMALL -1
#define TINICONV_CONVERT_OUT_TOO_SMALL -2
#define TINICONV_CONVERT_IN_ILSEQ -3
#define TINICONV_CONVERT_OUT_ILSEQ -4

EXPORT(int) tiniconv_convert(struct tiniconv_ctx_s *ctx,
  unsigned char const *in_buf, int in_size, int *p_in_size_consumed,
  unsigned char *out_buf, int out_size, int *p_out_size_consumed);
#ifdef __cplusplus
}
#endif

#endif /*TINICONV_H_*/
