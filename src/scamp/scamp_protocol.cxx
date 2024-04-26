/*  scamp_protocol.cxx*/

/*
 * Copyright (c) 2024 Daniel Marks

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
 */

#include <string.h>
#include <config.h>
#include <iostream>
#include <fstream>

#include "scamp_protocol.h"

static const uint16_t golay_matrix[12] =
{
    0b110111000101,
    0b101110001011,
    0b011100010111,
    0b111000101101,
    0b110001011011,
    0b100010110111,
    0b000101101111,
    0b001011011101,
    0b010110111001,
    0b101101110001,
    0b011011100011,
    0b111111111110
};

static const uint8_t scamp_6bit_codesymbols[60] = 
{'\0',   '\b',   '\r',    ' ',    '!',   0x22,   0x27,    '(',
  ')',    '*',    '+',    ',',    '-',    '.',    '/',    '0',
  '1',    '2',    '3',    '4',    '5',    '6',    '7',    '8',
  '9',    ':',    ';',    '=',    '?',    '@',    'A',    'B',
  'C',    'D',    'E',    'F',    'G',    'H',    'I',    'J',
  'K',    'L',    'M',    'N',    'O',    'P',    'Q',    'R',
  'S',    'T',    'U',    'V',    'W',    'X',    'Y',    'Z',
 0x5C,   '^',    '`',    '~' };
   /* Last 4 symbols can be interpreted as diacritical marks, 0x5C is diaeresis/umlaut */
   /* 0x27 can be interpreted as acute diacritical mark */

static uint8_t hamming_weight_16(uint16_t n)
{
  uint8_t s = 0, v;
  v = (n >> 8) & 0xFF;
  while (v)
  {
      v &= (v - 1);
      s++;
  }
  v = n & 0xFF;
  while (v)
  {
      v &= (v - 1);
      s++;
  }
  return s;
}

/* calculate the hamming weight a 30 bit number for
   to find hamming distance with sync word */
static uint8_t hamming_weight_30(uint32_t n)
{
  uint8_t s = 0, v;
  v = (n >> 24) & 0x3F;

  while (v)
  {
      v &= (v - 1);
      s++;
  }
  v = (n >> 16) & 0xFF;
  while (v)
  {
      v &= (v - 1);
      s++;
  }
  v = (n >> 8) & 0xFF;
  while (v)
  {
      v &= (v - 1);
      s++;
  }
  v = n & 0xFF;
  while (v)
  {
      v &= (v - 1);
      s++;
  }
  return s;
}

static uint16_t golay_mult(uint16_t wd_enc)
{
   uint16_t enc = 0;
   uint8_t i;
   for (i=12;i>0;)
   {
      i--;
      if (wd_enc & 1) enc ^= golay_matrix[i];
      wd_enc >>= 1;
   }
   return enc;
}

static uint32_t golay_encode(uint16_t wd_enc)
{
  uint16_t enc = golay_mult(wd_enc);
  return (((uint32_t)enc) << 12) | wd_enc;
}

static uint16_t golay_decode(uint32_t codeword, uint8_t *biterrs)
{
  uint16_t enc = codeword & 0xFFF;
  uint16_t parity = codeword >> 12;
  uint8_t i, biterr;
  uint16_t syndrome, parity_syndrome;

  /* if there are three or fewer errors in the parity bits, then
     we hope that there are no errors in the data bits, otherwise
     the error is uncorrected */
  syndrome = golay_mult(enc) ^ parity;
  biterr = hamming_weight_16(syndrome);
  if (biterr <= 3)
  {
     *biterrs = biterr;
     return enc;
  }

  /* check to see if the parity bits have no errors */
  parity_syndrome = golay_mult(parity) ^ enc;
  biterr = hamming_weight_16(parity_syndrome);
  if (biterr <= 3)
  {
     *biterrs = biterr;
     return enc ^ parity_syndrome;
  }

  /* we flip each bit of the data to see if we have two or fewer errors */
  for (i=12;i>0;)
  {
      i--;
      biterr = hamming_weight_16(syndrome ^ golay_matrix[i]);
      if (biterr <= 2)
      {
          *biterrs = biterr+1;
          return enc ^ (((uint16_t)0x800) >> i);
      }
  }

  /* we flip each bit of the parity to see if we have two or fewer errors */
  for (i=12;i>0;)
  {
      i--;
      uint16_t par_bit_synd = parity_syndrome ^ golay_matrix[i];
      biterr = hamming_weight_16(par_bit_synd);
      if (biterr <= 2)
      {
          *biterrs = biterr+1;
          return enc ^ par_bit_synd;
      }
  }
  return 0xFFFF;   /* uncorrectable error */
}

/*  The input word has 24 bits.  The output word has 30 bits, with bits
    1, 5, 9, 13, 17, 21 preceded by its complement bit inserted into
    the word */

static uint32_t scamp_add_reversal_bits(uint32_t codeword)
{
  uint32_t outword = 0;
  uint8_t i;

  for (i=0;i<6;i++)
  {
      if (i>0)
          outword = (outword << 5);
      codeword = (codeword << 4);
      uint8_t temp = (codeword >> 24) & 0x0F;
      outword |= (temp | (((temp & 0x08) ^ 0x08) << 1));
  }
  return outword;
}

/* remove complement bits from inside 30 bit word to yield 24 bit Golay
   code word */

static uint32_t scamp_remove_reversal_bits(uint32_t outword)
{
  uint32_t codeword = 0;
  uint8_t i;

  for (i=0;i<6;i++)
  {
      if (i>0)
      {
          outword = (outword << 5);
          codeword = (codeword << 4);
      }
      uint8_t temp = (outword >> 25) & 0x0F;
      codeword |= temp;
  }
  return codeword;
}

#define SCAMP_IS_DATA_CODE(x) (((uint8_t)((x) >> 8)) == 0x0F)
#define SCAMP_RES_CODE_BITS_SET(x) ((((uint8_t)(x)) & 0x3C) == 0x3C)
#define SCAMP_IS_RES_CODE(x) (SCAMP_RES_CODE_BITS_SET(x) && (!SCAMP_IS_DATA_CODE(x)))

/* decode a 12 bit code word to 8 bit bytes */
static uint8_t scamp_code_word_to_bytes(uint16_t code, uint8_t bytes[])
{
  uint8_t c;
  if (SCAMP_IS_DATA_CODE(code))
  {
    bytes[0] = code & 0xFF;
    return 1;
  }
  c = (code & 0x3F);
  if ((c == 0) || (c >= (sizeof(scamp_6bit_codesymbols)/sizeof(uint8_t))))
    return 0;
  bytes[0] = scamp_6bit_codesymbols[c];
  c = ((code >> 6) & 0x3F);
  if ((c == 0) || (c >= (sizeof(scamp_6bit_codesymbols)/sizeof(uint8_t))))
     return 1;
  bytes[1] = scamp_6bit_codesymbols[c];
  return 2;
}

/* encode 8 bit bytes as 12 bit code words, always encoding as 8-bit raw */
static void scamp_raw_bytes_to_code_words(uint8_t *bytes, uint8_t num_bytes, scamp_code_word_put ecwp, void *st)
{
  uint8_t cur_byte = 0;
  uint8_t frame = 0;
  while (cur_byte < num_bytes)
  {
     uint8_t b = bytes[cur_byte++];
     ecwp( ((uint16_t)(0xF00)) | b, st, cur_byte, frame++ );
  }
}

/* convert character to 6-bit code if it exists */
/* should probably replace this with an inverse look up table later */
static uint8_t scamp_find_code_in_table(uint8_t c)
{
    uint8_t i;
    if ((c >= 'a') && (c <= 'z')) c -= 32;
	if (c =='\n') c = '\r';
	if (c == 127) c = '\b';
    for (i=1;i<(sizeof(scamp_6bit_codesymbols)/sizeof(uint8_t));i++)
        if (scamp_6bit_codesymbols[i] == c) return i;
    return 0xFF;
}

#if 0
/* encode 8 bit bytes to 12 bit code words.
   code words that correspond to a 6-bit symbols are encoded as 6 bits.
   otherwise they are encoded as an 8-bit binary raw data word.
   If a byte that can be encoded as a 6 bit symbol precedes one that can
   not be encoded as a 6 bit symbol, and there is an extra symbol slot
   in the current word, fill it with a zero. */
static void scamp_bytes_to_code_words(uint8_t *bytes, uint8_t num_bytes, scamp_code_word_put ecwp, void *st)
{
  uint8_t cur_byte = 0;
  uint8_t frame = 0;
  uint16_t last_code_word = 0;
  uint16_t code_word;
  while (cur_byte < num_bytes)
  {
     uint8_t b = bytes[cur_byte++];
     uint8_t code1 = scamp_find_code_in_table(b);
     if (code1 == 0xFF)
     {
         code_word = ((uint16_t)(0xF00)) | b;
     } else
     {
         code_word = (uint16_t)code1;
         if (cur_byte < num_bytes)
         {
            b = bytes[cur_byte];
            uint8_t code2 = scamp_find_code_in_table(b);
            if (code2 != 0xFF)
            {
               code_word |= (((uint16_t)code2) << 6);
               cur_byte++;
            }
         }
         if (code_word == last_code_word)
         {
           if (ecwp(0, st, cur_byte, frame++)) break;
         }
     }
     if (ecwp(code_word, st, cur_byte, frame++)) break;
     last_code_word = code_word;
  }
}
#endif

#if 0
typedef struct _scamp_code_word_transmit_data
{
  dsp_dispatch_callback ddc;
  dsp_txmit_message_state *dtms;
} scamp_code_word_transmit_data;

static uint8_t scamp_code_word_transmit(uint16_t code, void *st, uint8_t pos, uint8_t frame)
{
  uint8_t data_code = SCAMP_IS_DATA_CODE(code);
  uint32_t code_30;
  scamp_code_word_transmit_data *scwtd = (scamp_code_word_transmit_data *) st;
  
  code_30 = golay_encode(code);
  code_30 = scamp_add_reversal_bits(code_30);
  scamp_send_frame_rep(code_30, data_code ? 1 : rc.scamp_resend_frames);
  if ((rc.scamp_resync_frames != 0) && (((frame+1) % rc.scamp_resync_frames) == 0))
  {
    scamp_send_frame(SCAMP_INIT_CODEWORD);
    scamp_send_frame(SCAMP_SYNC_CODEWORD);
  }
  while (scwtd->dtms->current_symbol <= pos)
  {
    scwtd->ddc(scwtd->dtms);
    scwtd->dtms->current_symbol++;
  }
  return scwtd->dtms->aborted;
}
#endif

static void scamp_reset_codeword(scamp_state *sc_st)
{
   sc_st->current_bit_no = 0;
   sc_st->current_word = SCAMP_BLANK_CODEWORD;
   sc_st->bitflips_in_phase = 0;
   sc_st->bitflips_lag = 0;
   sc_st->bitflips_lead = 0;
   sc_st->bitflips_ctr = 0;
}

static void scamp_retrain(scamp_state *sc_st)
{
  sc_st->cur_demod_edge_window = sc_st->demod_samples_per_bit;
  scamp_reset_codeword(sc_st);
  sc_st->last_code = 0;
  sc_st->resync = 0;  
}

static void scamp_decode_process(scamp_state *sc_st, uint32_t fr)
{
  uint8_t biterrs, bytes[2], nb;
  uint16_t gf;
  fr = scamp_remove_reversal_bits(fr);
  gf = golay_decode(fr,&biterrs);
  if (gf == 0xFFFF)
  {
    // decode_insert_into_fifo('#');
    return;
  }
  if (!SCAMP_IS_DATA_CODE(gf))
  {
    if (SCAMP_RES_CODE_BITS_SET(gf))
    {
      switch (gf)
      {
        case SCAMP_RES_CODE_END_TRANSMISSION:
            sc_st->reset_protocol = 1;
            sc_st->recv_chars[0] = '\\';
            sc_st->recv_chars[1] = '\r';
            break;
      }
      return;
    }
    if (gf == sc_st->last_code) return;
  }
  sc_st->last_code = gf;
  nb = scamp_code_word_to_bytes(gf, bytes);
  if (nb > 0) sc_st->recv_chars[0] = bytes[0];
  if (nb > 1) sc_st->recv_chars[1] = bytes[1];
}

/* this is called by the interrupt handle to decode the SCAMP frames from
   the spectral channels.  it figures out the magnitude of
   the signal given the current modulation type (OOK/FSK).  it tries to detect
   an edge to determine when the current bit has finished and when it should
   expect the next bit.  it resets the bit counter when the sync signal has
   been received, and when 30 bits of a frame have been received, stores the
   frame in the frame FIFO */
static void scamp_new_sample(scamp_state *sc_st, uint16_t channel_1, uint16_t channel_2)
{
    int16_t demod_sample, temp;
    uint8_t received_bit;
    uint8_t hamming_weight;
    uint16_t max_val;

    /* demodulate a sample either based on the amplitude in one frequency
       channel for OOK, or the difference in amplitude between two frequency
       channels for FSK */

    switch (sc_st->protocol)
    {
#ifdef SCAMP_VERY_SLOW_MODES
        case PROTOCOL_SCAMP_OOK_SLOW:
#endif
        case PROTOCOL_SCAMP_OOK:        demod_sample = channel_1 - sc_st->power_thr;
                                        max_val = channel_1;
                                        break;
#ifdef SCAMP_VERY_SLOW_MODES
        case PROTOCOL_SCAMP_FSK_VSLW:
        case PROTOCOL_SCAMP_FSK_SLOW:
#endif
	    case PROTOCOL_SCAMP_FSK:
	    case PROTOCOL_SCAMP_FSK_FAST:   demod_sample = channel_1 - channel_2;
                                        max_val = channel_2 > channel_1 ? channel_2 : channel_1;
                                        break;
    }
    sc_st->ct_sum += max_val;
    /* This is the automatic "gain" control (threshold level control).
       find the average of a certain number of samples (power of two for calculation
       speed) so that we can determine what to set the thresholds at for the bit on/off
       threshold for ook and the edge threshold for ook/fsk */
    if (sc_st->fsk)
    {
      if ((++sc_st->ct_average) >= (((uint16_t)1) << SCAMP_AVG_CT_PWR2_FSK))
      {
       temp = (sc_st->ct_sum) >> (SCAMP_AVG_CT_PWR2_FSK);
       /* don't allow threshold to get too low, or we'll be having bit edges constantly */
       sc_st->edge_thr = temp > sc_st->power_thr_min ? temp : sc_st->power_thr_min;
       sc_st->power_thr = sc_st->edge_thr >> 1;
       if (sc_st->protocol != PROTOCOL_SCAMP_FSK_FAST) sc_st->edge_thr = (sc_st->edge_thr * 3) / 2;
       sc_st->squelch_thr = sc_st->power_thr; 
       sc_st->ct_average = 0;
       sc_st->ct_sum = 0;
      }
    } else
    {
      if ((++sc_st->ct_average) >= (((uint16_t)1) << SCAMP_AVG_CT_PWR2_OOK))
      {
        temp = (sc_st->ct_sum) >> (SCAMP_AVG_CT_PWR2_OOK);
        /* don't allow threshold to get too low, or we'll be having bit edges constantly */
        sc_st->edge_thr = temp > sc_st->power_thr_min ? temp : sc_st->power_thr_min;
        sc_st->edge_thr >>= 1;
        sc_st->power_thr = sc_st->edge_thr;
        sc_st->squelch_thr = sc_st->power_thr >> 1; 
        sc_st->ct_average = 0;
        sc_st->ct_sum = 0;
      }
    }

    /* calculate the difference between the modulated signal between now and one bit period ago to see
       if there is a bit edge */
    temp = sc_st->demod_buffer[sc_st->demod_sample_no];
    sc_st->bit_edge_val = demod_sample > temp ? (demod_sample - temp) : (temp - demod_sample);
    /* store the demodulated sample into a circular buffer to calculate the edge */
    sc_st->demod_buffer[sc_st->demod_sample_no] = demod_sample;
    if ((++sc_st->demod_sample_no) >= sc_st->demod_samples_per_bit)
        sc_st->demod_sample_no = 0;

    if (sc_st->edge_ctr > 0)  /* count down the edge counter */
        --sc_st->edge_ctr;

    received_bit = 0;

    if (sc_st->edge_ctr < sc_st->cur_demod_edge_window)  /* if it below the edge window, start looking for the bit edge */
    {
        if (sc_st->bit_edge_val > sc_st->edge_thr)  /* the edge should come around now, does it exceed threshold */
        {
            if (sc_st->bit_edge_val > sc_st->max_bit_edge_val)  /* if so, have we reached the peak of the edge */
            {
                sc_st->max_bit_edge_val = sc_st->bit_edge_val;  /* if so, reset the edge center counter */
                sc_st->next_edge_ctr = 1;
                sc_st->cur_bit = demod_sample;              /* save the bit occurring at the edge */
            } else
                sc_st->next_edge_ctr++;                     /* otherwise count that we have passed the edge peak */
        } else
        {
            if (sc_st->max_bit_edge_val != 0)   /* if we have detected an edge */
            {
                sc_st->edge_ctr = sc_st->demod_samples_per_bit > sc_st->next_edge_ctr ? sc_st->demod_samples_per_bit - sc_st->next_edge_ctr : 0;
                                /* reset edge ctr to look for next edge */
                sc_st->max_bit_edge_val = 0;    /* reset max_bit_edge_val for next edge peak */
                received_bit = 1;
                sc_st->cur_demod_edge_window = sc_st->demod_edge_window;
            } else /* we haven't detected an edge */
            {
                if (sc_st->edge_ctr == 0)
                {
                    sc_st->cur_bit = demod_sample;              /* save the bit */
                    sc_st->edge_ctr = sc_st->demod_samples_per_bit;  /* reset and wait for the next bit edge to come along */
                    received_bit = 1;                       /* an edge hasn't been detected but a bit interval happened */
                }
            }
        }
    }
    if (!received_bit)   /* no bit available yet, wait */
        return;

    /* add the bit to the current word */
    sc_st->current_word = (sc_st->current_word << 1) | (sc_st->polarity ^ (sc_st->cur_bit > 0));
    sc_st->bitflips_ctr++;
    /* this is done on the bit before it is needed to reduce worst case latency */
    if (sc_st->bitflips_ctr == 4)
    {
        /* keep track of the number of complement bits in the word. should be 6 */
        uint8_t maskbits, lowerbyte = ((uint8_t)sc_st->current_word);
        maskbits = lowerbyte & 0x0C;
        sc_st->bitflips_in_phase += (maskbits == 0x08) || (maskbits == 0x04);
        maskbits = lowerbyte & 0x18;
        sc_st->bitflips_lag += (maskbits == 0x10) || (maskbits == 0x08);
        maskbits = lowerbyte & 0x06;
        sc_st->bitflips_lead += (maskbits == 0x04) || (maskbits == 0x02);
    } else if (sc_st->bitflips_ctr >= 5)
        sc_st->bitflips_ctr = 0;

    uint8_t thr = max_val >= sc_st->squelch_thr;
    
    if (thr)
      sc_st->threshold_counter = SCAMP_THRESHOLD_COUNTER_MAX;
    else
    {
      if (sc_st->threshold_counter > 0)
          sc_st->threshold_counter--;
    }
    if (sc_st->reset_protocol != 0)
    {
      scamp_retrain(sc_st);
      sc_st->reset_protocol = 0;
    }
    if ((!sc_st->fsk) || thr)  /* if there is a bit to sync to */
    {
      hamming_weight = hamming_weight_30(sc_st->current_word ^ SCAMP_SYNC_CODEWORD);
      if (hamming_weight < 4)  /* 30-bit resync word has occurred! */
      {
        scamp_reset_codeword(sc_st);
        sc_st->resync = 1;
        return;
      }
      /* if we 15 of the last 16 bits zeros with fsk, that means we have a reversed polarity start */
      hamming_weight = hamming_weight_16(sc_st->current_word);
      if ((hamming_weight < 2) && (sc_st->fsk))
      {
        sc_st->polarity = !sc_st->polarity; /* reverse the polarity of FSK frequencies and 0/1 */
        scamp_retrain(sc_st);
        return;
      }
      /* if we have 15 of the last 16 bits ones, that is a ook key down start, or a fsk start */
      if (hamming_weight >= 15)
      {
        scamp_retrain(sc_st);
        return;
      }
    } 
    /* if we have synced, and we have 30 bits, we have a frame */
    sc_st->current_bit_no++;
    if ((sc_st->current_bit_no >= 30) && (sc_st->resync))  /* we have a complete frame */
    {
       if ((sc_st->bitflips_lead > sc_st->bitflips_in_phase) && (sc_st->bitflips_lead >= 5))
        /* we are at least one bit flip ahead, we probably registered a spurious bit flip */
       {
         /* back up and get one more bit.*/
         sc_st->current_bit_no--;
         sc_st->bitflips_in_phase = sc_st->bitflips_lead;  /*lead now becomes in phase */
         sc_st->bitflips_lag = 0;  /* lag is now two behind, so we can't use it */
         sc_st->bitflips_lead = 0; /* clear bit_lead so we don't try a second time */
       } else
       {
          if ((sc_st->bitflips_lag > sc_st->bitflips_in_phase) && (sc_st->bitflips_lag >= 5))
          /* we are at least one bit flip short, we probably fell a bit behind */
          {
            if (sc_st->threshold_counter > 0)
            {
			   scamp_decode_process(sc_st, sc_st->current_word >> 1);
             // scamp_insert_into_frame_fifo(&sc_st->scamp_output_fifo, sc_st->current_word >> 1);
		    }
            /* start with the next word with one flip */
            sc_st->current_bit_no = 1;
            sc_st->bitflips_ctr = 1;
         } else
         {
             if (sc_st->bitflips_in_phase >= 4)
             {
               /* otherwise we just place in buffer and the code word is probably aligned */
               if (sc_st->threshold_counter > 0)
               {
  			     scamp_decode_process(sc_st, sc_st->current_word);
                 //scamp_insert_into_frame_fifo(&sc_st->scamp_output_fifo, sc_st->current_word);
               }
             }
             sc_st->current_bit_no = 0;
             sc_st->bitflips_ctr = 0;
          }
          sc_st->bitflips_in_phase = 0;
          sc_st->bitflips_lag = 0;
          sc_st->bitflips_lead = 0;
       }
    }
}

void scamp_add_frame_rep(scamp_state *sc_st, uint32_t frame, uint8_t rep)
{
	while (rep > 0)
	{
		if (sc_st->frames_num >= SCAMP_FRAMES_MAX)
			return;
		sc_st->frames[sc_st->frames_num++] = frame;
		rep--;
	}
}

void scamp_add_code_rep(scamp_state *sc_st, uint16_t code, uint8_t rep)
{
  uint32_t frame = golay_encode(code);
  frame = scamp_add_reversal_bits(frame);
  scamp_add_frame_rep(sc_st, frame, rep);
}

void scamp_send_start(scamp_state *sc_st)
{
  if (sc_st->fsk) 
  {
    scamp_add_frame_rep(sc_st, SCAMP_SOLID_CODEWORD, 1);
  } else
  {
    for (uint8_t i=0;i<4;i++)
      scamp_add_frame_rep(sc_st, SCAMP_DOTTED_CODEWORD, 1);
  }
  scamp_add_frame_rep(sc_st, SCAMP_INIT_CODEWORD, sc_st->repeat_frames);
  scamp_add_frame_rep(sc_st, SCAMP_SYNC_CODEWORD, sc_st->repeat_frames);
  sc_st->trans_state = SCAMP_TRANS_STATE_TRANS;
}

void scamp_send_code(scamp_state *sc_st, uint16_t code, uint8_t rep)
{
	scamp_add_code_rep(sc_st, code, rep ? sc_st->repeat_frames : 1);
	if (sc_st->resync_frames > 0)
	{
		if ((++sc_st->resync_frames_count) >= sc_st->resync_frames)
		{
			sc_st->resync_frames_count = 0;
			scamp_add_frame_rep(sc_st, SCAMP_INIT_CODEWORD, 1);
			scamp_add_frame_rep(sc_st, SCAMP_SYNC_CODEWORD, 1);
			sc_st->trans_state = SCAMP_TRANS_STATE_TRANS;
		}
	}
}

uint8_t scamp_txmit(scamp_state *sc_st, int16_t c)
{ 
  sc_st->frames_num = 0;
  if (c == -1) /* end of transmission */
  {
	  if (sc_st->trans_state == SCAMP_TRANS_STATE_WAIT_CHAR)
		 scamp_send_code(sc_st, sc_st->code_word, 1);
	  if (sc_st->trans_state != SCAMP_TRANS_STATE_IDLE)
		 scamp_add_frame_rep(sc_st, SCAMP_RES_CODE_END_TRANSMISSION_FRAME, sc_st->repeat_frames);	
      sc_st->trans_state = SCAMP_TRANS_STATE_IDLE;
      sc_st->resync_frames_count = 0;
      sc_st->reset_protocol = 1;
      return sc_st->frames_num;
  }
  if (c == -2) /* idle state */
  {
	  if (sc_st->trans_state == SCAMP_TRANS_STATE_WAIT_CHAR)
		 scamp_send_code(sc_st, sc_st->code_word, 1);
	  if (sc_st->trans_state != SCAMP_TRANS_STATE_IDLE)
		 scamp_send_code(sc_st, sc_st->code_word, 0);
      return sc_st->frames_num;
  }
  if ((c >= 'a') && (c <= 'z')) c -= 32;
  if ((sc_st->trans_state == SCAMP_TRANS_STATE_WAIT_CHAR))
  {
	 if ((c > 0) && (c <= 0xFF))
	 {
		uint8_t code = scamp_find_code_in_table(c);
		if (code == 0xFF)
		{
			scamp_send_code(sc_st, sc_st->code_word, 1);
			sc_st->code_word = ((uint16_t)(0xF00)) | c;     
			scamp_send_code(sc_st, sc_st->code_word, 1);
		} else
		{
			sc_st->code_word |= (((uint16_t)code) << 6);
			scamp_send_code(sc_st, sc_st->code_word, 1);
        }
     }
#if 0 
	 if ((c >= 0x100) && (c <= 0xFFF))
	 {
		scamp_send_code(sc_st, sc_st->code_word, 1);
		sc_st->code_word = c;
		scamp_send_code(sc_st, sc_st->code_word, 1);
	 }
#endif
	 sc_st->trans_state = SCAMP_TRANS_STATE_TRANS;
     return sc_st->frames_num;   
  }	      			
  if ((sc_st->trans_state == SCAMP_TRANS_STATE_IDLE) || (sc_st->trans_state == SCAMP_TRANS_STATE_TRANS))
  {
	 if (sc_st->trans_state == SCAMP_TRANS_STATE_IDLE)
		scamp_send_start(sc_st);
	 if ((c > 0) && (c <= 0xFF))
	 {
		uint8_t code = scamp_find_code_in_table(c);
		if (code == 0xFF)
		{
			sc_st->code_word = ((uint16_t)(0xF00)) | c;     
			scamp_send_code(sc_st, sc_st->code_word, 1);
		}
		else
		{
			sc_st->code_word = code;  
			sc_st->trans_state = SCAMP_TRANS_STATE_WAIT_CHAR;
		}
	 }
#if 0 
	 if ((c >= 0x100) && (c <= 0xFFF))
	 {
		sc_st->code_word = c;
		scamp_send_code(sc_st, sc_st->code_word, 1);
	 }
#endif
     return sc_st->frames_num;   
   }
   return sc_st->frames_num;
}
     
void SCAMP_protocol::init(uint8_t protocol)
{
	memset(&sc,'\000',sizeof(sc));
	sc.protocol = protocol;
    switch (sc.protocol)
    {
        case PROTOCOL_SCAMP_OOK:        sc.demod_samples_per_bit = 64 / 4;
                                        sc.fsk = 0;
                                        sc.demod_edge_window = 4;
                                        break;
        case PROTOCOL_SCAMP_FSK:        sc.demod_samples_per_bit = 60 / 4;
                                        sc.fsk = 1;
                                        sc.demod_edge_window = 4;
                                        break;
        case PROTOCOL_SCAMP_FSK_FAST:   sc.demod_samples_per_bit = 24 / 4;
                                        sc.fsk = 1;
                                        sc.demod_edge_window = 2;
                                        break;
#ifdef SCAMP_VERY_SLOW_MODES
        case PROTOCOL_SCAMP_OOK_SLOW:   sc.demod_samples_per_bit =  144 / 4;
                                        sc.demod_edge_window = 4;
                                        sc.fsk = 0;
                                        break;
        case PROTOCOL_SCAMP_FSK_VSLW:   
        case PROTOCOL_SCAMP_FSK_SLOW:   sc.demod_samples_per_bit =  144 / 4;
                                        sc.demod_edge_window = 4;
                                        sc.fsk = 1;
                                        break;
#endif // SCAMP_VERY_SLOW_MODES
    }

    sc.power_thr_min = ((uint16_t)sc.demod_samples_per_bit * 4) * (sc.fsk ? SCAMP_PWR_THR_DEF_FSK : SCAMP_PWR_THR_DEF_OOK);
    if (sc.fsk)
       sc.edge_thr = sc.power_thr;
    else
       sc.edge_thr = sc.power_thr << 1;
    scamp_retrain(&sc);
}

void SCAMP_protocol::decode_process(double mag1, double mag2, int recv_chars[2])
{
	sc.recv_chars[0] = sc.recv_chars[1] = -1;
	uint16_t imag1 = (uint16_t)(mag1 * 16384.0);
	uint16_t imag2 = (uint16_t)(mag2 * 16384.0);
	scamp_new_sample(&sc, imag1, imag2); 
    recv_chars[0] = sc.recv_chars[0];
    recv_chars[1] = sc.recv_chars[1];
}

void SCAMP_protocol::set_resync_repeat_frames(int resync_frames, int repeat_frames)
{ 
	if ((resync_frames < 0) || (resync_frames > 9))
		resync_frames = 0;
	if ((repeat_frames < 1) || (repeat_frames > 9))
		repeat_frames = 1;
	sc.resync_frames = resync_frames;
	sc.repeat_frames = repeat_frames;
}

int SCAMP_protocol::send_char(int c, uint32_t *fr[])
{
	*fr = sc.frames;
	return scamp_txmit(&sc, c);
}
