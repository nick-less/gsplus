/************************************************************************/
/*			KEGS: Apple //gs Emulator			*/
/*			Copyright 2002-2021 by Kent Dickey		*/
/*									*/
/*	This code is covered by the GNU GPL v3				*/
/*	See the file COPYING.txt or https://www.gnu.org/licenses/	*/
/*	This program is provided with no warranty			*/
/*									*/
/*	The KEGS web page is kegs.sourceforge.net			*/
/*	You may contact the author at: kadickey@alumni.princeton.edu	*/
/************************************************************************/

#include "defc.h"
#include "sound.h"

#include "mockingboard.h"
// Mockingboard contains two pairs of a 6522/AY-8913, where the 6522 interfaces
//  the AY-8913 (which makes the sounds) to the Apple II.  Each AY-8913
//  contains 3 channels of sound: A,B,C.  Model each pair separately.
// The AY-8913 has 16 registers.  The documentation numbers them using octal!
// The AY8913 is accessed using ORB of the 6522 as control: 0 = Reset, 4=Idle
//	5=Reg read; 6=Write reg; 7=Latch reg address

extern Mockingboard g_mockingboard;
dword64 g_mockingboard_last_int_dcycs =  0;
dword64 g_mockingboard_event_int_dcycs = 0; // 0 -> no event pending

extern int g_irq_pending;
extern double g_dsamps_per_dcyc;
extern double g_cur_dcycs;
extern double	g_dcycs_per_samp;


Mockingboard g_mockingboard;

// The AY8913 chip has non-linear amplitudes (it has 16 levels) and the
//  documentation does not match measured results.  But all the measurements
//  should really be done at the final speaker/jack since all the stuff in
//  the path affects it.  But: no one's done this for Mockingboard that I
//  have found, so I'm taking measurements from the AY8913 chip itself.
// AY8913 amplitudes from https://groups.google.com/forum/#!original/
//                              comp.sys.sinclair/-zCR2kxMryY/XgvaDICaldUJ
// by Matthew Westcott on December 21, 2001.
double g_ay8913_ampl_factor_westcott[16] = {            // NOT USED
        0.000,  // level[0]
        0.010,  // level[1]
        0.015,  // level[2]
        0.022,  // level[3]
        0.031,  // level[4]
        0.046,  // level[5]
        0.064,  // level[6]
        0.106,  // level[7]
        0.132,  // level[8]
        0.216,  // level[9]
        0.297,  // level[10]
        0.391,  // level[11]
        0.513,  // level[12]
        0.637,  // level[13]
        0.819,  // level[14]
        1.000,  // level[15]
};

// https://sourceforge.net/p/fuse-emulator/mailman/message/34065660/
//  refers to some Russian-language measurements at:
//  http://forum.tslabs.info/viewtopic.php?f=6&t=539 (translate from
//  Russian), they give:
// 0000,028F,03B3,0564, 07DC,0BA9,1083,1B7C,
// 2068,347A,4ACE,5F72, 7E16,A2A4,CE3A,FFFF
double g_ay8913_ampl_factor[16] = {
        0.000,  // level[0]
        0.010,  // level[1]
        0.014,  // level[2]
        0.021,  // level[3]
        0.031,  // level[4]
        0.046,  // level[5]
        0.064,  // level[6]
        0.107,  // level[7]
        0.127,  // level[8]
        0.205,  // level[9]
        0.292,  // level[10]
        0.373,  // level[11]
        0.493,  // level[12]
        0.635,  // level[13]
        0.806,  // level[14]
        1.000,  // level[15]
};
#define VAL_MOCK_RANGE		(39000)

// MAME also appears to try to figure out how the channels get "summed"
//  together.  KEGS code adds them in a completely independent way, and due
//  to the circuit used on the s, this is certainly incorrect.
#define MAX_MOCK_ENV_SAMPLES    2000
int     g_mock_env_vol[MAX_MOCK_ENV_SAMPLES];
byte    g_mock_noise_bytes[MAX_MOCK_ENV_SAMPLES];
int     g_mock_volume[16];              // Sample for each of the 16 amplitudes

word32	g_last_mock_vbl_count = 0;
extern word32	g_last_c030_vbl_count;

void remove_event_mockingboard(void) {
  (void)remove_event_entry(EV_MOCKINGBOARD);
}

void add_event_mockingboard(double dcycs) {
  if (dcycs < g_cur_dcycs) {
    dcycs = g_cur_dcycs;
  }
  add_event_entry(dcycs, EV_MOCKINGBOARD);
}


// from sound loop 
/*
	num_pairs = 0;
	// Do Mockinboard channels
	for(i = 0; i < 2; i++) {			// Pair: 0 or 1
		ay8913ptr = &(g_mockingboard.pair[i].ay8913);
		for(j = 0; j < 3; j++) {		// Channels: A, B, or C
			if((ay8913ptr->regs[8 + j] & 0x1f) == 0) {
				continue;
			}
			num_pairs = 2;
			g_last_mock_vbl_count = g_vbl_count;
			break;
		}
	}
	if((g_vbl_count - g_last_mock_vbl_count) < 120) {
		// Keep playing for 2 seconds, to avoid some static issues
		num_pairs = 2;
	}
	if(num_pairs) {
		sound_mask = -1;
		if(snd_buf_init == 0) {
			sound_mask = 0;
			snd_buf_init++;
		}
		outptr = outptr_start;
		ivol = -((VAL_MOCK_RANGE * 3 / (8 * 15)) * g_doc_vol);
			// Do 3/8 of range below 0, leaving 5/8 above 0
		for(i = 0; i < num_samps; i++) {
			outptr[0] = (outptr[0] & sound_mask) + ivol;
			outptr[1] = (outptr[1] & sound_mask) + ivol;
			outptr += 2;
		}
		for(i = 0; i < 16; i++) {
			dvolume = (g_doc_vol * VAL_MOCK_RANGE) / (15.0 * 3.0);
			ivol = (int)(g_ay8913_ampl_factor[i] * dvolume);
			g_mock_volume[i] = ivol;
		}
	}
	for(i = 0; i < num_pairs; i++) {
		if(g_mockingboard.disable_mask) {
			printf("dsamp:%lf\n", dsamps);
		}

		sound_mock_envelope(i, &(g_mock_env_vol[0]), num_samps,
							&(g_mock_volume[0]));
		sound_mock_noise(i, &(g_mock_noise_bytes[0]), num_samps);
		for(j = 0; j < 3; j++) {
			sound_mock_play(i, j, outptr_start,
				&(g_mock_env_vol[0]), &(g_mock_noise_bytes[0]),
				&(g_mock_volume[0]), num_samps);
		}
	}
*/

void
sound_mock_envelope(int pair, int *env_ptr, int num_samps, int *vol_ptr)
{
	Ay8913	*ay8913ptr;
	double	dmul, denv_period;
	dword64	env_dsamp, dsamp_inc;
	word32	ampl, eff_ampl, reg13, env_val, env_period;
	int	i;

	// This routine calculates a fixed-point increment to apply
	//  to env_dsamp, where the envelope value is in bits 44:40 (bit
	//  44 is to track the alternating waveform, 43:40 is the env_ampl).
	// This algorithm does not properly handle dynamically changing the
	//  envelope period in the middle of a step.  In the AY8913, the
	//  part counts up to the env_period, and if the period is changed
	//  to a value smaller than the current count, it steps immediately
	//  to the next step.  This routine will wait for enough fraction
	//  to accumulate before stepping.  At most, this can delay the step
	//  by almost the new count time (if the new period is smaller), but
	//  no more.  I suspect this is not noticeable.
	if(num_samps > MAX_MOCK_ENV_SAMPLES) {
		halt_printf("envelope overflow!: %d\n", num_samps);
		return;
	}

	ay8913ptr = &(g_mockingboard.pair[pair].ay8913);
	ampl = ay8913ptr->regs[8] | ay8913ptr->regs[9] | ay8913ptr->regs[10];
	if((ampl & 0x10) == 0) {
		// No one uses the envelope
		return;
	}

	env_dsamp = ay8913ptr->env_dsamp;
	env_period = ay8913ptr->regs[11] + (256 * ay8913ptr->regs[12]);
	if(env_period == 0) {
		denv_period = 0.5;		// To match MAME
	} else {
		denv_period = (double)env_period;
	}
	dmul = (1.0 / 16.0) * (1 << 20) * (1 << 20);	// (1ULL << 40) / 16.0
	// Calculate amount counter will count in one sample.
	// inc_per_tick 62.5KHz tick: (1/env_period)
	// inc_per_dcyc: (1/(16*env_period))
	// inc_per_samp = inc_per_dcyc * g_dcycs_per_samp
	dsamp_inc = (dword64)((dmul * g_dcycs_per_samp / denv_period));
			// Amount to inc per sample, fixed point, 40 bit frac

	reg13 = ay8913ptr->regs[13];			// "reg15", env ctrl
	eff_ampl = 0;
	for(i = 0; i < num_samps; i++) {
		env_dsamp = (env_dsamp + dsamp_inc) & 0x9fffffffffffULL;
		env_val = (env_dsamp >> 40) & 0xff;
		eff_ampl = env_val & 0xf;
		if((reg13 & 4) == 0) {
			eff_ampl = 15 - eff_ampl;	// not attack
		}
		if((reg13 & 8) && (reg13 & 2)) {
			// continue and alternate
			if(env_val & 0x10) {
				eff_ampl = 15 - eff_ampl;
			}
		}
		if(((reg13 & 8) == 0) && (env_val >= 0x10)) {
			eff_ampl = 0;
			ampl = 0;		// Turn off envelope
			env_dsamp |= (0x80ULL << 40);
		} else if((reg13 & 1) && (env_val >= 0x10)) {
			eff_ampl = ((reg13 >> 1) ^ (reg13 >> 2)) & 1;
			eff_ampl = eff_ampl * 15;
			ampl = eff_ampl;	// Turn off envelope
			env_dsamp |= (0x80ULL << 40);
		}
		*env_ptr++ = vol_ptr[eff_ampl & 0xf];
	}
	ay8913ptr->env_dsamp = env_dsamp;
}

void
sound_mock_noise(int pair, byte *noise_ptr, int num_samps)
{
	Ay8913	*ay8913ptr;
	word32	ampl, mix, noise_val, noise_samp, noise_period, xor, samp_inc;
	int	doit;
	int	i;

	if(num_samps > MAX_MOCK_ENV_SAMPLES) {
		halt_printf("noise overflow!: %d\n", num_samps);
		return;
	}

	ay8913ptr = &(g_mockingboard.pair[pair].ay8913);
	doit = 0;
	for(i = 0; i < 3; i++) {
		ampl = ay8913ptr->regs[8 + i];
		mix = ay8913ptr->regs[7] >> i;
		if((ampl != 0) && ((mix & 8) == 0)) {
			doit = 1;
			break;
		}
	}
	if(!doit) {
		// No channel looks at noise, don't bother
		return;
	}

	noise_val = ay8913ptr->noise_val;
	noise_samp = ay8913ptr->noise_samp;
	noise_period = (ay8913ptr->regs[6] & 0x1f);
	noise_period = noise_period << 16;
	samp_inc = (word32)(65536 * g_dcycs_per_samp / 16.0);
			// Amount to inc per sample
	if(noise_samp >= noise_period) {
		// Period changed during sound, reset
		noise_samp = noise_period;
	}
	for(i = 0; i < num_samps; i++) {
		noise_samp += samp_inc;
		if(noise_samp >= noise_period) {
			// HACK: handle fraction
			// 17-bit LFSR, algorithm from MAME:sound/ay8910.cpp
			// val = val ^ (((val & 1) ^ ((val >> 3) & 1)) << 17)
			xor = 0;
			xor = (noise_val ^ (noise_val >> 3)) & 1;
			noise_val = (noise_val ^ (xor << 17)) >> 1;
			noise_samp -= noise_period;
		}
		noise_ptr[i] = noise_val & 1;
	}
	ay8913ptr->noise_samp = noise_samp;
	ay8913ptr->noise_val = noise_val;
}

int g_did_mock_print = 100;

void
sound_mock_play(int pair, int channel, int *outptr, int *env_ptr,
				byte *noise_ptr, int *vol_ptr, int num_samps)
{
	Ay8913	*ay8913ptr;
	word32	ampl, mix, tone_samp, tone_period, toggle_tone;
	word32	samp_inc, noise_val;
	int	out, ival, do_print;
	int	i;

	if((g_mockingboard.disable_mask >> ((pair * 3) + channel)) & 1) {
		// This channel is disabled
		return;
	}

	ay8913ptr = &(g_mockingboard.pair[pair].ay8913);
	ampl = ay8913ptr->regs[8 + channel] & 0x1f;
	if(ampl == 0) {
		return;
	}
	toggle_tone = ay8913ptr->toggle_tone[channel];		// 0 or 1
	mix = (ay8913ptr->regs[7] >> channel) & 9;
	if(mix == 9) {
		// constant tone: output will be ampl for this channel.
		if(ampl & 0x10) {		// Envelope!
			// The envelope can make the tone, so must calculate it
		} else {
			// HACK: do nothing for now
			return;
		}
	}
	outptr += pair;			// pair[1] is right
	tone_samp = ay8913ptr->tone_samp[channel];
	tone_period = ay8913ptr->regs[2*channel] +
					(256 * ay8913ptr->regs[2*channel + 1]);
	tone_period = tone_period << 16;
	samp_inc = (word32)(65536 * g_dcycs_per_samp / 8.0);
			// Amount to inc per sample
	do_print = 0;
	if(g_mockingboard.disable_mask) {
		printf("Doing %d samps, mix:%d, ampl:%02x\n", num_samps, mix,
									ampl);
		do_print = 1;
		g_did_mock_print = 0;
	}
	if((num_samps > 500) && (g_did_mock_print == 0)) {
		do_print = 1;
		g_did_mock_print = 1;
		printf("Start of %d sample, channel %d mix:%02x ampl:%02x "
			"toggle_tone:%02x\n", num_samps, channel, mix, ampl,
			toggle_tone);
		printf(" tone_period:%08x, tone_samp:%08x, samp_inc:%08x\n",
			tone_period, tone_samp, samp_inc);
	}
	if(tone_samp >= tone_period) {
		// Period changed during sound, reset it
		tone_samp = tone_period;
	}
	for(i = 0; i < num_samps; i++) {
		tone_samp += samp_inc;
		if(tone_samp >= tone_period) {
			// HACK: handle toggling mid-sample...
			toggle_tone ^= 1;
			if(do_print) {
				printf("i:%d tone_toggled to %d, tone_period:"
					"%04x, pre tone_samp:%08x\n", i,
					toggle_tone, tone_period, tone_samp);
			}
			tone_samp -= tone_period;
			if(do_print) {
				printf("post tone_samp:%08x\n", tone_samp);
			}
		}
		noise_val = noise_ptr[i] & 1;
		out = (toggle_tone || (mix & 1)) &
						((noise_val & 1) || (mix & 8));
			// Careful mix of || and & above...
		ival = vol_ptr[ampl & 0xf];
		if(ampl & 0x10) {			// Envelope
			ival = env_ptr[i];
		}
		*outptr += ival*out;
		outptr += 2;
	}
	ay8913ptr->tone_samp[channel] = tone_samp;
	ay8913ptr->toggle_tone[channel] = toggle_tone;
}


void mock_ay8913_reset(int pair_num, double dcycs) {
  Ay8913 *ay8913ptr;
  int i;

  if (dcycs) {
    // Avoid unused parameter warning
  }
  ay8913ptr = &(g_mockingboard.pair[pair_num].ay8913);
  ay8913ptr->reg_addr_latch = 0;
  for (i = 0; i < 16; i++) {
    ay8913ptr->regs[i] = 0;
  }
  for (i = 0; i < 3; i++) {
    ay8913ptr->toggle_tone[i] = 0;
    ay8913ptr->tone_samp[i] = 0;
  }
  ay8913ptr->noise_val = 0x12345678;
  ay8913ptr->noise_samp = 0;
  ay8913ptr->env_dsamp = 0;
}

void mockingboard_reset(double dcycs) {
  word32 timer1_latch;

  timer1_latch = g_mockingboard.pair[0].mos6522.timer1_latch;
  memset(&g_mockingboard, 0, sizeof(g_mockingboard));

  g_mockingboard_last_int_dcycs = (dword64)dcycs;
  if (g_mockingboard_event_int_dcycs != 0) {
    (void)remove_event_mockingboard();
  }
  g_mockingboard_event_int_dcycs = 0;
  printf("At reset, timer1_latch: %08x\n", timer1_latch);
  if ((timer1_latch & 0xffff) == 0x234) { // MB-audit
    g_mockingboard.pair[0].mos6522.timer1_latch = timer1_latch;
  } else {
    g_mockingboard.pair[0].mos6522.timer1_latch = 0xff00;
  }
  g_mockingboard.pair[0].mos6522.timer1_counter = 0x2ff00;
  g_mockingboard.pair[0].mos6522.timer2_counter = 0x2ff00;
  g_mockingboard.pair[1].mos6522.timer1_latch = 0xff00;
  g_mockingboard.pair[1].mos6522.timer1_counter = 0x2ff00;
  g_mockingboard.pair[1].mos6522.timer2_counter = 0x2ff00;

  mock_ay8913_reset(0, dcycs);
  mock_ay8913_reset(1, dcycs);
}

void mock_show_pair(int pair_num, double dcycs, const char *str) {
  Mos6522 *mos6522ptr;

  mos6522ptr = &(g_mockingboard.pair[pair_num].mos6522);
  printf("Mock %d %s, t1_lat:%05x, t1_c:%05x, t2_l:%05x t2_c:%05x, ifr:"
         "%02x, acr:%02x, ier:%02x\n",
         pair_num, str, mos6522ptr->timer1_latch, mos6522ptr->timer1_counter,
         mos6522ptr->timer2_latch, mos6522ptr->timer2_counter, mos6522ptr->ifr,
         mos6522ptr->acr, mos6522ptr->ier);
  printf("  dcycs:%lf, event_int:%lld\n", dcycs,
         g_mockingboard_event_int_dcycs);
}

// Timers work as follows: if written with '10' in cycle 0, on cycle 1 the
//  counters loads with '10', and counts down to 9 on cycle 2...down to 1
//  on cycle 10, then 0 on cycle 11, and then signals an interrupt on cycle 12
// To handle this, timers are "read value + 1" so that timer==0 means the
//  timer should interrupt right now.  So to write '10' to a timer, the code
//  needs to write 10+2=12 to the variable, so that on the next cycle, it is
//  11, which will be returned as 10 to software reading the reg.

void mock_update_timers(int doit, double dcycs) {
  Mos6522 *mos6522ptr;
  dword64 dcycs_int, ddiff, dleft, timer1_int_dcycs, timer2_int_dcycs;
  dword64 closest_int_dcycs, event_int_dcycs;
  word32 timer_val, ier, timer_eff, timer_latch, log_ifr;
  int i;


  dcycs_int = (dword64)dcycs;
  ddiff = dcycs_int - g_mockingboard_last_int_dcycs;
  if (!doit && (dcycs_int <= g_mockingboard_last_int_dcycs)) {
    return; // Nothing more to do
  }

  // printf("mock_update_timers at %lf %016llx, ddiff:%llx\n", dcycs,
  //						dcycs_int, ddiff);
  // Update timers by ddiff integer cycles, calculate next event time
  g_mockingboard_last_int_dcycs = dcycs_int;
  closest_int_dcycs = 0;
  for (i = 0; i < 2; i++) { // pair_num
    mos6522ptr = &(g_mockingboard.pair[i].mos6522);
    timer1_int_dcycs = 0;
    timer_val = mos6522ptr->timer1_counter;
    ier = mos6522ptr->ier;
    timer_eff = (timer_val & 0x1ffff);
    dleft = ddiff;
    timer_latch = mos6522ptr->timer1_latch + 2;
    if (dleft < timer_eff) {
      // Move ahead only a little, no triggering
      timer_val = timer_val - (word32)dleft;
      if (ddiff) {
        // printf("New timer1_val:%05x, dleft:%08llx\n",
        //		timer_val, dleft);
      }
      if (((mos6522ptr->ifr & 0x40) == 0) && (ier & 0x40)) {
        // IFR not set yet, prepare an event
        timer1_int_dcycs = dcycs_int + (timer_val & 0x1ffff);
        // printf("t1_int_dcycs: %016llx\n",
        //			timer1_int_dcycs);
      }
    } else {
      // Timer1 has triggered now (maybe rolled over more
      //  than once).
      log_ifr = 0;
      if ((timer_val & 0x20000) == 0) {
        // Either free-running, or not one-shot already
        //  triggered
        // Set interrupt: Ensure IFR | 0x40 is set
        mos6522ptr->ifr =
            mock_6522_new_ifr(dcycs, i, mos6522ptr->ifr | 0x40, ier);
        log_ifr = 1;
      }
      dleft -= timer_eff;
      if (dleft >= timer_latch) {
        // It's rolled over several times, remove those
        dleft = dleft % timer_latch;
      }
      if (dleft == 0) {
        dleft = timer_latch;
      }
      timer_val = (timer_latch - dleft) & 0x1ffff;
      if ((mos6522ptr->acr & 0x40) == 0) {
        // ACR6=0: One shot mode, mark it as triggered
        timer_val |= 0x20000;
      }
    }

#if 0
		printf("%lf ch%d timer1 was %05x, now %05x\n", dcycs, i,
					mos6522ptr->timer1_counter, timer_val);
#endif

    mos6522ptr->timer1_counter = timer_val;

    // Handle timer2
    timer2_int_dcycs = 0;
    timer_val = mos6522ptr->timer2_counter;
    timer_eff = timer_val & 0x1ffff;
    dleft = ddiff;
    if (mos6522ptr->acr & 0x20) {
      // Count pulses mode.  Just don't count
      dleft = 0;
    }
    if (dleft < timer_eff) {
      // Move ahead only a little, no triggering
      timer_val = timer_val - (word32)dleft;
      if (((mos6522ptr->ifr & 0x20) == 0) && (ier & 0x20)) {
        // IFR not set yet, prepare an event
        timer2_int_dcycs = dcycs_int + (timer_val & 0x1ffff);
        // printf("t2_int_dcycs: %016llx\n",
        //			timer1_int_dcycs);
      }
    } else if (timer_val & 0x20000) {
      // And already triggered once, just update count
      timer_val = ((timer_eff - dleft) & 0xffff) | 0x20000;
    } else {
      // Has not triggered once yet, but it will now
      mos6522ptr->ifr =
          mock_6522_new_ifr(dcycs, i, mos6522ptr->ifr | 0x20, ier);
      timer_val = ((timer_val - dleft) & 0xffff) | 0x20000;
    }

    // printf("ch%d timer2 was %05x, now %05x\n", i,
    //			mos6522ptr->timer2_counter, timer_val);

    mos6522ptr->timer2_counter = timer_val;

    if (timer1_int_dcycs && timer2_int_dcycs) {
      timer1_int_dcycs = MY_MIN(timer1_int_dcycs, timer2_int_dcycs);
    }
    if (timer1_int_dcycs) {
      if (closest_int_dcycs) {
        closest_int_dcycs = MY_MIN(closest_int_dcycs, timer1_int_dcycs);
      } else {
        closest_int_dcycs = timer1_int_dcycs;
      }
    }
  }

  event_int_dcycs = g_mockingboard_event_int_dcycs;
  if (closest_int_dcycs) {
    // See if this is sooner than the current pending event
    // printf("closest_int_dcycs: %016llx, event_int:%016llx\n",
    //			closest_int_dcycs, event_int_dcycs);
    doit = 0;
    if (event_int_dcycs && (closest_int_dcycs < event_int_dcycs)) {
      // There was a pending event.  Discard it
      // printf("Call remove_event_mockingboard\n");
      remove_event_mockingboard();
      doit = 1;
    }
    if (!event_int_dcycs || doit) {
      // printf("Call add_event_mockingboard: %016llx %lld\n",
      //	closest_int_dcycs, closest_int_dcycs);
      add_event_mockingboard((double)closest_int_dcycs);
      g_mockingboard_event_int_dcycs = closest_int_dcycs;
    }
  }
}

void mockingboard_event(double dcycs) {
  // Received an event--we believe we may need to set an IRQ now.
  // Event was already removed from the event queue
  // printf("Mockingboard_event received at %lf\n", dcycs);
  g_mockingboard_event_int_dcycs = 0;
  mock_update_timers(1, dcycs);
}

word32 mockingboard_read(word32 loc, double dcycs) {
  int pair_num;

  // printf("mockingboard read: %04x\n", loc);
  pair_num = (loc >> 7) & 1; // 0 or 1
  return mock_6522_read(pair_num, loc & 0xf, dcycs);
}

void mockingboard_write(word32 loc, word32 val, double dcycs) {
  int pair_num;

  // printf("mockingboard write: %04x=%02x\n", loc, val);
  pair_num = (loc >> 7) & 1; // 0 or 1
  mock_6522_write(pair_num, loc & 0xf, val, dcycs);
}

word32 mock_6522_read(int pair_num, word32 loc, double dcycs) {
  Mos6522 *mos6522ptr;
  word32 val;

  // Read from 6522 #pair_num loc (0-15)
  mos6522ptr = &(g_mockingboard.pair[pair_num].mos6522);
  val = 0;
  switch (loc) {
  case 0x0: // ORB/IRB
    // Connected to AY8913 { RESET, BDIR, BC1 }
    val = mos6522ptr->orb;
    // There are no outputs from AY8913 to the 6522 B Port
    break;
  case 0x1: // ORA
  case 0xf: // ORA, no handshake
    val = mos6522ptr->ora;
    break;
  case 0x2: // DDRB
    val = mos6522ptr->ddrb;
    break;
  case 0x3: // DDRA
    val = mos6522ptr->ddra;
    break;
  case 0x4: // T1C-L (timer[0])
    mock_update_timers(1, dcycs);
    val = (mos6522ptr->timer1_counter - 1) & 0xff;
    mos6522ptr->ifr = mock_6522_new_ifr(
        dcycs, pair_num, mos6522ptr->ifr & (~0x40), mos6522ptr->ier);
    // Clear Bit 6
    mock_update_timers(1, dcycs); // Prepare another int (maybe)
    break;
  case 0x5: // T1C-H
    mock_update_timers(1, dcycs);
    val = ((mos6522ptr->timer1_counter - 1) >> 8) & 0xff;
    break;
  case 0x6: // T1L-L
    val = mos6522ptr->timer1_latch & 0xff;
    break;
  case 0x7: // T1L-H
    val = (mos6522ptr->timer1_latch >> 8) & 0xff;
    break;
  case 0x8: // T2C-L
    mock_update_timers(1, dcycs);
    val = (mos6522ptr->timer2_counter - 1) & 0xff;
    mos6522ptr->ifr = mock_6522_new_ifr(
        dcycs, pair_num, mos6522ptr->ifr & (~0x20), mos6522ptr->ier);
    // Clear Bit 5
    mock_update_timers(1, dcycs); // Prepare another int (maybe)
    break;
  case 0x9: // T2C-H
    mock_update_timers(1, dcycs);
    val = ((mos6522ptr->timer2_counter - 1) >> 8) & 0xff;
    break;
  case 0xa: // SR
    val = mos6522ptr->sr;
    halt_printf("Reading SR %d %02x\n", pair_num, val);
    break;
  case 0xb: // ACR
    val = mos6522ptr->acr;
    break;
  case 0xc: // PCR
    val = mos6522ptr->pcr;
    break;
  case 0xd: // IFR
    mock_update_timers(1, dcycs);
    val = mos6522ptr->ifr;
    break;
  case 0xe:                       // IER
    val = mos6522ptr->ier | 0x80; // Despite MOS6522
    break;                        //  datasheet, bit 7 = 1
  }
  // printf("6522 %d loc:%x ret:%02x\n", pair_num, loc, val);
  return val;
}

void mock_6522_write(int pair_num, word32 loc, word32 val, double dcycs) {
  Mos6522 *mos6522ptr;
  word32 ora, orb, new_val, mask;

  // Write to 6522 #num6522 loc (0-15)

  // printf("6522 %d loc:%x write:%02x\n", pair_num, loc, val);

  mos6522ptr = &(g_mockingboard.pair[pair_num].mos6522);
  switch (loc) {
  case 0x0: // ORB
    mask = mos6522ptr->ddrb;
    orb = mos6522ptr->orb;
    new_val = (val & mask) | (orb & (~mask));
    if (orb != new_val) {
      mock_ay8913_control_update(pair_num, new_val, orb, dcycs);
    }
    mos6522ptr->orb = new_val;
    break;
  case 0x1: // ORA
  case 0xf: // ORA, no handshake
    mask = mos6522ptr->ddra;
    ora = mos6522ptr->ora;
    new_val = (val & mask) | (ora & (~mask));
    mos6522ptr->ora = new_val;
    break;
  case 0x2: // DDRB
    orb = mos6522ptr->orb;
    new_val = (orb & val) | (orb & (~val));
    if (orb != new_val) {
      mock_ay8913_control_update(pair_num, new_val, orb, dcycs);
    }
    mos6522ptr->orb = new_val;
    mos6522ptr->ddrb = val;
    return;
  case 0x3: // DDRA
    ora = mos6522ptr->ora;
    mos6522ptr->ora = (ora & val) | (ora & (~val));
    mos6522ptr->ddra = val;
    return;
  case 0x4: // T1C-L
    mock_update_timers(0, dcycs);
    mos6522ptr->timer1_latch = (mos6522ptr->timer1_latch & 0x1ff00) | val;
    // printf("Set T1C-L, timer1_latch=%05x\n",
    //				mos6522ptr->timer1_latch);
    break;
  case 0x5: // T1C-H
    mock_update_timers(1, dcycs);
    val = (mos6522ptr->timer1_latch & 0xff) | (val << 8);
    mos6522ptr->timer1_latch = val;
    mos6522ptr->timer1_counter = val + 2;
    // The actual timer1_counter update happens next cycle,
    //  so we want val+1, plus another 1
    // printf("Set T1C-H, timer1_latch=%05x\n",
    //				mos6522ptr->timer1_latch);
    mos6522ptr->ifr = mock_6522_new_ifr(
        dcycs, pair_num, mos6522ptr->ifr & (~0x40), mos6522ptr->ier);
    // Clear Bit 6
    mock_update_timers(1, dcycs);
    break;
  case 0x6: // T1L-L
    mock_update_timers(0, dcycs);
    mos6522ptr->timer1_latch = (mos6522ptr->timer1_latch & 0x1ff00) | val;
    break;
  case 0x7: // T1L-H
    mock_update_timers(1, dcycs);
    val = (mos6522ptr->timer1_latch & 0xff) | (val << 8);
    mos6522ptr->timer1_latch = val;
    mos6522ptr->ifr = mock_6522_new_ifr(
        dcycs, pair_num, mos6522ptr->ifr & (~0x40), mos6522ptr->ier);
    // Clear Bit 6
    mock_update_timers(1, dcycs);
    // mock_show_pair(pair_num, dcycs, "Wrote T1L-H");
    break;
  case 0x8: // T2C-L
    mos6522ptr->timer2_latch = (mos6522ptr->timer2_latch & 0xff00) | val;
    break;
  case 0x9: // T2C-H
    mock_update_timers(1, dcycs);
    val = (mos6522ptr->timer2_latch & 0xff) | (val << 8);
    mos6522ptr->timer2_latch = val;
    mos6522ptr->timer2_counter = val + 2;
    mos6522ptr->ifr = mock_6522_new_ifr(
        dcycs, pair_num, mos6522ptr->ifr & (~0x20), mos6522ptr->ier);
    // Clear bit 5
    mock_update_timers(1, dcycs);
    break;
  case 0xa: // SR
    mos6522ptr->sr = val;
    halt_printf("Wrote SR reg: %d %02x\n", pair_num, val);
    break;
  case 0xb: // ACR
    mock_update_timers(0, dcycs);
    mos6522ptr->acr = val;
    mock_update_timers(1, dcycs);
    break;
  case 0xc: // PCR
    mos6522ptr->pcr = val;
    break;
  case 0xd: // IFR
    mock_update_timers(1, dcycs);
    mos6522ptr->ifr = mock_6522_new_ifr(
        dcycs, pair_num, mos6522ptr->ifr & (~val), mos6522ptr->ier);
    mock_update_timers(1, dcycs);
    break;
  case 0xe: // IER
    // Recalculate effective IFR with new IER
    mock_update_timers(1, dcycs);
    if (val & 0x80) { // Set EIR bits
      val = mos6522ptr->ier | val;
    } else { // Clear EIR bits
      val = mos6522ptr->ier & (~val);
    }
    val = val & 0x7f;
    mos6522ptr->ier = val;
    mos6522ptr->ifr = mock_6522_new_ifr(dcycs, pair_num, mos6522ptr->ifr, val);
    mock_update_timers(1, dcycs);
    // mock_show_pair(pair_num, dcycs, "Wrote IER");
    break;
  }
}

word32 mock_6522_new_ifr(double dcycs, int pair_num, word32 ifr, word32 ier) {
  word32 irq_mask;

  // Determine if there are any interrupts pending now
  irq_mask = IRQ_PENDING_MOCKINGBOARDA << pair_num;
  if ((ifr & ier & 0x7f) == 0) {
    // No IRQ pending anymore
    ifr = ifr & 0x7f; // Clear bit 7
    if (g_irq_pending & irq_mask) {
      // printf("MOCK INT OFF\n");
    }
    remove_irq(irq_mask);
  } else {
    // IRQ is pending
    ifr = ifr | 0x80; // Set bit 7
    if (!(g_irq_pending & irq_mask)) {
      // printf("MOCK INT ON\n");
    }
    add_irq(irq_mask);
  }
  return ifr;
}

void mock_ay8913_reg_read(int pair_num) {
  Mos6522 *mos6522ptr;
  Ay8913 *ay8913ptr;
  word32 reg_addr_latch, mask, val, ora;

  mos6522ptr = &(g_mockingboard.pair[pair_num].mos6522);
  ay8913ptr = &(g_mockingboard.pair[pair_num].ay8913);
  reg_addr_latch = ay8913ptr->reg_addr_latch;
  val = 0;
  if (reg_addr_latch < 16) {
    val = ay8913ptr->regs[reg_addr_latch];
  }
  // ORA at 6522 is merge of ORA using DDRA
  mask = mos6522ptr->ddra;
  ora = mos6522ptr->ora;
  mos6522ptr->ora = (ora & mask) | (val & (~mask));
}

word32 g_mock_channel_regs[3] = {
    0x39c3, // channel A: regs 0,1,6,7,8,11,12,13
    0x3acc, // channel B: regs 2,3,6,7,9,11,12,13
    0x3cf0  // channel C: regs 4,5,6,7,10,11,12,13
};

void mock_ay8913_reg_write(int pair_num, double dcycs) {
  Mos6522 *mos6522ptr;
  Ay8913 *ay8913ptr;
  double dsamps;
  word32 reg_addr_latch, ora, mask, rmask, do_print;
  int i;

  mos6522ptr = &(g_mockingboard.pair[pair_num].mos6522);
  ay8913ptr = &(g_mockingboard.pair[pair_num].ay8913);
  reg_addr_latch = ay8913ptr->reg_addr_latch;
  ora = mos6522ptr->ora;
  dsamps = dcycs * g_dsamps_per_dcyc;
  if (reg_addr_latch < 16) {
    mask = (g_mockingboard.disable_mask >> (3 * pair_num)) & 7;
    rmask = 0;
    do_print = 0;
    for (i = 0; i < 3; i++) {
      if (((mask >> i) & 1) == 0) {
        rmask |= g_mock_channel_regs[i];
      }
    }
    do_print = (rmask >> reg_addr_latch) & 1;
    if ((ora != ay8913ptr->regs[reg_addr_latch]) || (reg_addr_latch == 13)) {
      // New value, or writing to Envelope control
      do_print = 0;
      if (do_print) {
        printf("%.2lf %.2lf mock pair%d reg[%d]=%02x. "
               "[2,3]=%02x_%02x [67]=%02x,%02x, [9]="
               "%02x, [12,11]=%02x_%02x [13]=%02x\n",
               dsamps, dcycs, pair_num, reg_addr_latch, ora, ay8913ptr->regs[3],
               ay8913ptr->regs[2], ay8913ptr->regs[6], ay8913ptr->regs[7],
               ay8913ptr->regs[9], ay8913ptr->regs[12], ay8913ptr->regs[11],
               ay8913ptr->regs[13]);
      }
      sound_play(dsamps);
    }
    ay8913ptr->regs[reg_addr_latch] = ora;
    if (reg_addr_latch == 13) { // Envelope control
      ay8913ptr->env_dsamp &= 0x1fffffffffffULL;
      // Clear "hold" in (env_val & (0x80 << 40))
    }
  }
}

void mock_ay8913_control_update(int pair_num, word32 new_val, word32 prev_val,
                                double dcycs) {
  Mos6522 *mos6522ptr;
  Ay8913 *ay8913ptr;

  mos6522ptr = &(g_mockingboard.pair[pair_num].mos6522);
  ay8913ptr = &(g_mockingboard.pair[pair_num].ay8913);
  // printf("ay8913 %d control now:%02x\n", pair_num, new_val);

  // new_val and prev_val are { reset_l, BDIR, BC1 }
  // 4=Idle; 5=Read; 6=Write; 7=Latch_addr
  // Latch new address and write data at the time the ctl changes to Idle
  // Do read as soon as the ctl indicates to do a read.

  if ((new_val & 4) == 0) {
    mock_ay8913_reset(pair_num, dcycs);
    return;
  }
  new_val = new_val & 7;
  prev_val = prev_val & 7;
  if (prev_val == 7) { // Latch new address, latch it now
    ay8913ptr->reg_addr_latch = mos6522ptr->ora;
  } else if (prev_val == 6) { // Write data, do it now
    mock_ay8913_reg_write(pair_num, dcycs);
  }
  if (new_val == 5) {
    mock_ay8913_reg_read(pair_num);
  }
}

void mockingboard_show(int got_num, word32 disable_mask) {
  int i, j;

  if (got_num) {
    g_mockingboard.disable_mask = disable_mask;
  }
  printf("g_mockingboard.disable_mask:%02x\n", g_mockingboard.disable_mask);
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 14; j++) {
      printf("Mockingboard pair[%d].reg[%d]=%02x\n", i, j,
             g_mockingboard.pair[i].ay8913.regs[j]);
    }
  }
}
