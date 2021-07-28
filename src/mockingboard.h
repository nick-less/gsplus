#ifndef MOCKINGBOARD_H
#define MOCKINGBOARD_H


#define MY_MIN(a,b)	(((a) < (b)) ? (a) : (b))

typedef unsigned long long dword64;

#define EV_MOCKINGBOARD         8

#define IRQ_PENDING_MOCKINGBOARDA	0x10000
#define IRQ_PENDING_MOCKINGBOARDB	0x20000		/* must be BOARDA*2 */

// Mockingboard contains two pairs.  Each pair is a 6522 interfacing
//  to an AY-8913 to generate sounds.  Eacho AY-8913 contains 3 channels of
//  sound.  Model each pair separately.

STRUCT(Mos6522) {
	byte	orb;
	byte	ora;
	byte	ddrb;
	byte	ddra;
	word32	timer1_latch;
	word32	timer1_counter;
	word32	timer2_latch;
	word32	timer2_counter;
	byte	sr;
	byte	acr;
	byte	pcr;
	byte	ifr;
	byte	ier;
};

STRUCT(Ay8913) {
	byte	regs[16];
	byte	reg_addr_latch;
	byte	toggle_tone[3];		// Channel A,B,C: 0 = low, 1 = high
	word32	tone_samp[3];
	word32	noise_val;
	word32	noise_samp;
	dword64	env_dsamp;
};

STRUCT(Mock_pair) {
	Mos6522	mos6522;
	Ay8913	ay8913;
};

STRUCT(Mockingboard) {
	Mock_pair pair[2];
	word32	disable_mask;
};

void add_event_mockingboard(double dcycs);
void remove_event_mockingboard(void);


void mock_ay8913_reset(int pair_num, double dcycs);
void mockingboard_reset(double dcycs);
void mock_show_pair(int pair_num, double dcycs, const char *str);
void mock_update_timers(int doit, double dcycs);
void mockingboard_event(double dcycs);
word32 mockingboard_read(word32 loc, double dcycs);
void mockingboard_write(word32 loc, word32 val, double dcycs);
word32 mock_6522_read(int pair_num, word32 loc, double dcycs);
void mock_6522_write(int pair_num, word32 loc, word32 val, double dcycs);
word32 mock_6522_new_ifr(double dcycs, int pair_num, word32 ifr, word32 ier);
void mock_ay8913_reg_read(int pair_num);
void mock_ay8913_reg_write(int pair_num, double dcycs);
void mock_ay8913_control_update(int pair_num, word32 new_val, word32 prev_val, double dcycs);
void mockingboard_show(int got_num, word32 disable_mask);


#endif