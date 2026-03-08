#ifdef _MSC_VER
#	ifndef _CRT_SECURE_NO_WARNINGS
#		define _CRT_SECURE_NO_WARNINGS
#	endif
#elif defined __linux__ && !defined _GNU_SOURCE
#	define _GNU_SOURCE
#	include<stddef.h>//ptrdiff_t
#endif
#define _USE_MATH_DEFINES
#include<math.h>
#include<sys/stat.h>


	#define PROFILE_TIME		//should be on

#ifdef _MSC_VER
	#define PROFILE_SIZE
	#define LOUD			//size & time

	#define EXTRA_CHECKS
//	#define ESTIMATE_SIZE		//DEBUG		checks for zero frequency, visualizes context usage

//	#define TEST_INTERLEAVE
	#define ENABLE_GUIDE		//DEBUG		checks interleaved pixels
//	#define ANS_VAL			//DEBUG
#endif

	#define ENABLE_RCT_EXTENSION
	#define INTERLEAVESIMD


enum
{
	XCODERS=4,
	YCODERS=4,
	NLANES=XCODERS*YCODERS,

	ANALYSIS_XSTRIDE=4,
	ANALYSIS_YSTRIDE=4,

	DEFAULT_EFFORT_LEVEL=1,
	L1_NPREDS1=5,
	L1_NPREDS2=8,
	L1_NPREDS3=20,
	L1_SH1=16,
	L1_SH2=17,
	L1_SH3=19,

	GRBITS=3,
	NCTX=18,	//18*3+3 = 57 total

	PROBBITS=12,	//12 bit max	James Bonfield's CDF2sym: {freq<<20 | bias<<8 | sym}
	RANS_STATE_BITS=31,
	RANS_RENORM_BITS=16,

	XPAD=8,
	NROWS=4,
	NVAL=2*3*NLANES,//pY pU pV eY eU eV    for easy iteration, unlike avx2
};

#define SOFT_LG2
#define COMMON_rANS
#include"common.h"

INLINE void dec_yuv(
	uint32_t *mstate,
	const uint16_t *ctx,
	const uint32_t *CDF2syms,
	uint8_t **pstreamptr,
	const uint8_t *streamend,
	uint16_t *syms
)
{
	const uint8_t *streamptr=*pstreamptr;
	uint32_t decctx[NLANES];

	decctx[0x0]=ctx[0x0]<<PROBBITS|(mstate[0x0]&((1<<PROBBITS)-1));
	decctx[0x1]=ctx[0x1]<<PROBBITS|(mstate[0x1]&((1<<PROBBITS)-1));
	decctx[0x2]=ctx[0x2]<<PROBBITS|(mstate[0x2]&((1<<PROBBITS)-1));
	decctx[0x3]=ctx[0x3]<<PROBBITS|(mstate[0x3]&((1<<PROBBITS)-1));
	decctx[0x4]=ctx[0x4]<<PROBBITS|(mstate[0x4]&((1<<PROBBITS)-1));
	decctx[0x5]=ctx[0x5]<<PROBBITS|(mstate[0x5]&((1<<PROBBITS)-1));
	decctx[0x6]=ctx[0x6]<<PROBBITS|(mstate[0x6]&((1<<PROBBITS)-1));
	decctx[0x7]=ctx[0x7]<<PROBBITS|(mstate[0x7]&((1<<PROBBITS)-1));
	decctx[0x8]=ctx[0x8]<<PROBBITS|(mstate[0x8]&((1<<PROBBITS)-1));
	decctx[0x9]=ctx[0x9]<<PROBBITS|(mstate[0x9]&((1<<PROBBITS)-1));
	decctx[0xA]=ctx[0xA]<<PROBBITS|(mstate[0xA]&((1<<PROBBITS)-1));
	decctx[0xB]=ctx[0xB]<<PROBBITS|(mstate[0xB]&((1<<PROBBITS)-1));
	decctx[0xC]=ctx[0xC]<<PROBBITS|(mstate[0xC]&((1<<PROBBITS)-1));
	decctx[0xD]=ctx[0xD]<<PROBBITS|(mstate[0xD]&((1<<PROBBITS)-1));
	decctx[0xE]=ctx[0xE]<<PROBBITS|(mstate[0xE]&((1<<PROBBITS)-1));
	decctx[0xF]=ctx[0xF]<<PROBBITS|(mstate[0xF]&((1<<PROBBITS)-1));
	decctx[0x0]=CDF2syms[decctx[0x0]];
	decctx[0x1]=CDF2syms[decctx[0x1]];
	decctx[0x2]=CDF2syms[decctx[0x2]];
	decctx[0x3]=CDF2syms[decctx[0x3]];
	decctx[0x4]=CDF2syms[decctx[0x4]];
	decctx[0x5]=CDF2syms[decctx[0x5]];
	decctx[0x6]=CDF2syms[decctx[0x6]];
	decctx[0x7]=CDF2syms[decctx[0x7]];
	decctx[0x8]=CDF2syms[decctx[0x8]];
	decctx[0x9]=CDF2syms[decctx[0x9]];
	decctx[0xA]=CDF2syms[decctx[0xA]];
	decctx[0xB]=CDF2syms[decctx[0xB]];
	decctx[0xC]=CDF2syms[decctx[0xC]];
	decctx[0xD]=CDF2syms[decctx[0xD]];
	decctx[0xE]=CDF2syms[decctx[0xE]];
	decctx[0xF]=CDF2syms[decctx[0xF]];
	syms[0x0]=(uint8_t)decctx[0x0];
	syms[0x1]=(uint8_t)decctx[0x1];
	syms[0x2]=(uint8_t)decctx[0x2];
	syms[0x3]=(uint8_t)decctx[0x3];
	syms[0x4]=(uint8_t)decctx[0x4];
	syms[0x5]=(uint8_t)decctx[0x5];
	syms[0x6]=(uint8_t)decctx[0x6];
	syms[0x7]=(uint8_t)decctx[0x7];
	syms[0x8]=(uint8_t)decctx[0x8];
	syms[0x9]=(uint8_t)decctx[0x9];
	syms[0xA]=(uint8_t)decctx[0xA];
	syms[0xB]=(uint8_t)decctx[0xB];
	syms[0xC]=(uint8_t)decctx[0xC];
	syms[0xD]=(uint8_t)decctx[0xD];
	syms[0xE]=(uint8_t)decctx[0xE];
	syms[0xF]=(uint8_t)decctx[0xF];
#ifdef ANS_VAL
	ansval_check(mstate, sizeof(uint32_t), NLANES);
#endif

	//update		state = (state>>12)*freq+(rem-cdf)	rem-cdf is prebaked
	mstate[0x0]=(mstate[0x0]>>PROBBITS)*(decctx[0x0]>>(PROBBITS+8))+(decctx[0x0]<<PROBBITS>>(32-PROBBITS));
	mstate[0x1]=(mstate[0x1]>>PROBBITS)*(decctx[0x1]>>(PROBBITS+8))+(decctx[0x1]<<PROBBITS>>(32-PROBBITS));
	mstate[0x2]=(mstate[0x2]>>PROBBITS)*(decctx[0x2]>>(PROBBITS+8))+(decctx[0x2]<<PROBBITS>>(32-PROBBITS));
	mstate[0x3]=(mstate[0x3]>>PROBBITS)*(decctx[0x3]>>(PROBBITS+8))+(decctx[0x3]<<PROBBITS>>(32-PROBBITS));
	mstate[0x4]=(mstate[0x4]>>PROBBITS)*(decctx[0x4]>>(PROBBITS+8))+(decctx[0x4]<<PROBBITS>>(32-PROBBITS));
	mstate[0x5]=(mstate[0x5]>>PROBBITS)*(decctx[0x5]>>(PROBBITS+8))+(decctx[0x5]<<PROBBITS>>(32-PROBBITS));
	mstate[0x6]=(mstate[0x6]>>PROBBITS)*(decctx[0x6]>>(PROBBITS+8))+(decctx[0x6]<<PROBBITS>>(32-PROBBITS));
	mstate[0x7]=(mstate[0x7]>>PROBBITS)*(decctx[0x7]>>(PROBBITS+8))+(decctx[0x7]<<PROBBITS>>(32-PROBBITS));
	mstate[0x8]=(mstate[0x8]>>PROBBITS)*(decctx[0x8]>>(PROBBITS+8))+(decctx[0x8]<<PROBBITS>>(32-PROBBITS));
	mstate[0x9]=(mstate[0x9]>>PROBBITS)*(decctx[0x9]>>(PROBBITS+8))+(decctx[0x9]<<PROBBITS>>(32-PROBBITS));
	mstate[0xA]=(mstate[0xA]>>PROBBITS)*(decctx[0xA]>>(PROBBITS+8))+(decctx[0xA]<<PROBBITS>>(32-PROBBITS));
	mstate[0xB]=(mstate[0xB]>>PROBBITS)*(decctx[0xB]>>(PROBBITS+8))+(decctx[0xB]<<PROBBITS>>(32-PROBBITS));
	mstate[0xC]=(mstate[0xC]>>PROBBITS)*(decctx[0xC]>>(PROBBITS+8))+(decctx[0xC]<<PROBBITS>>(32-PROBBITS));
	mstate[0xD]=(mstate[0xD]>>PROBBITS)*(decctx[0xD]>>(PROBBITS+8))+(decctx[0xD]<<PROBBITS>>(32-PROBBITS));
	mstate[0xE]=(mstate[0xE]>>PROBBITS)*(decctx[0xE]>>(PROBBITS+8))+(decctx[0xE]<<PROBBITS>>(32-PROBBITS));
	mstate[0xF]=(mstate[0xF]>>PROBBITS)*(decctx[0xF]>>(PROBBITS+8))+(decctx[0xF]<<PROBBITS>>(32-PROBBITS));
#ifdef ANS_VAL
	ansval_check(mstate, sizeof(uint32_t), NLANES);
#endif
	{
		const uint8_t *p2;
		uint32_t s2;

		//renorm	if(state<(1<<(31-16)))state=state<<16|read16();
		s2=mstate[0x0]<<16|*(uint16_t*)streamptr; p2=streamptr+sizeof(uint16_t); if(mstate[0x0]<1<<(RANS_STATE_BITS-RANS_RENORM_BITS))mstate[0x0]=s2, streamptr=p2;
		s2=mstate[0x1]<<16|*(uint16_t*)streamptr; p2=streamptr+sizeof(uint16_t); if(mstate[0x1]<1<<(RANS_STATE_BITS-RANS_RENORM_BITS))mstate[0x1]=s2, streamptr=p2;
		s2=mstate[0x2]<<16|*(uint16_t*)streamptr; p2=streamptr+sizeof(uint16_t); if(mstate[0x2]<1<<(RANS_STATE_BITS-RANS_RENORM_BITS))mstate[0x2]=s2, streamptr=p2;
		s2=mstate[0x3]<<16|*(uint16_t*)streamptr; p2=streamptr+sizeof(uint16_t); if(mstate[0x3]<1<<(RANS_STATE_BITS-RANS_RENORM_BITS))mstate[0x3]=s2, streamptr=p2;
		s2=mstate[0x4]<<16|*(uint16_t*)streamptr; p2=streamptr+sizeof(uint16_t); if(mstate[0x4]<1<<(RANS_STATE_BITS-RANS_RENORM_BITS))mstate[0x4]=s2, streamptr=p2;
		s2=mstate[0x5]<<16|*(uint16_t*)streamptr; p2=streamptr+sizeof(uint16_t); if(mstate[0x5]<1<<(RANS_STATE_BITS-RANS_RENORM_BITS))mstate[0x5]=s2, streamptr=p2;
		s2=mstate[0x6]<<16|*(uint16_t*)streamptr; p2=streamptr+sizeof(uint16_t); if(mstate[0x6]<1<<(RANS_STATE_BITS-RANS_RENORM_BITS))mstate[0x6]=s2, streamptr=p2;
		s2=mstate[0x7]<<16|*(uint16_t*)streamptr; p2=streamptr+sizeof(uint16_t); if(mstate[0x7]<1<<(RANS_STATE_BITS-RANS_RENORM_BITS))mstate[0x7]=s2, streamptr=p2;
		s2=mstate[0x8]<<16|*(uint16_t*)streamptr; p2=streamptr+sizeof(uint16_t); if(mstate[0x8]<1<<(RANS_STATE_BITS-RANS_RENORM_BITS))mstate[0x8]=s2, streamptr=p2;
		s2=mstate[0x9]<<16|*(uint16_t*)streamptr; p2=streamptr+sizeof(uint16_t); if(mstate[0x9]<1<<(RANS_STATE_BITS-RANS_RENORM_BITS))mstate[0x9]=s2, streamptr=p2;
		s2=mstate[0xA]<<16|*(uint16_t*)streamptr; p2=streamptr+sizeof(uint16_t); if(mstate[0xA]<1<<(RANS_STATE_BITS-RANS_RENORM_BITS))mstate[0xA]=s2, streamptr=p2;
		s2=mstate[0xB]<<16|*(uint16_t*)streamptr; p2=streamptr+sizeof(uint16_t); if(mstate[0xB]<1<<(RANS_STATE_BITS-RANS_RENORM_BITS))mstate[0xB]=s2, streamptr=p2;
		s2=mstate[0xC]<<16|*(uint16_t*)streamptr; p2=streamptr+sizeof(uint16_t); if(mstate[0xC]<1<<(RANS_STATE_BITS-RANS_RENORM_BITS))mstate[0xC]=s2, streamptr=p2;
		s2=mstate[0xD]<<16|*(uint16_t*)streamptr; p2=streamptr+sizeof(uint16_t); if(mstate[0xD]<1<<(RANS_STATE_BITS-RANS_RENORM_BITS))mstate[0xD]=s2, streamptr=p2;
		s2=mstate[0xE]<<16|*(uint16_t*)streamptr; p2=streamptr+sizeof(uint16_t); if(mstate[0xE]<1<<(RANS_STATE_BITS-RANS_RENORM_BITS))mstate[0xE]=s2, streamptr=p2;
		s2=mstate[0xF]<<16|*(uint16_t*)streamptr; p2=streamptr+sizeof(uint16_t); if(mstate[0xF]<1<<(RANS_STATE_BITS-RANS_RENORM_BITS))mstate[0xF]=s2, streamptr=p2;
	}
	*pstreamptr=(uint8_t*)(size_t)streamptr;
}
INLINE void transpose16(uint8_t *data)
{
	uint64_t *x=(uint64_t*)data;
	/*
	Hacker's delight
	SWAP(A, B)  A^=B, B^=A, A^=B
	A	B
	1100	1010	start
	0110	1010	A^=B
	0110	1100	B^=A
	1010	1100	A^=B


	SWAR_SWAP(A, B, T)  T=(A^B<<k)&MASKHI, A^=T, B^=T>>k
	3  2  1  0
	00 01 10 11	A
	00 01 10 11	B

	00 00 01 10	A>>k
	00 01 11 01	A>>k^B
	00 01 00 01	T=(A>>k^B)&MASKLO
	00 00 10 10	B^=T
	01 00 01 00	T<<=k
	01 01 11 11	A^=T<<k

	01 01 11 11	A
	00 00 10 10	B

	
	00 01 10 11	A
	00 01 10 11	B

	01 10 11 00	B<<k
	01 11 01 11	A^B<<k
	01 00 01 00	T=(A^B<<k)&MASKHI
	01 01 11 11	A^=T
	00 01 00 01	T>>k
	00 00 10 10	B^=T>>k

	01 01 11 11	A
	00 00 10 10	B
	*/
#define SWAR_TRANSPOSE(SA0, SA1, SB0, SB1, DA0, DA1, DB0, DB1, K, MASK)\
	do\
	{\
		uint64_t a0, a1, b0, b1, t0, t1;\
		a0=SA0;\
		a1=SA1;\
		b0=SB0;\
		b1=SB1;\
		t0=(a0>>K^b0)&MASK,\
		t1=(a1>>K^b1)&MASK,\
		b0^=t0;\
		b1^=t1;\
		a0^=t0<<K;\
		a1^=t1<<K;\
		DA0=a0;\
		DA1=a1;\
		DB0=b0;\
		DB1=b1;\
	}while(0)
	/*
	0x00	0x01		0x00	0x10
	0x02	0x03		0x02	0x12
	0x04	0x05		0x04	0x14
	0x06	0x07		0x06	0x16
	0x08	0x09		0x08	0x18
	0x0A	0x0B		0x0A	0x1A
	0x0C	0x0D		0x0C	0x1C
	0x0E	0x0F		0x0E	0x1E

	0x10	0x11		0x01	0x11
	0x12	0x13		0x03	0x13
	0x14	0x15		0x05	0x15
	0x16	0x17		0x07	0x17
	0x18	0x19		0x09	0x19
	0x1A	0x1B		0x0B	0x1B
	0x1C	0x1D		0x0D	0x1D
	0x1E	0x1F		0x0F	0x1F
	*/
	SWAR_TRANSPOSE(x[0x00], x[0x11], x[0x02], x[0x13],    x[0x00], x[0x11], x[0x02], x[0x13], 0x08, 0x00FF00FF00FF00FFULL);
	SWAR_TRANSPOSE(x[0x04], x[0x15], x[0x06], x[0x17],    x[0x04], x[0x15], x[0x06], x[0x17], 0x08, 0x00FF00FF00FF00FFULL);
	SWAR_TRANSPOSE(x[0x08], x[0x19], x[0x0A], x[0x1B],    x[0x08], x[0x19], x[0x0A], x[0x1B], 0x08, 0x00FF00FF00FF00FFULL);
	SWAR_TRANSPOSE(x[0x0C], x[0x1D], x[0x0E], x[0x1F],    x[0x0C], x[0x1D], x[0x0E], x[0x1F], 0x08, 0x00FF00FF00FF00FFULL);

	SWAR_TRANSPOSE(x[0x00], x[0x11], x[0x04], x[0x15],    x[0x00], x[0x11], x[0x04], x[0x15], 0x10, 0x0000FFFF0000FFFFULL);
	SWAR_TRANSPOSE(x[0x02], x[0x13], x[0x06], x[0x17],    x[0x02], x[0x13], x[0x06], x[0x17], 0x10, 0x0000FFFF0000FFFFULL);
	SWAR_TRANSPOSE(x[0x08], x[0x19], x[0x0C], x[0x1D],    x[0x08], x[0x19], x[0x0C], x[0x1D], 0x10, 0x0000FFFF0000FFFFULL);
	SWAR_TRANSPOSE(x[0x0A], x[0x1B], x[0x0E], x[0x1F],    x[0x0A], x[0x1B], x[0x0E], x[0x1F], 0x10, 0x0000FFFF0000FFFFULL);

	SWAR_TRANSPOSE(x[0x00], x[0x11], x[0x08], x[0x19],    x[0x00], x[0x11], x[0x08], x[0x19], 0x20, 0x00000000FFFFFFFFULL);
	SWAR_TRANSPOSE(x[0x02], x[0x13], x[0x0A], x[0x1B],    x[0x02], x[0x13], x[0x0A], x[0x1B], 0x20, 0x00000000FFFFFFFFULL);
	SWAR_TRANSPOSE(x[0x04], x[0x15], x[0x0C], x[0x1D],    x[0x04], x[0x15], x[0x0C], x[0x1D], 0x20, 0x00000000FFFFFFFFULL);
	SWAR_TRANSPOSE(x[0x06], x[0x17], x[0x0E], x[0x1F],    x[0x06], x[0x17], x[0x0E], x[0x1F], 0x20, 0x00000000FFFFFFFFULL);
	

	SWAR_TRANSPOSE(x[0x10], x[0x01], x[0x12], x[0x03],    x[0x10], x[0x01], x[0x12], x[0x03], 0x08, 0x00FF00FF00FF00FFULL);
	SWAR_TRANSPOSE(x[0x14], x[0x05], x[0x16], x[0x07],    x[0x14], x[0x05], x[0x16], x[0x07], 0x08, 0x00FF00FF00FF00FFULL);
	SWAR_TRANSPOSE(x[0x18], x[0x09], x[0x1A], x[0x0B],    x[0x18], x[0x09], x[0x1A], x[0x0B], 0x08, 0x00FF00FF00FF00FFULL);
	SWAR_TRANSPOSE(x[0x1C], x[0x0D], x[0x1E], x[0x0F],    x[0x1C], x[0x0D], x[0x1E], x[0x0F], 0x08, 0x00FF00FF00FF00FFULL);

	SWAR_TRANSPOSE(x[0x10], x[0x01], x[0x14], x[0x05],    x[0x10], x[0x01], x[0x14], x[0x05], 0x10, 0x0000FFFF0000FFFFULL);
	SWAR_TRANSPOSE(x[0x12], x[0x03], x[0x16], x[0x07],    x[0x12], x[0x03], x[0x16], x[0x07], 0x10, 0x0000FFFF0000FFFFULL);
	SWAR_TRANSPOSE(x[0x18], x[0x09], x[0x1C], x[0x0D],    x[0x18], x[0x09], x[0x1C], x[0x0D], 0x10, 0x0000FFFF0000FFFFULL);
	SWAR_TRANSPOSE(x[0x1A], x[0x0B], x[0x1E], x[0x0F],    x[0x1A], x[0x0B], x[0x1E], x[0x0F], 0x10, 0x0000FFFF0000FFFFULL);

	SWAR_TRANSPOSE(x[0x10], x[0x01], x[0x18], x[0x09],    x[0x01], x[0x10], x[0x09], x[0x18], 0x20, 0x00000000FFFFFFFFULL);
	SWAR_TRANSPOSE(x[0x12], x[0x03], x[0x1A], x[0x0B],    x[0x03], x[0x12], x[0x0B], x[0x1A], 0x20, 0x00000000FFFFFFFFULL);
	SWAR_TRANSPOSE(x[0x14], x[0x05], x[0x1C], x[0x0D],    x[0x05], x[0x14], x[0x0D], x[0x1C], 0x20, 0x00000000FFFFFFFFULL);
	SWAR_TRANSPOSE(x[0x16], x[0x07], x[0x1E], x[0x0F],    x[0x07], x[0x16], x[0x0F], x[0x1E], 0x20, 0x00000000FFFFFFFFULL);

#undef  SWAR_TRANSPOSE
}
static void interleave_blocks_fwd(const uint8_t *original, int iw, int ih, uint8_t *interleaved)
{
	//original[ih][iw][3]
	//interleaved[ih/YCODERS][iw/XCODERS][3][NLANES]	contiguous & aligned
	//xrem[ih%YCODERS][iw][3]
	//yrem[ih-iw%YCODERS][iw%XCODERS][3]

	//only difference between fwd and inv:		swap assignments (const slow->fast)
	int rowstride=3*iw;
	int ixyblockw=iw/XCODERS;
	int ixyblockh=ih/YCODERS;
	int blockxbytes=3*NLANES*ixyblockw;
#ifdef INTERLEAVESIMD
	int SIMDxcount=blockxbytes&~((int)sizeof(uint8_t[16][NLANES])-1);
	int slowinc=sizeof(uint8_t[16]);
#endif
	uint8_t *fastptr=interleaved;
	ALIGN(32) const uint8_t *slowptrs[NLANES]={0}, *slowptrs0[NLANES]={0};
	for(int ky=0;ky<YCODERS;++ky)//spread slow pointers
	{
		for(int kx=0;kx<XCODERS;++kx)
			slowptrs0[XCODERS*ky+kx]=original+3*(iw*ixyblockh*ky+ixyblockw*kx);
	}
	for(int ky=0;ky<ixyblockh;++ky)//interleave
	{
		int kx=0;
		memcpy((void*)slowptrs, slowptrs0, sizeof(slowptrs));
#ifdef INTERLEAVESIMD
		for(;kx<SIMDxcount;kx+=(int)sizeof(uint8_t[16][NLANES]))
		{
			/*
			SIMD interleave
			|000111222333444555666777888999AAABBBCCCDDDEEEFFFGGGHHHIIIJJJKKKLLLMMMNNN	24 pixels * 3 channels
			|^        ^        ^        ^        ^        ^        ^        ^		8-lane serial interleave
			|^^^^     ^^^^     ^^^^     ^^^^     ^^^^     ^^^^     ^^^^     ^^^^		8-lane 4xSIMD interleave step 1
			|    ^^^^     ^^^^     ^^^^     ^^^^     ^^^^     ^^^^     ^^^^     ^^^^	8-lane 4xSIMD interleave step 2
			|        ^        ^        ^        ^        ^        ^        ^        ^	8-lane serial remainder
			|0369CFIL 0369CFIL 0369CFIL 147ADGJM 147ADGJM 147ADGJM 258BEHKN 258BEHKN 258BEHKN	serial behavior
			|0369CFIL0369CFIL0369CFIL147ADGJM  147ADGJM147ADGJM258BEHKN258BEHKN  2 5 8 B E H K N	SIMD quotients & serial remainder

			slowptrs[NLANES]		fastptr  (aligned because NLANES == sizeof(__m256i))
			A16x16 B16x16		->	A16x16T C16x16T B16x16T D16x16T
			C16x16 D16x16

			load, transpose, store blocks A, B
			load, transpose, store blocks C, D

			speed:
			SIMD maximum estimate	(3/load + 3/store)*32+80*0.5*2 = 272 cycles for 32*32 bytes  ->  3GHz*32*32/272 ~= 10770 MB/s without remainders, cache misses
			actual SIMD		MSVC 4900 MB/s		GCC 3900/3000 MB/s
			actual serial		1900 MB/s
			*/
			uint64_t block[32];
			block[0x0*2+0]=((uint64_t*)slowptrs[0x00])[0]; block[0x0*2+1]=((uint64_t*)slowptrs[0x00])[1]; slowptrs[0x00]+=slowinc;
			block[0x1*2+0]=((uint64_t*)slowptrs[0x01])[0]; block[0x1*2+1]=((uint64_t*)slowptrs[0x01])[1]; slowptrs[0x01]+=slowinc;
			block[0x2*2+0]=((uint64_t*)slowptrs[0x02])[0]; block[0x2*2+1]=((uint64_t*)slowptrs[0x02])[1]; slowptrs[0x02]+=slowinc;
			block[0x3*2+0]=((uint64_t*)slowptrs[0x03])[0]; block[0x3*2+1]=((uint64_t*)slowptrs[0x03])[1]; slowptrs[0x03]+=slowinc;
			block[0x4*2+0]=((uint64_t*)slowptrs[0x04])[0]; block[0x4*2+1]=((uint64_t*)slowptrs[0x04])[1]; slowptrs[0x04]+=slowinc;
			block[0x5*2+0]=((uint64_t*)slowptrs[0x05])[0]; block[0x5*2+1]=((uint64_t*)slowptrs[0x05])[1]; slowptrs[0x05]+=slowinc;
			block[0x6*2+0]=((uint64_t*)slowptrs[0x06])[0]; block[0x6*2+1]=((uint64_t*)slowptrs[0x06])[1]; slowptrs[0x06]+=slowinc;
			block[0x7*2+0]=((uint64_t*)slowptrs[0x07])[0]; block[0x7*2+1]=((uint64_t*)slowptrs[0x07])[1]; slowptrs[0x07]+=slowinc;
			block[0x8*2+0]=((uint64_t*)slowptrs[0x08])[0]; block[0x8*2+1]=((uint64_t*)slowptrs[0x08])[1]; slowptrs[0x08]+=slowinc;
			block[0x9*2+0]=((uint64_t*)slowptrs[0x09])[0]; block[0x9*2+1]=((uint64_t*)slowptrs[0x09])[1]; slowptrs[0x09]+=slowinc;
			block[0xA*2+0]=((uint64_t*)slowptrs[0x0A])[0]; block[0xA*2+1]=((uint64_t*)slowptrs[0x0A])[1]; slowptrs[0x0A]+=slowinc;
			block[0xB*2+0]=((uint64_t*)slowptrs[0x0B])[0]; block[0xB*2+1]=((uint64_t*)slowptrs[0x0B])[1]; slowptrs[0x0B]+=slowinc;
			block[0xC*2+0]=((uint64_t*)slowptrs[0x0C])[0]; block[0xC*2+1]=((uint64_t*)slowptrs[0x0C])[1]; slowptrs[0x0C]+=slowinc;
			block[0xD*2+0]=((uint64_t*)slowptrs[0x0D])[0]; block[0xD*2+1]=((uint64_t*)slowptrs[0x0D])[1]; slowptrs[0x0D]+=slowinc;
			block[0xE*2+0]=((uint64_t*)slowptrs[0x0E])[0]; block[0xE*2+1]=((uint64_t*)slowptrs[0x0E])[1]; slowptrs[0x0E]+=slowinc;
			block[0xF*2+0]=((uint64_t*)slowptrs[0x0F])[0]; block[0xF*2+1]=((uint64_t*)slowptrs[0x0F])[1]; slowptrs[0x0F]+=slowinc;
			transpose16((uint8_t*)block);
			((uint64_t*)fastptr)[0x00]=block[0x00]; ((uint64_t*)fastptr)[0x01]=block[0x01]; ((uint64_t*)fastptr)[0x02]=block[0x02]; ((uint64_t*)fastptr)[0x03]=block[0x03];
			((uint64_t*)fastptr)[0x04]=block[0x04]; ((uint64_t*)fastptr)[0x05]=block[0x05]; ((uint64_t*)fastptr)[0x06]=block[0x06]; ((uint64_t*)fastptr)[0x07]=block[0x07];
			((uint64_t*)fastptr)[0x08]=block[0x08]; ((uint64_t*)fastptr)[0x09]=block[0x09]; ((uint64_t*)fastptr)[0x0A]=block[0x0A]; ((uint64_t*)fastptr)[0x0B]=block[0x0B];
			((uint64_t*)fastptr)[0x0C]=block[0x0C]; ((uint64_t*)fastptr)[0x0D]=block[0x0D]; ((uint64_t*)fastptr)[0x0E]=block[0x0E]; ((uint64_t*)fastptr)[0x0F]=block[0x0F];
			((uint64_t*)fastptr)[0x10]=block[0x10]; ((uint64_t*)fastptr)[0x11]=block[0x11]; ((uint64_t*)fastptr)[0x12]=block[0x12]; ((uint64_t*)fastptr)[0x13]=block[0x13];
			((uint64_t*)fastptr)[0x14]=block[0x14]; ((uint64_t*)fastptr)[0x15]=block[0x15]; ((uint64_t*)fastptr)[0x16]=block[0x16]; ((uint64_t*)fastptr)[0x17]=block[0x17];
			((uint64_t*)fastptr)[0x18]=block[0x18]; ((uint64_t*)fastptr)[0x19]=block[0x19]; ((uint64_t*)fastptr)[0x1A]=block[0x1A]; ((uint64_t*)fastptr)[0x1B]=block[0x1B];
			((uint64_t*)fastptr)[0x1C]=block[0x1C]; ((uint64_t*)fastptr)[0x1D]=block[0x1D]; ((uint64_t*)fastptr)[0x1E]=block[0x1E]; ((uint64_t*)fastptr)[0x1F]=block[0x1F];

			fastptr+=sizeof(uint8_t[16][NLANES]);
		}
#endif
#if 1
		for(;kx<blockxbytes;kx+=NLANES)
		{
			/*
			toy example
			CWH 3*4*5 original
			4 coders
			2*2 blocks
			4 Y's then 4 U's then 4 V's

			original	channel[block][pixel]
			r00 g00 b00 r01 g01 b01  r10 g10 b10 r11 g11 b11
			r02 g02 b02 r03 g03 b03  r12 g12 b12 r13 g13 b13
			r20 g20 b20 r21 g21 b21  r30 g30 b30 r31 g31 b31
			r22 g22 b22 r23 g23 b23  r32 g32 b32 r33 g33 b33

			interleaved
			r00 r10 r20 r30 g00 g10 g20 g30 b00 b10 b20 b30  r01 r11 r21 r31 g01 g11 g21 g31 b01 b11 b21 b31
			r02 r12 r22 r32 g02 g12 g22 g32 b02 b12 b22 b32  r03 r13 r23 r33 g03 g13 g23 g33 b03 b13 b23 b33
			*/
			fastptr[0x00]=*slowptrs[0x00]++;
			fastptr[0x01]=*slowptrs[0x01]++;
			fastptr[0x02]=*slowptrs[0x02]++;
			fastptr[0x03]=*slowptrs[0x03]++;
			fastptr[0x04]=*slowptrs[0x04]++;
			fastptr[0x05]=*slowptrs[0x05]++;
			fastptr[0x06]=*slowptrs[0x06]++;
			fastptr[0x07]=*slowptrs[0x07]++;
			fastptr[0x08]=*slowptrs[0x08]++;
			fastptr[0x09]=*slowptrs[0x09]++;
			fastptr[0x0A]=*slowptrs[0x0A]++;
			fastptr[0x0B]=*slowptrs[0x0B]++;
			fastptr[0x0C]=*slowptrs[0x0C]++;
			fastptr[0x0D]=*slowptrs[0x0D]++;
			fastptr[0x0E]=*slowptrs[0x0E]++;
			fastptr[0x0F]=*slowptrs[0x0F]++;
			fastptr+=NLANES;
//#if defined __GNUC__ && !defined PROFILER
//#pragma GCC unroll 32
//#endif
//			for(int k=0;k<NLANES;++k)
//				*fastptr++=*slowptrs[k]++;
		}
#endif
		for(int k=0;k<NLANES;++k)
			slowptrs0[k]+=rowstride;
	}
}
static void interleave_blocks_inv(const uint8_t *interleaved, int iw, int ih, uint8_t *original)
{
	//original[ih][iw][3]
	//interleaved[ih/YCODERS][iw/XCODERS][3][NLANES]	contiguous & aligned
	//xrem[ih%YCODERS][iw][3]
	//yrem[ih-iw%YCODERS][iw%XCODERS][3]

	//only difference between fwd and inv:		swap assignments (const slow->fast)
	int rowstride=3*iw;
	int ixyblockw=iw/XCODERS;
	int ixyblockh=ih/YCODERS;
	int blockxbytes=3*NLANES*ixyblockw;
#ifdef INTERLEAVESIMD
	int SIMDxcount=blockxbytes&~((int)sizeof(uint8_t[16][NLANES])-1);
	int slowinc=sizeof(uint8_t[16]);
#endif
	const uint8_t *fastptr=interleaved;
	ALIGN(32) uint8_t *slowptrs[NLANES]={0}, *slowptrs0[NLANES]={0};
	for(int ky=0;ky<YCODERS;++ky)//spread slow pointers
	{
		for(int kx=0;kx<XCODERS;++kx)
			slowptrs0[XCODERS*ky+kx]=original+3*(iw*ixyblockh*ky+ixyblockw*kx);
	}
	for(int ky=0;ky<ixyblockh;++ky)//interleave
	{
		int kx=0;
		memcpy((void*)slowptrs, slowptrs0, sizeof(slowptrs));
#ifdef INTERLEAVESIMD
		for(;kx<SIMDxcount;kx+=(int)sizeof(uint8_t[16][NLANES]))
		{
			/*
			SIMD interleave
			|000111222333444555666777888999AAABBBCCCDDDEEEFFFGGGHHHIIIJJJKKKLLLMMMNNN	24 pixels * 3 channels
			|^        ^        ^        ^        ^        ^        ^        ^		8-lane serial interleave
			|^^^^     ^^^^     ^^^^     ^^^^     ^^^^     ^^^^     ^^^^     ^^^^		8-lane 4xSIMD interleave step 1
			|    ^^^^     ^^^^     ^^^^     ^^^^     ^^^^     ^^^^     ^^^^     ^^^^	8-lane 4xSIMD interleave step 2
			|        ^        ^        ^        ^        ^        ^        ^        ^	8-lane serial remainder
			|0369CFIL 0369CFIL 0369CFIL 147ADGJM 147ADGJM 147ADGJM 258BEHKN 258BEHKN 258BEHKN	serial behavior
			|0369CFIL0369CFIL0369CFIL147ADGJM  147ADGJM147ADGJM258BEHKN258BEHKN  2 5 8 B E H K N	SIMD quotients & serial remainder

			slowptrs[NLANES]		fastptr  (aligned because NLANES == sizeof(__m256i))
			A16x16 B16x16		->	A16x16T C16x16T B16x16T D16x16T
			C16x16 D16x16

			load, transpose, store blocks A, B
			load, transpose, store blocks C, D

			speed:
			SIMD maximum estimate	(3/load + 3/store)*32+80*0.5*2 = 272 cycles for 32*32 bytes  ->  3GHz*32*32/272 ~= 10770 MB/s without remainders, cache misses
			actual SIMD		MSVC 4900 MB/s		GCC 3900/3000 MB/s
			actual serial		1900 MB/s
			*/
			uint64_t block[32];
			block[0x00]=((uint64_t*)fastptr)[0x00]; block[0x01]=((uint64_t*)fastptr)[0x01]; block[0x02]=((uint64_t*)fastptr)[0x02]; block[0x03]=((uint64_t*)fastptr)[0x03];
			block[0x04]=((uint64_t*)fastptr)[0x04]; block[0x05]=((uint64_t*)fastptr)[0x05]; block[0x06]=((uint64_t*)fastptr)[0x06]; block[0x07]=((uint64_t*)fastptr)[0x07];
			block[0x08]=((uint64_t*)fastptr)[0x08]; block[0x09]=((uint64_t*)fastptr)[0x09]; block[0x0A]=((uint64_t*)fastptr)[0x0A]; block[0x0B]=((uint64_t*)fastptr)[0x0B];
			block[0x0C]=((uint64_t*)fastptr)[0x0C]; block[0x0D]=((uint64_t*)fastptr)[0x0D]; block[0x0E]=((uint64_t*)fastptr)[0x0E]; block[0x0F]=((uint64_t*)fastptr)[0x0F];
			block[0x10]=((uint64_t*)fastptr)[0x10]; block[0x11]=((uint64_t*)fastptr)[0x11]; block[0x12]=((uint64_t*)fastptr)[0x12]; block[0x13]=((uint64_t*)fastptr)[0x13];
			block[0x14]=((uint64_t*)fastptr)[0x14]; block[0x15]=((uint64_t*)fastptr)[0x15]; block[0x16]=((uint64_t*)fastptr)[0x16]; block[0x17]=((uint64_t*)fastptr)[0x17];
			block[0x18]=((uint64_t*)fastptr)[0x18]; block[0x19]=((uint64_t*)fastptr)[0x19]; block[0x1A]=((uint64_t*)fastptr)[0x1A]; block[0x1B]=((uint64_t*)fastptr)[0x1B];
			block[0x1C]=((uint64_t*)fastptr)[0x1C]; block[0x1D]=((uint64_t*)fastptr)[0x1D]; block[0x1E]=((uint64_t*)fastptr)[0x1E]; block[0x1F]=((uint64_t*)fastptr)[0x1F];
			transpose16((uint8_t*)block);
			((uint64_t*)slowptrs[0x00])[0]=block[0x0*2+0]; ((uint64_t*)slowptrs[0x00])[1]=block[0x0*2+1]; slowptrs[0x00]+=slowinc;
			((uint64_t*)slowptrs[0x01])[0]=block[0x1*2+0]; ((uint64_t*)slowptrs[0x01])[1]=block[0x1*2+1]; slowptrs[0x01]+=slowinc;
			((uint64_t*)slowptrs[0x02])[0]=block[0x2*2+0]; ((uint64_t*)slowptrs[0x02])[1]=block[0x2*2+1]; slowptrs[0x02]+=slowinc;
			((uint64_t*)slowptrs[0x03])[0]=block[0x3*2+0]; ((uint64_t*)slowptrs[0x03])[1]=block[0x3*2+1]; slowptrs[0x03]+=slowinc;
			((uint64_t*)slowptrs[0x04])[0]=block[0x4*2+0]; ((uint64_t*)slowptrs[0x04])[1]=block[0x4*2+1]; slowptrs[0x04]+=slowinc;
			((uint64_t*)slowptrs[0x05])[0]=block[0x5*2+0]; ((uint64_t*)slowptrs[0x05])[1]=block[0x5*2+1]; slowptrs[0x05]+=slowinc;
			((uint64_t*)slowptrs[0x06])[0]=block[0x6*2+0]; ((uint64_t*)slowptrs[0x06])[1]=block[0x6*2+1]; slowptrs[0x06]+=slowinc;
			((uint64_t*)slowptrs[0x07])[0]=block[0x7*2+0]; ((uint64_t*)slowptrs[0x07])[1]=block[0x7*2+1]; slowptrs[0x07]+=slowinc;
			((uint64_t*)slowptrs[0x08])[0]=block[0x8*2+0]; ((uint64_t*)slowptrs[0x08])[1]=block[0x8*2+1]; slowptrs[0x08]+=slowinc;
			((uint64_t*)slowptrs[0x09])[0]=block[0x9*2+0]; ((uint64_t*)slowptrs[0x09])[1]=block[0x9*2+1]; slowptrs[0x09]+=slowinc;
			((uint64_t*)slowptrs[0x0A])[0]=block[0xA*2+0]; ((uint64_t*)slowptrs[0x0A])[1]=block[0xA*2+1]; slowptrs[0x0A]+=slowinc;
			((uint64_t*)slowptrs[0x0B])[0]=block[0xB*2+0]; ((uint64_t*)slowptrs[0x0B])[1]=block[0xB*2+1]; slowptrs[0x0B]+=slowinc;
			((uint64_t*)slowptrs[0x0C])[0]=block[0xC*2+0]; ((uint64_t*)slowptrs[0x0C])[1]=block[0xC*2+1]; slowptrs[0x0C]+=slowinc;
			((uint64_t*)slowptrs[0x0D])[0]=block[0xD*2+0]; ((uint64_t*)slowptrs[0x0D])[1]=block[0xD*2+1]; slowptrs[0x0D]+=slowinc;
			((uint64_t*)slowptrs[0x0E])[0]=block[0xE*2+0]; ((uint64_t*)slowptrs[0x0E])[1]=block[0xE*2+1]; slowptrs[0x0E]+=slowinc;
			((uint64_t*)slowptrs[0x0F])[0]=block[0xF*2+0]; ((uint64_t*)slowptrs[0x0F])[1]=block[0xF*2+1]; slowptrs[0x0F]+=slowinc;

			fastptr+=sizeof(uint8_t[16][NLANES]);
		}
#endif
#if 1
		for(;kx<blockxbytes;kx+=NLANES)
		{
			/*
			toy example
			CWH 3*4*5 original
			4 coders
			2*2 blocks
			4 Y's then 4 U's then 4 V's

			original	channel[block][pixel]
			r00 g00 b00 r01 g01 b01  r10 g10 b10 r11 g11 b11
			r02 g02 b02 r03 g03 b03  r12 g12 b12 r13 g13 b13
			r20 g20 b20 r21 g21 b21  r30 g30 b30 r31 g31 b31
			r22 g22 b22 r23 g23 b23  r32 g32 b32 r33 g33 b33

			interleaved
			r00 r10 r20 r30 g00 g10 g20 g30 b00 b10 b20 b30  r01 r11 r21 r31 g01 g11 g21 g31 b01 b11 b21 b31
			r02 r12 r22 r32 g02 g12 g22 g32 b02 b12 b22 b32  r03 r13 r23 r33 g03 g13 g23 g33 b03 b13 b23 b33
			*/
			*slowptrs[0x00]++=fastptr[0x00];
			*slowptrs[0x01]++=fastptr[0x01];
			*slowptrs[0x02]++=fastptr[0x02];
			*slowptrs[0x03]++=fastptr[0x03];
			*slowptrs[0x04]++=fastptr[0x04];
			*slowptrs[0x05]++=fastptr[0x05];
			*slowptrs[0x06]++=fastptr[0x06];
			*slowptrs[0x07]++=fastptr[0x07];
			*slowptrs[0x08]++=fastptr[0x08];
			*slowptrs[0x09]++=fastptr[0x09];
			*slowptrs[0x0A]++=fastptr[0x0A];
			*slowptrs[0x0B]++=fastptr[0x0B];
			*slowptrs[0x0C]++=fastptr[0x0C];
			*slowptrs[0x0D]++=fastptr[0x0D];
			*slowptrs[0x0E]++=fastptr[0x0E];
			*slowptrs[0x0F]++=fastptr[0x0F];
			fastptr+=NLANES;
//#if defined __GNUC__ && !defined PROFILER
//#pragma GCC unroll 32
//#endif
//			for(int k=0;k<NLANES;++k)
//				*slowptrs[k]++=*fastptr++;
		}
#endif
		for(int k=0;k<NLANES;++k)
			slowptrs0[k]+=rowstride;
	}
}
int codec_l1_port(int argc, char **argv)
{
	if(argc!=3&&argc!=4&&argc!=5)
	{
		printf(
			"Usage:  \"%s\"  input  output  [Effort]  [Dist]    To encode/decode.\n"
			"  Effort  =  0 CG / 1~3 L1 | 4 Profiler.\n"
			"  Dist    =  lossy distortion. 4 <= Dist <= 31.\n"
			"Built on %s %s\n"
			, argv[0]
			, __DATE__, __TIME__
		);
		return 1;
	}
	const char *srcfn=argv[1], *dstfn=argv[2];
	int param1=argc<4?DEFAULT_EFFORT_LEVEL:atoi(argv[3]), dist=argc<5?1:atoi(argv[4]);
	int effort=param1&3, profile=param1>>2;
	if(dist>1)
		CLAMP2(dist, 3, 31);
#ifdef ESTIMATE_SIZE
	double esize[3*NLANES]={0};
#endif
#ifdef LOUD
	double t=time_sec();
	ptrdiff_t usize2=0;
	{
		struct stat info={0};
		stat(srcfn, &info);
		usize2=info.st_size;
	}
#endif
	prof_timestamp=time_sec();
	//prof_checkpoint(0, 0);
	if(!srcfn||!dstfn)
	{
		CRASH("Codec requires both source and destination filenames");
		return 1;
	}
	int fwd=0, iw=0, ih=0, rowstride=0;
	int bestrct=0, npreds=0, sh=0;
	uint64_t bypassmask=0;//0: emit stats  1: rare context (bypass)
	ptrdiff_t usize=0, cap=0;
	uint8_t *image=0, *imptr=0, *streamptr=0, *streamstart=0, *streamend=0;
	int psize=0;
	int16_t *pixels=0;
	ptrdiff_t cheadersize=0, csize=0;
	{
		FILE *fsrc=fopen(srcfn, "rb");
		if(!fsrc)
		{
			CRASH("Cannot open \"%s\"", srcfn);
			return 1;
		}
		int c=0;
		fread(&c, 1, 2, fsrc);
		fwd=c==('P'|'6'<<8);
		if(!fwd&&c!=('L'|'1'<<8))
		{
			CRASH("Unsupported file \"%s\"", srcfn);
			return 1;
		}
		if(fwd)
		{
			int nread=0, vmax=0;
#ifdef LOUD
			print_timestamp("%Y-%m-%d_%H%M%S\n");
#endif
			c=fgetc(fsrc);
			if(c!='\n')
			{
				CRASH("Invalid PPM file");
				return 1;
			}
			nread=fscanf(fsrc, "%d %d", &iw, &ih);
			if(nread!=2)
			{
				CRASH("Unsupported PPM file");
				return 1;
			}
			nread=fscanf(fsrc, "%d", &vmax);
			if(nread!=1||vmax!=255)
			{
				CRASH("Unsupported PPM file");
				return 1;
			}
			c=fgetc(fsrc);
			if(c!='\n')
			{
				CRASH("Invalid PPM file");
				return 1;
			}
		}
		else
		{
			int flags=0;

			iw=0;
			ih=0;
			dist=0;
			fread(&iw, 1, 3, fsrc);
			fread(&ih, 1, 3, fsrc);
			fread(&flags, 1, 1, fsrc);
			bestrct=flags>>2;
			effort=flags&3;
			fread(&dist, 1, 1, fsrc);
			fread(&bypassmask, 1, 8, fsrc);
			cheadersize=ftell(fsrc);
		}
		if(iw<1||ih<1)
		{
			CRASH("Unsupported source file");
			return 1;
		}
		rowstride=3*iw;
		usize=(ptrdiff_t)3*iw*ih;
		cap=(ptrdiff_t)4*iw*ih;
		image=(uint8_t*)malloc(cap+sizeof(uint8_t[32]));
		if(!image)
		{
			CRASH("Alloc error");
			return 1;
		}
		if(fwd)
		{
			fread(image, 1, usize, fsrc);//read image
			streamptr=streamstart=image+cap;//bwd-bwd ANS encoding
			profile_size(streamptr, "start");
		}
		else
		{
			struct stat info={0};
			stat(srcfn, &info);
			csize=info.st_size;
			streamptr=streamstart=image+cap-(csize-cheadersize)-sizeof(uint8_t[32]);
			streamend=image+cap-sizeof(uint8_t[32]);
			fread(streamstart, 1, csize-cheadersize, fsrc);//read stream
		}
		fclose(fsrc);
	}
	prof_checkpoint(fwd?usize:csize, "fread");
	int blockw=iw/XCODERS;
	int blockh=ih/YCODERS;
	int qxbytes=blockw*XCODERS*3;//iw/XCODERS*XCODERS*3
	int ixcount=blockw*NLANES, ixbytes=3*ixcount;//ix = interleaved circular buffer width		iw/XCODERS*NLANES
	int xremw=iw-blockw*XCODERS, yremh=ih-blockh*YCODERS;
	int xrembytes=3*xremw;
	int nctx=3*NCTX+3*(xremw||yremh);
	ptrdiff_t isize=(ptrdiff_t)ixbytes*blockh;
	ptrdiff_t interleavedsize=isize<<fwd;//fwd ? interleave residuals & context : pack residuals
	uint8_t *interleaved=(uint8_t*)malloc(interleavedsize);
	if(!interleaved)
	{
		CRASH("Alloc error");
		return 1;
	}
	(void)xrembytes;
	const int hsize=nctx*(int)sizeof(int[256]);
	int *hists=fwd?(int*)malloc(hsize):0;//fwd-only

	int CDF2syms_size=nctx*(int)sizeof(int[1<<PROBBITS]);
	if(fwd)//DIV-free rANS encoder reuses this as SIMD symbol info
		CDF2syms_size=nctx*(int)sizeof(rANS_SIMD_SymInfo[256]);
	uint32_t *CDF2syms=(uint32_t*)malloc(CDF2syms_size);

	psize=(blockw+2*XPAD)*(int)sizeof(int16_t[NROWS*NVAL]);//int16_t[blockw+2*XPAD][NROWS*NVAL]
	pixels=(int16_t*)malloc(psize);

	int ctsize=(int)sizeof(int16_t[512<<GRBITS]);
	int16_t *ctable=(int16_t*)malloc(ctsize);
	if((fwd&&!hists)||!CDF2syms||!pixels||!ctable)
	{
		CRASH("Alloc error");
		return 1;
	}
	for(int k=0;k<512<<GRBITS;++k)
	{
		int ctx=FLOOR_LOG2(k*k+1);
		if(ctx>NCTX-1)
			ctx=NCTX-1;
		ctable[k]=ctx;
	}
	if(fwd)
	{
		memset(hists, 0, hsize);
#ifdef TEST_INTERLEAVE
		{
			uint8_t x[256];
			for(int k=0;k<256;++k)
				x[k]=k;
			transpose16(x);
			for(int k=0;k<256;++k)
			{
				if(x[k]!=(uint8_t)(k<<4|k>>4))
				{
					CRASH("");
					return 1;
				}
			}
			transpose16(x);
			for(int k=0;k<256;++k)
			{
				if(x[k]!=k)
				{
					CRASH("");
					return 1;
				}
			}
		}
		system("cd");
		guide_save(image, iw, ih);
		save_ppm("20250227_1225AM_original.PPM", image, iw, ih);
		interleave_blocks_fwd(image, iw, ih, interleaved);
		save_ppm("20250226_1153PM_interleaved.PPM", interleaved, ixcount, blockh);
		interleave_blocks_inv(interleaved, iw, ih, image);
		save_ppm("20250227_1244AM_deinterleaved.PPM", image, iw, ih);
		if(memcmp(image, g_image, usize))
			CRASH("ERROR");
		printf("SUCCESS\n");
		exit(0);
#endif
		interleave_blocks_fwd(image, iw, ih, interleaved+isize);//reuse memory: read 8-bit pixel, write 16-bit context<<8|residual
		guide_save(interleaved+isize, ixcount, blockh);
		prof_checkpoint(usize, "interleave");
		{//analysis
			ALIGN(32) int64_t counters[OCH_COUNT]={0};
			imptr=interleaved+isize+ixbytes+3*NLANES;
			for(int ky=1;ky<blockh;ky+=ANALYSIS_YSTRIDE)
			{
				for(int kx=1;kx<blockw-(ANALYSIS_XSTRIDE-1);kx+=ANALYSIS_XSTRIDE)
				{
					for(int k=0;k<16;++k)
					{
						int r, g, b, rg, gb, br;

						r=(imptr[k+0*NLANES]-imptr[k+0*NLANES-3*NLANES]-imptr[k-ixbytes+0*NLANES]+imptr[k-ixbytes+0*NLANES-3*NLANES])<<2;
						g=(imptr[k+1*NLANES]-imptr[k+1*NLANES-3*NLANES]-imptr[k-ixbytes+1*NLANES]+imptr[k-ixbytes+1*NLANES-3*NLANES])<<2;
						b=(imptr[k+2*NLANES]-imptr[k+2*NLANES-3*NLANES]-imptr[k-ixbytes+2*NLANES]+imptr[k-ixbytes+2*NLANES-3*NLANES])<<2;
						rg=r-g;
						gb=g-b;
						br=b-r;
						counters[OCH_YX00]+=abs(r);
						counters[OCH_Y0X0]+=abs(g);
						counters[OCH_Y00X]+=abs(b);
						counters[OCH_CX40]+=abs(rg);
						counters[OCH_C0X4]+=abs(gb);
						counters[OCH_C40X]+=abs(br);
#ifdef ENABLE_RCT_EXTENSION
						counters[OCH_CX31]+=abs(rg+(gb>>2));//r-(3*g+b)/4 = r-g-(b-g)/4
						counters[OCH_C3X1]+=abs(rg+(br>>2));//g-(3*r+b)/4 = g-r-(b-r)/4
						counters[OCH_C31X]+=abs(br+(rg>>2));//b-(3*r+g)/4 = b-r-(g-r)/4
						counters[OCH_CX13]+=abs(br+(gb>>2));//r-(g+3*b)/4 = r-b-(g-b)/4
						counters[OCH_C1X3]+=abs(gb+(br>>2));//g-(r+3*b)/4 = g-b-(r-b)/4
						counters[OCH_C13X]+=abs(gb+(rg>>2));//b-(r+3*g)/4 = b-g-(r-g)/4
						counters[OCH_CX22]+=abs((rg-br)>>1);//r-(g+b)/2 = (r-g + r-b)/2
						counters[OCH_C2X2]+=abs((gb-rg)>>1);//g-(r+b)/2 = (g-r + g-b)/2
						counters[OCH_C22X]+=abs((br-gb)>>1);//b-(r+g)/2 = (b-r + b-g)/2
#endif
					}
					imptr+=3*NLANES*ANALYSIS_XSTRIDE;
				}
				imptr+=ixbytes*(ANALYSIS_YSTRIDE-1);
			}
			{
				int64_t minerr;
				int kt;

				for(kt=0, minerr=0;kt<RCT_COUNT;++kt)
				{
					const uint8_t *rct=rct_combinations[kt];
					int64_t currerr=
						+counters[rct[0]]
						+counters[rct[1]]
						+counters[rct[2]]
					;
#ifdef LOUD
					printf("%2d  %-14s %12lld + %12lld + %12lld = %12lld%s\n"
						, kt
						, rct_names[kt]
						, counters[rct[0]]
						, counters[rct[1]]
						, counters[rct[2]]
						, currerr
						, !kt||minerr>currerr?" <-":""
					);
#endif
					if(!kt||minerr>currerr)
					{
						minerr=currerr;
						bestrct=kt;
					}
				}
			}
			prof_checkpoint(usize, "analysis");
		}
	}
	else
	{
		//decode stats
		BitPackerLIFO ec;
		bitpacker_dec_init(&ec, streamptr, streamend);
		for(int kc=0;kc<nctx;++kc)
			dec_unpackhist(&ec, CDF2syms+((ptrdiff_t)kc<<PROBBITS), bypassmask, kc);
		streamptr=(uint8_t*)(size_t)ec.srcfwdptr;
		prof_checkpoint(CDF2syms_size, "unpack histograms");
	}
	int L1statesize=0;
	int *L1state=0;
	switch(effort)
	{
	case 0://use CG
		npreds=0;
		break;
	case 1://use Manhattan
		npreds=L1_NPREDS1;
		sh=L1_SH1;
		break;
	case 2:
		npreds=L1_NPREDS2;
		sh=L1_SH2;
		break;
	case 3:
		npreds=L1_NPREDS3;
		sh=L1_SH3;
		break;
	}
	if(fwd)
	{
#ifndef LOUD
		if(profile)
#endif
			printf("%s  NPREDS=%d  %td bytes\n", rct_names[bestrct], npreds, usize);
	}
	if(effort)
	{
		L1statesize=(int)sizeof(int[2*NLANES*3*(L1_NPREDS3+1)]);//{preds, coeffs} * (NPREDS+{bias}) * 3 channels * NLANES
		L1state=(int*)malloc(L1statesize);
		if(!L1state)
		{
			CRASH("Alloc error");
			return 1;
		}
		memset(L1state, 0, L1statesize);
	}
	const uint8_t *combination=rct_combinations[bestrct];
	int
		yidx=combination[II_PERM_Y]*NLANES,
		uidx=combination[II_PERM_U]*NLANES,
		vidx=combination[II_PERM_V]*NLANES;
	int uhelpmask=-(combination[II_COEFF_U_SUB_Y]!=0);
	int vc0=combination[II_COEFF_V_SUB_Y];
	int vc1=combination[II_COEFF_V_SUB_U];
	//int paddedwidth=blockw+2*XPAD;
	memset(pixels, 0, psize);
	int16_t myuv[3*NLANES];
	int dist_rcp=0x10000;
	if(dist>1)
		dist_rcp=((1<<16)+dist-1)/dist;//x/dist  ->  {x*=inv; x=(x>>16)+((uint32_t)x>>31);}
	memset(myuv, 0, sizeof(myuv));
	uint8_t *ctxptr=interleaved;
	imptr=interleaved+(fwd?isize:0);
	uint32_t mstate[NLANES];
	int16_t *L1preds=effort?(int16_t*)L1state:0;
	int *L1weights=effort?(int*)(L1state+1*(ptrdiff_t)NLANES*3*(L1_NPREDS3+1)):0;
	if(effort)
		FILLMEM(L1weights, (1<<sh)/npreds, (npreds+1)*sizeof(int[6*8]), sizeof(int));
	//{
	//	static const int weights0[]=
	//	{
	//		100000,//0	N
	//		100000,//1	W
	//		 80000,//2	3*(N-NN)+NNN
	//		 80000,//3	3*(W-WW)+WWW
	//		 50000,//4	W+NE-N
	//		 50000,//5	(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4
	//		150000,//6	N+W-NW
	//		 50000,//7	N+NE-NNE
	//		0,
	//	};
	//	int *L1coeffs=(int*)L1weights;
	//	for(int k=0;k<L1_NPREDS2+1;++k)
	//		FILLMEM(L1coeffs+6*8*k, weights0[k], sizeof(int[6*8]), sizeof(int));
	//}
	if(!fwd)
	{
#ifdef _DEBUG
		if(streamptr>streamend)
			CRASH("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
		memcpy(mstate, streamptr, sizeof(mstate));
		streamptr+=sizeof(mstate);
	}
	for(int ky=0;ky<blockh;++ky)//main coding loop
	{
		int16_t *rows[]=
		{
			pixels+(XPAD*NROWS+(ky-0LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NROWS+(ky-2LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NROWS+(ky-3LL+NROWS)%NROWS)*NVAL,
		};
		uint16_t syms[3*NLANES]={0};
		int16_t cW[3*NLANES]={0};
		int16_t pred[3*NLANES], pred0[3*NLANES], ctx[3*NLANES];
		int16_t msyms[3*NLANES], moffset[2*NLANES];
		for(int kx=0;kx<blockw;++kx)
		{
			//rows[Y][V+X*NROWS*NVAL]  add 3*NLANES for energy
			ctx[0x00]=ctable[rows[0][0x30 - 1*NROWS*NVAL]];
			ctx[0x01]=ctable[rows[0][0x31 - 1*NROWS*NVAL]];
			ctx[0x02]=ctable[rows[0][0x32 - 1*NROWS*NVAL]];
			ctx[0x03]=ctable[rows[0][0x33 - 1*NROWS*NVAL]];
			ctx[0x04]=ctable[rows[0][0x34 - 1*NROWS*NVAL]];
			ctx[0x05]=ctable[rows[0][0x35 - 1*NROWS*NVAL]];
			ctx[0x06]=ctable[rows[0][0x36 - 1*NROWS*NVAL]];
			ctx[0x07]=ctable[rows[0][0x37 - 1*NROWS*NVAL]];
			ctx[0x08]=ctable[rows[0][0x38 - 1*NROWS*NVAL]];
			ctx[0x09]=ctable[rows[0][0x39 - 1*NROWS*NVAL]];
			ctx[0x0A]=ctable[rows[0][0x3A - 1*NROWS*NVAL]];
			ctx[0x0B]=ctable[rows[0][0x3B - 1*NROWS*NVAL]];
			ctx[0x0C]=ctable[rows[0][0x3C - 1*NROWS*NVAL]];
			ctx[0x0D]=ctable[rows[0][0x3D - 1*NROWS*NVAL]];
			ctx[0x0E]=ctable[rows[0][0x3E - 1*NROWS*NVAL]];//not scientific notation!
			ctx[0x0F]=ctable[rows[0][0x3F - 1*NROWS*NVAL]];

			ctx[0x10]=ctable[rows[0][0x40 - 1*NROWS*NVAL]]+NCTX;
			ctx[0x11]=ctable[rows[0][0x41 - 1*NROWS*NVAL]]+NCTX;
			ctx[0x12]=ctable[rows[0][0x42 - 1*NROWS*NVAL]]+NCTX;
			ctx[0x13]=ctable[rows[0][0x43 - 1*NROWS*NVAL]]+NCTX;
			ctx[0x14]=ctable[rows[0][0x44 - 1*NROWS*NVAL]]+NCTX;
			ctx[0x15]=ctable[rows[0][0x45 - 1*NROWS*NVAL]]+NCTX;
			ctx[0x16]=ctable[rows[0][0x46 - 1*NROWS*NVAL]]+NCTX;
			ctx[0x17]=ctable[rows[0][0x47 - 1*NROWS*NVAL]]+NCTX;
			ctx[0x18]=ctable[rows[0][0x48 - 1*NROWS*NVAL]]+NCTX;
			ctx[0x19]=ctable[rows[0][0x49 - 1*NROWS*NVAL]]+NCTX;
			ctx[0x1A]=ctable[rows[0][0x4A - 1*NROWS*NVAL]]+NCTX;
			ctx[0x1B]=ctable[rows[0][0x4B - 1*NROWS*NVAL]]+NCTX;
			ctx[0x1C]=ctable[rows[0][0x4C - 1*NROWS*NVAL]]+NCTX;
			ctx[0x1D]=ctable[rows[0][0x4D - 1*NROWS*NVAL]]+NCTX;
			ctx[0x1E]=ctable[rows[0][0x4E - 1*NROWS*NVAL]]+NCTX;
			ctx[0x1F]=ctable[rows[0][0x4F - 1*NROWS*NVAL]]+NCTX;

			ctx[0x20]=ctable[rows[0][0x50 - 1*NROWS*NVAL]]+NCTX*2;
			ctx[0x21]=ctable[rows[0][0x51 - 1*NROWS*NVAL]]+NCTX*2;
			ctx[0x22]=ctable[rows[0][0x52 - 1*NROWS*NVAL]]+NCTX*2;
			ctx[0x23]=ctable[rows[0][0x53 - 1*NROWS*NVAL]]+NCTX*2;
			ctx[0x24]=ctable[rows[0][0x54 - 1*NROWS*NVAL]]+NCTX*2;
			ctx[0x25]=ctable[rows[0][0x55 - 1*NROWS*NVAL]]+NCTX*2;
			ctx[0x26]=ctable[rows[0][0x56 - 1*NROWS*NVAL]]+NCTX*2;
			ctx[0x27]=ctable[rows[0][0x57 - 1*NROWS*NVAL]]+NCTX*2;
			ctx[0x28]=ctable[rows[0][0x58 - 1*NROWS*NVAL]]+NCTX*2;
			ctx[0x29]=ctable[rows[0][0x59 - 1*NROWS*NVAL]]+NCTX*2;
			ctx[0x2A]=ctable[rows[0][0x5A - 1*NROWS*NVAL]]+NCTX*2;
			ctx[0x2B]=ctable[rows[0][0x5B - 1*NROWS*NVAL]]+NCTX*2;
			ctx[0x2C]=ctable[rows[0][0x5C - 1*NROWS*NVAL]]+NCTX*2;
			ctx[0x2D]=ctable[rows[0][0x5D - 1*NROWS*NVAL]]+NCTX*2;
			ctx[0x2E]=ctable[rows[0][0x5E - 1*NROWS*NVAL]]+NCTX*2;
			ctx[0x2F]=ctable[rows[0][0x5F - 1*NROWS*NVAL]]+NCTX*2;
#ifdef EXTRA_CHECKS
			for(int k=0;k<3*NLANES;++k)
			{
				uint32_t c=(uint32_t)(ctx[k]-k/NLANES*NCTX);
				int eW=rows[0][k+0x30-1*NROWS*NVAL];
				if(c>=NCTX||c!=FLOOR_LOG2(eW*eW+1))
				{
					CRASH("");
				}
			}
#endif
			{
				const int borderW=3;
				const int borderN=3;
				const int borderE=3;
				int cond_cg=
					(uint32_t)(kx-borderW)>=(uint32_t)(blockw-(borderW+borderE))||
					(uint32_t)(ky-borderN)>=(uint32_t)(blockh-borderN);
				int16_t vmin[3*NLANES], vmax[3*NLANES], grad[3*NLANES];

				for(int k=0;k<3*NLANES;++k)
				{
					int N, W, t0, t1;

					//rows[Y][V+X*NROWS*NVAL]
					N=rows[1][k+0*NROWS*NVAL];
					W=rows[0][k-1*NROWS*NVAL];
					t0=N, t1=W;
					if(N<W)t1=N, t0=W;
					vmin[k]=t1;
					vmax[k]=t0;
					pred[k]=grad[k]=N+W-rows[1][k-1*NROWS*NVAL];
				}
				if(!effort)
					goto effort0;
				if(effort==1)//predict
				{
					for(int k=0;k<3*NLANES;++k)
					{
						int N=rows[1][k+0*NROWS*NVAL];
						
						/*
						effort 1	5 preds
						0	W
						1	N+W-NW
						2	2*N-NN
						3	NE
						4	eW
						*/
						L1preds[0*3*NLANES+k]=rows[0][k-1*NROWS*NVAL];
						L1preds[1*3*NLANES+k]=pred[k];
						L1preds[2*3*NLANES+k]=2*N-rows[2][k+0*NROWS*NVAL];
						L1preds[3*3*NLANES+k]=rows[1][k+1*NROWS*NVAL];
						L1preds[4*3*NLANES+k]=cW[k];
					}
					for(int k=0;k<3*NLANES;++k)
					{
						pred[k]=((1<<L1_SH1>>1)
							+L1weights[0*3*NLANES+k]*L1preds[0*3*NLANES+k]
							+L1weights[1*3*NLANES+k]*L1preds[1*3*NLANES+k]
							+L1weights[2*3*NLANES+k]*L1preds[2*3*NLANES+k]
							+L1weights[3*3*NLANES+k]*L1preds[3*3*NLANES+k]
							+L1weights[4*3*NLANES+k]*L1preds[4*3*NLANES+k]
						)>>L1_SH1;
					}
					if(!cond_cg)//loosen pred range
					{
						for(int k=0;k<3*NLANES;++k)
						{
							int NE=rows[1][k+1*NROWS*NVAL], NEEE=rows[1][k+3*NROWS*NVAL];
							if(vmin[k]>NE)vmin[k]=NE;
							if(vmax[k]<NE)vmax[k]=NE;
							if(vmin[k]>NEEE)vmin[k]=NEEE;
							if(vmax[k]<NEEE)vmax[k]=NEEE;
						}
					}
					memcpy(pred0, pred, sizeof(pred0));
				}
				else if(effort==2)
				{
					for(int k=0;k<3*NLANES;++k)
					{
						int N=rows[1][k+0*NROWS*NVAL], W=rows[0][k-1*NROWS*NVAL];
						
						/*
						effort 2	8 preds
						0	N
						1	W
						2	3*(N-NN)+NNN
						3	3*(W-WW)+WWW
						4	W+NE-N
						5	(WWWW+WWW+NNN+NEE+NEEE+NEEEE-(NW+N))>>2
						6	N+W-NW
						7	N+NE-NNE
						*/
						L1preds[0*3*NLANES+k]=N;
						L1preds[1*3*NLANES+k]=W;
						L1preds[2*3*NLANES+k]=3*(N-rows[2][k+0*NROWS*NVAL])+rows[3][k+0*NROWS*NVAL];
						L1preds[3*3*NLANES+k]=3*(W-rows[0][k-2*NROWS*NVAL])+rows[0][k-3*NROWS*NVAL];
						L1preds[4*3*NLANES+k]=W+rows[1][k+1*NROWS*NVAL]-N;
						L1preds[5*3*NLANES+k]=(
							+rows[0][k-4*NROWS*NVAL]//WWWW
							+rows[0][k-3*NROWS*NVAL]//WWW
							+rows[3][k+0*NROWS*NVAL]//NNN
							+rows[1][k+2*NROWS*NVAL]//NEE
							+rows[1][k+3*NROWS*NVAL]//NEEE
							+rows[1][k+4*NROWS*NVAL]//NEEEE
							-rows[1][k-1*NROWS*NVAL]//NW
							-N
						)>>2;
						L1preds[6*3*NLANES+k]=pred[k];
						L1preds[7*3*NLANES+k]=N+rows[1][k+1*NROWS*NVAL]-rows[2][k+1*NROWS*NVAL];
					}
					for(int k=0;k<3*NLANES;++k)
					{
						pred[k]=((1<<L1_SH2>>1)
							+L1weights[0*3*NLANES+k]*L1preds[0*3*NLANES+k]
							+L1weights[1*3*NLANES+k]*L1preds[1*3*NLANES+k]
							+L1weights[2*3*NLANES+k]*L1preds[2*3*NLANES+k]
							+L1weights[3*3*NLANES+k]*L1preds[3*3*NLANES+k]
							+L1weights[4*3*NLANES+k]*L1preds[4*3*NLANES+k]
							+L1weights[5*3*NLANES+k]*L1preds[5*3*NLANES+k]
							+L1weights[6*3*NLANES+k]*L1preds[6*3*NLANES+k]
							+L1weights[7*3*NLANES+k]*L1preds[7*3*NLANES+k]
						)>>L1_SH2;
					}
					if(!cond_cg)//loosen pred range
					{
						for(int k=0;k<3*NLANES;++k)
						{
							int NE=rows[1][k+1*NROWS*NVAL], NEEE=rows[1][k+3*NROWS*NVAL];
							if(vmin[k]>NE)vmin[k]=NE;
							if(vmax[k]<NE)vmax[k]=NE;
							if(vmin[k]>NEEE)vmin[k]=NEEE;
							if(vmax[k]<NEEE)vmax[k]=NEEE;
						}
					}
					memcpy(pred0, pred, sizeof(pred0));
				}
				else if(effort==3)
				{
					for(int k=0;k<3*NLANES;++k)
					{
						int N=rows[1][k+0*NROWS*NVAL], W=rows[0][k-1*NROWS*NVAL];
						
						/*
						effort 3	20 preds
						0x00	N+W-NW
						0x01	N
						0x02	W
						0x03	W+NE-N
						0x04	3*(N-NN)+NNN
						0x05	3*(W-WW)+WWW
						0x06	N+NE-NNE
						0x07	NEE
						0x08	NN
						0x09	WW
						0x0A	2*N-NN
						0x0B	2*W-WW
						0x0C	NEEE
						0x0D	NEEEE
						0x0E	NNWW
						0x0F	NNEE
						0x10	N+NW-NNW
						0x11	W+NW-NWW
						0x12	(WWWW+NEEEE)>>1
						0x13	(WWW+NNN+NEEE-NW)>>1
						*/
						L1preds[0x00*3*NLANES+k]=pred[k];
						L1preds[0x01*3*NLANES+k]=N;
						L1preds[0x02*3*NLANES+k]=W;
						L1preds[0x03*3*NLANES+k]=W+rows[1][k+1*NROWS*NVAL]-N;
						L1preds[0x04*3*NLANES+k]=3*(N-rows[2][k+0*NROWS*NVAL])+rows[3][k+0*NROWS*NVAL];
						L1preds[0x05*3*NLANES+k]=3*(W-rows[0][k-2*NROWS*NVAL])+rows[0][k-3*NROWS*NVAL];
						L1preds[0x06*3*NLANES+k]=N+rows[1][k+1*NROWS*NVAL]-rows[2][k+1*NROWS*NVAL];
						L1preds[0x07*3*NLANES+k]=rows[1][k+2*NROWS*NVAL];
						L1preds[0x08*3*NLANES+k]=rows[2][k+0*NROWS*NVAL];
						L1preds[0x09*3*NLANES+k]=rows[0][k-2*NROWS*NVAL];
						L1preds[0x0A*3*NLANES+k]=2*N-rows[2][k+0*NROWS*NVAL];
						L1preds[0x0B*3*NLANES+k]=2*W-rows[0][k-2*NROWS*NVAL];
						L1preds[0x0C*3*NLANES+k]=rows[1][k+3*NROWS*NVAL];
						L1preds[0x0D*3*NLANES+k]=rows[1][k+4*NROWS*NVAL];
						L1preds[0x0E*3*NLANES+k]=rows[2][k-2*NROWS*NVAL];
						L1preds[0x0F*3*NLANES+k]=rows[2][k+2*NROWS*NVAL];
						L1preds[0x10*3*NLANES+k]=N+rows[1][k-1*NROWS*NVAL]-rows[2][k-1*NROWS*NVAL];
						L1preds[0x11*3*NLANES+k]=W+rows[1][k-1*NROWS*NVAL]-rows[1][k-2*NROWS*NVAL];
						L1preds[0x12*3*NLANES+k]=(rows[0][k-4*NROWS*NVAL]+rows[1][k+4*NROWS*NVAL])>>1;
						L1preds[0x13*3*NLANES+k]=(rows[0][k-3*NROWS*NVAL]+rows[3][k+0*NROWS*NVAL]+rows[1][k+3*NROWS*NVAL]-rows[1][k-1*NROWS*NVAL])>>1;
					}
					for(int k=0;k<3*NLANES;++k)
					{
						pred[k]=((1<<L1_SH3>>1)
							+L1weights[0x00*3*NLANES+k]*L1preds[0x00*3*NLANES+k]
							+L1weights[0x01*3*NLANES+k]*L1preds[0x01*3*NLANES+k]
							+L1weights[0x02*3*NLANES+k]*L1preds[0x02*3*NLANES+k]
							+L1weights[0x03*3*NLANES+k]*L1preds[0x03*3*NLANES+k]
							+L1weights[0x04*3*NLANES+k]*L1preds[0x04*3*NLANES+k]
							+L1weights[0x05*3*NLANES+k]*L1preds[0x05*3*NLANES+k]
							+L1weights[0x06*3*NLANES+k]*L1preds[0x06*3*NLANES+k]
							+L1weights[0x07*3*NLANES+k]*L1preds[0x07*3*NLANES+k]
							+L1weights[0x08*3*NLANES+k]*L1preds[0x08*3*NLANES+k]
							+L1weights[0x09*3*NLANES+k]*L1preds[0x09*3*NLANES+k]
							+L1weights[0x0A*3*NLANES+k]*L1preds[0x0A*3*NLANES+k]
							+L1weights[0x0B*3*NLANES+k]*L1preds[0x0B*3*NLANES+k]
							+L1weights[0x0C*3*NLANES+k]*L1preds[0x0C*3*NLANES+k]
							+L1weights[0x0D*3*NLANES+k]*L1preds[0x0D*3*NLANES+k]
							+L1weights[0x0E*3*NLANES+k]*L1preds[0x0E*3*NLANES+k]
							+L1weights[0x0F*3*NLANES+k]*L1preds[0x0F*3*NLANES+k]
							+L1weights[0x10*3*NLANES+k]*L1preds[0x10*3*NLANES+k]
							+L1weights[0x11*3*NLANES+k]*L1preds[0x11*3*NLANES+k]
							+L1weights[0x12*3*NLANES+k]*L1preds[0x12*3*NLANES+k]
							+L1weights[0x13*3*NLANES+k]*L1preds[0x13*3*NLANES+k]
						)>>L1_SH3;
					}
					if(!cond_cg)//loosen pred range
					{
						for(int k=0;k<3*NLANES;++k)
						{
							int NE=rows[1][k+1*NROWS*NVAL], NEEE=rows[1][k+3*NROWS*NVAL];
							if(vmin[k]>NE)vmin[k]=NE;
							if(vmax[k]<NE)vmax[k]=NE;
							if(vmin[k]>NEEE)vmin[k]=NEEE;
							if(vmax[k]<NEEE)vmax[k]=NEEE;
						}
					}
					memcpy(pred0, pred, sizeof(pred0));
				}
				if(cond_cg)
					memcpy(pred, grad, sizeof(pred0));
			effort0:
				for(int k=0;k<3*NLANES;++k)
					CLAMP2(pred[k], vmin[k], vmax[k]);
			}
			if(fwd)
			{
				for(int k=0;k<NLANES;++k)
				{
					myuv[k+0*NLANES]=imptr[k+yidx]-128;
					myuv[k+1*NLANES]=imptr[k+uidx]-128;
					myuv[k+2*NLANES]=imptr[k+vidx]-128;
				}

				//decorrelate Y
				for(int k=0;k<NLANES;++k)
					msyms[k]=myuv[k]-pred[k];
				if(dist>1)
				{
#ifdef ENABLE_GUIDE
					int16_t yuv0[NLANES];
					memcpy(yuv0, myuv+0*NLANES, sizeof(yuv0));
#endif
					for(int k=0;k<NLANES;++k)
					{
						int sym, recon;

						sym=msyms[k+0*NLANES];
						msyms[k+0*NLANES]=sym=(sym*dist_rcp>>16)-(sym>>15);
						recon=sym*dist+pred[k+0*NLANES];
						CLAMP2(recon, -128, 127);
						myuv[k+0*NLANES]=recon;
#ifdef ENABLE_GUIDE
						g_image[k+yidx]=myuv[k+0*NLANES]+128;
						g_sqe[0]+=abs(yuv0[k]-recon);
#endif
					}
				}
				for(int k=0;k<NLANES;++k)
					((uint16_t*)ctxptr)[k+0*NLANES]=syms[k+0*NLANES]=ctx[k+0*NLANES]<<8|(uint8_t)(msyms[k+0*NLANES]+128);

				//decorrelate U
				for(int k=0;k<NLANES;++k)
				{
					moffset[k+0*NLANES]=myuv[k+0*NLANES]&uhelpmask;
					pred[k+1*NLANES]+=moffset[k];
					CLAMP2(pred[k+1*NLANES], -128, 127);
					msyms[k+1*NLANES]=myuv[k+1*NLANES]-pred[k+1*NLANES];
				}
				if(dist>1)
				{
#ifdef ENABLE_GUIDE
					int16_t yuv0[NLANES];
					memcpy(yuv0, myuv+1*NLANES, sizeof(yuv0));
#endif
					for(int k=0;k<NLANES;++k)
					{
						int sym, recon;

						sym=msyms[k+1*NLANES];
						msyms[k+1*NLANES]=sym=(sym*dist_rcp>>16)-(sym>>15);
						recon=sym*dist+pred[k+1*NLANES];
						CLAMP2(recon, -128, 127);
						myuv[k+1*NLANES]=recon;
#ifdef ENABLE_GUIDE
						g_image[k+uidx]=myuv[k+1*NLANES]+128;
						g_sqe[1]+=abs(yuv0[k]-recon);
#endif
					}
				}
				for(int k=0;k<NLANES;++k)
					((uint16_t*)ctxptr)[k+1*NLANES]=syms[k+1*NLANES]=ctx[k+1*NLANES]<<8|(uint8_t)(msyms[k+1*NLANES]+128);
				
				//decorrelate V
				for(int k=0;k<NLANES;++k)
				{
					moffset[k+1*NLANES]=(vc0*myuv[k+0*NLANES]+vc1*myuv[k+1*NLANES])>>2;
					pred[k+2*NLANES]+=moffset[k+1*NLANES];
					CLAMP2(pred[k+2*NLANES], -128, 127);
					msyms[k+2*NLANES]=myuv[k+2*NLANES]-pred[k+2*NLANES];
				}
				if(dist>1)
				{
#ifdef ENABLE_GUIDE
					int16_t yuv0[NLANES];
					memcpy(yuv0, myuv+2*NLANES, sizeof(yuv0));
#endif
					for(int k=0;k<NLANES;++k)
					{
						int sym, recon;

						sym=msyms[k+2*NLANES];
						msyms[k+2*NLANES]=sym=(sym*dist_rcp>>16)-(sym>>15);
						recon=sym*dist+pred[k+1*NLANES];
						CLAMP2(recon, -128, 127);
						myuv[k+2*NLANES]=recon;
#ifdef ENABLE_GUIDE
						g_image[k+vidx]=myuv[k+2*NLANES]+128;
						g_sqe[2]+=abs(yuv0[k]-recon);
#endif
					}
				}
				for(int k=0;k<NLANES;++k)
					((uint16_t*)ctxptr)[k+2*NLANES]=syms[k+2*NLANES]=ctx[k+2*NLANES]<<8|(uint8_t)(msyms[k+2*NLANES]+128);
				{
					int *pa, *pb, *pc, va, vb, vc;
					pa=hists+syms[0*NLANES+0x0]; pb=hists+syms[1*NLANES+0x0]; pc=hists+syms[2*NLANES+0x0]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NLANES+0x1]; pb=hists+syms[1*NLANES+0x1]; pc=hists+syms[2*NLANES+0x1]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NLANES+0x2]; pb=hists+syms[1*NLANES+0x2]; pc=hists+syms[2*NLANES+0x2]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NLANES+0x3]; pb=hists+syms[1*NLANES+0x3]; pc=hists+syms[2*NLANES+0x3]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NLANES+0x4]; pb=hists+syms[1*NLANES+0x4]; pc=hists+syms[2*NLANES+0x4]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NLANES+0x5]; pb=hists+syms[1*NLANES+0x5]; pc=hists+syms[2*NLANES+0x5]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NLANES+0x6]; pb=hists+syms[1*NLANES+0x6]; pc=hists+syms[2*NLANES+0x6]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NLANES+0x7]; pb=hists+syms[1*NLANES+0x7]; pc=hists+syms[2*NLANES+0x7]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NLANES+0x8]; pb=hists+syms[1*NLANES+0x8]; pc=hists+syms[2*NLANES+0x8]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NLANES+0x9]; pb=hists+syms[1*NLANES+0x9]; pc=hists+syms[2*NLANES+0x9]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NLANES+0xA]; pb=hists+syms[1*NLANES+0xA]; pc=hists+syms[2*NLANES+0xA]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NLANES+0xB]; pb=hists+syms[1*NLANES+0xB]; pc=hists+syms[2*NLANES+0xB]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NLANES+0xC]; pb=hists+syms[1*NLANES+0xC]; pc=hists+syms[2*NLANES+0xC]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NLANES+0xD]; pb=hists+syms[1*NLANES+0xD]; pc=hists+syms[2*NLANES+0xD]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NLANES+0xE]; pb=hists+syms[1*NLANES+0xE]; pc=hists+syms[2*NLANES+0xE]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NLANES+0xF]; pb=hists+syms[1*NLANES+0xF]; pc=hists+syms[2*NLANES+0xF]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
				}
				ctxptr+=sizeof(int16_t[3][NLANES]);
			}
			else
			{
				//decode main

				//yuv = (char)(sym+pred-128)	= (uint8_t)(sym+pred)-128
				dec_yuv(mstate, (uint16_t*)ctx+0*NLANES, CDF2syms, &streamptr, streamend, (uint16_t*)myuv+0*NLANES);//residuals from [0 ~ 255]
				dec_yuv(mstate, (uint16_t*)ctx+1*NLANES, CDF2syms, &streamptr, streamend, (uint16_t*)myuv+1*NLANES);
				dec_yuv(mstate, (uint16_t*)ctx+2*NLANES, CDF2syms, &streamptr, streamend, (uint16_t*)myuv+2*NLANES);
				
				//reconstruct Y
				if(dist>1)
				{
					for(int k=0;k<NLANES;++k)
					{
						int p=(myuv[k+0*NLANES]-128)*dist+pred[k+0*NLANES];
						CLAMP2(p, -128, 127);
						myuv[k+0*NLANES]=p;
					}
				}
				else
				{
					for(int k=0;k<NLANES;++k)
						myuv[k+0*NLANES]=(uint8_t)(myuv[k+0*NLANES]+pred[k+0*NLANES])-128;
				}
				for(int k=0;k<NLANES;++k)
					imptr[k+yidx]=myuv[k+0*NLANES]+128;
#ifdef ENABLE_GUIDE
				if(memcmp(imptr+yidx, g_image+(imptr-interleaved)+yidx, NLANES))
				{
					printf("original  decoded  original-decoded  XYC0 %d %d %d\n", kx, ky, yidx);
					for(int k=0;k<NLANES;++k)
						printf("0x%02X  0x%02X  %4d\n",
							g_image[imptr-interleaved+yidx+k],
							imptr[yidx+k],
							g_image[imptr-interleaved+yidx+k]-imptr[yidx+k]
						);
					CRASH("guide error XYC0 %d %d %d/%d", kx, ky, yidx, NLANES);
				}
#endif

				//reconstruct U
				for(int k=0;k<NLANES;++k)
				{
					moffset[k+0*NLANES]=myuv[k+0*NLANES]&uhelpmask;
					pred[k+1*NLANES]+=moffset[k+0*NLANES];
					CLAMP2(pred[k+1*NLANES], -128, 127);
				}
				if(dist>1)
				{
					for(int k=0;k<NLANES;++k)
					{
						int p=(myuv[k+1*NLANES]-128)*dist+pred[k+1*NLANES];
						CLAMP2(p, -128, 127);
						myuv[k+1*NLANES]=p;
					}
				}
				else
				{
					for(int k=0;k<NLANES;++k)
						myuv[k+1*NLANES]=(uint8_t)(myuv[k+1*NLANES]+pred[k+1*NLANES])-128;
				}
				for(int k=0;k<NLANES;++k)
					imptr[k+uidx]=myuv[k+1*NLANES]+128;
#ifdef ENABLE_GUIDE
				if(memcmp(imptr+uidx, g_image+(imptr-interleaved)+uidx, NLANES))
				{
					printf("original  decoded  original-decoded  XYC1 %d %d %d\n", kx, ky, uidx);
					for(int k=0;k<NLANES;++k)
						printf("0x%02X  0x%02X  %4d\n",
							g_image[imptr-interleaved+uidx+k],
							imptr[uidx+k],
							g_image[imptr-interleaved+uidx+k]-imptr[uidx+k]
						);
					CRASH("guide error XYC1 %d %d %d/%d", kx, ky, uidx, NLANES);
				}
#endif
				
				//reconstruct V
				for(int k=0;k<NLANES;++k)
				{
					moffset[k+1*NLANES]=(vc0*myuv[k+0*NLANES]+vc1*myuv[k+1*NLANES])>>2;
					pred[k+2*NLANES]+=moffset[k+1*NLANES];
					CLAMP2(pred[k+2*NLANES], -128, 127);
				}
				if(dist>1)
				{
					for(int k=0;k<NLANES;++k)
					{
						int p=(myuv[k+2*NLANES]-128)*dist+pred[k+2*NLANES];
						CLAMP2(p, -128, 127);
						myuv[k+2*NLANES]=p;
					}
				}
				else
				{
					for(int k=0;k<NLANES;++k)
						myuv[k+2*NLANES]=(uint8_t)(myuv[k+2*NLANES]+pred[k+2*NLANES])-128;
				}
				for(int k=0;k<NLANES;++k)
					imptr[k+vidx]=myuv[k+2*NLANES]+128;
#ifdef ENABLE_GUIDE
				if(memcmp(imptr+vidx, g_image+(imptr-interleaved)+vidx, NLANES))
				{
					printf("original  decoded  original-decoded  XYC2 %d %d %d\n", kx, ky, vidx);
					for(int k=0;k<NLANES;++k)
						printf("0x%02X  0x%02X  %4d\n",
							g_image[imptr-interleaved+vidx+k],
							imptr[vidx+k],
							g_image[imptr-interleaved+vidx+k]-imptr[vidx+k]
						);
					CRASH("guide error XYC2 %d %d %d/%d", kx, ky, vidx, NLANES);
				}
#endif
			}
			for(int k=0;k<NLANES;++k)
			{
				int s0, s1, s2;

				s0=myuv[k+0*NLANES]-pred[k+0*NLANES];
				s1=myuv[k+1*NLANES]-pred[k+1*NLANES];
				s2=myuv[k+2*NLANES]-pred[k+2*NLANES];
				rows[0][k+0*NROWS*NVAL+0*NLANES]=myuv[k+0*NLANES];
				rows[0][k+0*NROWS*NVAL+1*NLANES]=myuv[k+1*NLANES]-moffset[k+0*NLANES];
				rows[0][k+0*NROWS*NVAL+2*NLANES]=myuv[k+2*NLANES]-moffset[k+1*NLANES];
				s0=(((s0<<1^s0>>31)<<GRBITS)+(rows[1][k+2*NROWS*NVAL+3*NLANES]>rows[1][k+3*NROWS*NVAL+3*NLANES]?rows[1][k+2*NROWS*NVAL+3*NLANES]:rows[1][k+3*NROWS*NVAL+3*NLANES])+1)>>1;
				s1=(((s1<<1^s1>>31)<<GRBITS)+(rows[1][k+2*NROWS*NVAL+4*NLANES]>rows[1][k+3*NROWS*NVAL+4*NLANES]?rows[1][k+2*NROWS*NVAL+4*NLANES]:rows[1][k+3*NROWS*NVAL+4*NLANES])+1)>>1;
				s2=(((s2<<1^s2>>31)<<GRBITS)+(rows[1][k+2*NROWS*NVAL+5*NLANES]>rows[1][k+3*NROWS*NVAL+5*NLANES]?rows[1][k+2*NROWS*NVAL+5*NLANES]:rows[1][k+3*NROWS*NVAL+5*NLANES])+1)>>1;
				rows[0][k+0*NROWS*NVAL+3*NLANES]=(rows[0][k-1*NROWS*NVAL+3*NLANES]+s0+1)>>1;
				rows[0][k+0*NROWS*NVAL+4*NLANES]=(rows[0][k-1*NROWS*NVAL+4*NLANES]+s1+1)>>1;
				rows[0][k+0*NROWS*NVAL+5*NLANES]=(rows[0][k-1*NROWS*NVAL+5*NLANES]+s2+1)>>1;
			}
			if(effort)
			{
				for(int k=0;k<3*NLANES;++k)
				{
					int t0, t1;

					t0=rows[0][k+0*NROWS*NVAL];
					t1=pred0[k];
					msyms[k]=(t0>t1)-(t0<t1);
				}
				if(effort==1)//update
				{
					for(int kp=0;kp<L1_NPREDS1;++kp)//update
					{
						for(int k=0;k<3*NLANES;++k)
							L1weights[3*NLANES*kp+k]+=msyms[k]*L1preds[3*NLANES*kp+k];
					}
					for(int k=0;k<3*NLANES;++k)
						cW[k]=rows[0][k+0*NROWS*NVAL]-pred0[k];
				}
				else if(effort==2)//update
				{
					for(int kp=0;kp<L1_NPREDS2;++kp)//update
					{
						for(int k=0;k<3*NLANES;++k)
							L1weights[3*NLANES*kp+k]+=msyms[k]*L1preds[3*NLANES*kp+k];
					}
				}
				else if(effort==3)//update
				{
					for(int kp=0;kp<L1_NPREDS3;++kp)//update
					{
						for(int k=0;k<3*NLANES;++k)
							L1weights[3*NLANES*kp+k]+=msyms[k]*L1preds[3*NLANES*kp+k];
					}
				}
			}
			rows[0]+=NROWS*NVAL;
			rows[1]+=NROWS*NVAL;
			rows[2]+=NROWS*NVAL;
			rows[3]+=NROWS*NVAL;
			imptr+=3*NLANES;
		}
		//printf("%8d/%8d\r", ky+1, blockh);
	}
	prof_checkpoint(isize, "main");
#ifdef ENABLE_GUIDE
	if(dist>1&&fwd)
	{
		double rmse[]=
		{
			sqrt(g_sqe[0]*3/isize),
			sqrt(g_sqe[1]*3/isize),
			sqrt(g_sqe[2]*3/isize),
			sqrt((g_sqe[0]+g_sqe[1]+g_sqe[2])/isize),
		};
		double psnr[]=
		{
			20*log10(255/rmse[0]),
			20*log10(255/rmse[1]),
			20*log10(255/rmse[2]),
			20*log10(255/rmse[3]),
		};
		printf("RMSE PSNR\n");
		printf("T %12.6lf %12.6lf\n", rmse[3], psnr[3]);
		printf("Y %12.6lf %12.6lf\n", rmse[0], psnr[0]);
		printf("U %12.6lf %12.6lf\n", rmse[1], psnr[1]);
		printf("V %12.6lf %12.6lf\n", rmse[2], psnr[2]);
	}
#endif

	if(effort)
		free(L1state);
	if(fwd)//all rANS encoding is bwd-bwd
	{
		rANS_SIMD_SymInfo *syminfo=(rANS_SIMD_SymInfo*)CDF2syms;
		rANS_SIMD_SymInfo *rsyminfo=(rANS_SIMD_SymInfo*)CDF2syms+(ptrdiff_t)3*NCTX*256;
		int32_t *rhist=hists+(ptrdiff_t)3*NCTX*256;
		
		if(xremw||yremh)
		{
			for(int ky=0;ky<yremh;++ky)
				decorr1d(image+rowstride*(blockh*YCODERS+ky), iw, 3, bestrct, rhist);
			for(int kx=0;kx<xremw;++kx)
				decorr1d(image+qxbytes+3*kx, blockh*YCODERS, rowstride, bestrct, rhist);
		}

		//normalize/integrate hists
		for(int kc=0;kc<nctx;++kc)
			enc_hist2stats(hists+(ptrdiff_t)256*kc, syminfo+(ptrdiff_t)256*kc, &bypassmask, kc, 0, 1);
			
		if(xremw||yremh)//encode remainder
		{
			uint32_t state=1<<(RANS_STATE_BITS-RANS_RENORM_BITS);
			for(int kx=xremw-1;kx>=0;--kx)
				encode1d_port(image+qxbytes+3*kx, blockh*YCODERS, rowstride, &state, &streamptr, image, rsyminfo);
			for(int ky=yremh-1;ky>=0;--ky)
				encode1d_port(image+rowstride*(blockh*YCODERS+ky), iw, 3, &state, &streamptr, image, rsyminfo);
			//flush
			streamptr-=4;
#ifdef _DEBUG
			if(streamptr<=image)
				CRASH("OOB ptr %016zX <= %016zX", streamptr, image);
#endif
			*(uint32_t*)streamptr=state;
			prof_checkpoint(usize-isize, "encode remainder");
			profile_size(streamptr, "/ %9td bytes remainders", usize-isize);
		}

		//encode main
		uint16_t *ctxptr2=(uint16_t*)(interleaved+(isize<<1));
		for(int k=0;k<NLANES;++k)
			mstate[k]=1<<(RANS_STATE_BITS-RANS_RENORM_BITS);
		for(int ky=blockh-1;ky>=0;--ky)
		{
#ifdef ESTIMATE_SIZE
			int kc=2;
#endif
			for(int kx=3*blockw-1;kx>=0;--kx)//blockw = iw/XCODERS
			{
				rANS_SIMD_SymInfo p[NLANES];

				ctxptr2-=NLANES;
#ifdef ANS_VAL
				for(int k=0;k<NLANES;++k)
				{
					if((uint32_t)((ctxptr2[k]>>8)-kc*NCTX)>=NCTX)
						CRASH("XYC %d %d %d", kx, ky, kc);
				}
#endif
				p[0x0]=syminfo[ctxptr2[0x0]];
				p[0x1]=syminfo[ctxptr2[0x1]];
				p[0x2]=syminfo[ctxptr2[0x2]];
				p[0x3]=syminfo[ctxptr2[0x3]];
				p[0x4]=syminfo[ctxptr2[0x4]];
				p[0x5]=syminfo[ctxptr2[0x5]];
				p[0x6]=syminfo[ctxptr2[0x6]];
				p[0x7]=syminfo[ctxptr2[0x7]];
				p[0x8]=syminfo[ctxptr2[0x8]];
				p[0x9]=syminfo[ctxptr2[0x9]];
				p[0xA]=syminfo[ctxptr2[0xA]];
				p[0xB]=syminfo[ctxptr2[0xB]];
				p[0xC]=syminfo[ctxptr2[0xC]];
				p[0xD]=syminfo[ctxptr2[0xD]];
				p[0xE]=syminfo[ctxptr2[0xE]];
				p[0xF]=syminfo[ctxptr2[0xF]];
#ifdef ESTIMATE_SIZE
				{
					const double norm=1./(1<<PROBBITS);
					for(int k=0;k<NLANES;++k)
					{
						if((uint32_t)((ctxptr2[k]>>8)-kc*NCTX)>=NCTX)
							CRASH("XYC %d %d %d", kx, ky, kc);
						int freq=(1<<PROBBITS)-(p[k].negf&0xFFFF);
						if((uint32_t)(freq-1)>=(uint32_t)((1<<PROBBITS)-1))
							CRASH("freq = %d", freq);
						esize[kc*NLANES+k]-=log2(freq*norm)*0.125;
					}
				}
				--kc;
				if(kc<0)
					kc=2;
				if(streamptr-sizeof(uint16_t[NLANES])<image)
				{
					double ctotal=0;
					for(int k=0;k<3*NLANES;++k)
						ctotal+=esize[k];
					CRASH("inflation XYC %d %d %d  %8.4lf%%  %12.2lf", kx, ky, kc, 100.*blockh/(blockh-1-(ky+1)), ctotal);
				}
#endif
				for(int k=NLANES-1;k>=0;--k)//enc renorm		if(state>(freq<<(31-12))-1){write16(state); state>>=16;}
				{
					if(mstate[k]>p[k].smax)
					{
						streamptr-=sizeof(uint16_t);
						*(uint16_t*)streamptr=(uint16_t)mstate[k];
						mstate[k]>>=16;
					}
				}
#ifdef ANS_VAL
				ansval_push(mstate, sizeof(uint32_t), NLANES);
#endif
				for(int k=0;k<NLANES;++k)//enc update		state += (state*invf>>sh)*negf+cdf		state = state/freq<<12|(cdf+state%freq)
					mstate[k]+=(uint32_t)((uint64_t)mstate[k]*p[k].invf>>p[k].sh)*p[k].negf+p[k].cdf;
#ifdef ANS_VAL
				ansval_push(mstate, sizeof(uint32_t), NLANES);
#endif
			}
		}
		//flush
		streamptr-=sizeof(mstate);
#ifdef _DEBUG
		if(streamptr<=image)
			CRASH("OOB ptr %016zX <= %016zX", streamptr, image);
#endif
		memcpy(streamptr, mstate, sizeof(mstate));
		prof_checkpoint(isize, "encode main");
		profile_size(streamptr, "/ %9td bytes main", isize);

		//pack hists
		{
			BitPackerLIFO ec;
			bitpacker_enc_init(&ec, image, streamptr);
			for(int kc=nctx-1;kc>=0;--kc)
				enc_packhist(&ec, hists+(ptrdiff_t)256*kc, bypassmask, kc);
			bitpacker_enc_flush(&ec);
			streamptr=ec.dstbwdptr;
		}
		prof_checkpoint((ptrdiff_t)nctx*256, "pack histograms");
		profile_size(streamptr, "/ %9d bytes overhead", nctx*PROBBITS<<8>>3);

		//save compressed file
		{
			ptrdiff_t csize2=0;
			FILE *fdst=fopen(dstfn, "wb");
			if(!fdst)
			{
				CRASH("Cannot open \"%s\" for writing", fdst);
				return 1;
			}
			csize2+=fwrite("L1", 1, 2, fdst);
			csize2+=fwrite(&iw, 1, 3, fdst);
			csize2+=fwrite(&ih, 1, 3, fdst);
			{
				int flags=bestrct<<2|(effort&3);
				csize2+=fwrite(&flags, 1, 1, fdst);
			}
			csize2+=fwrite(&dist, 1, 1, fdst);
			csize2+=fwrite(&bypassmask, 1, 8, fdst);
#ifdef _DEBUG
			if(streamptr>streamstart)
				CRASH("OOB ptr %016zX > %016zX", streamptr, streamstart);
			if(streamptr<image)
				CRASH("OOB ptr %016zX < %016zX", streamptr, image);
#endif
			csize2+=fwrite(streamptr, 1, streamstart-streamptr, fdst);
			fclose(fdst);
			
#ifdef ESTIMATE_SIZE
			double etotal=0;
			for(int k=0;k<3*NLANES;++k)
			{
				etotal+=esize[k];
				printf("E %12.2lf\n", esize[k]);//
			}
			printf("Total estimate  %12.2lf bytes\n", etotal);
#endif
#ifdef LOUD
			printf("L1C Port WH %d*%d  RCT %2d %s  effort %d  dist %3d  \"%s\"\n"
				, iw
				, ih
				, bestrct
				, rct_names[bestrct]
				, effort
				, dist
				, srcfn
			);
			printf("%8td/%8td bytes\n", csize2, usize2);
#endif
			prof_checkpoint(csize2, "fwrite");
			(void)csize2;
		}
		free(hists);
	}
	else
	{
		//deinterleave
		interleave_blocks_inv(interleaved, iw, ih, image);
		prof_checkpoint(usize, "deinterleave");

		if(xremw||yremh)
		{
			uint32_t *rCDF2syms=(uint32_t*)CDF2syms+((ptrdiff_t)3*NCTX<<PROBBITS);
#ifdef _DEBUG
			if(streamptr>streamend)
				CRASH("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
			uint32_t state=*(uint32_t*)streamptr;
			streamptr+=4;
			for(int ky=0;ky<yremh;++ky)
				decode1d(image+rowstride*(blockh*YCODERS+ky), iw, 3, bestrct, &state, (const uint8_t**)&streamptr, streamend, rCDF2syms);
			for(int kx=0;kx<xremw;++kx)
				decode1d(image+qxbytes+3*kx, blockh*YCODERS, rowstride, bestrct, &state, (const uint8_t**)&streamptr, streamend, rCDF2syms);
			prof_checkpoint(usize-isize, "remainder");
		}

		//save PPM file
		save_ppm(dstfn, image, iw, ih);
		prof_checkpoint(usize, "fwrite");
	}
	free(pixels);
	free(CDF2syms);
	free(interleaved);
	free(image);
	free(ctable);

#ifdef LOUD
	t=time_sec()-t;
	printf("%c  %12.6lf sec  %12.6lf MB/s\n", 'D'+fwd, t, usize/(t*1024*1024));
#endif
#ifdef PROFILE_TIME
#ifdef __GNUC__
	if(profile)
#endif
		prof_print(usize);
#endif
	(void)och_names;
	(void)rct_names;
	(void)print_timestamp;
	(void)encode1d_sse41;
	(void)encode1d;
	return 0;
}
