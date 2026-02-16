#ifdef _MSC_VER
#	ifndef _CRT_SECURE_NO_WARNINGS
#		define _CRT_SECURE_NO_WARNINGS
#	endif
#elif defined __linux__ && !defined _GNU_SOURCE
#	define _GNU_SOURCE
#	include<stddef.h>//ptrdiff_t
#endif
#include<stdint.h>
#include<math.h>
#include<immintrin.h>
#include<sys/stat.h>


	#define PROFILE_TIME		//should be on

#ifdef _MSC_VER
	#define PROFILE_SIZE
	#define LOUD			//size & time

//	#define ESTIMATE_SIZE		//DEBUG		checks for zero frequency, visualizes context usage

//	#define SAVE_RESIDUALS
//	#define TEST_INTERLEAVE
#endif

//	#define ANALYSIS_GRAD
//	#define USE_L2			//bad
	#define ENABLE_RCT_EXTENSION
	#define INTERLEAVESIMD		//tailored for 64 lanes		2.5x faster interleave
//	#define EMULATE_GATHER		//gather is a little faster


enum
{
	XCODERS=8,
	YCODERS=8,
	NCODERS=XCODERS*YCODERS,

	ANALYSIS_XSTRIDE=2,
	ANALYSIS_YSTRIDE=2,

	DEFAULT_EFFORT_LEVEL=2,
	L1_NPREDS1=4,
	L1_NPREDS2=8,
	L1_NPREDS3=20,
	L1_SH1=15,	//L1SH1 <= 16
	L1_SH2=17,	//L1SH2 >= 16
	L1_SH3=19,	//L1SH3 >= 16

	GRBITS=3,
	NCTX=18,	//18*3+3 = 57 total

	PROBBITS=12,	//12 bit max	James Bonfield's CDF2sym: {freq<<20 | bias<<8 | sym}
	RANS_STATE_BITS=31,
	RANS_RENORM_BITS=16,

	XPAD=8,
	NROWS=4,
	NCH=3,
	NVAL=2,
	NREG=sizeof(int16_t[NCODERS])/sizeof(__m512i),
};

#define COMMON_rANS
#include"common.h"
AWM_INLINE void gather32(int *dst, const int *src, const int *offsets)
{
#ifdef EMULATE_GATHER
	volatile int *ptr=dst;
	ptr[0x0]=src[offsets[0x0]];
	ptr[0x1]=src[offsets[0x1]];
	ptr[0x2]=src[offsets[0x2]];
	ptr[0x3]=src[offsets[0x3]];
	ptr[0x4]=src[offsets[0x4]];
	ptr[0x5]=src[offsets[0x5]];
	ptr[0x6]=src[offsets[0x6]];
	ptr[0x7]=src[offsets[0x7]];
	ptr[0x8]=src[offsets[0x8]];
	ptr[0x9]=src[offsets[0x9]];
	ptr[0xA]=src[offsets[0xA]];
	ptr[0xB]=src[offsets[0xB]];
	ptr[0xC]=src[offsets[0xC]];
	ptr[0xD]=src[offsets[0xD]];
	ptr[0xE]=src[offsets[0xE]];
	ptr[0xF]=src[offsets[0xF]];
#else
	_mm512_store_si512((__m512i*)dst, _mm512_i32gather_epi32(_mm512_load_si512((__m512i*)offsets), src, sizeof(int)));
#endif
}
AWM_INLINE void dec_yuv(
	__m512i *mstate,
	const __m512i *ctx0,
	const uint32_t *CDF2syms,
	uint8_t **pstreamptr,
	const uint8_t *streamend,
	__m512i *syms
)
{
	const uint8_t *streamptr=*pstreamptr;
	__m512i decctx[4];
	decctx[0]=_mm512_slli_epi32(ctx0[0], 16);
	decctx[1]=_mm512_slli_epi32(ctx0[1], 16);
	decctx[2]=_mm512_srli_epi32(ctx0[0], 16);
	decctx[3]=_mm512_srli_epi32(ctx0[1], 16);
	decctx[0]=_mm512_srli_epi32(decctx[0], 16-PROBBITS);
	decctx[1]=_mm512_srli_epi32(decctx[1], 16-PROBBITS);
	decctx[2]=_mm512_slli_epi32(decctx[2], PROBBITS);
	decctx[3]=_mm512_slli_epi32(decctx[3], PROBBITS);
	{
		__m512i mprobmask=_mm512_set1_epi32((1<<PROBBITS)-1);
		__m512i rem0=_mm512_and_si512(mstate[0], mprobmask);
		__m512i rem1=_mm512_and_si512(mstate[1], mprobmask);
		__m512i rem2=_mm512_and_si512(mstate[2], mprobmask);
		__m512i rem3=_mm512_and_si512(mstate[3], mprobmask);
		decctx[0]=_mm512_or_si512(decctx[0], rem0);
		decctx[1]=_mm512_or_si512(decctx[1], rem1);
		decctx[2]=_mm512_or_si512(decctx[2], rem2);
		decctx[3]=_mm512_or_si512(decctx[3], rem3);
	}
#ifdef ANS_VAL
	ALIGN(64) int debugctx[NCODERS];
	memcpy(debugctx, decctx, sizeof(int[NCODERS]));
#endif
	gather32((int*)(decctx+0), (const int*)CDF2syms, (int*)(decctx+0));
	gather32((int*)(decctx+1), (const int*)CDF2syms, (int*)(decctx+1));
	gather32((int*)(decctx+2), (const int*)CDF2syms, (int*)(decctx+2));
	gather32((int*)(decctx+3), (const int*)CDF2syms, (int*)(decctx+3));

	//update		state = (state>>12)*freq+(rem-cdf)	rem-cdf is prebaked		{freq<<20 | bias<<8 | sym}
	{
		__m512i mfreq0=_mm512_srli_epi32(decctx[0], PROBBITS+8);//1 <= freq <= 0xF01
		__m512i mfreq1=_mm512_srli_epi32(decctx[1], PROBBITS+8);
		__m512i mfreq2=_mm512_srli_epi32(decctx[2], PROBBITS+8);
		__m512i mfreq3=_mm512_srli_epi32(decctx[3], PROBBITS+8);
#ifdef ANS_VAL
		__m512i mdebugfreq[2];
		mdebugfreq[0]=_mm512_or_si512(mfreq0, _mm512_slli_epi32(mfreq2, 16));
		mdebugfreq[1]=_mm512_or_si512(mfreq1, _mm512_slli_epi32(mfreq3, 16));
		ansval_check(mdebugfreq, sizeof(int16_t), NCODERS);
#endif
		mstate[0]=_mm512_srli_epi32(mstate[0], PROBBITS);
		mstate[1]=_mm512_srli_epi32(mstate[1], PROBBITS);
		mstate[2]=_mm512_srli_epi32(mstate[2], PROBBITS);
		mstate[3]=_mm512_srli_epi32(mstate[3], PROBBITS);
		mstate[0]=_mm512_mullo_epi32(mstate[0], mfreq0);//10 cycles
		mstate[1]=_mm512_mullo_epi32(mstate[1], mfreq1);
		mstate[2]=_mm512_mullo_epi32(mstate[2], mfreq2);
		mstate[3]=_mm512_mullo_epi32(mstate[3], mfreq3);
	}
	{
		__m512i mbias0=_mm512_slli_epi32(decctx[0], PROBBITS);
		__m512i mbias1=_mm512_slli_epi32(decctx[1], PROBBITS);
		__m512i mbias2=_mm512_slli_epi32(decctx[2], PROBBITS);
		__m512i mbias3=_mm512_slli_epi32(decctx[3], PROBBITS);
		mbias0=_mm512_srli_epi32(mbias0, 32-PROBBITS);
		mbias1=_mm512_srli_epi32(mbias1, 32-PROBBITS);
		mbias2=_mm512_srli_epi32(mbias2, 32-PROBBITS);
		mbias3=_mm512_srli_epi32(mbias3, 32-PROBBITS);
		mstate[0]=_mm512_add_epi32(mstate[0], mbias0);
		mstate[1]=_mm512_add_epi32(mstate[1], mbias1);
		mstate[2]=_mm512_add_epi32(mstate[2], mbias2);
		mstate[3]=_mm512_add_epi32(mstate[3], mbias3);
	}
	syms[0]=_mm512_mask_blend_epi16(0xAAAAAAAA, decctx[0], _mm512_slli_epi32(decctx[2], 16));
	syms[1]=_mm512_mask_blend_epi16(0xAAAAAAAA, decctx[1], _mm512_slli_epi32(decctx[3], 16));
	syms[0]=_mm512_maskz_mov_epi8(0x5555555555555555, syms[0]);
	syms[1]=_mm512_maskz_mov_epi8(0x5555555555555555, syms[1]);
#ifdef ANS_VAL
	ansval_check(mstate, sizeof(int), NCODERS);
#endif
	//renorm
	{
		__m512i smin=_mm512_set1_epi32(1<<(RANS_STATE_BITS-RANS_RENORM_BITS));
#ifdef _DEBUG
		if(streamptr>streamend)
			CRASH("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
		__mmask16 mask0=_mm512_cmpgt_epu32_mask(smin, mstate[0]);
		__mmask16 mask1=_mm512_cmpgt_epu32_mask(smin, mstate[1]);
		__mmask16 mask2=_mm512_cmpgt_epu32_mask(smin, mstate[2]);
		__mmask16 mask3=_mm512_cmpgt_epu32_mask(smin, mstate[3]);
		int step0=_mm_popcnt_u32(mask0);
		int step1=_mm_popcnt_u32(mask1);
		int step2=_mm_popcnt_u32(mask2);
		int step3=_mm_popcnt_u32(mask3);
		__m512i lo0, lo1, lo2, lo3;
		__m256i tlo0=_mm256_loadu_si256((const __m256i*)streamptr); streamptr+=step0*sizeof(int16_t);
		__m256i tlo1=_mm256_loadu_si256((const __m256i*)streamptr); streamptr+=step1*sizeof(int16_t);
		__m256i tlo2=_mm256_loadu_si256((const __m256i*)streamptr); streamptr+=step2*sizeof(int16_t);
		__m256i tlo3=_mm256_loadu_si256((const __m256i*)streamptr); streamptr+=step3*sizeof(int16_t);
#ifdef _DEBUG
		if(streamptr>streamend)
			CRASH("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
		tlo0=_mm256_maskz_expand_epi16(mask0, tlo0);
		tlo1=_mm256_maskz_expand_epi16(mask1, tlo1);
		tlo2=_mm256_maskz_expand_epi16(mask2, tlo2);
		tlo3=_mm256_maskz_expand_epi16(mask3, tlo3);
		lo0=_mm512_cvtepu16_epi32(tlo0);
		lo1=_mm512_cvtepu16_epi32(tlo1);
		lo2=_mm512_cvtepu16_epi32(tlo2);
		lo3=_mm512_cvtepu16_epi32(tlo3);
		mstate[0]=_mm512_mask_slli_epi32(mstate[0], mask0, mstate[0], 16);
		mstate[1]=_mm512_mask_slli_epi32(mstate[1], mask1, mstate[1], 16);
		mstate[2]=_mm512_mask_slli_epi32(mstate[2], mask2, mstate[2], 16);
		mstate[3]=_mm512_mask_slli_epi32(mstate[3], mask3, mstate[3], 16);
		mstate[0]=_mm512_add_epi32(mstate[0], lo0);
		mstate[1]=_mm512_add_epi32(mstate[1], lo1);
		mstate[2]=_mm512_add_epi32(mstate[2], lo2);
		mstate[3]=_mm512_add_epi32(mstate[3], lo3);
#ifdef ANS_VAL
		ansval_check(mstate, sizeof(int), NCODERS);
#endif
	}
	*pstreamptr=(uint8_t*)(size_t)streamptr;
}

AWM_INLINE void transpose8(__m512i *data)
{
#if 1
	__m512i swap1=_mm512_set_epi8(
	//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0,
		14, 15, 12, 13, 10, 11,  8,  9,  6,  7,  4,  5,  2,  3,  0,  1,
		14, 15, 12, 13, 10, 11,  8,  9,  6,  7,  4,  5,  2,  3,  0,  1,
		14, 15, 12, 13, 10, 11,  8,  9,  6,  7,  4,  5,  2,  3,  0,  1,
		14, 15, 12, 13, 10, 11,  8,  9,  6,  7,  4,  5,  2,  3,  0,  1
	);
	__m512i swap2=_mm512_set_epi8(
	//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0,
		13, 12, 15, 14,  9,  8, 11, 10,  5,  4,  7,  6,  1,  0,  3,  2,
		13, 12, 15, 14,  9,  8, 11, 10,  5,  4,  7,  6,  1,  0,  3,  2,
		13, 12, 15, 14,  9,  8, 11, 10,  5,  4,  7,  6,  1,  0,  3,  2,
		13, 12, 15, 14,  9,  8, 11, 10,  5,  4,  7,  6,  1,  0,  3,  2
	);
	__m512i swap4=_mm512_set_epi8(
	//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0,
		11, 10,  9,  8, 15, 14, 13, 12,  3,  2,  1,  0,  7,  6,  5,  4,
		11, 10,  9,  8, 15, 14, 13, 12,  3,  2,  1,  0,  7,  6,  5,  4,
		11, 10,  9,  8, 15, 14, 13, 12,  3,  2,  1,  0,  7,  6,  5,  4,
		11, 10,  9,  8, 15, 14, 13, 12,  3,  2,  1,  0,  7,  6,  5,  4
	);
	//[01 23 45 67]
	__m512i t0=data[0];
	__m512i t2=data[2];
	__m512i t4=data[4];
	__m512i t6=data[6];
	__m512i t1=_mm512_shuffle_epi8(data[1], swap1);
	__m512i t3=_mm512_shuffle_epi8(data[3], swap1);
	__m512i t5=_mm512_shuffle_epi8(data[5], swap1);
	__m512i t7=_mm512_shuffle_epi8(data[7], swap1);
	__m512i s0=_mm512_mask_mov_epi8(t0, 0xAAAAAAAAAAAAAAAA, t1);
	__m512i s1=_mm512_mask_mov_epi8(t0, 0x5555555555555555, t1);
	__m512i s2=_mm512_mask_mov_epi8(t2, 0xAAAAAAAAAAAAAAAA, t3);
	__m512i s3=_mm512_mask_mov_epi8(t2, 0x5555555555555555, t3);
	__m512i s4=_mm512_mask_mov_epi8(t4, 0xAAAAAAAAAAAAAAAA, t5);
	__m512i s5=_mm512_mask_mov_epi8(t4, 0x5555555555555555, t5);
	__m512i s6=_mm512_mask_mov_epi8(t6, 0xAAAAAAAAAAAAAAAA, t7);
	__m512i s7=_mm512_mask_mov_epi8(t6, 0x5555555555555555, t7);
	s1=_mm512_shuffle_epi8(s1, swap1);
	s3=_mm512_shuffle_epi8(s3, swap1);
	s5=_mm512_shuffle_epi8(s5, swap1);
	s7=_mm512_shuffle_epi8(s7, swap1);

	//[02 46] [13 57]
	s2=_mm512_shuffle_epi8(s2, swap2);
	s6=_mm512_shuffle_epi8(s6, swap2);
	s3=_mm512_shuffle_epi8(s3, swap2);
	s7=_mm512_shuffle_epi8(s7, swap2);
	t0=_mm512_mask_mov_epi16(s0, 0xAAAAAAAA, s2);
	t2=_mm512_mask_mov_epi16(s0, 0x55555555, s2);
	t4=_mm512_mask_mov_epi16(s4, 0xAAAAAAAA, s6);
	t6=_mm512_mask_mov_epi16(s4, 0x55555555, s6);
	t1=_mm512_mask_mov_epi16(s1, 0xAAAAAAAA, s3);
	t3=_mm512_mask_mov_epi16(s1, 0x55555555, s3);
	t5=_mm512_mask_mov_epi16(s5, 0xAAAAAAAA, s7);
	t7=_mm512_mask_mov_epi16(s5, 0x55555555, s7);
	t2=_mm512_shuffle_epi8(t2, swap2);
	t6=_mm512_shuffle_epi8(t6, swap2);
	t3=_mm512_shuffle_epi8(t3, swap2);
	t7=_mm512_shuffle_epi8(t7, swap2);

	//[04] [15] [26] [37]
	t4=_mm512_shuffle_epi8(t4, swap4);
	t5=_mm512_shuffle_epi8(t5, swap4);
	t6=_mm512_shuffle_epi8(t6, swap4);
	t7=_mm512_shuffle_epi8(t7, swap4);
	s0=_mm512_mask_mov_epi32(t0, 0xAAAA, t4);
	s4=_mm512_mask_mov_epi32(t0, 0x5555, t4);
	s1=_mm512_mask_mov_epi32(t1, 0xAAAA, t5);
	s5=_mm512_mask_mov_epi32(t1, 0x5555, t5);
	s2=_mm512_mask_mov_epi32(t2, 0xAAAA, t6);
	s6=_mm512_mask_mov_epi32(t2, 0x5555, t6);
	s3=_mm512_mask_mov_epi32(t3, 0xAAAA, t7);
	s7=_mm512_mask_mov_epi32(t3, 0x5555, t7);
	s4=_mm512_shuffle_epi8(s4, swap4);
	s5=_mm512_shuffle_epi8(s5, swap4);
	s6=_mm512_shuffle_epi8(s6, swap4);
	s7=_mm512_shuffle_epi8(s7, swap4);
	data[0]=s0;
	data[1]=s1;
	data[2]=s2;
	data[3]=s3;
	data[4]=s4;
	data[5]=s5;
	data[6]=s6;
	data[7]=s7;
#else
	uint8_t d2[8][8][8];
	for(int kr=0;kr<8;++kr)
	{
		for(int kl=0;kl<8;++kl)
		{
			for(int kx=0;kx<8;++kx)
				d2[kx][kl][kr]=((uint8_t*)data)[kr<<6|kl<<3|kx];
		}
	}
	memcpy(data, d2, sizeof(d2));
#endif
}
static void interleave_blocks_fwd(const uint8_t *original, int iw, int ih, uint8_t *interleaved)
{
#if 0
	{
		ALIGN(64) uint8_t d1[64*8], d1t[64*8], d0[64*8];
		for(int ky=0;ky<8;++ky)
		{
			for(int kl=0;kl<8;++kl)
			{
				for(int kx=0;kx<8;++kx)
				{
					int val=rand();
					d1[ky<<6|kl<<3|kx]=(uint8_t)val;
					d1t[kx<<6|ky<<3|kl]=(uint8_t)val;
				}
			}
		}
		memcpy(d0, d1, sizeof(d0));
		transpose_bytes((__m512i*)d1);
		transpose_lanes((__m512i*)d1);
		transpose_qwords((__m512i*)d1);
		if(memcmp(d1, d1t, sizeof(d1)))
			CRASH("");
		transpose_qwords((__m512i*)d1);
		transpose_lanes((__m512i*)d1);
		transpose_bytes((__m512i*)d1);
		//permute_fwd((__m512i*)d1);
		//permute_inv((__m512i*)d1);
		if(memcmp(d1, d0, sizeof(d0)))
			CRASH("");
	}
#endif
	//original[ih][iw][3]
	//interleaved[ih/YCODERS][iw/XCODERS][3][NCODERS]	contiguous & aligned
	//xrem[ih%YCODERS][iw][3]
	//yrem[ih-iw%YCODERS][iw%XCODERS][3]

	//only difference between fwd and inv:		swap assignments (const slow->fast)
	int rowstride=3*iw;
	int ixyblockw=iw/XCODERS;
	int ixyblockh=ih/YCODERS;
	int blockxbytes=3*NCODERS*ixyblockw;
#ifdef INTERLEAVESIMD
	int SIMDxcount=blockxbytes&~((int)sizeof(char[8][NCODERS])-1);
	__m512i slowinc=_mm512_set1_epi64(sizeof(char[8]));
#endif
	uint8_t *fastptr=interleaved;
	ALIGN(64) const uint8_t *slowptrs[NCODERS]={0}, *slowptrs0[NCODERS]={0};
	for(int ky=0;ky<YCODERS;++ky)//spread slow pointers
	{
		for(int kx=0;kx<XCODERS;++kx)
		//	slowptrs0[XCODERS*ky+kx]=original+3*(iw*ixyblockh*ky+ixyblockw*kx);
			slowptrs0[XCODERS*kx+ky]=original+3*(iw*ixyblockh*ky+ixyblockw*kx);//transposed pointers
	}
	for(int ky=0;ky<ixyblockh;++ky)//interleave
	{
		int kx=0;
		memcpy((void*)slowptrs, slowptrs0, sizeof(slowptrs));
#ifdef INTERLEAVESIMD
		__m512i mslow[8];
		mslow[0]=_mm512_load_si512((__m512i*)slowptrs+0);
		mslow[1]=_mm512_load_si512((__m512i*)slowptrs+1);
		mslow[2]=_mm512_load_si512((__m512i*)slowptrs+2);
		mslow[3]=_mm512_load_si512((__m512i*)slowptrs+3);
		mslow[4]=_mm512_load_si512((__m512i*)slowptrs+4);
		mslow[5]=_mm512_load_si512((__m512i*)slowptrs+5);
		mslow[6]=_mm512_load_si512((__m512i*)slowptrs+6);
		mslow[7]=_mm512_load_si512((__m512i*)slowptrs+7);
		for(;kx<SIMDxcount;kx+=(int)sizeof(char[8][NCODERS]))
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

			slowptrs[NCODERS]		fastptr  (aligned because NCODERS == sizeof(__m256i))
			A16x16 B16x16		->	A16x16T C16x16T B16x16T D16x16T
			C16x16 D16x16

			load, transpose, store blocks A, B
			load, transpose, store blocks C, D

			speed:
			SIMD maximum estimate	(3/load + 3/store)*32+80*0.5*2 = 272 cycles for 32*32 bytes  ->  3GHz*32*32/272 ~= 10770 MB/s without remainders, cache misses
			actual SIMD		MSVC 4900 MB/s		GCC 3900/3000 MB/s
			actual serial		1900 MB/s
			*/
			__m512i mdata[8];
			mdata[0]=_mm512_i64gather_epi64(mslow[0], 0, 1);
			mdata[1]=_mm512_i64gather_epi64(mslow[1], 0, 1);
			mdata[2]=_mm512_i64gather_epi64(mslow[2], 0, 1);
			mdata[3]=_mm512_i64gather_epi64(mslow[3], 0, 1);
			mdata[4]=_mm512_i64gather_epi64(mslow[4], 0, 1);
			mdata[5]=_mm512_i64gather_epi64(mslow[5], 0, 1);
			mdata[6]=_mm512_i64gather_epi64(mslow[6], 0, 1);
			mdata[7]=_mm512_i64gather_epi64(mslow[7], 0, 1);
			//transpose_lanes(mdata);
			//transpose_qwords(mdata);
			//transpose_bytes(mdata);
			transpose8(mdata);
			_mm512_store_si512((__m512i*)fastptr+0, mdata[0]);
			_mm512_store_si512((__m512i*)fastptr+1, mdata[1]);
			_mm512_store_si512((__m512i*)fastptr+2, mdata[2]);
			_mm512_store_si512((__m512i*)fastptr+3, mdata[3]);
			_mm512_store_si512((__m512i*)fastptr+4, mdata[4]);
			_mm512_store_si512((__m512i*)fastptr+5, mdata[5]);
			_mm512_store_si512((__m512i*)fastptr+6, mdata[6]);
			_mm512_store_si512((__m512i*)fastptr+7, mdata[7]);

			fastptr+=sizeof(char[8][NCODERS]);
			mslow[0]=_mm512_add_epi64(mslow[0], slowinc);
			mslow[1]=_mm512_add_epi64(mslow[1], slowinc);
			mslow[2]=_mm512_add_epi64(mslow[2], slowinc);
			mslow[3]=_mm512_add_epi64(mslow[3], slowinc);
			mslow[4]=_mm512_add_epi64(mslow[4], slowinc);
			mslow[5]=_mm512_add_epi64(mslow[5], slowinc);
			mslow[6]=_mm512_add_epi64(mslow[6], slowinc);
			mslow[7]=_mm512_add_epi64(mslow[7], slowinc);
		}
		_mm512_store_si512((__m512i*)slowptrs+0, mslow[0]);
		_mm512_store_si512((__m512i*)slowptrs+1, mslow[1]);
		_mm512_store_si512((__m512i*)slowptrs+2, mslow[2]);
		_mm512_store_si512((__m512i*)slowptrs+3, mslow[3]);
		_mm512_store_si512((__m512i*)slowptrs+4, mslow[4]);
		_mm512_store_si512((__m512i*)slowptrs+5, mslow[5]);
		_mm512_store_si512((__m512i*)slowptrs+6, mslow[6]);
		_mm512_store_si512((__m512i*)slowptrs+7, mslow[7]);
#else
		(void)transpose8;
#endif
		for(;kx<blockxbytes;kx+=NCODERS)
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
			//for(int k=0;k<NCODERS;++k)
			//	*fastptr++=*slowptrs[k]++;
			for(int k1=0;k1<8;++k1)//transposed pointers
				for(int k2=0;k2<8;++k2)
					*fastptr++=*slowptrs[k2<<3|k1]++;
		}
		for(int k=0;k<NCODERS;++k)
			slowptrs0[k]+=rowstride;
	}
}
static void interleave_blocks_inv(const uint8_t *interleaved, int iw, int ih, uint8_t *original)
{
	//original[ih][iw][3]
	//interleaved[ih/YCODERS][iw/XCODERS][3][NCODERS]	contiguous & aligned
	//xrem[ih%YCODERS][iw][3]
	//yrem[ih-iw%YCODERS][iw%XCODERS][3]

	//only difference between fwd and inv:		swap assignments (const slow->fast)
	int rowstride=3*iw;
	int ixyblockw=iw/XCODERS;
	int ixyblockh=ih/YCODERS;
	int blockxbytes=3*NCODERS*ixyblockw;
#ifdef INTERLEAVESIMD
	int SIMDxcount=blockxbytes&~((int)sizeof(char[8][NCODERS])-1);
	__m512i slowinc=_mm512_set1_epi64(sizeof(char[8]));
#endif
	const uint8_t *fastptr=interleaved;
	ALIGN(64) uint8_t *slowptrs[NCODERS]={0}, *slowptrs0[NCODERS]={0};
	for(int ky=0;ky<YCODERS;++ky)//spread slow pointers
	{
		for(int kx=0;kx<XCODERS;++kx)
		//	slowptrs0[XCODERS*ky+kx]=original+3*(iw*ixyblockh*ky+ixyblockw*kx);
			slowptrs0[XCODERS*kx+ky]=original+3*(iw*ixyblockh*ky+ixyblockw*kx);//transposed pointers
	}
	for(int ky=0;ky<ixyblockh;++ky)//interleave
	{
		int kx=0;
		memcpy((void*)slowptrs, slowptrs0, sizeof(slowptrs));
#ifdef INTERLEAVESIMD
		__m512i mslow[8];
		mslow[0]=_mm512_load_si512((__m512i*)slowptrs+0);
		mslow[1]=_mm512_load_si512((__m512i*)slowptrs+1);
		mslow[2]=_mm512_load_si512((__m512i*)slowptrs+2);
		mslow[3]=_mm512_load_si512((__m512i*)slowptrs+3);
		mslow[4]=_mm512_load_si512((__m512i*)slowptrs+4);
		mslow[5]=_mm512_load_si512((__m512i*)slowptrs+5);
		mslow[6]=_mm512_load_si512((__m512i*)slowptrs+6);
		mslow[7]=_mm512_load_si512((__m512i*)slowptrs+7);
		for(;kx<SIMDxcount;kx+=(int)sizeof(char[8][NCODERS]))
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

			slowptrs[NCODERS]		fastptr  (aligned because NCODERS == sizeof(__m256i))
			A16x16 B16x16		->	A16x16T C16x16T B16x16T D16x16T
			C16x16 D16x16

			load, transpose, store blocks A, B
			load, transpose, store blocks C, D

			speed:
			SIMD maximum estimate	(3/load + 3/store)*32+80*0.5*2 = 272 cycles for 32*32 bytes  ->  3GHz*32*32/272 ~= 10770 MB/s without remainders, cache misses
			actual SIMD		MSVC 4900 MB/s		GCC 3900/3000 MB/s
			actual serial		1900 MB/s
			*/
			__m512i mdata[8];
			mdata[0]=_mm512_load_si512((__m512i*)fastptr+0);
			mdata[1]=_mm512_load_si512((__m512i*)fastptr+1);
			mdata[2]=_mm512_load_si512((__m512i*)fastptr+2);
			mdata[3]=_mm512_load_si512((__m512i*)fastptr+3);
			mdata[4]=_mm512_load_si512((__m512i*)fastptr+4);
			mdata[5]=_mm512_load_si512((__m512i*)fastptr+5);
			mdata[6]=_mm512_load_si512((__m512i*)fastptr+6);
			mdata[7]=_mm512_load_si512((__m512i*)fastptr+7);
			//transpose_bytes(mdata);
			//transpose_qwords(mdata);
			//transpose_lanes(mdata);
			transpose8(mdata);
			_mm512_i64scatter_epi64(0, mslow[0], mdata[0], 1);
			_mm512_i64scatter_epi64(0, mslow[1], mdata[1], 1);
			_mm512_i64scatter_epi64(0, mslow[2], mdata[2], 1);
			_mm512_i64scatter_epi64(0, mslow[3], mdata[3], 1);
			_mm512_i64scatter_epi64(0, mslow[4], mdata[4], 1);
			_mm512_i64scatter_epi64(0, mslow[5], mdata[5], 1);
			_mm512_i64scatter_epi64(0, mslow[6], mdata[6], 1);
			_mm512_i64scatter_epi64(0, mslow[7], mdata[7], 1);

			fastptr+=sizeof(char[8][NCODERS]);
			mslow[0]=_mm512_add_epi64(mslow[0], slowinc);
			mslow[1]=_mm512_add_epi64(mslow[1], slowinc);
			mslow[2]=_mm512_add_epi64(mslow[2], slowinc);
			mslow[3]=_mm512_add_epi64(mslow[3], slowinc);
			mslow[4]=_mm512_add_epi64(mslow[4], slowinc);
			mslow[5]=_mm512_add_epi64(mslow[5], slowinc);
			mslow[6]=_mm512_add_epi64(mslow[6], slowinc);
			mslow[7]=_mm512_add_epi64(mslow[7], slowinc);
		}
		_mm512_store_si512((__m512i*)slowptrs+0, mslow[0]);
		_mm512_store_si512((__m512i*)slowptrs+1, mslow[1]);
		_mm512_store_si512((__m512i*)slowptrs+2, mslow[2]);
		_mm512_store_si512((__m512i*)slowptrs+3, mslow[3]);
		_mm512_store_si512((__m512i*)slowptrs+4, mslow[4]);
		_mm512_store_si512((__m512i*)slowptrs+5, mslow[5]);
		_mm512_store_si512((__m512i*)slowptrs+6, mslow[6]);
		_mm512_store_si512((__m512i*)slowptrs+7, mslow[7]);
#else
		(void)transpose8;
#endif
#if 1
		for(;kx<blockxbytes;kx+=NCODERS)
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
			//for(int k=0;k<NCODERS;++k)
			//	*slowptrs[k]++=*fastptr++;
			for(int k1=0;k1<8;++k1)//transposed pointers
				for(int k2=0;k2<8;++k2)
					*slowptrs[k2<<3|k1]++=*fastptr++;
		}
#endif
		for(int k=0;k<NCODERS;++k)
			slowptrs0[k]+=rowstride;
	}
}
int codec_l1_avx512(int argc, char **argv)
{
	if(argc!=3&&argc!=4&&argc!=5)
	{
		printf(
			"Usage: \"%s\"  input  output  [Effort]  [Dist]    Encode/decode.\n"
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
	double esize[3*NCODERS]={0};
#endif
#ifdef LOUD
	ptrdiff_t usize2=0;
	double t=time_sec();
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
	uint64_t bypassmask=0;//0: emit stats	1: rare context (bypass)
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
		int tag=0;
		fread(&tag, 1, 2, fsrc);
		fwd=tag==('P'|'6'<<8);
		if(!fwd&&tag!=('V'|'1'<<8))
		{
			CRASH("Unsupported file \"%s\"", srcfn);
			return 1;
		}
		if(fwd)
		{
#ifdef LOUD
			print_timestamp("%Y-%m-%d_%H%M%S\n");
#endif
			int temp=fgetc(fsrc);
			if(temp!='\n')
			{
				CRASH("Invalid PPM file");
				return 1;
			}
			int nread=fscanf(fsrc, "%d %d", &iw, &ih);
			if(nread!=2)
			{
				CRASH("Unsupported PPM file");
				return 1;
			}
			int vmax=0;
			nread=fscanf(fsrc, "%d", &vmax);
			if(nread!=1||vmax!=255)
			{
				CRASH("Unsupported PPM file");
				return 1;
			}
			temp=fgetc(fsrc);
			if(temp!='\n')
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
	//	image=(uint8_t*)_mm_malloc(cap+sizeof(__m512i), 0x1000);
		image=(uint8_t*)malloc(cap+sizeof(__m512i));
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
			streamptr=streamstart=image+cap-(csize-cheadersize)-sizeof(__m512i);
			streamend=image+cap-sizeof(__m512i);
			fread(streamstart, 1, csize-cheadersize, fsrc);//read stream
		}
		fclose(fsrc);
	}
	prof_checkpoint(fwd?usize:csize, "fread");
	int blockw=iw/XCODERS;
	int blockh=ih/YCODERS;
	int qxbytes=blockw*XCODERS*3;//iw/XCODERS*XCODERS*3
	int ixcount=blockw*NCODERS, ixbytes=3*ixcount;//ix = interleaved circular buffer width		iw/XCODERS*NCODERS
	int xremw=iw-blockw*XCODERS, yremh=ih-blockh*YCODERS;
	int xrembytes=3*xremw;
	int nctx=3*NCTX+3*(xremw||yremh);
	ptrdiff_t isize=(ptrdiff_t)ixbytes*blockh;
	ptrdiff_t interleavedsize=isize<<fwd;//fwd ? interleave residuals & context : pack residuals
	uint8_t *interleaved=(uint8_t*)_mm_malloc(interleavedsize, sizeof(__m512i));
	if(!interleaved)
	{
		CRASH("Alloc error");
		return 1;
	}
	(void)xrembytes;
	int hsize=nctx*(int)sizeof(int[256]);
	int *hists=fwd?(int*)malloc(hsize):0;//fwd-only

	int CDF2syms_size=nctx*(int)sizeof(int32_t[1<<PROBBITS]);
	if(fwd)//DIV-free rANS encoder reuses this as SIMD symbol info
		CDF2syms_size=nctx*(int)sizeof(rANS_SIMD_SymInfo[256]);
	uint32_t *CDF2syms=(uint32_t*)_mm_malloc(CDF2syms_size, sizeof(__m512i));

	psize=(blockw+2*XPAD)*(int)sizeof(int16_t[NCH*NROWS*NVAL*NCODERS]);//int16_t[blockw+2*XPAD][NCH*NROWS*NVAL*NCODERS]
	pixels=(int16_t*)_mm_malloc(psize, sizeof(__m512i));
	if((fwd&&!hists)||!CDF2syms||!pixels)
	{
		CRASH("Alloc error");
		return 1;
	}
	if(fwd)
	{
		memset(hists, 0, hsize);
#ifdef TEST_INTERLEAVE
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
#ifdef ANALYSIS_GRAD
			const int ystart=1;
#else
			const int ystart=0;
#endif
			int64_t counters[OCH_COUNT]={0};
			__m512i mcounters[OCH_COUNT];//64-bit
			__m256i half8=_mm256_set1_epi8(-128);
			__m512i wordmask=_mm512_set1_epi64(0xFFFF);
			memset(mcounters, 0, sizeof(mcounters));
			imptr=interleaved+isize;
			for(int ky=ystart;ky<blockh;ky+=ANALYSIS_YSTRIDE)
			{
				__m512i prev[OCH_COUNT][2];//16-bit
				memset(prev, 0, sizeof(prev));
				for(int kx=0;kx<blockw;kx+=ANALYSIS_XSTRIDE)
				{
					__m512i r0=_mm512_cvtepi8_epi16(_mm256_add_epi8(_mm256_load_si256((__m256i*)imptr+0), half8));
					__m512i r1=_mm512_cvtepi8_epi16(_mm256_add_epi8(_mm256_load_si256((__m256i*)imptr+1), half8));
					__m512i g0=_mm512_cvtepi8_epi16(_mm256_add_epi8(_mm256_load_si256((__m256i*)imptr+2), half8));
					__m512i g1=_mm512_cvtepi8_epi16(_mm256_add_epi8(_mm256_load_si256((__m256i*)imptr+3), half8));
					__m512i b0=_mm512_cvtepi8_epi16(_mm256_add_epi8(_mm256_load_si256((__m256i*)imptr+4), half8));
					__m512i b1=_mm512_cvtepi8_epi16(_mm256_add_epi8(_mm256_load_si256((__m256i*)imptr+5), half8));
#ifdef ANALYSIS_GRAD
					__m512i rN0=_mm512_cvtepi8_epi16(_mm256_add_epi8(_mm256_load_si256((__m256i*)(imptr-ixbytes)+0), half8));
					__m512i rN1=_mm512_cvtepi8_epi16(_mm256_add_epi8(_mm256_load_si256((__m256i*)(imptr-ixbytes)+1), half8));
					__m512i gN0=_mm512_cvtepi8_epi16(_mm256_add_epi8(_mm256_load_si256((__m256i*)(imptr-ixbytes)+2), half8));
					__m512i gN1=_mm512_cvtepi8_epi16(_mm256_add_epi8(_mm256_load_si256((__m256i*)(imptr-ixbytes)+3), half8));
					__m512i bN0=_mm512_cvtepi8_epi16(_mm256_add_epi8(_mm256_load_si256((__m256i*)(imptr-ixbytes)+4), half8));
					__m512i bN1=_mm512_cvtepi8_epi16(_mm256_add_epi8(_mm256_load_si256((__m256i*)(imptr-ixbytes)+5), half8));
					r0=_mm512_sub_epi16(r0, rN0);
					g0=_mm512_sub_epi16(g0, gN0);
					b0=_mm512_sub_epi16(b0, bN0);
					r1=_mm512_sub_epi16(r1, rN1);
					g1=_mm512_sub_epi16(g1, gN1);
					b1=_mm512_sub_epi16(b1, bN1);
#endif
					imptr+=3*NCODERS*ANALYSIS_XSTRIDE;
					r0=_mm512_slli_epi16(r0, 2);
					g0=_mm512_slli_epi16(g0, 2);
					b0=_mm512_slli_epi16(b0, 2);
					r1=_mm512_slli_epi16(r1, 2);
					g1=_mm512_slli_epi16(g1, 2);
					b1=_mm512_slli_epi16(b1, 2);
					__m512i rg0=_mm512_sub_epi16(r0, g0);
					__m512i gb0=_mm512_sub_epi16(g0, b0);
					__m512i br0=_mm512_sub_epi16(b0, r0);
					__m512i rg1=_mm512_sub_epi16(r1, g1);
					__m512i gb1=_mm512_sub_epi16(g1, b1);
					__m512i br1=_mm512_sub_epi16(b1, r1);
#define UPDATE(LANE, IDXA, A0, IDXB, B0, IDXC, C0)\
	do\
	{\
		__m512i t0=A0, t1=B0, t2=C0;\
		__m512i ta=_mm512_sub_epi16(t0, prev[IDXA][LANE]);\
		__m512i tb=_mm512_sub_epi16(t1, prev[IDXB][LANE]);\
		__m512i tc=_mm512_sub_epi16(t2, prev[IDXC][LANE]);\
		prev[IDXA][LANE]=t0;\
		prev[IDXB][LANE]=t1;\
		prev[IDXC][LANE]=t2;\
		ta=_mm512_abs_epi16(ta);\
		tb=_mm512_abs_epi16(tb);\
		tc=_mm512_abs_epi16(tc);\
		ta=_mm512_add_epi16(ta, _mm512_srli_epi64(ta, 32));\
		tb=_mm512_add_epi16(tb, _mm512_srli_epi64(tb, 32));\
		tc=_mm512_add_epi16(tc, _mm512_srli_epi64(tc, 32));\
		ta=_mm512_add_epi16(ta, _mm512_srli_epi64(ta, 16));\
		tb=_mm512_add_epi16(tb, _mm512_srli_epi64(tb, 16));\
		tc=_mm512_add_epi16(tc, _mm512_srli_epi64(tc, 16));\
		mcounters[IDXA]=_mm512_add_epi64(mcounters[IDXA], _mm512_and_si512(ta, wordmask));\
		mcounters[IDXB]=_mm512_add_epi64(mcounters[IDXB], _mm512_and_si512(tb, wordmask));\
		mcounters[IDXC]=_mm512_add_epi64(mcounters[IDXC], _mm512_and_si512(tc, wordmask));\
	}while(0)
					UPDATE(0, OCH_YX00, r0, OCH_Y0X0, g0, OCH_Y00X, b0);
					UPDATE(1, OCH_YX00, r1, OCH_Y0X0, g1, OCH_Y00X, b1);
					UPDATE(0, OCH_CX40, rg0, OCH_C0X4, gb0, OCH_C40X, br0);
					UPDATE(1, OCH_CX40, rg1, OCH_C0X4, gb1, OCH_C40X, br1);
#ifdef ENABLE_RCT_EXTENSION
					UPDATE(0
						, OCH_CX31, _mm512_add_epi16(rg0, _mm512_srai_epi16(gb0, 2))//r-(3*g+b)/4 = r-g-(b-g)/4
						, OCH_C3X1, _mm512_add_epi16(rg0, _mm512_srai_epi16(br0, 2))//g-(3*r+b)/4 = g-r-(b-r)/4
						, OCH_C31X, _mm512_add_epi16(br0, _mm512_srai_epi16(rg0, 2))//b-(3*r+g)/4 = b-r-(g-r)/4
					);
					UPDATE(1
						, OCH_CX31, _mm512_add_epi16(rg1, _mm512_srai_epi16(gb1, 2))
						, OCH_C3X1, _mm512_add_epi16(rg1, _mm512_srai_epi16(br1, 2))
						, OCH_C31X, _mm512_add_epi16(br1, _mm512_srai_epi16(rg1, 2))
					);
					UPDATE(0
						, OCH_CX13, _mm512_add_epi16(br0, _mm512_srai_epi16(gb0, 2))//r-(g+3*b)/4 = r-b-(g-b)/4
						, OCH_C1X3, _mm512_add_epi16(gb0, _mm512_srai_epi16(br0, 2))//g-(r+3*b)/4 = g-b-(r-b)/4
						, OCH_C13X, _mm512_add_epi16(gb0, _mm512_srai_epi16(rg0, 2))//b-(r+3*g)/4 = b-g-(r-g)/4
					);
					UPDATE(1
						, OCH_CX13, _mm512_add_epi16(br1, _mm512_srai_epi16(gb1, 2))
						, OCH_C1X3, _mm512_add_epi16(gb1, _mm512_srai_epi16(br1, 2))
						, OCH_C13X, _mm512_add_epi16(gb1, _mm512_srai_epi16(rg1, 2))
					);
					UPDATE(0
						, OCH_CX22, _mm512_srai_epi16(_mm512_sub_epi16(rg0, br0), 1)//r-(g+b)/2 = (r-g + r-b)/2
						, OCH_C2X2, _mm512_srai_epi16(_mm512_sub_epi16(gb0, rg0), 1)//g-(r+b)/2 = (g-r + g-b)/2
						, OCH_C22X, _mm512_srai_epi16(_mm512_sub_epi16(br0, gb0), 1)//b-(r+g)/2 = (b-r + b-g)/2
					);
					UPDATE(1
						, OCH_CX22, _mm512_srai_epi16(_mm512_sub_epi16(rg1, br1), 1)
						, OCH_C2X2, _mm512_srai_epi16(_mm512_sub_epi16(gb1, rg1), 1)
						, OCH_C22X, _mm512_srai_epi16(_mm512_sub_epi16(br1, gb1), 1)
					);
#endif
#undef  UPDATE
				}
				imptr+=ixbytes*(ANALYSIS_YSTRIDE-1);
			}
			for(int k=0;k<OCH_COUNT;++k)
			{
				ALIGN(64) int64_t temp[8]={0};
				_mm512_store_si512((__m512i*)temp, mcounters[k]);
				counters[k]=temp[0]+temp[1]+temp[2]+temp[3]+temp[4]+temp[5]+temp[6]+temp[7];
			}
			int64_t minerr=0;
			for(int kt=0;kt<RCT_COUNT;++kt)
			{
				const uint8_t *rct=rct_combinations[kt];
				int64_t currerr=
					+counters[rct[0]]
					+counters[rct[1]]
					+counters[rct[2]]
				;
#ifdef LOUD
				printf("%-14s %12lld + %12lld + %12lld = %12lld%s\n",
					rct_names[kt],
					counters[rct[0]],
					counters[rct[1]],
					counters[rct[2]],
					currerr,
					!kt||minerr>currerr?" <-":""
				);
#endif
				if(!kt||minerr>currerr)
				{
					minerr=currerr;
					bestrct=kt;
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
	case 1://use L1
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
		L1statesize=(int)sizeof(int[2*NCODERS*3*(L1_NPREDS3+1)]);//{preds, coeffs} * (NPREDS+{bias}) * 3 channels * NCODERS
		L1state=(int*)_mm_malloc(L1statesize, sizeof(__m512i));
		if(!L1state)
		{
			CRASH("Alloc error");
			return 1;
		}
		memset(L1state, 0, L1statesize);
	}
	const uint8_t *combination=rct_combinations[bestrct];
	int
		yidx=combination[II_PERM_Y]*NCODERS,
		uidx=combination[II_PERM_U]*NCODERS,
		vidx=combination[II_PERM_V]*NCODERS;
	__m512i uhelpmask=_mm512_set1_epi16(-(combination[II_COEFF_U_SUB_Y]!=0));
	__m512i vc0=_mm512_set1_epi16(combination[II_COEFF_V_SUB_Y]);
	__m512i vc1=_mm512_set1_epi16(combination[II_COEFF_V_SUB_U]);
	memset(pixels, 0, psize);
	__m512i mctxmax=_mm512_set1_epi16(NCTX-1);
	__m512i mctxuoffset=_mm512_set1_epi16(NCTX);
	__m512i mctxvoffset=_mm512_set1_epi16(NCTX*2);
	__m512i amin=_mm512_set1_epi16(-128);
	__m512i amax=_mm512_set1_epi16(127);
	__m256i half8=_mm256_set1_epi8(-128);
	__m512i wordmask=_mm512_set1_epi32(0xFFFF);
	__m512i myuv[6];
	__m512i dist_rcp=_mm512_set1_epi16(0x7FFF), mdist=_mm512_set1_epi16(1);
#ifdef SAVE_RESIDUALS
	uint8_t *residuals=0;
	if(fwd)
	{
		residuals=(uint8_t*)malloc(isize);
		if(!residuals)
		{
			CRASH("Alloc error");
			return 1;
		}
		memset(residuals, 0, isize);
	}
#endif
	if(dist>1)
	{
		dist_rcp=_mm512_set1_epi16(((1<<16)+dist-1)/dist);//x/dist  ->  {x*=inv; x=(x>>16)+((uint32_t)x>>31);}
		mdist=_mm512_set1_epi16(dist);
	}
	memset(myuv, 0, sizeof(myuv));
	uint8_t *ctxptr=interleaved;
	imptr=interleaved+(fwd?isize:0);
	__m512i mstate[4];
	__m512i *L1preds=effort?(__m512i*)L1state:0;
	int *L1weights=effort?(int*)(L1state+1*(ptrdiff_t)NCODERS*3*(L1_NPREDS3+1)):0;
	if(effort)
		FILLMEM(L1weights, (1<<sh)/npreds, (npreds+1)*sizeof(int[3*NCODERS]), sizeof(int));
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
		ALIGN(64) int16_t *rows[]=
		{
			//X=0: [row0: pY pY pU pU pV pV eY eY eU eU eV eV  row1: ...]  X=1: ...
			pixels+(XPAD*NROWS+(ky-0LL+NROWS)%NROWS)*NCODERS*NCH*NVAL,
			pixels+(XPAD*NROWS+(ky-1LL+NROWS)%NROWS)*NCODERS*NCH*NVAL,
			pixels+(XPAD*NROWS+(ky-2LL+NROWS)%NROWS)*NCODERS*NCH*NVAL,
			pixels+(XPAD*NROWS+(ky-3LL+NROWS)%NROWS)*NCODERS*NCH*NVAL,

			//X=0: [row0: pY pY eY eY pU pU eU eU pV pV eV eV  row1: ...]  X=1: ...
			//pixels+(XPAD*NROWS*NCH+(ky-0LL+NROWS)%NROWS)*NVAL*NCODERS,
			//pixels+(XPAD*NROWS*NCH+(ky-1LL+NROWS)%NROWS)*NVAL*NCODERS,
			//pixels+(XPAD*NROWS*NCH+(ky-2LL+NROWS)%NROWS)*NVAL*NCODERS,
			//pixels+(XPAD*NROWS*NCH+(ky-3LL+NROWS)%NROWS)*NVAL*NCODERS,
		};
		ALIGN(64) uint16_t syms[3*NCODERS]={0};
		__m512i NW[6], N[6], W[6];
		__m512i eW[6], ecurr[6], eNEE[6], eNEEE[6];
		memset(NW, 0, sizeof(NW));
		memset(N, 0, sizeof(N));
		memset(W, 0, sizeof(W));
		memset(eW, 0, sizeof(eW));
		memset(ecurr, 0, sizeof(ecurr));
		memset(eNEE, 0, sizeof(eNEE));
		memset(eNEEE, 0, sizeof(eNEEE));
		//(__m512i*)rows[-Y]+C+(E+X*NROWS*NVAL)*NREG*NCH
		eNEE[0]=_mm512_load_si512((__m512i*)rows[1]+0+(1+2*NROWS*NVAL)*NREG*NCH);
		eNEE[1]=_mm512_load_si512((__m512i*)rows[1]+1+(1+2*NROWS*NVAL)*NREG*NCH);
		eNEE[2]=_mm512_load_si512((__m512i*)rows[1]+2+(1+2*NROWS*NVAL)*NREG*NCH);
		eNEE[3]=_mm512_load_si512((__m512i*)rows[1]+3+(1+2*NROWS*NVAL)*NREG*NCH);
		eNEE[4]=_mm512_load_si512((__m512i*)rows[1]+4+(1+2*NROWS*NVAL)*NREG*NCH);
		eNEE[5]=_mm512_load_si512((__m512i*)rows[1]+5+(1+2*NROWS*NVAL)*NREG*NCH);
		for(int kx=0;kx<blockw;++kx)
		{
			__m512i
				pred[6], ctx[6];
			__m512i predYUV0[6];
			N[0]=_mm512_load_si512((__m512i*)rows[1]+0+(0+0*NROWS*NVAL)*NREG*NCH);//y neighbors
			N[1]=_mm512_load_si512((__m512i*)rows[1]+1+(0+0*NROWS*NVAL)*NREG*NCH);//y
			N[2]=_mm512_load_si512((__m512i*)rows[1]+2+(0+0*NROWS*NVAL)*NREG*NCH);//u
			N[3]=_mm512_load_si512((__m512i*)rows[1]+3+(0+0*NROWS*NVAL)*NREG*NCH);//u
			N[4]=_mm512_load_si512((__m512i*)rows[1]+4+(0+0*NROWS*NVAL)*NREG*NCH);//v
			N[5]=_mm512_load_si512((__m512i*)rows[1]+5+(0+0*NROWS*NVAL)*NREG*NCH);//v
			{
				//context = FLOOR_LOG2(eW*eW+1)
				__m512i one=_mm512_set1_epi32(1);
				__m512i cy0=_mm512_and_si512(eW[0], wordmask), cy1=_mm512_srli_epi32(eW[0], 16);
				__m512i cy2=_mm512_and_si512(eW[1], wordmask), cy3=_mm512_srli_epi32(eW[1], 16);
				__m512i cu0=_mm512_and_si512(eW[2], wordmask), cu1=_mm512_srli_epi32(eW[2], 16);
				__m512i cu2=_mm512_and_si512(eW[3], wordmask), cu3=_mm512_srli_epi32(eW[3], 16);
				__m512i cv0=_mm512_and_si512(eW[4], wordmask), cv1=_mm512_srli_epi32(eW[4], 16);
				__m512i cv2=_mm512_and_si512(eW[5], wordmask), cv3=_mm512_srli_epi32(eW[5], 16);
				cy0=_mm512_mullo_epi32(cy0, cy0);
				cy1=_mm512_mullo_epi32(cy1, cy1);
				cy2=_mm512_mullo_epi32(cy2, cy2);
				cy3=_mm512_mullo_epi32(cy3, cy3);
				cu0=_mm512_mullo_epi32(cu0, cu0);
				cu1=_mm512_mullo_epi32(cu1, cu1);
				cu2=_mm512_mullo_epi32(cu2, cu2);
				cu3=_mm512_mullo_epi32(cu3, cu3);
				cv0=_mm512_mullo_epi32(cv0, cv0);
				cv1=_mm512_mullo_epi32(cv1, cv1);
				cv2=_mm512_mullo_epi32(cv2, cv2);
				cv3=_mm512_mullo_epi32(cv3, cv3);
				cy0=_mm512_add_epi32(cy0, one);
				cy1=_mm512_add_epi32(cy1, one);
				cy2=_mm512_add_epi32(cy2, one);
				cy3=_mm512_add_epi32(cy3, one);
				cu0=_mm512_add_epi32(cu0, one);
				cu1=_mm512_add_epi32(cu1, one);
				cu2=_mm512_add_epi32(cu2, one);
				cu3=_mm512_add_epi32(cu3, one);
				cv0=_mm512_add_epi32(cv0, one);
				cv1=_mm512_add_epi32(cv1, one);
				cv2=_mm512_add_epi32(cv2, one);
				cv3=_mm512_add_epi32(cv3, one);
				cy0=_mm512_lzcnt_epi32(cy0);
				cy1=_mm512_lzcnt_epi32(cy1);
				cy2=_mm512_lzcnt_epi32(cy2);
				cy3=_mm512_lzcnt_epi32(cy3);
				cu0=_mm512_lzcnt_epi32(cu0);
				cu1=_mm512_lzcnt_epi32(cu1);
				cu2=_mm512_lzcnt_epi32(cu2);
				cu3=_mm512_lzcnt_epi32(cu3);
				cv0=_mm512_lzcnt_epi32(cv0);
				cv1=_mm512_lzcnt_epi32(cv1);
				cv2=_mm512_lzcnt_epi32(cv2);
				cv3=_mm512_lzcnt_epi32(cv3);
				cy1=_mm512_slli_epi32(cy1, 16);
				cu1=_mm512_slli_epi32(cu1, 16);
				cv1=_mm512_slli_epi32(cv1, 16);
				cy3=_mm512_slli_epi32(cy3, 16);
				cu3=_mm512_slli_epi32(cu3, 16);
				cv3=_mm512_slli_epi32(cv3, 16);
				ctx[0]=_mm512_or_si512(cy0, cy1);
				ctx[1]=_mm512_or_si512(cy2, cy3);
				ctx[2]=_mm512_or_si512(cu0, cu1);
				ctx[3]=_mm512_or_si512(cu2, cu3);
				ctx[4]=_mm512_or_si512(cv0, cv1);
				ctx[5]=_mm512_or_si512(cv2, cv3);
				{
					__m512i mask=_mm512_set1_epi16(31);
					ctx[0]=_mm512_sub_epi16(mask, ctx[0]);
					ctx[1]=_mm512_sub_epi16(mask, ctx[1]);
					ctx[2]=_mm512_sub_epi16(mask, ctx[2]);
					ctx[3]=_mm512_sub_epi16(mask, ctx[3]);
					ctx[4]=_mm512_sub_epi16(mask, ctx[4]);
					ctx[5]=_mm512_sub_epi16(mask, ctx[5]);
				}
				ctx[0]=_mm512_min_epi16(ctx[0], mctxmax);
				ctx[1]=_mm512_min_epi16(ctx[1], mctxmax);
				ctx[2]=_mm512_min_epi16(ctx[2], mctxmax);
				ctx[3]=_mm512_min_epi16(ctx[3], mctxmax);
				ctx[4]=_mm512_min_epi16(ctx[4], mctxmax);
				ctx[5]=_mm512_min_epi16(ctx[5], mctxmax);
			}
			{
				const int borderW=3;
				const int borderN=3;
				const int borderE=3;
				int cond_cg=(uint32_t)(kx-borderW)>=(uint32_t)(ixbytes-(borderW+borderE))
					||(uint32_t)(ky-borderN)>=(uint32_t)(blockh-borderN);
				__m512i mcg[6];
				__m512i xmin[6], xmax[6];
				xmin[0]=_mm512_min_epi16(N[0], W[0]); xmax[0]=_mm512_max_epi16(N[0], W[0]);
				xmin[1]=_mm512_min_epi16(N[1], W[1]); xmax[1]=_mm512_max_epi16(N[1], W[1]);
				xmin[2]=_mm512_min_epi16(N[2], W[2]); xmax[2]=_mm512_max_epi16(N[2], W[2]);
				xmin[3]=_mm512_min_epi16(N[3], W[3]); xmax[3]=_mm512_max_epi16(N[3], W[3]);
				xmin[4]=_mm512_min_epi16(N[4], W[4]); xmax[4]=_mm512_max_epi16(N[4], W[4]);
				xmin[5]=_mm512_min_epi16(N[5], W[5]); xmax[5]=_mm512_max_epi16(N[5], W[5]);
				pred[0]=_mm512_add_epi16(N[0], W[0]);//N+W-NW
				pred[1]=_mm512_add_epi16(N[1], W[1]);
				pred[2]=_mm512_add_epi16(N[2], W[2]);
				pred[3]=_mm512_add_epi16(N[3], W[3]);
				pred[4]=_mm512_add_epi16(N[4], W[4]);
				pred[5]=_mm512_add_epi16(N[5], W[5]);
				pred[0]=_mm512_sub_epi16(pred[0], NW[0]);
				pred[1]=_mm512_sub_epi16(pred[1], NW[1]);
				pred[2]=_mm512_sub_epi16(pred[2], NW[2]);
				pred[3]=_mm512_sub_epi16(pred[3], NW[3]);
				pred[4]=_mm512_sub_epi16(pred[4], NW[4]);
				pred[5]=_mm512_sub_epi16(pred[5], NW[5]);
				mcg[0]=pred[0];
				mcg[1]=pred[1];
				mcg[2]=pred[2];
				mcg[3]=pred[3];
				mcg[4]=pred[4];
				mcg[5]=pred[5];
				if(effort==1)//predict
				{
					/*
					effort 1
					0	N+W-NW
					1	2*N-NN
					2	W
					3	NE
					*/

					//N+W-NW
					L1preds[0*6+0]=pred[0];
					L1preds[0*6+1]=pred[1];
					L1preds[0*6+2]=pred[2];
					L1preds[0*6+3]=pred[3];
					L1preds[0*6+4]=pred[4];
					L1preds[0*6+5]=pred[5];

					//2*N-NN
					L1preds[1*6+0]=_mm512_sub_epi16(_mm512_add_epi16(N[0], N[0]), _mm512_load_si512((__m512i*)rows[2]+0+(0+0*NCH)*NROWS*NVAL));
					L1preds[1*6+1]=_mm512_sub_epi16(_mm512_add_epi16(N[1], N[1]), _mm512_load_si512((__m512i*)rows[2]+1+(0+0*NCH)*NROWS*NVAL));
					L1preds[1*6+2]=_mm512_sub_epi16(_mm512_add_epi16(N[2], N[2]), _mm512_load_si512((__m512i*)rows[2]+2+(0+0*NCH)*NROWS*NVAL));
					L1preds[1*6+3]=_mm512_sub_epi16(_mm512_add_epi16(N[3], N[3]), _mm512_load_si512((__m512i*)rows[2]+3+(0+0*NCH)*NROWS*NVAL));
					L1preds[1*6+4]=_mm512_sub_epi16(_mm512_add_epi16(N[4], N[4]), _mm512_load_si512((__m512i*)rows[2]+4+(0+0*NCH)*NROWS*NVAL));
					L1preds[1*6+5]=_mm512_sub_epi16(_mm512_add_epi16(N[5], N[5]), _mm512_load_si512((__m512i*)rows[2]+5+(0+0*NCH)*NROWS*NVAL));

					//W
					L1preds[2*6+0]=W[0];
					L1preds[2*6+1]=W[1];
					L1preds[2*6+2]=W[2];
					L1preds[2*6+3]=W[3];
					L1preds[2*6+4]=W[4];
					L1preds[2*6+5]=W[5];

					//NE
					L1preds[3*6+0]=_mm512_load_si512((__m512i*)rows[1]+0+(0+1*NROWS*NVAL)*NREG*NCH);
					L1preds[3*6+1]=_mm512_load_si512((__m512i*)rows[1]+1+(0+1*NROWS*NVAL)*NREG*NCH);
					L1preds[3*6+2]=_mm512_load_si512((__m512i*)rows[1]+2+(0+1*NROWS*NVAL)*NREG*NCH);
					L1preds[3*6+3]=_mm512_load_si512((__m512i*)rows[1]+3+(0+1*NROWS*NVAL)*NREG*NCH);
					L1preds[3*6+4]=_mm512_load_si512((__m512i*)rows[1]+4+(0+1*NROWS*NVAL)*NREG*NCH);
					L1preds[3*6+5]=_mm512_load_si512((__m512i*)rows[1]+5+(0+1*NROWS*NVAL)*NREG*NCH);
					

					//mix
					__m512i mp[12], t[12];
					mp[0x0]=_mm512_setzero_si512();
					mp[0x1]=_mm512_setzero_si512();
					mp[0x2]=_mm512_setzero_si512();
					mp[0x3]=_mm512_setzero_si512();
					mp[0x4]=_mm512_setzero_si512();
					mp[0x5]=_mm512_setzero_si512();
					mp[0x6]=_mm512_setzero_si512();
					mp[0x7]=_mm512_setzero_si512();
					mp[0x8]=_mm512_setzero_si512();
					mp[0x9]=_mm512_setzero_si512();
					mp[0xA]=_mm512_setzero_si512();
					mp[0xB]=_mm512_setzero_si512();
					for(int k=0;k<L1_NPREDS1;++k)
					{
						//16 -> 32		3 lo 3 hi registers
						t[0x0]=_mm512_slli_epi32(L1preds[k*6+0], 16);
						t[0x1]=_mm512_slli_epi32(L1preds[k*6+1], 16);
						t[0x2]=_mm512_slli_epi32(L1preds[k*6+2], 16);
						t[0x3]=_mm512_slli_epi32(L1preds[k*6+3], 16);
						t[0x4]=_mm512_slli_epi32(L1preds[k*6+4], 16);
						t[0x5]=_mm512_slli_epi32(L1preds[k*6+5], 16);
						t[0x6]=_mm512_srai_epi32(L1preds[k*6+0], 16);
						t[0x7]=_mm512_srai_epi32(L1preds[k*6+1], 16);
						t[0x8]=_mm512_srai_epi32(L1preds[k*6+2], 16);
						t[0x9]=_mm512_srai_epi32(L1preds[k*6+3], 16);
						t[0xA]=_mm512_srai_epi32(L1preds[k*6+4], 16);
						t[0xB]=_mm512_srai_epi32(L1preds[k*6+5], 16);
						t[0x0]=_mm512_srai_epi32(t[0x0], 16);
						t[0x1]=_mm512_srai_epi32(t[0x1], 16);
						t[0x2]=_mm512_srai_epi32(t[0x2], 16);
						t[0x3]=_mm512_srai_epi32(t[0x3], 16);
						t[0x4]=_mm512_srai_epi32(t[0x4], 16);
						t[0x5]=_mm512_srai_epi32(t[0x5], 16);
						t[0x0]=_mm512_mullo_epi32(t[0x0], _mm512_load_si512((__m512i*)L1weights+k*12+0x0));
						t[0x1]=_mm512_mullo_epi32(t[0x1], _mm512_load_si512((__m512i*)L1weights+k*12+0x1));
						t[0x2]=_mm512_mullo_epi32(t[0x2], _mm512_load_si512((__m512i*)L1weights+k*12+0x2));
						t[0x3]=_mm512_mullo_epi32(t[0x3], _mm512_load_si512((__m512i*)L1weights+k*12+0x3));
						t[0x4]=_mm512_mullo_epi32(t[0x4], _mm512_load_si512((__m512i*)L1weights+k*12+0x4));
						t[0x5]=_mm512_mullo_epi32(t[0x5], _mm512_load_si512((__m512i*)L1weights+k*12+0x5));
						t[0x6]=_mm512_mullo_epi32(t[0x6], _mm512_load_si512((__m512i*)L1weights+k*12+0x6));
						t[0x7]=_mm512_mullo_epi32(t[0x7], _mm512_load_si512((__m512i*)L1weights+k*12+0x7));
						t[0x8]=_mm512_mullo_epi32(t[0x8], _mm512_load_si512((__m512i*)L1weights+k*12+0x8));
						t[0x9]=_mm512_mullo_epi32(t[0x9], _mm512_load_si512((__m512i*)L1weights+k*12+0x9));
						t[0xA]=_mm512_mullo_epi32(t[0xA], _mm512_load_si512((__m512i*)L1weights+k*12+0xA));
						t[0xB]=_mm512_mullo_epi32(t[0xB], _mm512_load_si512((__m512i*)L1weights+k*12+0xB));
						mp[0x0]=_mm512_add_epi32(mp[0x0], t[0x0]);
						mp[0x1]=_mm512_add_epi32(mp[0x1], t[0x1]);
						mp[0x2]=_mm512_add_epi32(mp[0x2], t[0x2]);
						mp[0x3]=_mm512_add_epi32(mp[0x3], t[0x3]);
						mp[0x4]=_mm512_add_epi32(mp[0x4], t[0x4]);
						mp[0x5]=_mm512_add_epi32(mp[0x5], t[0x5]);
						mp[0x6]=_mm512_add_epi32(mp[0x6], t[0x6]);
						mp[0x7]=_mm512_add_epi32(mp[0x7], t[0x7]);
						mp[0x8]=_mm512_add_epi32(mp[0x8], t[0x8]);
						mp[0x9]=_mm512_add_epi32(mp[0x9], t[0x9]);
						mp[0xA]=_mm512_add_epi32(mp[0xA], t[0xA]);
						mp[0xB]=_mm512_add_epi32(mp[0xB], t[0xB]);
					}

					__m512i rcon=_mm512_set1_epi32(1<<L1_SH1>>1);
					mp[0x0]=_mm512_add_epi32(mp[0x0], rcon);//rounding to nearest
					mp[0x1]=_mm512_add_epi32(mp[0x1], rcon);
					mp[0x2]=_mm512_add_epi32(mp[0x2], rcon);
					mp[0x3]=_mm512_add_epi32(mp[0x3], rcon);
					mp[0x4]=_mm512_add_epi32(mp[0x4], rcon);
					mp[0x5]=_mm512_add_epi32(mp[0x5], rcon);
					mp[0x6]=_mm512_add_epi32(mp[0x6], rcon);
					mp[0x7]=_mm512_add_epi32(mp[0x7], rcon);
					mp[0x8]=_mm512_add_epi32(mp[0x8], rcon);
					mp[0x9]=_mm512_add_epi32(mp[0x9], rcon);
					mp[0xA]=_mm512_add_epi32(mp[0xA], rcon);
					mp[0xB]=_mm512_add_epi32(mp[0xB], rcon);

					mp[0x0]=_mm512_srai_epi32(mp[0x0], L1_SH1);
					mp[0x1]=_mm512_srai_epi32(mp[0x1], L1_SH1);
					mp[0x2]=_mm512_srai_epi32(mp[0x2], L1_SH1);
					mp[0x3]=_mm512_srai_epi32(mp[0x3], L1_SH1);
					mp[0x4]=_mm512_srai_epi32(mp[0x4], L1_SH1);
					mp[0x5]=_mm512_srai_epi32(mp[0x5], L1_SH1);
					mp[0x6]=_mm512_slli_epi32(mp[0x6], 16-L1_SH1);
					mp[0x7]=_mm512_slli_epi32(mp[0x7], 16-L1_SH1);
					mp[0x8]=_mm512_slli_epi32(mp[0x8], 16-L1_SH1);
					mp[0x9]=_mm512_slli_epi32(mp[0x9], 16-L1_SH1);
					mp[0xA]=_mm512_slli_epi32(mp[0xA], 16-L1_SH1);
					mp[0xB]=_mm512_slli_epi32(mp[0xB], 16-L1_SH1);
					//32 -> 16
					pred[0]=_mm512_mask_blend_epi16(0xAAAAAAAA, mp[0x0], mp[0x6]);
					pred[1]=_mm512_mask_blend_epi16(0xAAAAAAAA, mp[0x1], mp[0x7]);
					pred[2]=_mm512_mask_blend_epi16(0xAAAAAAAA, mp[0x2], mp[0x8]);
					pred[3]=_mm512_mask_blend_epi16(0xAAAAAAAA, mp[0x3], mp[0x9]);
					pred[4]=_mm512_mask_blend_epi16(0xAAAAAAAA, mp[0x4], mp[0xA]);
					pred[5]=_mm512_mask_blend_epi16(0xAAAAAAAA, mp[0x5], mp[0xB]);


					//loosen pred range
					if(!cond_cg)
					{
						t[0]=_mm512_load_si512((__m512i*)rows[1]+0+(0+1*NROWS*NVAL)*NREG*NCH);//NE
						t[1]=_mm512_load_si512((__m512i*)rows[1]+1+(0+1*NROWS*NVAL)*NREG*NCH);
						t[2]=_mm512_load_si512((__m512i*)rows[1]+2+(0+1*NROWS*NVAL)*NREG*NCH);
						t[3]=_mm512_load_si512((__m512i*)rows[1]+3+(0+1*NROWS*NVAL)*NREG*NCH);
						t[4]=_mm512_load_si512((__m512i*)rows[1]+4+(0+1*NROWS*NVAL)*NREG*NCH);
						t[5]=_mm512_load_si512((__m512i*)rows[1]+5+(0+1*NROWS*NVAL)*NREG*NCH);
						xmin[0]=_mm512_min_epi16(xmin[0], t[0]); xmax[0]=_mm512_max_epi16(xmax[0], t[0]);
						xmin[1]=_mm512_min_epi16(xmin[1], t[1]); xmax[1]=_mm512_max_epi16(xmax[1], t[1]);
						xmin[2]=_mm512_min_epi16(xmin[2], t[2]); xmax[2]=_mm512_max_epi16(xmax[2], t[2]);
						xmin[3]=_mm512_min_epi16(xmin[3], t[3]); xmax[3]=_mm512_max_epi16(xmax[3], t[3]);
						xmin[4]=_mm512_min_epi16(xmin[4], t[4]); xmax[4]=_mm512_max_epi16(xmax[4], t[4]);
						xmin[5]=_mm512_min_epi16(xmin[5], t[5]); xmax[5]=_mm512_max_epi16(xmax[5], t[5]);
						t[0]=_mm512_load_si512((__m512i*)rows[1]+0+(0+3*NROWS*NVAL)*NREG*NCH);//NEEE
						t[1]=_mm512_load_si512((__m512i*)rows[1]+1+(0+3*NROWS*NVAL)*NREG*NCH);
						t[2]=_mm512_load_si512((__m512i*)rows[1]+2+(0+3*NROWS*NVAL)*NREG*NCH);
						t[3]=_mm512_load_si512((__m512i*)rows[1]+3+(0+3*NROWS*NVAL)*NREG*NCH);
						t[4]=_mm512_load_si512((__m512i*)rows[1]+4+(0+3*NROWS*NVAL)*NREG*NCH);
						t[5]=_mm512_load_si512((__m512i*)rows[1]+5+(0+3*NROWS*NVAL)*NREG*NCH);
						xmin[0]=_mm512_min_epi16(xmin[0], t[0]); xmax[0]=_mm512_max_epi16(xmax[0], t[0]);
						xmin[1]=_mm512_min_epi16(xmin[1], t[1]); xmax[1]=_mm512_max_epi16(xmax[1], t[1]);
						xmin[2]=_mm512_min_epi16(xmin[2], t[2]); xmax[2]=_mm512_max_epi16(xmax[2], t[2]);
						xmin[3]=_mm512_min_epi16(xmin[3], t[3]); xmax[3]=_mm512_max_epi16(xmax[3], t[3]);
						xmin[4]=_mm512_min_epi16(xmin[4], t[4]); xmax[4]=_mm512_max_epi16(xmax[4], t[4]);
						xmin[5]=_mm512_min_epi16(xmin[5], t[5]); xmax[5]=_mm512_max_epi16(xmax[5], t[5]);
					}
				}
				else if(effort==2)
				{
					__m512i cache[6];
					/*
					effort 2
					0	N
					1	W
					2	3*(N-NN)+NNN
					3	3*(W-WW)+WWW
					4	W+NE-N
					5	(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4
					6	N+W-NW
					7	N+NE-NNE
					*/

					//N
					L1preds[0*6+0]=N[0];
					L1preds[0*6+1]=N[1];
					L1preds[0*6+2]=N[2];
					L1preds[0*6+3]=N[3];
					L1preds[0*6+4]=N[4];
					L1preds[0*6+5]=N[5];

					//W
					L1preds[1*6+0]=W[0];
					L1preds[1*6+1]=W[1];
					L1preds[1*6+2]=W[2];
					L1preds[1*6+3]=W[3];
					L1preds[1*6+4]=W[4];
					L1preds[1*6+5]=W[5];

					//3*(N-NN)+NNN
					cache[0]=_mm512_sub_epi16(N[0], _mm512_load_si512((__m512i*)rows[2]+0+(0+0*NROWS*NVAL)*NREG*NCH));//N-NN
					cache[1]=_mm512_sub_epi16(N[1], _mm512_load_si512((__m512i*)rows[2]+1+(0+0*NROWS*NVAL)*NREG*NCH));
					cache[2]=_mm512_sub_epi16(N[2], _mm512_load_si512((__m512i*)rows[2]+2+(0+0*NROWS*NVAL)*NREG*NCH));
					cache[3]=_mm512_sub_epi16(N[3], _mm512_load_si512((__m512i*)rows[2]+3+(0+0*NROWS*NVAL)*NREG*NCH));
					cache[4]=_mm512_sub_epi16(N[4], _mm512_load_si512((__m512i*)rows[2]+4+(0+0*NROWS*NVAL)*NREG*NCH));
					cache[5]=_mm512_sub_epi16(N[5], _mm512_load_si512((__m512i*)rows[2]+5+(0+0*NROWS*NVAL)*NREG*NCH));
					cache[0]=_mm512_add_epi16(cache[0], _mm512_slli_epi16(cache[0], 1));//*3
					cache[1]=_mm512_add_epi16(cache[1], _mm512_slli_epi16(cache[1], 1));
					cache[2]=_mm512_add_epi16(cache[2], _mm512_slli_epi16(cache[2], 1));
					cache[3]=_mm512_add_epi16(cache[3], _mm512_slli_epi16(cache[3], 1));
					cache[4]=_mm512_add_epi16(cache[4], _mm512_slli_epi16(cache[4], 1));
					cache[5]=_mm512_add_epi16(cache[5], _mm512_slli_epi16(cache[5], 1));
					L1preds[2*6+0]=_mm512_add_epi16(cache[0], _mm512_load_si512((__m512i*)rows[3]+0+(0+0*NROWS*NVAL)*NREG*NCH));//+NNN
					L1preds[2*6+1]=_mm512_add_epi16(cache[1], _mm512_load_si512((__m512i*)rows[3]+1+(0+0*NROWS*NVAL)*NREG*NCH));
					L1preds[2*6+2]=_mm512_add_epi16(cache[2], _mm512_load_si512((__m512i*)rows[3]+2+(0+0*NROWS*NVAL)*NREG*NCH));
					L1preds[2*6+3]=_mm512_add_epi16(cache[3], _mm512_load_si512((__m512i*)rows[3]+3+(0+0*NROWS*NVAL)*NREG*NCH));
					L1preds[2*6+4]=_mm512_add_epi16(cache[4], _mm512_load_si512((__m512i*)rows[3]+4+(0+0*NROWS*NVAL)*NREG*NCH));
					L1preds[2*6+5]=_mm512_add_epi16(cache[5], _mm512_load_si512((__m512i*)rows[3]+5+(0+0*NROWS*NVAL)*NREG*NCH));

					//3*(W-WW)+WWW
					cache[0]=_mm512_sub_epi16(W[0], _mm512_load_si512((__m512i*)rows[0]+0+(0-2*NROWS*NVAL)*NREG*NCH));//W-WW
					cache[1]=_mm512_sub_epi16(W[1], _mm512_load_si512((__m512i*)rows[0]+1+(0-2*NROWS*NVAL)*NREG*NCH));
					cache[2]=_mm512_sub_epi16(W[2], _mm512_load_si512((__m512i*)rows[0]+2+(0-2*NROWS*NVAL)*NREG*NCH));
					cache[3]=_mm512_sub_epi16(W[3], _mm512_load_si512((__m512i*)rows[0]+3+(0-2*NROWS*NVAL)*NREG*NCH));
					cache[4]=_mm512_sub_epi16(W[4], _mm512_load_si512((__m512i*)rows[0]+4+(0-2*NROWS*NVAL)*NREG*NCH));
					cache[5]=_mm512_sub_epi16(W[5], _mm512_load_si512((__m512i*)rows[0]+5+(0-2*NROWS*NVAL)*NREG*NCH));
					cache[0]=_mm512_add_epi16(cache[0], _mm512_slli_epi16(cache[0], 1));//*3
					cache[1]=_mm512_add_epi16(cache[1], _mm512_slli_epi16(cache[1], 1));
					cache[2]=_mm512_add_epi16(cache[2], _mm512_slli_epi16(cache[2], 1));
					cache[3]=_mm512_add_epi16(cache[3], _mm512_slli_epi16(cache[3], 1));
					cache[4]=_mm512_add_epi16(cache[4], _mm512_slli_epi16(cache[4], 1));
					cache[5]=_mm512_add_epi16(cache[5], _mm512_slli_epi16(cache[5], 1));
					L1preds[3*6+0]=_mm512_add_epi16(cache[0], _mm512_load_si512((__m512i*)rows[0]+0+(0-3*NROWS*NVAL)*NREG*NCH));//+WWW
					L1preds[3*6+1]=_mm512_add_epi16(cache[1], _mm512_load_si512((__m512i*)rows[0]+1+(0-3*NROWS*NVAL)*NREG*NCH));
					L1preds[3*6+2]=_mm512_add_epi16(cache[2], _mm512_load_si512((__m512i*)rows[0]+2+(0-3*NROWS*NVAL)*NREG*NCH));
					L1preds[3*6+3]=_mm512_add_epi16(cache[3], _mm512_load_si512((__m512i*)rows[0]+3+(0-3*NROWS*NVAL)*NREG*NCH));
					L1preds[3*6+4]=_mm512_add_epi16(cache[4], _mm512_load_si512((__m512i*)rows[0]+4+(0-3*NROWS*NVAL)*NREG*NCH));
					L1preds[3*6+5]=_mm512_add_epi16(cache[5], _mm512_load_si512((__m512i*)rows[0]+5+(0-3*NROWS*NVAL)*NREG*NCH));

					//W+NE-N
					cache[0]=_mm512_sub_epi16(W[0], N[0]);
					cache[1]=_mm512_sub_epi16(W[1], N[1]);
					cache[2]=_mm512_sub_epi16(W[2], N[2]);
					cache[3]=_mm512_sub_epi16(W[3], N[3]);
					cache[4]=_mm512_sub_epi16(W[4], N[4]);
					cache[5]=_mm512_sub_epi16(W[5], N[5]);
					L1preds[4*3+0]=_mm512_add_epi16(cache[0], _mm512_load_si512((__m512i*)rows[1]+0+(0+1*NROWS*NVAL)*NREG*NCH));//+NE
					L1preds[4*3+1]=_mm512_add_epi16(cache[1], _mm512_load_si512((__m512i*)rows[1]+1+(0+1*NROWS*NVAL)*NREG*NCH));
					L1preds[4*3+2]=_mm512_add_epi16(cache[2], _mm512_load_si512((__m512i*)rows[1]+2+(0+1*NROWS*NVAL)*NREG*NCH));
					L1preds[4*3+3]=_mm512_add_epi16(cache[3], _mm512_load_si512((__m512i*)rows[1]+3+(0+1*NROWS*NVAL)*NREG*NCH));
					L1preds[4*3+4]=_mm512_add_epi16(cache[4], _mm512_load_si512((__m512i*)rows[1]+4+(0+1*NROWS*NVAL)*NREG*NCH));
					L1preds[4*3+5]=_mm512_add_epi16(cache[5], _mm512_load_si512((__m512i*)rows[1]+5+(0+1*NROWS*NVAL)*NREG*NCH));

					//N+W-NW
					L1preds[5*6+0]=pred[0];
					L1preds[5*6+1]=pred[1];
					L1preds[5*6+2]=pred[2];
					L1preds[5*6+3]=pred[3];
					L1preds[5*6+4]=pred[4];
					L1preds[5*6+5]=pred[5];

					//N+NE-NNE
					cache[0]=_mm512_add_epi16(N[0], _mm512_load_si512((__m512i*)rows[1]+0+(0+1*NROWS*NVAL)*NREG*NCH));//N+NE
					cache[1]=_mm512_add_epi16(N[1], _mm512_load_si512((__m512i*)rows[1]+1+(0+1*NROWS*NVAL)*NREG*NCH));
					cache[2]=_mm512_add_epi16(N[2], _mm512_load_si512((__m512i*)rows[1]+2+(0+1*NROWS*NVAL)*NREG*NCH));
					cache[3]=_mm512_add_epi16(N[3], _mm512_load_si512((__m512i*)rows[1]+3+(0+1*NROWS*NVAL)*NREG*NCH));
					cache[4]=_mm512_add_epi16(N[4], _mm512_load_si512((__m512i*)rows[1]+4+(0+1*NROWS*NVAL)*NREG*NCH));
					cache[5]=_mm512_add_epi16(N[5], _mm512_load_si512((__m512i*)rows[1]+5+(0+1*NROWS*NVAL)*NREG*NCH));
					L1preds[6*6+0]=_mm512_sub_epi16(cache[0], _mm512_load_si512((__m512i*)rows[2]+0+(0+1*NROWS*NVAL)*NREG*NCH));//NNE
					L1preds[6*6+1]=_mm512_sub_epi16(cache[1], _mm512_load_si512((__m512i*)rows[2]+1+(0+1*NROWS*NVAL)*NREG*NCH));
					L1preds[6*6+2]=_mm512_sub_epi16(cache[2], _mm512_load_si512((__m512i*)rows[2]+2+(0+1*NROWS*NVAL)*NREG*NCH));
					L1preds[6*6+3]=_mm512_sub_epi16(cache[3], _mm512_load_si512((__m512i*)rows[2]+3+(0+1*NROWS*NVAL)*NREG*NCH));
					L1preds[6*6+4]=_mm512_sub_epi16(cache[4], _mm512_load_si512((__m512i*)rows[2]+4+(0+1*NROWS*NVAL)*NREG*NCH));
					L1preds[6*6+5]=_mm512_sub_epi16(cache[5], _mm512_load_si512((__m512i*)rows[2]+5+(0+1*NROWS*NVAL)*NREG*NCH));

					//(WWWW+WWW+NNN+NNEE+NEEE+NEEEE-(N+W))>>2
					cache[0]=_mm512_load_si512((__m512i*)rows[0]+0+(0-4*NROWS*NVAL)*NREG*NCH);//WWWW
					cache[1]=_mm512_load_si512((__m512i*)rows[0]+1+(0-4*NROWS*NVAL)*NREG*NCH);
					cache[2]=_mm512_load_si512((__m512i*)rows[0]+2+(0-4*NROWS*NVAL)*NREG*NCH);
					cache[3]=_mm512_load_si512((__m512i*)rows[0]+3+(0-4*NROWS*NVAL)*NREG*NCH);
					cache[4]=_mm512_load_si512((__m512i*)rows[0]+4+(0-4*NROWS*NVAL)*NREG*NCH);
					cache[5]=_mm512_load_si512((__m512i*)rows[0]+5+(0-4*NROWS*NVAL)*NREG*NCH);
					cache[0]=_mm512_add_epi16(cache[0], _mm512_load_si512((__m512i*)rows[0]+0+(0-3*NROWS*NVAL)*NREG*NCH));//+WWW
					cache[1]=_mm512_add_epi16(cache[1], _mm512_load_si512((__m512i*)rows[0]+1+(0-3*NROWS*NVAL)*NREG*NCH));
					cache[2]=_mm512_add_epi16(cache[2], _mm512_load_si512((__m512i*)rows[0]+2+(0-3*NROWS*NVAL)*NREG*NCH));
					cache[3]=_mm512_add_epi16(cache[3], _mm512_load_si512((__m512i*)rows[0]+3+(0-3*NROWS*NVAL)*NREG*NCH));
					cache[4]=_mm512_add_epi16(cache[4], _mm512_load_si512((__m512i*)rows[0]+4+(0-3*NROWS*NVAL)*NREG*NCH));
					cache[5]=_mm512_add_epi16(cache[5], _mm512_load_si512((__m512i*)rows[0]+5+(0-3*NROWS*NVAL)*NREG*NCH));
					cache[0]=_mm512_add_epi16(cache[0], _mm512_load_si512((__m512i*)rows[3]+0+(0+0*NROWS*NVAL)*NREG*NCH));//+NNN
					cache[1]=_mm512_add_epi16(cache[1], _mm512_load_si512((__m512i*)rows[3]+1+(0+0*NROWS*NVAL)*NREG*NCH));
					cache[2]=_mm512_add_epi16(cache[2], _mm512_load_si512((__m512i*)rows[3]+2+(0+0*NROWS*NVAL)*NREG*NCH));
					cache[3]=_mm512_add_epi16(cache[3], _mm512_load_si512((__m512i*)rows[3]+3+(0+0*NROWS*NVAL)*NREG*NCH));
					cache[4]=_mm512_add_epi16(cache[4], _mm512_load_si512((__m512i*)rows[3]+4+(0+0*NROWS*NVAL)*NREG*NCH));
					cache[5]=_mm512_add_epi16(cache[5], _mm512_load_si512((__m512i*)rows[3]+5+(0+0*NROWS*NVAL)*NREG*NCH));
					cache[0]=_mm512_add_epi16(cache[0], _mm512_load_si512((__m512i*)rows[2]+0+(0+2*NROWS*NVAL)*NREG*NCH));//+NNEE
					cache[1]=_mm512_add_epi16(cache[1], _mm512_load_si512((__m512i*)rows[2]+1+(0+2*NROWS*NVAL)*NREG*NCH));
					cache[2]=_mm512_add_epi16(cache[2], _mm512_load_si512((__m512i*)rows[2]+2+(0+2*NROWS*NVAL)*NREG*NCH));
					cache[3]=_mm512_add_epi16(cache[3], _mm512_load_si512((__m512i*)rows[2]+3+(0+2*NROWS*NVAL)*NREG*NCH));
					cache[4]=_mm512_add_epi16(cache[4], _mm512_load_si512((__m512i*)rows[2]+4+(0+2*NROWS*NVAL)*NREG*NCH));
					cache[5]=_mm512_add_epi16(cache[5], _mm512_load_si512((__m512i*)rows[2]+5+(0+2*NROWS*NVAL)*NREG*NCH));
					cache[0]=_mm512_add_epi16(cache[0], _mm512_load_si512((__m512i*)rows[1]+0+(0+3*NROWS*NVAL)*NREG*NCH));//+NEEE
					cache[1]=_mm512_add_epi16(cache[1], _mm512_load_si512((__m512i*)rows[1]+1+(0+3*NROWS*NVAL)*NREG*NCH));
					cache[2]=_mm512_add_epi16(cache[2], _mm512_load_si512((__m512i*)rows[1]+2+(0+3*NROWS*NVAL)*NREG*NCH));
					cache[3]=_mm512_add_epi16(cache[3], _mm512_load_si512((__m512i*)rows[1]+3+(0+3*NROWS*NVAL)*NREG*NCH));
					cache[4]=_mm512_add_epi16(cache[4], _mm512_load_si512((__m512i*)rows[1]+4+(0+3*NROWS*NVAL)*NREG*NCH));
					cache[5]=_mm512_add_epi16(cache[5], _mm512_load_si512((__m512i*)rows[1]+5+(0+3*NROWS*NVAL)*NREG*NCH));
					cache[0]=_mm512_add_epi16(cache[0], _mm512_load_si512((__m512i*)rows[1]+0+(0+4*NROWS*NVAL)*NREG*NCH));//+NEEEE
					cache[1]=_mm512_add_epi16(cache[1], _mm512_load_si512((__m512i*)rows[1]+1+(0+4*NROWS*NVAL)*NREG*NCH));
					cache[2]=_mm512_add_epi16(cache[2], _mm512_load_si512((__m512i*)rows[1]+2+(0+4*NROWS*NVAL)*NREG*NCH));
					cache[3]=_mm512_add_epi16(cache[3], _mm512_load_si512((__m512i*)rows[1]+3+(0+4*NROWS*NVAL)*NREG*NCH));
					cache[4]=_mm512_add_epi16(cache[4], _mm512_load_si512((__m512i*)rows[1]+4+(0+4*NROWS*NVAL)*NREG*NCH));
					cache[5]=_mm512_add_epi16(cache[5], _mm512_load_si512((__m512i*)rows[1]+5+(0+4*NROWS*NVAL)*NREG*NCH));
					cache[0]=_mm512_sub_epi16(cache[0], _mm512_add_epi16(N[0], W[0]));
					cache[1]=_mm512_sub_epi16(cache[1], _mm512_add_epi16(N[1], W[1]));
					cache[2]=_mm512_sub_epi16(cache[2], _mm512_add_epi16(N[2], W[2]));
					cache[3]=_mm512_sub_epi16(cache[3], _mm512_add_epi16(N[3], W[3]));
					cache[4]=_mm512_sub_epi16(cache[4], _mm512_add_epi16(N[4], W[4]));
					cache[5]=_mm512_sub_epi16(cache[5], _mm512_add_epi16(N[5], W[5]));
					L1preds[7*6+0]=_mm512_srai_epi16(cache[0], 2);
					L1preds[7*6+1]=_mm512_srai_epi16(cache[1], 2);
					L1preds[7*6+2]=_mm512_srai_epi16(cache[2], 2);
					L1preds[7*6+3]=_mm512_srai_epi16(cache[3], 2);
					L1preds[7*6+4]=_mm512_srai_epi16(cache[4], 2);
					L1preds[7*6+5]=_mm512_srai_epi16(cache[5], 2);


					//mix
					__m512i mp[12], t[12];
					mp[0x0]=_mm512_setzero_si512();
					mp[0x1]=_mm512_setzero_si512();
					mp[0x2]=_mm512_setzero_si512();
					mp[0x3]=_mm512_setzero_si512();
					mp[0x4]=_mm512_setzero_si512();
					mp[0x5]=_mm512_setzero_si512();
					mp[0x6]=_mm512_setzero_si512();
					mp[0x7]=_mm512_setzero_si512();
					mp[0x8]=_mm512_setzero_si512();
					mp[0x9]=_mm512_setzero_si512();
					mp[0xA]=_mm512_setzero_si512();
					mp[0xB]=_mm512_setzero_si512();
					for(int k=0;k<L1_NPREDS2;++k)
					{
						//16 -> 32		3 lo 3 hi registers
						t[0x0]=_mm512_slli_epi32(L1preds[k*6+0], 16);
						t[0x1]=_mm512_slli_epi32(L1preds[k*6+1], 16);
						t[0x2]=_mm512_slli_epi32(L1preds[k*6+2], 16);
						t[0x3]=_mm512_slli_epi32(L1preds[k*6+3], 16);
						t[0x4]=_mm512_slli_epi32(L1preds[k*6+4], 16);
						t[0x5]=_mm512_slli_epi32(L1preds[k*6+5], 16);
						t[0x6]=_mm512_srai_epi32(L1preds[k*6+0], 16);
						t[0x7]=_mm512_srai_epi32(L1preds[k*6+1], 16);
						t[0x8]=_mm512_srai_epi32(L1preds[k*6+2], 16);
						t[0x9]=_mm512_srai_epi32(L1preds[k*6+3], 16);
						t[0xA]=_mm512_srai_epi32(L1preds[k*6+4], 16);
						t[0xB]=_mm512_srai_epi32(L1preds[k*6+5], 16);
						t[0x0]=_mm512_srai_epi32(t[0x0], 16);
						t[0x1]=_mm512_srai_epi32(t[0x1], 16);
						t[0x2]=_mm512_srai_epi32(t[0x2], 16);
						t[0x3]=_mm512_srai_epi32(t[0x3], 16);
						t[0x4]=_mm512_srai_epi32(t[0x4], 16);
						t[0x5]=_mm512_srai_epi32(t[0x5], 16);
						t[0x0]=_mm512_mullo_epi32(t[0x0], _mm512_load_si512((__m512i*)L1weights+k*12+0x0));
						t[0x1]=_mm512_mullo_epi32(t[0x1], _mm512_load_si512((__m512i*)L1weights+k*12+0x1));
						t[0x2]=_mm512_mullo_epi32(t[0x2], _mm512_load_si512((__m512i*)L1weights+k*12+0x2));
						t[0x3]=_mm512_mullo_epi32(t[0x3], _mm512_load_si512((__m512i*)L1weights+k*12+0x3));
						t[0x4]=_mm512_mullo_epi32(t[0x4], _mm512_load_si512((__m512i*)L1weights+k*12+0x4));
						t[0x5]=_mm512_mullo_epi32(t[0x5], _mm512_load_si512((__m512i*)L1weights+k*12+0x5));
						t[0x6]=_mm512_mullo_epi32(t[0x6], _mm512_load_si512((__m512i*)L1weights+k*12+0x6));
						t[0x7]=_mm512_mullo_epi32(t[0x7], _mm512_load_si512((__m512i*)L1weights+k*12+0x7));
						t[0x8]=_mm512_mullo_epi32(t[0x8], _mm512_load_si512((__m512i*)L1weights+k*12+0x8));
						t[0x9]=_mm512_mullo_epi32(t[0x9], _mm512_load_si512((__m512i*)L1weights+k*12+0x9));
						t[0xA]=_mm512_mullo_epi32(t[0xA], _mm512_load_si512((__m512i*)L1weights+k*12+0xA));
						t[0xB]=_mm512_mullo_epi32(t[0xB], _mm512_load_si512((__m512i*)L1weights+k*12+0xB));
						mp[0x0]=_mm512_add_epi32(mp[0x0], t[0x0]);
						mp[0x1]=_mm512_add_epi32(mp[0x1], t[0x1]);
						mp[0x2]=_mm512_add_epi32(mp[0x2], t[0x2]);
						mp[0x3]=_mm512_add_epi32(mp[0x3], t[0x3]);
						mp[0x4]=_mm512_add_epi32(mp[0x4], t[0x4]);
						mp[0x5]=_mm512_add_epi32(mp[0x5], t[0x5]);
						mp[0x6]=_mm512_add_epi32(mp[0x6], t[0x6]);
						mp[0x7]=_mm512_add_epi32(mp[0x7], t[0x7]);
						mp[0x8]=_mm512_add_epi32(mp[0x8], t[0x8]);
						mp[0x9]=_mm512_add_epi32(mp[0x9], t[0x9]);
						mp[0xA]=_mm512_add_epi32(mp[0xA], t[0xA]);
						mp[0xB]=_mm512_add_epi32(mp[0xB], t[0xB]);
					}
					__m512i rcon=_mm512_set1_epi32(1<<L1_SH2>>1);
					mp[0x0]=_mm512_add_epi32(mp[0x0], rcon);//rounding to nearest
					mp[0x1]=_mm512_add_epi32(mp[0x1], rcon);
					mp[0x2]=_mm512_add_epi32(mp[0x2], rcon);
					mp[0x3]=_mm512_add_epi32(mp[0x3], rcon);
					mp[0x4]=_mm512_add_epi32(mp[0x4], rcon);
					mp[0x5]=_mm512_add_epi32(mp[0x5], rcon);
					mp[0x6]=_mm512_add_epi32(mp[0x6], rcon);
					mp[0x7]=_mm512_add_epi32(mp[0x7], rcon);
					mp[0x8]=_mm512_add_epi32(mp[0x8], rcon);
					mp[0x9]=_mm512_add_epi32(mp[0x9], rcon);
					mp[0xA]=_mm512_add_epi32(mp[0xA], rcon);
					mp[0xB]=_mm512_add_epi32(mp[0xB], rcon);

					mp[0x0]=_mm512_srai_epi32(mp[0x0], L1_SH2);
					mp[0x1]=_mm512_srai_epi32(mp[0x1], L1_SH2);
					mp[0x2]=_mm512_srai_epi32(mp[0x2], L1_SH2);
					mp[0x3]=_mm512_srai_epi32(mp[0x3], L1_SH2);
					mp[0x4]=_mm512_srai_epi32(mp[0x4], L1_SH2);
					mp[0x5]=_mm512_srai_epi32(mp[0x5], L1_SH2);
					mp[0x6]=_mm512_srai_epi32(mp[0x6], L1_SH2-16);
					mp[0x7]=_mm512_srai_epi32(mp[0x7], L1_SH2-16);
					mp[0x8]=_mm512_srai_epi32(mp[0x8], L1_SH2-16);
					mp[0x9]=_mm512_srai_epi32(mp[0x9], L1_SH2-16);
					mp[0xA]=_mm512_srai_epi32(mp[0xA], L1_SH2-16);
					mp[0xB]=_mm512_srai_epi32(mp[0xB], L1_SH2-16);
					//32 -> 16
					pred[0]=_mm512_mask_blend_epi16(0xAAAAAAAA, mp[0x0], mp[0x6]);
					pred[1]=_mm512_mask_blend_epi16(0xAAAAAAAA, mp[0x1], mp[0x7]);
					pred[2]=_mm512_mask_blend_epi16(0xAAAAAAAA, mp[0x2], mp[0x8]);
					pred[3]=_mm512_mask_blend_epi16(0xAAAAAAAA, mp[0x3], mp[0x9]);
					pred[4]=_mm512_mask_blend_epi16(0xAAAAAAAA, mp[0x4], mp[0xA]);
					pred[5]=_mm512_mask_blend_epi16(0xAAAAAAAA, mp[0x5], mp[0xB]);


					//loosen pred range
					if(!cond_cg)
					{
						cache[0]=_mm512_load_si512((__m512i*)rows[1]+0+(0+1*NROWS*NVAL)*NREG*NCH);//NE
						cache[1]=_mm512_load_si512((__m512i*)rows[1]+1+(0+1*NROWS*NVAL)*NREG*NCH);
						cache[2]=_mm512_load_si512((__m512i*)rows[1]+2+(0+1*NROWS*NVAL)*NREG*NCH);
						cache[3]=_mm512_load_si512((__m512i*)rows[1]+3+(0+1*NROWS*NVAL)*NREG*NCH);
						cache[4]=_mm512_load_si512((__m512i*)rows[1]+4+(0+1*NROWS*NVAL)*NREG*NCH);
						cache[5]=_mm512_load_si512((__m512i*)rows[1]+5+(0+1*NROWS*NVAL)*NREG*NCH);
						xmin[0]=_mm512_min_epi16(xmin[0], cache[0]); xmax[0]=_mm512_max_epi16(xmax[0], cache[0]);
						xmin[1]=_mm512_min_epi16(xmin[1], cache[1]); xmax[1]=_mm512_max_epi16(xmax[1], cache[1]);
						xmin[2]=_mm512_min_epi16(xmin[2], cache[2]); xmax[2]=_mm512_max_epi16(xmax[2], cache[2]);
						xmin[3]=_mm512_min_epi16(xmin[3], cache[3]); xmax[3]=_mm512_max_epi16(xmax[3], cache[3]);
						xmin[4]=_mm512_min_epi16(xmin[4], cache[4]); xmax[4]=_mm512_max_epi16(xmax[4], cache[4]);
						xmin[5]=_mm512_min_epi16(xmin[5], cache[5]); xmax[5]=_mm512_max_epi16(xmax[5], cache[5]);
						cache[0]=_mm512_load_si512((__m512i*)rows[1]+0+(0+3*NROWS*NVAL)*NREG*NCH);//NEEE
						cache[1]=_mm512_load_si512((__m512i*)rows[1]+1+(0+3*NROWS*NVAL)*NREG*NCH);
						cache[2]=_mm512_load_si512((__m512i*)rows[1]+2+(0+3*NROWS*NVAL)*NREG*NCH);
						cache[3]=_mm512_load_si512((__m512i*)rows[1]+3+(0+3*NROWS*NVAL)*NREG*NCH);
						cache[4]=_mm512_load_si512((__m512i*)rows[1]+4+(0+3*NROWS*NVAL)*NREG*NCH);
						cache[5]=_mm512_load_si512((__m512i*)rows[1]+5+(0+3*NROWS*NVAL)*NREG*NCH);
						xmin[0]=_mm512_min_epi16(xmin[0], cache[0]); xmax[0]=_mm512_max_epi16(xmax[0], cache[0]);
						xmin[1]=_mm512_min_epi16(xmin[1], cache[1]); xmax[1]=_mm512_max_epi16(xmax[1], cache[1]);
						xmin[2]=_mm512_min_epi16(xmin[2], cache[2]); xmax[2]=_mm512_max_epi16(xmax[2], cache[2]);
						xmin[3]=_mm512_min_epi16(xmin[3], cache[3]); xmax[3]=_mm512_max_epi16(xmax[3], cache[3]);
						xmin[4]=_mm512_min_epi16(xmin[4], cache[4]); xmax[4]=_mm512_max_epi16(xmax[4], cache[4]);
						xmin[5]=_mm512_min_epi16(xmin[5], cache[5]); xmax[5]=_mm512_max_epi16(xmax[5], cache[5]);
					}
				}
				else if(effort==3)
				{
					__m512i cache[6];
					/*
					effort 3
					0	N+W-NW
					1	N
					2	W
					3	W+NE-N
					4	3*(N-NN)+NNN
					5	3*(W-WW)+WWW
					6	N+NE-NNE
					7	NEE
					8	NN
					9	WW
					10	2*N-NN
					11	2*W-WW
					12	NEEE
					13	NEEEE
					14	NNWW
					15	NNEE
					16	N+NW-NNW
					17	W+NW-NWW
					18	(WWWW+NEEEE)>>1
					19	(WWW+NNN+NEEE-NW)>>1

					(__m256i*)rows[-Y]+E*3+C+X*6
					*/

					//N+W-NW
					L1preds[0*6+0]=pred[0];
					L1preds[0*6+1]=pred[1];
					L1preds[0*6+2]=pred[2];
					L1preds[0*6+3]=pred[3];
					L1preds[0*6+4]=pred[4];
					L1preds[0*6+5]=pred[5];

					//N
					L1preds[1*6+0]=N[0];
					L1preds[1*6+1]=N[1];
					L1preds[1*6+2]=N[2];
					L1preds[1*6+3]=N[3];
					L1preds[1*6+4]=N[4];
					L1preds[1*6+5]=N[5];

					//W
					L1preds[2*6+0]=W[0];
					L1preds[2*6+1]=W[1];
					L1preds[2*6+2]=W[2];
					L1preds[2*6+3]=W[3];
					L1preds[2*6+4]=W[4];
					L1preds[2*6+5]=W[5];

					//W+NE-N
					cache[0]=_mm512_add_epi16(W[0], _mm512_load_si512((__m512i*)rows[1]+0+(0+1*NROWS*NVAL)*NREG*NCH));
					cache[1]=_mm512_add_epi16(W[1], _mm512_load_si512((__m512i*)rows[1]+1+(0+1*NROWS*NVAL)*NREG*NCH));
					cache[2]=_mm512_add_epi16(W[2], _mm512_load_si512((__m512i*)rows[1]+2+(0+1*NROWS*NVAL)*NREG*NCH));
					cache[3]=_mm512_add_epi16(W[3], _mm512_load_si512((__m512i*)rows[1]+3+(0+1*NROWS*NVAL)*NREG*NCH));
					cache[4]=_mm512_add_epi16(W[4], _mm512_load_si512((__m512i*)rows[1]+4+(0+1*NROWS*NVAL)*NREG*NCH));
					cache[5]=_mm512_add_epi16(W[5], _mm512_load_si512((__m512i*)rows[1]+5+(0+1*NROWS*NVAL)*NREG*NCH));
					L1preds[3*6+0]=_mm512_sub_epi16(cache[0], N[0]);
					L1preds[3*6+1]=_mm512_sub_epi16(cache[1], N[1]);
					L1preds[3*6+2]=_mm512_sub_epi16(cache[2], N[2]);
					L1preds[3*6+3]=_mm512_sub_epi16(cache[3], N[3]);
					L1preds[3*6+4]=_mm512_sub_epi16(cache[4], N[4]);
					L1preds[3*6+5]=_mm512_sub_epi16(cache[5], N[5]);

					//3*(N-NN)+NNN
					cache[0]=_mm512_sub_epi16(N[0], _mm512_load_si512((__m512i*)rows[2]+0+(0+0*NROWS*NVAL)*NREG*NCH));//N-NN
					cache[1]=_mm512_sub_epi16(N[1], _mm512_load_si512((__m512i*)rows[2]+1+(0+0*NROWS*NVAL)*NREG*NCH));
					cache[2]=_mm512_sub_epi16(N[2], _mm512_load_si512((__m512i*)rows[2]+2+(0+0*NROWS*NVAL)*NREG*NCH));
					cache[3]=_mm512_sub_epi16(N[3], _mm512_load_si512((__m512i*)rows[2]+3+(0+0*NROWS*NVAL)*NREG*NCH));
					cache[4]=_mm512_sub_epi16(N[4], _mm512_load_si512((__m512i*)rows[2]+4+(0+0*NROWS*NVAL)*NREG*NCH));
					cache[5]=_mm512_sub_epi16(N[5], _mm512_load_si512((__m512i*)rows[2]+5+(0+0*NROWS*NVAL)*NREG*NCH));
					cache[0]=_mm512_add_epi16(cache[0], _mm512_slli_epi16(cache[0], 1));//*3
					cache[1]=_mm512_add_epi16(cache[1], _mm512_slli_epi16(cache[1], 1));
					cache[2]=_mm512_add_epi16(cache[2], _mm512_slli_epi16(cache[2], 1));
					cache[3]=_mm512_add_epi16(cache[3], _mm512_slli_epi16(cache[3], 1));
					cache[4]=_mm512_add_epi16(cache[4], _mm512_slli_epi16(cache[4], 1));
					cache[5]=_mm512_add_epi16(cache[5], _mm512_slli_epi16(cache[5], 1));
					L1preds[4*6+0]=_mm512_add_epi16(cache[0], _mm512_load_si512((__m512i*)rows[3]+0+(0+0*NROWS*NVAL)*NREG*NCH));//+NNN
					L1preds[4*6+1]=_mm512_add_epi16(cache[1], _mm512_load_si512((__m512i*)rows[3]+1+(0+0*NROWS*NVAL)*NREG*NCH));
					L1preds[4*6+2]=_mm512_add_epi16(cache[2], _mm512_load_si512((__m512i*)rows[3]+2+(0+0*NROWS*NVAL)*NREG*NCH));
					L1preds[4*6+3]=_mm512_add_epi16(cache[3], _mm512_load_si512((__m512i*)rows[3]+3+(0+0*NROWS*NVAL)*NREG*NCH));
					L1preds[4*6+4]=_mm512_add_epi16(cache[4], _mm512_load_si512((__m512i*)rows[3]+4+(0+0*NROWS*NVAL)*NREG*NCH));
					L1preds[4*6+5]=_mm512_add_epi16(cache[5], _mm512_load_si512((__m512i*)rows[3]+5+(0+0*NROWS*NVAL)*NREG*NCH));

					//3*(W-WW)+WWW
					cache[0]=_mm512_sub_epi16(W[0], _mm512_load_si512((__m512i*)rows[0]+0+(0-2*NROWS*NVAL)*NREG*NCH));//W-WW
					cache[1]=_mm512_sub_epi16(W[1], _mm512_load_si512((__m512i*)rows[0]+1+(0-2*NROWS*NVAL)*NREG*NCH));
					cache[2]=_mm512_sub_epi16(W[2], _mm512_load_si512((__m512i*)rows[0]+2+(0-2*NROWS*NVAL)*NREG*NCH));
					cache[3]=_mm512_sub_epi16(W[3], _mm512_load_si512((__m512i*)rows[0]+3+(0-2*NROWS*NVAL)*NREG*NCH));
					cache[4]=_mm512_sub_epi16(W[4], _mm512_load_si512((__m512i*)rows[0]+4+(0-2*NROWS*NVAL)*NREG*NCH));
					cache[5]=_mm512_sub_epi16(W[5], _mm512_load_si512((__m512i*)rows[0]+5+(0-2*NROWS*NVAL)*NREG*NCH));
					cache[0]=_mm512_add_epi16(cache[0], _mm512_slli_epi16(cache[0], 1));//*3
					cache[1]=_mm512_add_epi16(cache[1], _mm512_slli_epi16(cache[1], 1));
					cache[2]=_mm512_add_epi16(cache[2], _mm512_slli_epi16(cache[2], 1));
					cache[3]=_mm512_add_epi16(cache[3], _mm512_slli_epi16(cache[3], 1));
					cache[4]=_mm512_add_epi16(cache[4], _mm512_slli_epi16(cache[4], 1));
					cache[5]=_mm512_add_epi16(cache[5], _mm512_slli_epi16(cache[5], 1));
					L1preds[5*6+0]=_mm512_add_epi16(cache[0], _mm512_load_si512((__m512i*)rows[0]+0+(0-3*NROWS*NVAL)*NREG*NCH));//+WWW
					L1preds[5*6+1]=_mm512_add_epi16(cache[1], _mm512_load_si512((__m512i*)rows[0]+1+(0-3*NROWS*NVAL)*NREG*NCH));
					L1preds[5*6+2]=_mm512_add_epi16(cache[2], _mm512_load_si512((__m512i*)rows[0]+2+(0-3*NROWS*NVAL)*NREG*NCH));
					L1preds[5*6+3]=_mm512_add_epi16(cache[3], _mm512_load_si512((__m512i*)rows[0]+3+(0-3*NROWS*NVAL)*NREG*NCH));
					L1preds[5*6+4]=_mm512_add_epi16(cache[4], _mm512_load_si512((__m512i*)rows[0]+4+(0-3*NROWS*NVAL)*NREG*NCH));
					L1preds[5*6+5]=_mm512_add_epi16(cache[5], _mm512_load_si512((__m512i*)rows[0]+5+(0-3*NROWS*NVAL)*NREG*NCH));

					//N+NE-NNE
					cache[0]=_mm512_add_epi16(N[0], _mm512_load_si512((__m512i*)rows[1]+0+(0+1*NROWS*NVAL)*NREG*NCH));//N+NE
					cache[1]=_mm512_add_epi16(N[1], _mm512_load_si512((__m512i*)rows[1]+1+(0+1*NROWS*NVAL)*NREG*NCH));
					cache[2]=_mm512_add_epi16(N[2], _mm512_load_si512((__m512i*)rows[1]+2+(0+1*NROWS*NVAL)*NREG*NCH));
					cache[3]=_mm512_add_epi16(N[3], _mm512_load_si512((__m512i*)rows[1]+3+(0+1*NROWS*NVAL)*NREG*NCH));
					cache[4]=_mm512_add_epi16(N[4], _mm512_load_si512((__m512i*)rows[1]+4+(0+1*NROWS*NVAL)*NREG*NCH));
					cache[5]=_mm512_add_epi16(N[5], _mm512_load_si512((__m512i*)rows[1]+5+(0+1*NROWS*NVAL)*NREG*NCH));
					L1preds[6*6+0]=_mm512_sub_epi16(cache[0], _mm512_load_si512((__m512i*)rows[2]+0+(0+1*NROWS*NVAL)*NREG*NCH));//NNE
					L1preds[6*6+1]=_mm512_sub_epi16(cache[1], _mm512_load_si512((__m512i*)rows[2]+1+(0+1*NROWS*NVAL)*NREG*NCH));
					L1preds[6*6+2]=_mm512_sub_epi16(cache[2], _mm512_load_si512((__m512i*)rows[2]+2+(0+1*NROWS*NVAL)*NREG*NCH));
					L1preds[6*6+3]=_mm512_sub_epi16(cache[3], _mm512_load_si512((__m512i*)rows[2]+3+(0+1*NROWS*NVAL)*NREG*NCH));
					L1preds[6*6+4]=_mm512_sub_epi16(cache[4], _mm512_load_si512((__m512i*)rows[2]+4+(0+1*NROWS*NVAL)*NREG*NCH));
					L1preds[6*6+5]=_mm512_sub_epi16(cache[5], _mm512_load_si512((__m512i*)rows[2]+5+(0+1*NROWS*NVAL)*NREG*NCH));

					//NEE
					L1preds[7*6+0]=_mm512_load_si512((__m512i*)rows[1]+0+(0+2*NROWS*NVAL)*NREG*NCH);
					L1preds[7*6+1]=_mm512_load_si512((__m512i*)rows[1]+1+(0+2*NROWS*NVAL)*NREG*NCH);
					L1preds[7*6+2]=_mm512_load_si512((__m512i*)rows[1]+2+(0+2*NROWS*NVAL)*NREG*NCH);
					L1preds[7*6+3]=_mm512_load_si512((__m512i*)rows[1]+3+(0+2*NROWS*NVAL)*NREG*NCH);
					L1preds[7*6+4]=_mm512_load_si512((__m512i*)rows[1]+4+(0+2*NROWS*NVAL)*NREG*NCH);
					L1preds[7*6+5]=_mm512_load_si512((__m512i*)rows[1]+5+(0+2*NROWS*NVAL)*NREG*NCH);

					//NN
					L1preds[8*6+0]=_mm512_load_si512((__m512i*)rows[2]+0+(0+0*NROWS*NVAL)*NREG*NCH);
					L1preds[8*6+1]=_mm512_load_si512((__m512i*)rows[2]+1+(0+0*NROWS*NVAL)*NREG*NCH);
					L1preds[8*6+2]=_mm512_load_si512((__m512i*)rows[2]+2+(0+0*NROWS*NVAL)*NREG*NCH);
					L1preds[8*6+3]=_mm512_load_si512((__m512i*)rows[2]+3+(0+0*NROWS*NVAL)*NREG*NCH);
					L1preds[8*6+4]=_mm512_load_si512((__m512i*)rows[2]+4+(0+0*NROWS*NVAL)*NREG*NCH);
					L1preds[8*6+5]=_mm512_load_si512((__m512i*)rows[2]+5+(0+0*NROWS*NVAL)*NREG*NCH);

					//WW
					L1preds[9*6+0]=_mm512_load_si512((__m512i*)rows[0]+0+(0-2*NROWS*NVAL)*NREG*NCH);
					L1preds[9*6+1]=_mm512_load_si512((__m512i*)rows[0]+1+(0-2*NROWS*NVAL)*NREG*NCH);
					L1preds[9*6+2]=_mm512_load_si512((__m512i*)rows[0]+2+(0-2*NROWS*NVAL)*NREG*NCH);
					L1preds[9*6+3]=_mm512_load_si512((__m512i*)rows[0]+3+(0-2*NROWS*NVAL)*NREG*NCH);
					L1preds[9*6+4]=_mm512_load_si512((__m512i*)rows[0]+4+(0-2*NROWS*NVAL)*NREG*NCH);
					L1preds[9*6+5]=_mm512_load_si512((__m512i*)rows[0]+5+(0-2*NROWS*NVAL)*NREG*NCH);

					//2*N-NN
					L1preds[10*6+0]=_mm512_sub_epi16(_mm512_slli_epi16(N[0], 1), _mm512_load_si512((__m512i*)rows[2]+0+(0+0*NROWS*NVAL)*NREG*NCH));
					L1preds[10*6+1]=_mm512_sub_epi16(_mm512_slli_epi16(N[1], 1), _mm512_load_si512((__m512i*)rows[2]+1+(0+0*NROWS*NVAL)*NREG*NCH));
					L1preds[10*6+2]=_mm512_sub_epi16(_mm512_slli_epi16(N[2], 1), _mm512_load_si512((__m512i*)rows[2]+2+(0+0*NROWS*NVAL)*NREG*NCH));
					L1preds[10*6+3]=_mm512_sub_epi16(_mm512_slli_epi16(N[3], 1), _mm512_load_si512((__m512i*)rows[2]+3+(0+0*NROWS*NVAL)*NREG*NCH));
					L1preds[10*6+4]=_mm512_sub_epi16(_mm512_slli_epi16(N[4], 1), _mm512_load_si512((__m512i*)rows[2]+4+(0+0*NROWS*NVAL)*NREG*NCH));
					L1preds[10*6+5]=_mm512_sub_epi16(_mm512_slli_epi16(N[5], 1), _mm512_load_si512((__m512i*)rows[2]+5+(0+0*NROWS*NVAL)*NREG*NCH));

					//2*W-WW
					L1preds[11*6+0]=_mm512_sub_epi16(_mm512_slli_epi16(W[0], 1), _mm512_load_si512((__m512i*)rows[0]+0+(0-2*NROWS*NVAL)*NREG*NCH));
					L1preds[11*6+1]=_mm512_sub_epi16(_mm512_slli_epi16(W[1], 1), _mm512_load_si512((__m512i*)rows[0]+1+(0-2*NROWS*NVAL)*NREG*NCH));
					L1preds[11*6+2]=_mm512_sub_epi16(_mm512_slli_epi16(W[2], 1), _mm512_load_si512((__m512i*)rows[0]+2+(0-2*NROWS*NVAL)*NREG*NCH));
					L1preds[11*6+3]=_mm512_sub_epi16(_mm512_slli_epi16(W[3], 1), _mm512_load_si512((__m512i*)rows[0]+3+(0-2*NROWS*NVAL)*NREG*NCH));
					L1preds[11*6+4]=_mm512_sub_epi16(_mm512_slli_epi16(W[4], 1), _mm512_load_si512((__m512i*)rows[0]+4+(0-2*NROWS*NVAL)*NREG*NCH));
					L1preds[11*6+5]=_mm512_sub_epi16(_mm512_slli_epi16(W[5], 1), _mm512_load_si512((__m512i*)rows[0]+5+(0-2*NROWS*NVAL)*NREG*NCH));

					//NEEE
					L1preds[12*6+0]=_mm512_load_si512((__m512i*)rows[1]+0+(0+3*NROWS*NVAL)*NREG*NCH);
					L1preds[12*6+1]=_mm512_load_si512((__m512i*)rows[1]+1+(0+3*NROWS*NVAL)*NREG*NCH);
					L1preds[12*6+2]=_mm512_load_si512((__m512i*)rows[1]+2+(0+3*NROWS*NVAL)*NREG*NCH);
					L1preds[12*6+3]=_mm512_load_si512((__m512i*)rows[1]+3+(0+3*NROWS*NVAL)*NREG*NCH);
					L1preds[12*6+4]=_mm512_load_si512((__m512i*)rows[1]+4+(0+3*NROWS*NVAL)*NREG*NCH);
					L1preds[12*6+5]=_mm512_load_si512((__m512i*)rows[1]+5+(0+3*NROWS*NVAL)*NREG*NCH);

					//NEEEE
					L1preds[13*6+0]=_mm512_load_si512((__m512i*)rows[1]+0+(0+4*NROWS*NVAL)*NREG*NCH);
					L1preds[13*6+1]=_mm512_load_si512((__m512i*)rows[1]+1+(0+4*NROWS*NVAL)*NREG*NCH);
					L1preds[13*6+2]=_mm512_load_si512((__m512i*)rows[1]+2+(0+4*NROWS*NVAL)*NREG*NCH);
					L1preds[13*6+3]=_mm512_load_si512((__m512i*)rows[1]+3+(0+4*NROWS*NVAL)*NREG*NCH);
					L1preds[13*6+4]=_mm512_load_si512((__m512i*)rows[1]+4+(0+4*NROWS*NVAL)*NREG*NCH);
					L1preds[13*6+5]=_mm512_load_si512((__m512i*)rows[1]+5+(0+4*NROWS*NVAL)*NREG*NCH);

					//NNWW
					L1preds[14*6+0]=_mm512_load_si512((__m512i*)rows[2]+0+(0-2*NROWS*NVAL)*NREG*NCH);
					L1preds[14*6+1]=_mm512_load_si512((__m512i*)rows[2]+1+(0-2*NROWS*NVAL)*NREG*NCH);
					L1preds[14*6+2]=_mm512_load_si512((__m512i*)rows[2]+2+(0-2*NROWS*NVAL)*NREG*NCH);
					L1preds[14*6+3]=_mm512_load_si512((__m512i*)rows[2]+3+(0-2*NROWS*NVAL)*NREG*NCH);
					L1preds[14*6+4]=_mm512_load_si512((__m512i*)rows[2]+4+(0-2*NROWS*NVAL)*NREG*NCH);
					L1preds[14*6+5]=_mm512_load_si512((__m512i*)rows[2]+5+(0-2*NROWS*NVAL)*NREG*NCH);

					//NNEE
					L1preds[15*6+0]=_mm512_load_si512((__m512i*)rows[2]+0+(0+2*NROWS*NVAL)*NREG*NCH);
					L1preds[15*6+1]=_mm512_load_si512((__m512i*)rows[2]+1+(0+2*NROWS*NVAL)*NREG*NCH);
					L1preds[15*6+2]=_mm512_load_si512((__m512i*)rows[2]+2+(0+2*NROWS*NVAL)*NREG*NCH);
					L1preds[15*6+3]=_mm512_load_si512((__m512i*)rows[2]+3+(0+2*NROWS*NVAL)*NREG*NCH);
					L1preds[15*6+4]=_mm512_load_si512((__m512i*)rows[2]+4+(0+2*NROWS*NVAL)*NREG*NCH);
					L1preds[15*6+5]=_mm512_load_si512((__m512i*)rows[2]+5+(0+2*NROWS*NVAL)*NREG*NCH);

					//N+NW-NNW
					L1preds[16*6+0]=_mm512_sub_epi16(_mm512_add_epi16(N[0], NW[0]), _mm512_load_si512((__m512i*)rows[2]+0+(0-1*NROWS*NVAL)*NREG*NCH));
					L1preds[16*6+1]=_mm512_sub_epi16(_mm512_add_epi16(N[1], NW[1]), _mm512_load_si512((__m512i*)rows[2]+1+(0-1*NROWS*NVAL)*NREG*NCH));
					L1preds[16*6+2]=_mm512_sub_epi16(_mm512_add_epi16(N[2], NW[2]), _mm512_load_si512((__m512i*)rows[2]+2+(0-1*NROWS*NVAL)*NREG*NCH));
					L1preds[16*6+3]=_mm512_sub_epi16(_mm512_add_epi16(N[3], NW[3]), _mm512_load_si512((__m512i*)rows[2]+3+(0-1*NROWS*NVAL)*NREG*NCH));
					L1preds[16*6+4]=_mm512_sub_epi16(_mm512_add_epi16(N[4], NW[4]), _mm512_load_si512((__m512i*)rows[2]+4+(0-1*NROWS*NVAL)*NREG*NCH));
					L1preds[16*6+5]=_mm512_sub_epi16(_mm512_add_epi16(N[5], NW[5]), _mm512_load_si512((__m512i*)rows[2]+5+(0-1*NROWS*NVAL)*NREG*NCH));

					//W+NW-NWW
					L1preds[17*6+0]=_mm512_sub_epi16(_mm512_add_epi16(W[0], NW[0]), _mm512_load_si512((__m512i*)rows[1]+0+(0-2*NROWS*NVAL)*NREG*NCH));
					L1preds[17*6+1]=_mm512_sub_epi16(_mm512_add_epi16(W[1], NW[1]), _mm512_load_si512((__m512i*)rows[1]+1+(0-2*NROWS*NVAL)*NREG*NCH));
					L1preds[17*6+2]=_mm512_sub_epi16(_mm512_add_epi16(W[2], NW[2]), _mm512_load_si512((__m512i*)rows[1]+2+(0-2*NROWS*NVAL)*NREG*NCH));
					L1preds[17*6+3]=_mm512_sub_epi16(_mm512_add_epi16(W[3], NW[3]), _mm512_load_si512((__m512i*)rows[1]+3+(0-2*NROWS*NVAL)*NREG*NCH));
					L1preds[17*6+4]=_mm512_sub_epi16(_mm512_add_epi16(W[4], NW[4]), _mm512_load_si512((__m512i*)rows[1]+4+(0-2*NROWS*NVAL)*NREG*NCH));
					L1preds[17*6+5]=_mm512_sub_epi16(_mm512_add_epi16(W[5], NW[5]), _mm512_load_si512((__m512i*)rows[1]+5+(0-2*NROWS*NVAL)*NREG*NCH));

					//(WWWW+NEEEE)>>1
					L1preds[18*6+0]=_mm512_srai_epi16(_mm512_add_epi16(_mm512_load_si512((__m512i*)rows[0]+0+(0-4*NROWS*NVAL)*NREG*NCH), _mm512_load_si512((__m512i*)rows[1]+0+(0+4*NROWS*NVAL)*NREG*NCH)), 1);
					L1preds[18*6+1]=_mm512_srai_epi16(_mm512_add_epi16(_mm512_load_si512((__m512i*)rows[0]+1+(0-4*NROWS*NVAL)*NREG*NCH), _mm512_load_si512((__m512i*)rows[1]+1+(0+4*NROWS*NVAL)*NREG*NCH)), 1);
					L1preds[18*6+2]=_mm512_srai_epi16(_mm512_add_epi16(_mm512_load_si512((__m512i*)rows[0]+2+(0-4*NROWS*NVAL)*NREG*NCH), _mm512_load_si512((__m512i*)rows[1]+2+(0+4*NROWS*NVAL)*NREG*NCH)), 1);
					L1preds[18*6+3]=_mm512_srai_epi16(_mm512_add_epi16(_mm512_load_si512((__m512i*)rows[0]+3+(0-4*NROWS*NVAL)*NREG*NCH), _mm512_load_si512((__m512i*)rows[1]+3+(0+4*NROWS*NVAL)*NREG*NCH)), 1);
					L1preds[18*6+4]=_mm512_srai_epi16(_mm512_add_epi16(_mm512_load_si512((__m512i*)rows[0]+4+(0-4*NROWS*NVAL)*NREG*NCH), _mm512_load_si512((__m512i*)rows[1]+4+(0+4*NROWS*NVAL)*NREG*NCH)), 1);
					L1preds[18*6+5]=_mm512_srai_epi16(_mm512_add_epi16(_mm512_load_si512((__m512i*)rows[0]+5+(0-4*NROWS*NVAL)*NREG*NCH), _mm512_load_si512((__m512i*)rows[1]+5+(0+4*NROWS*NVAL)*NREG*NCH)), 1);

					//(WWW+NNN+NEEE-NW)>>1
					cache[0]=_mm512_add_epi16(_mm512_load_si512((__m512i*)rows[0]+0+(0-3*NROWS*NVAL)*NREG*NCH), _mm512_load_si512((__m512i*)rows[3]+0+(0+0*NROWS*NVAL)*NREG*NCH));//WWW+NNN
					cache[1]=_mm512_add_epi16(_mm512_load_si512((__m512i*)rows[0]+1+(0-3*NROWS*NVAL)*NREG*NCH), _mm512_load_si512((__m512i*)rows[3]+1+(0+0*NROWS*NVAL)*NREG*NCH));
					cache[2]=_mm512_add_epi16(_mm512_load_si512((__m512i*)rows[0]+2+(0-3*NROWS*NVAL)*NREG*NCH), _mm512_load_si512((__m512i*)rows[3]+2+(0+0*NROWS*NVAL)*NREG*NCH));
					cache[3]=_mm512_add_epi16(_mm512_load_si512((__m512i*)rows[0]+3+(0-3*NROWS*NVAL)*NREG*NCH), _mm512_load_si512((__m512i*)rows[3]+3+(0+0*NROWS*NVAL)*NREG*NCH));
					cache[4]=_mm512_add_epi16(_mm512_load_si512((__m512i*)rows[0]+4+(0-3*NROWS*NVAL)*NREG*NCH), _mm512_load_si512((__m512i*)rows[3]+4+(0+0*NROWS*NVAL)*NREG*NCH));
					cache[5]=_mm512_add_epi16(_mm512_load_si512((__m512i*)rows[0]+5+(0-3*NROWS*NVAL)*NREG*NCH), _mm512_load_si512((__m512i*)rows[3]+5+(0+0*NROWS*NVAL)*NREG*NCH));
					cache[0]=_mm512_add_epi16(cache[0], _mm512_load_si512((__m512i*)rows[1]+0+(0+3*NROWS*NVAL)*NREG*NCH));//+NEEE
					cache[1]=_mm512_add_epi16(cache[1], _mm512_load_si512((__m512i*)rows[1]+1+(0+3*NROWS*NVAL)*NREG*NCH));
					cache[2]=_mm512_add_epi16(cache[2], _mm512_load_si512((__m512i*)rows[1]+2+(0+3*NROWS*NVAL)*NREG*NCH));
					cache[3]=_mm512_add_epi16(cache[3], _mm512_load_si512((__m512i*)rows[1]+3+(0+3*NROWS*NVAL)*NREG*NCH));
					cache[4]=_mm512_add_epi16(cache[4], _mm512_load_si512((__m512i*)rows[1]+4+(0+3*NROWS*NVAL)*NREG*NCH));
					cache[5]=_mm512_add_epi16(cache[5], _mm512_load_si512((__m512i*)rows[1]+5+(0+3*NROWS*NVAL)*NREG*NCH));
					cache[0]=_mm512_sub_epi16(cache[0], NW[0]);//-NW
					cache[1]=_mm512_sub_epi16(cache[1], NW[1]);
					cache[2]=_mm512_sub_epi16(cache[2], NW[2]);
					cache[3]=_mm512_sub_epi16(cache[3], NW[3]);
					cache[4]=_mm512_sub_epi16(cache[4], NW[4]);
					cache[5]=_mm512_sub_epi16(cache[5], NW[5]);
					L1preds[19*6+0]=_mm512_srai_epi16(cache[0], 1);
					L1preds[19*6+1]=_mm512_srai_epi16(cache[1], 1);
					L1preds[19*6+2]=_mm512_srai_epi16(cache[2], 1);
					L1preds[19*6+3]=_mm512_srai_epi16(cache[3], 1);
					L1preds[19*6+4]=_mm512_srai_epi16(cache[4], 1);
					L1preds[19*6+5]=_mm512_srai_epi16(cache[5], 1);


					//mix
					__m512i mp[12], t[12];
					mp[0x0]=_mm512_setzero_si512();
					mp[0x1]=_mm512_setzero_si512();
					mp[0x2]=_mm512_setzero_si512();
					mp[0x3]=_mm512_setzero_si512();
					mp[0x4]=_mm512_setzero_si512();
					mp[0x5]=_mm512_setzero_si512();
					mp[0x6]=_mm512_setzero_si512();
					mp[0x7]=_mm512_setzero_si512();
					mp[0x8]=_mm512_setzero_si512();
					mp[0x9]=_mm512_setzero_si512();
					mp[0xA]=_mm512_setzero_si512();
					mp[0xB]=_mm512_setzero_si512();
					for(int k=0;k<L1_NPREDS3;++k)
					{
						//16 -> 32		3 lo 3 hi registers
						t[0x0]=_mm512_slli_epi32(L1preds[k*6+0], 16);
						t[0x1]=_mm512_slli_epi32(L1preds[k*6+1], 16);
						t[0x2]=_mm512_slli_epi32(L1preds[k*6+2], 16);
						t[0x3]=_mm512_slli_epi32(L1preds[k*6+3], 16);
						t[0x4]=_mm512_slli_epi32(L1preds[k*6+4], 16);
						t[0x5]=_mm512_slli_epi32(L1preds[k*6+5], 16);
						t[0x6]=_mm512_srai_epi32(L1preds[k*6+0], 16);
						t[0x7]=_mm512_srai_epi32(L1preds[k*6+1], 16);
						t[0x8]=_mm512_srai_epi32(L1preds[k*6+2], 16);
						t[0x9]=_mm512_srai_epi32(L1preds[k*6+3], 16);
						t[0xA]=_mm512_srai_epi32(L1preds[k*6+4], 16);
						t[0xB]=_mm512_srai_epi32(L1preds[k*6+5], 16);
						t[0x0]=_mm512_srai_epi32(t[0x0], 16);
						t[0x1]=_mm512_srai_epi32(t[0x1], 16);
						t[0x2]=_mm512_srai_epi32(t[0x2], 16);
						t[0x3]=_mm512_srai_epi32(t[0x3], 16);
						t[0x4]=_mm512_srai_epi32(t[0x4], 16);
						t[0x5]=_mm512_srai_epi32(t[0x5], 16);
						t[0x0]=_mm512_mullo_epi32(t[0x0], _mm512_load_si512((__m512i*)L1weights+k*12+0x0));
						t[0x1]=_mm512_mullo_epi32(t[0x1], _mm512_load_si512((__m512i*)L1weights+k*12+0x1));
						t[0x2]=_mm512_mullo_epi32(t[0x2], _mm512_load_si512((__m512i*)L1weights+k*12+0x2));
						t[0x3]=_mm512_mullo_epi32(t[0x3], _mm512_load_si512((__m512i*)L1weights+k*12+0x3));
						t[0x4]=_mm512_mullo_epi32(t[0x4], _mm512_load_si512((__m512i*)L1weights+k*12+0x4));
						t[0x5]=_mm512_mullo_epi32(t[0x5], _mm512_load_si512((__m512i*)L1weights+k*12+0x5));
						t[0x6]=_mm512_mullo_epi32(t[0x6], _mm512_load_si512((__m512i*)L1weights+k*12+0x6));
						t[0x7]=_mm512_mullo_epi32(t[0x7], _mm512_load_si512((__m512i*)L1weights+k*12+0x7));
						t[0x8]=_mm512_mullo_epi32(t[0x8], _mm512_load_si512((__m512i*)L1weights+k*12+0x8));
						t[0x9]=_mm512_mullo_epi32(t[0x9], _mm512_load_si512((__m512i*)L1weights+k*12+0x9));
						t[0xA]=_mm512_mullo_epi32(t[0xA], _mm512_load_si512((__m512i*)L1weights+k*12+0xA));
						t[0xB]=_mm512_mullo_epi32(t[0xB], _mm512_load_si512((__m512i*)L1weights+k*12+0xB));
						mp[0x0]=_mm512_add_epi32(mp[0x0], t[0x0]);
						mp[0x1]=_mm512_add_epi32(mp[0x1], t[0x1]);
						mp[0x2]=_mm512_add_epi32(mp[0x2], t[0x2]);
						mp[0x3]=_mm512_add_epi32(mp[0x3], t[0x3]);
						mp[0x4]=_mm512_add_epi32(mp[0x4], t[0x4]);
						mp[0x5]=_mm512_add_epi32(mp[0x5], t[0x5]);
						mp[0x6]=_mm512_add_epi32(mp[0x6], t[0x6]);
						mp[0x7]=_mm512_add_epi32(mp[0x7], t[0x7]);
						mp[0x8]=_mm512_add_epi32(mp[0x8], t[0x8]);
						mp[0x9]=_mm512_add_epi32(mp[0x9], t[0x9]);
						mp[0xA]=_mm512_add_epi32(mp[0xA], t[0xA]);
						mp[0xB]=_mm512_add_epi32(mp[0xB], t[0xB]);
					}
					__m512i rcon=_mm512_set1_epi32(1<<L1_SH3>>1);
					mp[0x0]=_mm512_add_epi32(mp[0x0], rcon);//rounding to nearest
					mp[0x1]=_mm512_add_epi32(mp[0x1], rcon);
					mp[0x2]=_mm512_add_epi32(mp[0x2], rcon);
					mp[0x3]=_mm512_add_epi32(mp[0x3], rcon);
					mp[0x4]=_mm512_add_epi32(mp[0x4], rcon);
					mp[0x5]=_mm512_add_epi32(mp[0x5], rcon);
					mp[0x6]=_mm512_add_epi32(mp[0x6], rcon);
					mp[0x7]=_mm512_add_epi32(mp[0x7], rcon);
					mp[0x8]=_mm512_add_epi32(mp[0x8], rcon);
					mp[0x9]=_mm512_add_epi32(mp[0x9], rcon);
					mp[0xA]=_mm512_add_epi32(mp[0xA], rcon);
					mp[0xB]=_mm512_add_epi32(mp[0xB], rcon);

					mp[0x0]=_mm512_srai_epi32(mp[0x0], L1_SH3);
					mp[0x1]=_mm512_srai_epi32(mp[0x1], L1_SH3);
					mp[0x2]=_mm512_srai_epi32(mp[0x2], L1_SH3);
					mp[0x3]=_mm512_srai_epi32(mp[0x3], L1_SH3);
					mp[0x4]=_mm512_srai_epi32(mp[0x4], L1_SH3);
					mp[0x5]=_mm512_srai_epi32(mp[0x5], L1_SH3);
					mp[0x6]=_mm512_srai_epi32(mp[0x6], L1_SH3-16);
					mp[0x7]=_mm512_srai_epi32(mp[0x7], L1_SH3-16);
					mp[0x8]=_mm512_srai_epi32(mp[0x8], L1_SH3-16);
					mp[0x9]=_mm512_srai_epi32(mp[0x9], L1_SH3-16);
					mp[0xA]=_mm512_srai_epi32(mp[0xA], L1_SH3-16);
					mp[0xB]=_mm512_srai_epi32(mp[0xB], L1_SH3-16);
					//32 -> 16
					pred[0]=_mm512_mask_blend_epi16(0xAAAAAAAA, mp[0x0], mp[0x6]);
					pred[1]=_mm512_mask_blend_epi16(0xAAAAAAAA, mp[0x1], mp[0x7]);
					pred[2]=_mm512_mask_blend_epi16(0xAAAAAAAA, mp[0x2], mp[0x8]);
					pred[3]=_mm512_mask_blend_epi16(0xAAAAAAAA, mp[0x3], mp[0x9]);
					pred[4]=_mm512_mask_blend_epi16(0xAAAAAAAA, mp[0x4], mp[0xA]);
					pred[5]=_mm512_mask_blend_epi16(0xAAAAAAAA, mp[0x5], mp[0xB]);


					//loosen pred range
					if(!cond_cg)
					{
						cache[0]=_mm512_load_si512((__m512i*)rows[1]+0+(0+1*NROWS*NVAL)*NREG*NCH);//NE
						cache[1]=_mm512_load_si512((__m512i*)rows[1]+1+(0+1*NROWS*NVAL)*NREG*NCH);
						cache[2]=_mm512_load_si512((__m512i*)rows[1]+2+(0+1*NROWS*NVAL)*NREG*NCH);
						cache[3]=_mm512_load_si512((__m512i*)rows[1]+3+(0+1*NROWS*NVAL)*NREG*NCH);
						cache[4]=_mm512_load_si512((__m512i*)rows[1]+4+(0+1*NROWS*NVAL)*NREG*NCH);
						cache[5]=_mm512_load_si512((__m512i*)rows[1]+5+(0+1*NROWS*NVAL)*NREG*NCH);
						xmin[0]=_mm512_min_epi16(xmin[0], cache[0]); xmax[0]=_mm512_max_epi16(xmax[0], cache[0]);
						xmin[1]=_mm512_min_epi16(xmin[1], cache[1]); xmax[1]=_mm512_max_epi16(xmax[1], cache[1]);
						xmin[2]=_mm512_min_epi16(xmin[2], cache[2]); xmax[2]=_mm512_max_epi16(xmax[2], cache[2]);
						xmin[3]=_mm512_min_epi16(xmin[3], cache[3]); xmax[3]=_mm512_max_epi16(xmax[3], cache[3]);
						xmin[4]=_mm512_min_epi16(xmin[4], cache[4]); xmax[4]=_mm512_max_epi16(xmax[4], cache[4]);
						xmin[5]=_mm512_min_epi16(xmin[5], cache[5]); xmax[5]=_mm512_max_epi16(xmax[5], cache[5]);
						cache[0]=_mm512_load_si512((__m512i*)rows[1]+0+(0+3*NROWS*NVAL)*NREG*NCH);//NEEE
						cache[1]=_mm512_load_si512((__m512i*)rows[1]+1+(0+3*NROWS*NVAL)*NREG*NCH);
						cache[2]=_mm512_load_si512((__m512i*)rows[1]+2+(0+3*NROWS*NVAL)*NREG*NCH);
						cache[3]=_mm512_load_si512((__m512i*)rows[1]+3+(0+3*NROWS*NVAL)*NREG*NCH);
						cache[4]=_mm512_load_si512((__m512i*)rows[1]+4+(0+3*NROWS*NVAL)*NREG*NCH);
						cache[5]=_mm512_load_si512((__m512i*)rows[1]+5+(0+3*NROWS*NVAL)*NREG*NCH);
						xmin[0]=_mm512_min_epi16(xmin[0], cache[0]); xmax[0]=_mm512_max_epi16(xmax[0], cache[0]);
						xmin[1]=_mm512_min_epi16(xmin[1], cache[1]); xmax[1]=_mm512_max_epi16(xmax[1], cache[1]);
						xmin[2]=_mm512_min_epi16(xmin[2], cache[2]); xmax[2]=_mm512_max_epi16(xmax[2], cache[2]);
						xmin[3]=_mm512_min_epi16(xmin[3], cache[3]); xmax[3]=_mm512_max_epi16(xmax[3], cache[3]);
						xmin[4]=_mm512_min_epi16(xmin[4], cache[4]); xmax[4]=_mm512_max_epi16(xmax[4], cache[4]);
						xmin[5]=_mm512_min_epi16(xmin[5], cache[5]); xmax[5]=_mm512_max_epi16(xmax[5], cache[5]);
					}
				}
				predYUV0[0]=pred[0];
				predYUV0[1]=pred[1];
				predYUV0[2]=pred[2];
				predYUV0[3]=pred[3];
				predYUV0[4]=pred[4];
				predYUV0[5]=pred[5];
				if(cond_cg)
				{
					pred[0]=mcg[0];
					pred[1]=mcg[1];
					pred[2]=mcg[2];
					pred[3]=mcg[3];
					pred[4]=mcg[4];
					pred[5]=mcg[5];
				}

				pred[0]=_mm512_max_epi16(pred[0], xmin[0]); pred[0]=_mm512_min_epi16(pred[0], xmax[0]);
				pred[1]=_mm512_max_epi16(pred[1], xmin[1]); pred[1]=_mm512_min_epi16(pred[1], xmax[1]);
				pred[2]=_mm512_max_epi16(pred[2], xmin[2]); pred[2]=_mm512_min_epi16(pred[2], xmax[2]);
				pred[3]=_mm512_max_epi16(pred[3], xmin[3]); pred[3]=_mm512_min_epi16(pred[3], xmax[3]);
				pred[4]=_mm512_max_epi16(pred[4], xmin[4]); pred[4]=_mm512_min_epi16(pred[4], xmax[4]);
				pred[5]=_mm512_max_epi16(pred[5], xmin[5]); pred[5]=_mm512_min_epi16(pred[5], xmax[5]);
			}
			__m512i msyms[6], moffset[4];
			if(fwd)
			{
				myuv[0]=_mm512_cvtepi8_epi16(_mm256_add_epi8(_mm256_load_si256((__m256i*)(imptr+yidx)+0), half8));//load yuv
				myuv[1]=_mm512_cvtepi8_epi16(_mm256_add_epi8(_mm256_load_si256((__m256i*)(imptr+yidx)+1), half8));
				myuv[2]=_mm512_cvtepi8_epi16(_mm256_add_epi8(_mm256_load_si256((__m256i*)(imptr+uidx)+0), half8));
				myuv[3]=_mm512_cvtepi8_epi16(_mm256_add_epi8(_mm256_load_si256((__m256i*)(imptr+uidx)+1), half8));
				myuv[4]=_mm512_cvtepi8_epi16(_mm256_add_epi8(_mm256_load_si256((__m256i*)(imptr+vidx)+0), half8));
				myuv[5]=_mm512_cvtepi8_epi16(_mm256_add_epi8(_mm256_load_si256((__m256i*)(imptr+vidx)+1), half8));
				
				//decorrelate Y
				msyms[0]=_mm512_sub_epi16(myuv[0], pred[0]);
				msyms[1]=_mm512_sub_epi16(myuv[1], pred[1]);
				if(dist>1)
				{
#ifdef ENABLE_GUIDE
					ALIGN(64) uint16_t diff[32];
					__m512i modified, mdiff0, mdiff1;
					__m512i yuv0=myuv[0*2+0];
					__m512i yuv1=myuv[0*2+1];
#endif
					__m512i t0=_mm512_srai_epi16(msyms[0], 15);
					__m512i t1=_mm512_srai_epi16(msyms[1], 15);
					msyms[0]=_mm512_mulhi_epi16(msyms[0], dist_rcp);
					msyms[1]=_mm512_mulhi_epi16(msyms[1], dist_rcp);
					msyms[0]=_mm512_sub_epi16(msyms[0], t0);
					msyms[1]=_mm512_sub_epi16(msyms[1], t1);
					myuv[0*2+0]=_mm512_mullo_epi16(msyms[0], mdist);
					myuv[0*2+1]=_mm512_mullo_epi16(msyms[1], mdist);
					myuv[0*2+0]=_mm512_add_epi16(myuv[0*2+0], pred[0]);
					myuv[0*2+1]=_mm512_add_epi16(myuv[0*2+1], pred[1]);
					myuv[0*2+0]=_mm512_max_epi16(myuv[0*2+0], amin);
					myuv[0*2+1]=_mm512_max_epi16(myuv[0*2+1], amin);
					myuv[0*2+0]=_mm512_min_epi16(myuv[0*2+0], amax);
					myuv[0*2+1]=_mm512_min_epi16(myuv[0*2+1], amax);
#ifdef ENABLE_GUIDE
					modified=_mm512_inserti64x4(_mm512_castsi256_si512(_mm512_cvtepi16_epi8(myuv[0*2+0])), _mm512_cvtepi16_epi8(myuv[0*2+1]), 1);
					modified=_mm512_xor_si512(modified, _mm512_set1_epi8(-128));
					_mm512_storeu_si512((__m512i*)(g_image+(imptr-interleaved-isize)+yidx), modified);

					mdiff0=_mm512_sub_epi16(yuv0, myuv[0*2+0]);
					mdiff1=_mm512_sub_epi16(yuv1, myuv[0*2+1]);
					mdiff0=_mm512_abs_epi16(mdiff0);
					mdiff1=_mm512_abs_epi16(mdiff1);
					mdiff0=_mm512_mullo_epi16(mdiff0, mdiff0);
					mdiff1=_mm512_mullo_epi16(mdiff1, mdiff1);
					mdiff0=_mm512_add_epi16(mdiff0, mdiff1);
					_mm512_store_si512((__m512i*)diff, mdiff0);
					for(int k=0;k<_countof(diff);++k)
						g_sqe[0]+=diff[k];
#endif
				}
				
				//decorrelate U
				moffset[0]=_mm512_and_si512(myuv[0], uhelpmask);
				moffset[1]=_mm512_and_si512(myuv[1], uhelpmask);
				pred[2]=_mm512_add_epi16(pred[2], moffset[0]);
				pred[3]=_mm512_add_epi16(pred[3], moffset[1]);
				pred[2]=_mm512_max_epi16(pred[2], amin); pred[2]=_mm512_min_epi16(pred[2], amax);
				pred[3]=_mm512_max_epi16(pred[3], amin); pred[3]=_mm512_min_epi16(pred[3], amax);
				msyms[2]=_mm512_sub_epi16(myuv[2], pred[2]);
				msyms[3]=_mm512_sub_epi16(myuv[3], pred[3]);
				if(dist>1)
				{
#ifdef ENABLE_GUIDE
					ALIGN(64) uint16_t diff[32];
					__m512i modified, mdiff0, mdiff1;
					__m512i yuv0=myuv[1*2+0];
					__m512i yuv1=myuv[1*2+1];
#endif
					__m512i t0=_mm512_srai_epi16(msyms[2], 15);
					__m512i t1=_mm512_srai_epi16(msyms[3], 15);
					msyms[2]=_mm512_mulhi_epi16(msyms[2], dist_rcp);
					msyms[3]=_mm512_mulhi_epi16(msyms[3], dist_rcp);
					msyms[2]=_mm512_sub_epi16(msyms[2], t0);
					msyms[3]=_mm512_sub_epi16(msyms[3], t1);
					myuv[1*2+0]=_mm512_mullo_epi16(msyms[2], mdist);
					myuv[1*2+1]=_mm512_mullo_epi16(msyms[3], mdist);
					myuv[1*2+0]=_mm512_add_epi16(myuv[1*2+0], pred[2]);
					myuv[1*2+1]=_mm512_add_epi16(myuv[1*2+1], pred[3]);
					myuv[1*2+0]=_mm512_max_epi16(myuv[1*2+0], amin);
					myuv[1*2+1]=_mm512_max_epi16(myuv[1*2+1], amin);
					myuv[1*2+0]=_mm512_min_epi16(myuv[1*2+0], amax);
					myuv[1*2+1]=_mm512_min_epi16(myuv[1*2+1], amax);
#ifdef ENABLE_GUIDE
					modified=_mm512_inserti64x4(_mm512_castsi256_si512(_mm512_cvtepi16_epi8(myuv[1*2+0])), _mm512_cvtepi16_epi8(myuv[1*2+1]), 1);
					modified=_mm512_xor_si512(modified, _mm512_set1_epi8(-128));
					_mm512_storeu_si512((__m512i*)(g_image+(imptr-interleaved-isize)+uidx), modified);

					mdiff0=_mm512_sub_epi16(yuv0, myuv[1*2+0]);
					mdiff1=_mm512_sub_epi16(yuv1, myuv[1*2+1]);
					mdiff0=_mm512_abs_epi16(mdiff0);
					mdiff1=_mm512_abs_epi16(mdiff1);
					mdiff0=_mm512_mullo_epi16(mdiff0, mdiff0);
					mdiff1=_mm512_mullo_epi16(mdiff1, mdiff1);
					mdiff0=_mm512_add_epi16(mdiff0, mdiff1);
					_mm512_store_si512((__m512i*)diff, mdiff0);
					for(int k=0;k<_countof(diff);++k)
						g_sqe[1]+=diff[k];
#endif
				}

				//decorrelate V
				moffset[2]=_mm512_mullo_epi16(vc0, myuv[0]);
				moffset[3]=_mm512_mullo_epi16(vc0, myuv[1]);
				moffset[2]=_mm512_add_epi16(moffset[2], _mm512_mullo_epi16(vc1, myuv[2]));
				moffset[3]=_mm512_add_epi16(moffset[3], _mm512_mullo_epi16(vc1, myuv[3]));
				moffset[2]=_mm512_srai_epi16(moffset[2], 2);
				moffset[3]=_mm512_srai_epi16(moffset[3], 2);
				pred[4]=_mm512_add_epi16(pred[4], moffset[2]);
				pred[5]=_mm512_add_epi16(pred[5], moffset[3]);
				pred[4]=_mm512_max_epi16(pred[4], amin); pred[4]=_mm512_min_epi16(pred[4], amax);
				pred[5]=_mm512_max_epi16(pred[5], amin); pred[5]=_mm512_min_epi16(pred[5], amax);
				msyms[4]=_mm512_sub_epi16(myuv[4], pred[4]);
				msyms[5]=_mm512_sub_epi16(myuv[5], pred[5]);
				if(dist>1)
				{
#ifdef ENABLE_GUIDE
					ALIGN(64) uint16_t diff[32];
					__m512i modified, mdiff0, mdiff1;
					__m512i yuv0=myuv[2*2+0];
					__m512i yuv1=myuv[2*2+1];
#endif
					__m512i t0=_mm512_srai_epi16(msyms[4], 15);
					__m512i t1=_mm512_srai_epi16(msyms[5], 15);
					msyms[4]=_mm512_mulhi_epi16(msyms[4], dist_rcp);
					msyms[5]=_mm512_mulhi_epi16(msyms[5], dist_rcp);
					msyms[4]=_mm512_sub_epi16(msyms[4], t0);
					msyms[5]=_mm512_sub_epi16(msyms[5], t1);
					myuv[2*2+0]=_mm512_mullo_epi16(msyms[4], mdist);
					myuv[2*2+1]=_mm512_mullo_epi16(msyms[5], mdist);
					myuv[2*2+0]=_mm512_add_epi16(myuv[2*2+0], pred[4]);
					myuv[2*2+1]=_mm512_add_epi16(myuv[2*2+1], pred[5]);
					myuv[2*2+0]=_mm512_max_epi16(myuv[2*2+0], amin);
					myuv[2*2+1]=_mm512_max_epi16(myuv[2*2+1], amin);
					myuv[2*2+0]=_mm512_min_epi16(myuv[2*2+0], amax);
					myuv[2*2+1]=_mm512_min_epi16(myuv[2*2+1], amax);
#ifdef ENABLE_GUIDE
					modified=_mm512_inserti64x4(_mm512_castsi256_si512(_mm512_cvtepi16_epi8(myuv[2*2+0])), _mm512_cvtepi16_epi8(myuv[2*2+1]), 1);
					modified=_mm512_xor_si512(modified, _mm512_set1_epi8(-128));
					_mm512_storeu_si512((__m512i*)(g_image+(imptr-interleaved-isize)+vidx), modified);

					mdiff0=_mm512_sub_epi16(yuv0, myuv[2*2+0]);
					mdiff1=_mm512_sub_epi16(yuv1, myuv[2*2+1]);
					mdiff0=_mm512_abs_epi16(mdiff0);
					mdiff1=_mm512_abs_epi16(mdiff1);
					mdiff0=_mm512_mullo_epi16(mdiff0, mdiff0);
					mdiff1=_mm512_mullo_epi16(mdiff1, mdiff1);
					mdiff0=_mm512_add_epi16(mdiff0, mdiff1);
					_mm512_store_si512((__m512i*)diff, mdiff0);
					for(int k=0;k<_countof(diff);++k)
						g_sqe[2]+=diff[k];
#endif
				}
				
				ecurr[0]=msyms[0];
				ecurr[1]=msyms[1];
				ecurr[2]=msyms[2];
				ecurr[3]=msyms[3];
				ecurr[4]=msyms[4];
				ecurr[5]=msyms[5];
				msyms[0]=_mm512_sub_epi16(msyms[0], amin);
				msyms[1]=_mm512_sub_epi16(msyms[1], amin);
				msyms[2]=_mm512_sub_epi16(msyms[2], amin);
				msyms[3]=_mm512_sub_epi16(msyms[3], amin);
				msyms[4]=_mm512_sub_epi16(msyms[4], amin);
				msyms[5]=_mm512_sub_epi16(msyms[5], amin);
#ifdef SAVE_RESIDUALS
				{
					ptrdiff_t idx=imptr-interleaved-isize;
					for(int k=0;k<3*NCODERS;++k)
						residuals[idx+k]=(uint8_t)((uint16_t*)msyms)[k];
				}
#endif
				ctx[2]=_mm512_add_epi16(ctx[2], mctxuoffset);
				ctx[3]=_mm512_add_epi16(ctx[3], mctxuoffset);
				ctx[4]=_mm512_add_epi16(ctx[4], mctxvoffset);
				ctx[5]=_mm512_add_epi16(ctx[5], mctxvoffset);
				ctx[0]=_mm512_slli_epi16(ctx[0], 8);
				ctx[1]=_mm512_slli_epi16(ctx[1], 8);
				ctx[2]=_mm512_slli_epi16(ctx[2], 8);
				ctx[3]=_mm512_slli_epi16(ctx[3], 8);
				ctx[4]=_mm512_slli_epi16(ctx[4], 8);
				ctx[5]=_mm512_slli_epi16(ctx[5], 8);
				ctx[0]=_mm512_mask_mov_epi8(ctx[0], 0x5555555555555555, msyms[0]);
				ctx[1]=_mm512_mask_mov_epi8(ctx[1], 0x5555555555555555, msyms[1]);
				ctx[2]=_mm512_mask_mov_epi8(ctx[2], 0x5555555555555555, msyms[2]);
				ctx[3]=_mm512_mask_mov_epi8(ctx[3], 0x5555555555555555, msyms[3]);
				ctx[4]=_mm512_mask_mov_epi8(ctx[4], 0x5555555555555555, msyms[4]);
				ctx[5]=_mm512_mask_mov_epi8(ctx[5], 0x5555555555555555, msyms[5]);
				_mm512_store_si512((__m512i*)syms+0, ctx[0]);
				_mm512_store_si512((__m512i*)syms+1, ctx[1]);
				_mm512_store_si512((__m512i*)syms+2, ctx[2]);
				_mm512_store_si512((__m512i*)syms+3, ctx[3]);
				_mm512_store_si512((__m512i*)syms+4, ctx[4]);
				_mm512_store_si512((__m512i*)syms+5, ctx[5]);
				_mm512_store_si512((__m512i*)ctxptr+0, ctx[0]);//store ctx|residuals
				_mm512_store_si512((__m512i*)ctxptr+1, ctx[1]);
				_mm512_store_si512((__m512i*)ctxptr+2, ctx[2]);
				_mm512_store_si512((__m512i*)ctxptr+3, ctx[3]);
				_mm512_store_si512((__m512i*)ctxptr+4, ctx[4]);
				_mm512_store_si512((__m512i*)ctxptr+5, ctx[5]);
				
				{
					int *pa, *pb, *pc, va, vb, vc;
					pa=hists+syms[0*NCODERS+0x00]; pb=hists+syms[1*NCODERS+0x00]; pc=hists+syms[2*NCODERS+0x00]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x01]; pb=hists+syms[1*NCODERS+0x01]; pc=hists+syms[2*NCODERS+0x01]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x02]; pb=hists+syms[1*NCODERS+0x02]; pc=hists+syms[2*NCODERS+0x02]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x03]; pb=hists+syms[1*NCODERS+0x03]; pc=hists+syms[2*NCODERS+0x03]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x04]; pb=hists+syms[1*NCODERS+0x04]; pc=hists+syms[2*NCODERS+0x04]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x05]; pb=hists+syms[1*NCODERS+0x05]; pc=hists+syms[2*NCODERS+0x05]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x06]; pb=hists+syms[1*NCODERS+0x06]; pc=hists+syms[2*NCODERS+0x06]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x07]; pb=hists+syms[1*NCODERS+0x07]; pc=hists+syms[2*NCODERS+0x07]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x08]; pb=hists+syms[1*NCODERS+0x08]; pc=hists+syms[2*NCODERS+0x08]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x09]; pb=hists+syms[1*NCODERS+0x09]; pc=hists+syms[2*NCODERS+0x09]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x0A]; pb=hists+syms[1*NCODERS+0x0A]; pc=hists+syms[2*NCODERS+0x0A]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x0B]; pb=hists+syms[1*NCODERS+0x0B]; pc=hists+syms[2*NCODERS+0x0B]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x0C]; pb=hists+syms[1*NCODERS+0x0C]; pc=hists+syms[2*NCODERS+0x0C]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x0D]; pb=hists+syms[1*NCODERS+0x0D]; pc=hists+syms[2*NCODERS+0x0D]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x0E]; pb=hists+syms[1*NCODERS+0x0E]; pc=hists+syms[2*NCODERS+0x0E]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x0F]; pb=hists+syms[1*NCODERS+0x0F]; pc=hists+syms[2*NCODERS+0x0F]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x10]; pb=hists+syms[1*NCODERS+0x10]; pc=hists+syms[2*NCODERS+0x10]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x11]; pb=hists+syms[1*NCODERS+0x11]; pc=hists+syms[2*NCODERS+0x11]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x12]; pb=hists+syms[1*NCODERS+0x12]; pc=hists+syms[2*NCODERS+0x12]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x13]; pb=hists+syms[1*NCODERS+0x13]; pc=hists+syms[2*NCODERS+0x13]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x14]; pb=hists+syms[1*NCODERS+0x14]; pc=hists+syms[2*NCODERS+0x14]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x15]; pb=hists+syms[1*NCODERS+0x15]; pc=hists+syms[2*NCODERS+0x15]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x16]; pb=hists+syms[1*NCODERS+0x16]; pc=hists+syms[2*NCODERS+0x16]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x17]; pb=hists+syms[1*NCODERS+0x17]; pc=hists+syms[2*NCODERS+0x17]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x18]; pb=hists+syms[1*NCODERS+0x18]; pc=hists+syms[2*NCODERS+0x18]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x19]; pb=hists+syms[1*NCODERS+0x19]; pc=hists+syms[2*NCODERS+0x19]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x1A]; pb=hists+syms[1*NCODERS+0x1A]; pc=hists+syms[2*NCODERS+0x1A]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x1B]; pb=hists+syms[1*NCODERS+0x1B]; pc=hists+syms[2*NCODERS+0x1B]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x1C]; pb=hists+syms[1*NCODERS+0x1C]; pc=hists+syms[2*NCODERS+0x1C]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x1D]; pb=hists+syms[1*NCODERS+0x1D]; pc=hists+syms[2*NCODERS+0x1D]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x1E]; pb=hists+syms[1*NCODERS+0x1E]; pc=hists+syms[2*NCODERS+0x1E]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x1F]; pb=hists+syms[1*NCODERS+0x1F]; pc=hists+syms[2*NCODERS+0x1F]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x20]; pb=hists+syms[1*NCODERS+0x20]; pc=hists+syms[2*NCODERS+0x20]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x21]; pb=hists+syms[1*NCODERS+0x21]; pc=hists+syms[2*NCODERS+0x21]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x22]; pb=hists+syms[1*NCODERS+0x22]; pc=hists+syms[2*NCODERS+0x22]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x23]; pb=hists+syms[1*NCODERS+0x23]; pc=hists+syms[2*NCODERS+0x23]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x24]; pb=hists+syms[1*NCODERS+0x24]; pc=hists+syms[2*NCODERS+0x24]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x25]; pb=hists+syms[1*NCODERS+0x25]; pc=hists+syms[2*NCODERS+0x25]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x26]; pb=hists+syms[1*NCODERS+0x26]; pc=hists+syms[2*NCODERS+0x26]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x27]; pb=hists+syms[1*NCODERS+0x27]; pc=hists+syms[2*NCODERS+0x27]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x28]; pb=hists+syms[1*NCODERS+0x28]; pc=hists+syms[2*NCODERS+0x28]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x29]; pb=hists+syms[1*NCODERS+0x29]; pc=hists+syms[2*NCODERS+0x29]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x2A]; pb=hists+syms[1*NCODERS+0x2A]; pc=hists+syms[2*NCODERS+0x2A]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x2B]; pb=hists+syms[1*NCODERS+0x2B]; pc=hists+syms[2*NCODERS+0x2B]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x2C]; pb=hists+syms[1*NCODERS+0x2C]; pc=hists+syms[2*NCODERS+0x2C]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x2D]; pb=hists+syms[1*NCODERS+0x2D]; pc=hists+syms[2*NCODERS+0x2D]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x2E]; pb=hists+syms[1*NCODERS+0x2E]; pc=hists+syms[2*NCODERS+0x2E]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x2F]; pb=hists+syms[1*NCODERS+0x2F]; pc=hists+syms[2*NCODERS+0x2F]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x30]; pb=hists+syms[1*NCODERS+0x30]; pc=hists+syms[2*NCODERS+0x30]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x31]; pb=hists+syms[1*NCODERS+0x31]; pc=hists+syms[2*NCODERS+0x31]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x32]; pb=hists+syms[1*NCODERS+0x32]; pc=hists+syms[2*NCODERS+0x32]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x33]; pb=hists+syms[1*NCODERS+0x33]; pc=hists+syms[2*NCODERS+0x33]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x34]; pb=hists+syms[1*NCODERS+0x34]; pc=hists+syms[2*NCODERS+0x34]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x35]; pb=hists+syms[1*NCODERS+0x35]; pc=hists+syms[2*NCODERS+0x35]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x36]; pb=hists+syms[1*NCODERS+0x36]; pc=hists+syms[2*NCODERS+0x36]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x37]; pb=hists+syms[1*NCODERS+0x37]; pc=hists+syms[2*NCODERS+0x37]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x38]; pb=hists+syms[1*NCODERS+0x38]; pc=hists+syms[2*NCODERS+0x38]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x39]; pb=hists+syms[1*NCODERS+0x39]; pc=hists+syms[2*NCODERS+0x39]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x3A]; pb=hists+syms[1*NCODERS+0x3A]; pc=hists+syms[2*NCODERS+0x3A]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x3B]; pb=hists+syms[1*NCODERS+0x3B]; pc=hists+syms[2*NCODERS+0x3B]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x3C]; pb=hists+syms[1*NCODERS+0x3C]; pc=hists+syms[2*NCODERS+0x3C]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x3D]; pb=hists+syms[1*NCODERS+0x3D]; pc=hists+syms[2*NCODERS+0x3D]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x3E]; pb=hists+syms[1*NCODERS+0x3E]; pc=hists+syms[2*NCODERS+0x3E]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
					pa=hists+syms[0*NCODERS+0x3F]; pb=hists+syms[1*NCODERS+0x3F]; pc=hists+syms[2*NCODERS+0x3F]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
				}
				ctxptr+=sizeof(int16_t[3][NCODERS]);
			}
			else
			{
				//decode
				__m256i msyms8[6], tmp;
				
				//yuv = (char)(sym+pred-128)	= (uint8_t)(sym+pred)-128
				dec_yuv(mstate, ctx+0, CDF2syms+((ptrdiff_t)NCTX*0<<PROBBITS), &streamptr, streamend, myuv+0);//residuals from [0 ~ 255]
				dec_yuv(mstate, ctx+2, CDF2syms+((ptrdiff_t)NCTX*1<<PROBBITS), &streamptr, streamend, myuv+2);
				dec_yuv(mstate, ctx+4, CDF2syms+((ptrdiff_t)NCTX*2<<PROBBITS), &streamptr, streamend, myuv+4);
				tmp=_mm256_set1_epi8(-128);
				
				//reconstruct Y
				if(dist>1)
				{
					ecurr[0]=_mm512_add_epi16(myuv[0], amin);
					ecurr[1]=_mm512_add_epi16(myuv[1], amin);
					myuv[0]=_mm512_mullo_epi16(ecurr[0], mdist);
					myuv[1]=_mm512_mullo_epi16(ecurr[1], mdist);
					myuv[0]=_mm512_add_epi16(myuv[0], pred[0]);
					myuv[1]=_mm512_add_epi16(myuv[1], pred[1]);
					myuv[0]=_mm512_max_epi16(myuv[0], amin);
					myuv[1]=_mm512_max_epi16(myuv[1], amin);
					myuv[0]=_mm512_min_epi16(myuv[0], amax);
					myuv[1]=_mm512_min_epi16(myuv[1], amax);
				}
				else
				{
					myuv[0]=_mm512_add_epi16(myuv[0], pred[0]);
					myuv[1]=_mm512_add_epi16(myuv[1], pred[1]);
					myuv[0]=_mm512_maskz_mov_epi8(0x5555555555555555, myuv[0]);
					myuv[1]=_mm512_maskz_mov_epi8(0x5555555555555555, myuv[1]);
					myuv[0]=_mm512_add_epi16(myuv[0], amin);
					myuv[1]=_mm512_add_epi16(myuv[1], amin);
					ecurr[0]=_mm512_sub_epi16(myuv[0], pred[0]);//sub pred
					ecurr[1]=_mm512_sub_epi16(myuv[1], pred[1]);
				}
				msyms8[0]=_mm512_cvtepi16_epi8(myuv[0]);
				msyms8[1]=_mm512_cvtepi16_epi8(myuv[1]);
				msyms8[0]=_mm256_add_epi8(msyms8[0], tmp);
				msyms8[1]=_mm256_add_epi8(msyms8[1], tmp);
				_mm256_store_si256((__m256i*)(imptr+yidx)+0, msyms8[0]);//store bytes
				_mm256_store_si256((__m256i*)(imptr+yidx)+1, msyms8[1]);
#ifdef ENABLE_GUIDE
				if(memcmp(imptr+yidx, g_image+(imptr-interleaved)+yidx, NCODERS))
				{
					printf("original  decoded  original-decoded  XYC0 %d %d %d\n", kx, ky, yidx);
					for(int k=0;k<NCODERS;++k)
						printf("L%2d  0x%02X  0x%02X  %4d\n"
							, k
							, g_image[imptr-interleaved+yidx+k]
							, imptr[yidx+k]
							, g_image[imptr-interleaved+yidx+k]-imptr[yidx+k]
						);
					CRASH("guide error XYC0 %d %d %d/%d", kx, ky, yidx, NCODERS);
				}
#endif

				//reconstruct U
				moffset[0]=_mm512_and_si512(myuv[0], uhelpmask);
				moffset[1]=_mm512_and_si512(myuv[1], uhelpmask);
				pred[2]=_mm512_add_epi16(pred[2], moffset[0]);
				pred[3]=_mm512_add_epi16(pred[3], moffset[1]);
				pred[2]=_mm512_max_epi16(pred[2], amin); pred[2]=_mm512_min_epi16(pred[2], amax);
				pred[3]=_mm512_max_epi16(pred[3], amin); pred[3]=_mm512_min_epi16(pred[3], amax);
				if(dist>1)
				{
					ecurr[2]=_mm512_add_epi16(myuv[2], amin);
					ecurr[3]=_mm512_add_epi16(myuv[3], amin);
					myuv[2]=_mm512_mullo_epi16(ecurr[2], mdist);
					myuv[3]=_mm512_mullo_epi16(ecurr[3], mdist);
					myuv[2]=_mm512_add_epi16(myuv[2], pred[2]);
					myuv[3]=_mm512_add_epi16(myuv[3], pred[3]);
					myuv[2]=_mm512_max_epi16(myuv[2], amin);
					myuv[3]=_mm512_max_epi16(myuv[3], amin);
					myuv[2]=_mm512_min_epi16(myuv[2], amax);
					myuv[3]=_mm512_min_epi16(myuv[3], amax);
				}
				else
				{
					myuv[2]=_mm512_add_epi16(myuv[2], pred[2]);
					myuv[3]=_mm512_add_epi16(myuv[3], pred[3]);
					myuv[2]=_mm512_maskz_mov_epi8(0x5555555555555555, myuv[2]);
					myuv[3]=_mm512_maskz_mov_epi8(0x5555555555555555, myuv[3]);
					myuv[2]=_mm512_add_epi16(myuv[2], amin);
					myuv[3]=_mm512_add_epi16(myuv[3], amin);
					ecurr[2]=_mm512_sub_epi16(myuv[2], pred[2]);//sub pred
					ecurr[3]=_mm512_sub_epi16(myuv[3], pred[3]);
				}
				msyms8[2]=_mm512_cvtepi16_epi8(myuv[2]);
				msyms8[3]=_mm512_cvtepi16_epi8(myuv[3]);
				msyms8[2]=_mm256_add_epi8(msyms8[2], tmp);
				msyms8[3]=_mm256_add_epi8(msyms8[3], tmp);
				_mm256_store_si256((__m256i*)(imptr+uidx)+0, msyms8[2]);
				_mm256_store_si256((__m256i*)(imptr+uidx)+1, msyms8[3]);
#ifdef ENABLE_GUIDE
				if(memcmp(imptr+uidx, g_image+(imptr-interleaved)+uidx, NCODERS))
				{
					printf("original  decoded  original-decoded  XYC1 %d %d %d\n", kx, ky, uidx);
					for(int k=0;k<NCODERS;++k)
						printf("L%2d  0x%02X  0x%02X  %4d\n"
							, k
							, g_image[imptr-interleaved+uidx+k]
							, imptr[uidx+k]
							, g_image[imptr-interleaved+uidx+k]-imptr[uidx+k]
						);
					CRASH("guide error XYC1 %d %d %d/%d", kx, ky, uidx, NCODERS);
				}
#endif
				
				//reconstruct V
				moffset[2]=_mm512_mullo_epi16(vc0, myuv[0]);
				moffset[3]=_mm512_mullo_epi16(vc0, myuv[1]);
				moffset[2]=_mm512_add_epi16(moffset[2], _mm512_mullo_epi16(vc1, myuv[2]));
				moffset[3]=_mm512_add_epi16(moffset[3], _mm512_mullo_epi16(vc1, myuv[3]));
				moffset[2]=_mm512_srai_epi16(moffset[2], 2);
				moffset[3]=_mm512_srai_epi16(moffset[3], 2);
				pred[4]=_mm512_add_epi16(pred[4], moffset[2]);
				pred[5]=_mm512_add_epi16(pred[5], moffset[3]);
				pred[4]=_mm512_max_epi16(pred[4], amin); pred[4]=_mm512_min_epi16(pred[4], amax);
				pred[5]=_mm512_max_epi16(pred[5], amin); pred[5]=_mm512_min_epi16(pred[5], amax);
				if(dist>1)
				{
					ecurr[4]=_mm512_add_epi16(myuv[4], amin);
					ecurr[5]=_mm512_add_epi16(myuv[5], amin);
					myuv[4]=_mm512_mullo_epi16(ecurr[4], mdist);
					myuv[5]=_mm512_mullo_epi16(ecurr[5], mdist);
					myuv[4]=_mm512_add_epi16(myuv[4], pred[4]);
					myuv[5]=_mm512_add_epi16(myuv[5], pred[5]);
					myuv[4]=_mm512_max_epi16(myuv[4], amin);
					myuv[5]=_mm512_max_epi16(myuv[5], amin);
					myuv[4]=_mm512_min_epi16(myuv[4], amax);
					myuv[5]=_mm512_min_epi16(myuv[5], amax);
				}
				else
				{
					myuv[4]=_mm512_add_epi16(myuv[4], pred[4]);
					myuv[5]=_mm512_add_epi16(myuv[5], pred[5]);
					myuv[4]=_mm512_maskz_mov_epi8(0x5555555555555555, myuv[4]);
					myuv[5]=_mm512_maskz_mov_epi8(0x5555555555555555, myuv[5]);
					myuv[4]=_mm512_add_epi16(myuv[4], amin);
					myuv[5]=_mm512_add_epi16(myuv[5], amin);
					ecurr[4]=_mm512_sub_epi16(myuv[4], pred[4]);//sub pred
					ecurr[5]=_mm512_sub_epi16(myuv[5], pred[5]);
				}
				msyms8[4]=_mm512_cvtepi16_epi8(myuv[4]);
				msyms8[5]=_mm512_cvtepi16_epi8(myuv[5]);
				msyms8[4]=_mm256_add_epi8(msyms8[4], tmp);
				msyms8[5]=_mm256_add_epi8(msyms8[5], tmp);
				_mm256_store_si256((__m256i*)(imptr+vidx)+0, msyms8[4]);
				_mm256_store_si256((__m256i*)(imptr+vidx)+1, msyms8[5]);
#ifdef ENABLE_GUIDE
				if(memcmp(imptr+vidx, g_image+(imptr-interleaved)+vidx, NCODERS))
				{
					printf("original  decoded  original-decoded  XYC2 %d %d %d\n", kx, ky, vidx);
					for(int k=0;k<NCODERS;++k)
						printf("L%2d  0x%02X  0x%02X  %4d\n"
							, k
							, g_image[imptr-interleaved+vidx+k]
							, imptr[vidx+k]
							, g_image[imptr-interleaved+vidx+k]-imptr[vidx+k]
						);
					CRASH("guide error XYC2 %d %d %d/%d", kx, ky, vidx, NCODERS);
				}
#endif
			}
			W[0]=myuv[0];
			W[1]=myuv[1];
			W[2]=_mm512_sub_epi16(myuv[2], moffset[0]);
			W[3]=_mm512_sub_epi16(myuv[3], moffset[1]);
			W[4]=_mm512_sub_epi16(myuv[4], moffset[2]);
			W[5]=_mm512_sub_epi16(myuv[5], moffset[3]);
			_mm512_store_si512((__m512i*)rows[0]+0+(0+0*NROWS*NVAL)*NREG*NCH, W[0]);//store neighbors
			_mm512_store_si512((__m512i*)rows[0]+1+(0+0*NROWS*NVAL)*NREG*NCH, W[1]);
			_mm512_store_si512((__m512i*)rows[0]+2+(0+0*NROWS*NVAL)*NREG*NCH, W[2]);
			_mm512_store_si512((__m512i*)rows[0]+3+(0+0*NROWS*NVAL)*NREG*NCH, W[3]);
			_mm512_store_si512((__m512i*)rows[0]+4+(0+0*NROWS*NVAL)*NREG*NCH, W[4]);
			_mm512_store_si512((__m512i*)rows[0]+5+(0+0*NROWS*NVAL)*NREG*NCH, W[5]);
			ecurr[0]=_mm512_xor_si512(_mm512_slli_epi16(ecurr[0], 1), _mm512_srai_epi16(ecurr[0], 15));
			ecurr[1]=_mm512_xor_si512(_mm512_slli_epi16(ecurr[1], 1), _mm512_srai_epi16(ecurr[1], 15));
			ecurr[2]=_mm512_xor_si512(_mm512_slli_epi16(ecurr[2], 1), _mm512_srai_epi16(ecurr[2], 15));
			ecurr[3]=_mm512_xor_si512(_mm512_slli_epi16(ecurr[3], 1), _mm512_srai_epi16(ecurr[3], 15));
			ecurr[4]=_mm512_xor_si512(_mm512_slli_epi16(ecurr[4], 1), _mm512_srai_epi16(ecurr[4], 15));
			ecurr[5]=_mm512_xor_si512(_mm512_slli_epi16(ecurr[5], 1), _mm512_srai_epi16(ecurr[5], 15));
			if(effort==1)//update
			{
				__m512i mu[6];

				mu[0]=_mm512_sub_epi16(W[0], predYUV0[0]);
				mu[1]=_mm512_sub_epi16(W[1], predYUV0[1]);
				mu[2]=_mm512_sub_epi16(W[2], predYUV0[2]);
				mu[3]=_mm512_sub_epi16(W[3], predYUV0[3]);
				mu[4]=_mm512_sub_epi16(W[4], predYUV0[4]);
				mu[5]=_mm512_sub_epi16(W[5], predYUV0[5]);
				for(int k=0;k<L1_NPREDS1;++k)//update
				{
					__m512i mc[12];
#ifdef USE_L2
					mc[0x0]=_mm512_mullo_epi16(L1preds[k*6+0], mu[0]);//L2
					mc[0x1]=_mm512_mullo_epi16(L1preds[k*6+1], mu[1]);
					mc[0x2]=_mm512_mullo_epi16(L1preds[k*6+2], mu[2]);
					mc[0x3]=_mm512_mullo_epi16(L1preds[k*6+3], mu[3]);
					mc[0x4]=_mm512_mullo_epi16(L1preds[k*6+4], mu[4]);
					mc[0x5]=_mm512_mullo_epi16(L1preds[k*6+5], mu[5]);
#else
					__mmask32 mask0=_mm512_movepi16_mask(mu[0]);//L1
					__mmask32 mask1=_mm512_movepi16_mask(mu[1]);
					__mmask32 mask2=_mm512_movepi16_mask(mu[2]);
					__mmask32 mask3=_mm512_movepi16_mask(mu[3]);
					__mmask32 mask4=_mm512_movepi16_mask(mu[4]);
					__mmask32 mask5=_mm512_movepi16_mask(mu[5]);
					__m512i zero=_mm512_setzero_si512();
					mc[0x0]=_mm512_mask_sub_epi16(L1preds[k*6+0], mask0, zero, L1preds[k*6+0]);
					mc[0x1]=_mm512_mask_sub_epi16(L1preds[k*6+1], mask1, zero, L1preds[k*6+1]);
					mc[0x2]=_mm512_mask_sub_epi16(L1preds[k*6+2], mask2, zero, L1preds[k*6+2]);
					mc[0x3]=_mm512_mask_sub_epi16(L1preds[k*6+3], mask3, zero, L1preds[k*6+3]);
					mc[0x4]=_mm512_mask_sub_epi16(L1preds[k*6+4], mask4, zero, L1preds[k*6+4]);
					mc[0x5]=_mm512_mask_sub_epi16(L1preds[k*6+5], mask5, zero, L1preds[k*6+5]);
#endif
					//16 -> 32	3 lo 3 hi registers
					mc[0x6]=_mm512_srai_epi32(mc[0x0], 16);
					mc[0x7]=_mm512_srai_epi32(mc[0x1], 16);
					mc[0x8]=_mm512_srai_epi32(mc[0x2], 16);
					mc[0x9]=_mm512_srai_epi32(mc[0x3], 16);
					mc[0xA]=_mm512_srai_epi32(mc[0x4], 16);
					mc[0xB]=_mm512_srai_epi32(mc[0x5], 16);

					mc[0x0]=_mm512_slli_epi32(mc[0x0], 16);
					mc[0x1]=_mm512_slli_epi32(mc[0x1], 16);
					mc[0x2]=_mm512_slli_epi32(mc[0x2], 16);
					mc[0x3]=_mm512_slli_epi32(mc[0x3], 16);
					mc[0x4]=_mm512_slli_epi32(mc[0x4], 16);
					mc[0x5]=_mm512_slli_epi32(mc[0x5], 16);

					mc[0x0]=_mm512_srai_epi32(mc[0x0], 16);
					mc[0x1]=_mm512_srai_epi32(mc[0x1], 16);
					mc[0x2]=_mm512_srai_epi32(mc[0x2], 16);
					mc[0x3]=_mm512_srai_epi32(mc[0x3], 16);
					mc[0x4]=_mm512_srai_epi32(mc[0x4], 16);
					mc[0x5]=_mm512_srai_epi32(mc[0x5], 16);
					mc[0x0]=_mm512_add_epi32(mc[0x0], _mm512_load_si512((__m512i*)L1weights+k*12+0x0));//update coeffs
					mc[0x1]=_mm512_add_epi32(mc[0x1], _mm512_load_si512((__m512i*)L1weights+k*12+0x1));
					mc[0x2]=_mm512_add_epi32(mc[0x2], _mm512_load_si512((__m512i*)L1weights+k*12+0x2));
					mc[0x3]=_mm512_add_epi32(mc[0x3], _mm512_load_si512((__m512i*)L1weights+k*12+0x3));
					mc[0x4]=_mm512_add_epi32(mc[0x4], _mm512_load_si512((__m512i*)L1weights+k*12+0x4));
					mc[0x5]=_mm512_add_epi32(mc[0x5], _mm512_load_si512((__m512i*)L1weights+k*12+0x5));
					mc[0x6]=_mm512_add_epi32(mc[0x6], _mm512_load_si512((__m512i*)L1weights+k*12+0x6));
					mc[0x7]=_mm512_add_epi32(mc[0x7], _mm512_load_si512((__m512i*)L1weights+k*12+0x7));
					mc[0x8]=_mm512_add_epi32(mc[0x8], _mm512_load_si512((__m512i*)L1weights+k*12+0x8));
					mc[0x9]=_mm512_add_epi32(mc[0x9], _mm512_load_si512((__m512i*)L1weights+k*12+0x9));
					mc[0xA]=_mm512_add_epi32(mc[0xA], _mm512_load_si512((__m512i*)L1weights+k*12+0xA));
					mc[0xB]=_mm512_add_epi32(mc[0xB], _mm512_load_si512((__m512i*)L1weights+k*12+0xB));
					_mm512_store_si512((__m512i*)L1weights+k*12+0x0, mc[0x0]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0x1, mc[0x1]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0x2, mc[0x2]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0x3, mc[0x3]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0x4, mc[0x4]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0x5, mc[0x5]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0x6, mc[0x6]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0x7, mc[0x7]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0x8, mc[0x8]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0x9, mc[0x9]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0xA, mc[0xA]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0xB, mc[0xB]);
				}
			}
			else if(effort==2)//update
			{
				__m512i mu[6];

				mu[0]=_mm512_sub_epi16(W[0], predYUV0[0]);
				mu[1]=_mm512_sub_epi16(W[1], predYUV0[1]);
				mu[2]=_mm512_sub_epi16(W[2], predYUV0[2]);
				mu[3]=_mm512_sub_epi16(W[3], predYUV0[3]);
				mu[4]=_mm512_sub_epi16(W[4], predYUV0[4]);
				mu[5]=_mm512_sub_epi16(W[5], predYUV0[5]);
				for(int k=0;k<L1_NPREDS2;++k)//update
				{
					__m512i mc[12];
#ifdef USE_L2
					mc[0x0]=_mm512_mullo_epi16(L1preds[k*6+0], mu[0]);//L2
					mc[0x1]=_mm512_mullo_epi16(L1preds[k*6+1], mu[1]);
					mc[0x2]=_mm512_mullo_epi16(L1preds[k*6+2], mu[2]);
					mc[0x3]=_mm512_mullo_epi16(L1preds[k*6+3], mu[3]);
					mc[0x4]=_mm512_mullo_epi16(L1preds[k*6+4], mu[4]);
					mc[0x5]=_mm512_mullo_epi16(L1preds[k*6+5], mu[5]);
#else
					__mmask32 mask0=_mm512_movepi16_mask(mu[0]);//L1
					__mmask32 mask1=_mm512_movepi16_mask(mu[1]);
					__mmask32 mask2=_mm512_movepi16_mask(mu[2]);
					__mmask32 mask3=_mm512_movepi16_mask(mu[3]);
					__mmask32 mask4=_mm512_movepi16_mask(mu[4]);
					__mmask32 mask5=_mm512_movepi16_mask(mu[5]);
					__m512i zero=_mm512_setzero_si512();
					mc[0x0]=_mm512_mask_sub_epi16(L1preds[k*6+0], mask0, zero, L1preds[k*6+0]);
					mc[0x1]=_mm512_mask_sub_epi16(L1preds[k*6+1], mask1, zero, L1preds[k*6+1]);
					mc[0x2]=_mm512_mask_sub_epi16(L1preds[k*6+2], mask2, zero, L1preds[k*6+2]);
					mc[0x3]=_mm512_mask_sub_epi16(L1preds[k*6+3], mask3, zero, L1preds[k*6+3]);
					mc[0x4]=_mm512_mask_sub_epi16(L1preds[k*6+4], mask4, zero, L1preds[k*6+4]);
					mc[0x5]=_mm512_mask_sub_epi16(L1preds[k*6+5], mask5, zero, L1preds[k*6+5]);
#endif
					//16 -> 32	3 lo 3 hi registers
					mc[0x6]=_mm512_srai_epi32(mc[0x0], 16);
					mc[0x7]=_mm512_srai_epi32(mc[0x1], 16);
					mc[0x8]=_mm512_srai_epi32(mc[0x2], 16);
					mc[0x9]=_mm512_srai_epi32(mc[0x3], 16);
					mc[0xA]=_mm512_srai_epi32(mc[0x4], 16);
					mc[0xB]=_mm512_srai_epi32(mc[0x5], 16);

					mc[0x0]=_mm512_slli_epi32(mc[0x0], 16);
					mc[0x1]=_mm512_slli_epi32(mc[0x1], 16);
					mc[0x2]=_mm512_slli_epi32(mc[0x2], 16);
					mc[0x3]=_mm512_slli_epi32(mc[0x3], 16);
					mc[0x4]=_mm512_slli_epi32(mc[0x4], 16);
					mc[0x5]=_mm512_slli_epi32(mc[0x5], 16);

					mc[0x0]=_mm512_srai_epi32(mc[0x0], 16);
					mc[0x1]=_mm512_srai_epi32(mc[0x1], 16);
					mc[0x2]=_mm512_srai_epi32(mc[0x2], 16);
					mc[0x3]=_mm512_srai_epi32(mc[0x3], 16);
					mc[0x4]=_mm512_srai_epi32(mc[0x4], 16);
					mc[0x5]=_mm512_srai_epi32(mc[0x5], 16);
					mc[0x0]=_mm512_add_epi32(mc[0x0], _mm512_load_si512((__m512i*)L1weights+k*12+0x0));//update coeffs
					mc[0x1]=_mm512_add_epi32(mc[0x1], _mm512_load_si512((__m512i*)L1weights+k*12+0x1));
					mc[0x2]=_mm512_add_epi32(mc[0x2], _mm512_load_si512((__m512i*)L1weights+k*12+0x2));
					mc[0x3]=_mm512_add_epi32(mc[0x3], _mm512_load_si512((__m512i*)L1weights+k*12+0x3));
					mc[0x4]=_mm512_add_epi32(mc[0x4], _mm512_load_si512((__m512i*)L1weights+k*12+0x4));
					mc[0x5]=_mm512_add_epi32(mc[0x5], _mm512_load_si512((__m512i*)L1weights+k*12+0x5));
					mc[0x6]=_mm512_add_epi32(mc[0x6], _mm512_load_si512((__m512i*)L1weights+k*12+0x6));
					mc[0x7]=_mm512_add_epi32(mc[0x7], _mm512_load_si512((__m512i*)L1weights+k*12+0x7));
					mc[0x8]=_mm512_add_epi32(mc[0x8], _mm512_load_si512((__m512i*)L1weights+k*12+0x8));
					mc[0x9]=_mm512_add_epi32(mc[0x9], _mm512_load_si512((__m512i*)L1weights+k*12+0x9));
					mc[0xA]=_mm512_add_epi32(mc[0xA], _mm512_load_si512((__m512i*)L1weights+k*12+0xA));
					mc[0xB]=_mm512_add_epi32(mc[0xB], _mm512_load_si512((__m512i*)L1weights+k*12+0xB));
					_mm512_store_si512((__m512i*)L1weights+k*12+0x0, mc[0x0]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0x1, mc[0x1]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0x2, mc[0x2]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0x3, mc[0x3]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0x4, mc[0x4]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0x5, mc[0x5]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0x6, mc[0x6]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0x7, mc[0x7]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0x8, mc[0x8]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0x9, mc[0x9]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0xA, mc[0xA]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0xB, mc[0xB]);
				}
			}
			else if(effort==3)//update
			{
				__m512i mu[6];

				mu[0]=_mm512_sub_epi16(W[0], predYUV0[0]);
				mu[1]=_mm512_sub_epi16(W[1], predYUV0[1]);
				mu[2]=_mm512_sub_epi16(W[2], predYUV0[2]);
				mu[3]=_mm512_sub_epi16(W[3], predYUV0[3]);
				mu[4]=_mm512_sub_epi16(W[4], predYUV0[4]);
				mu[5]=_mm512_sub_epi16(W[5], predYUV0[5]);
				for(int k=0;k<L1_NPREDS3;++k)//update
				{
					__m512i mc[12];
#ifdef USE_L2
					mc[0x0]=_mm512_mullo_epi16(L1preds[k*6+0], mu[0]);//L2
					mc[0x1]=_mm512_mullo_epi16(L1preds[k*6+1], mu[1]);
					mc[0x2]=_mm512_mullo_epi16(L1preds[k*6+2], mu[2]);
					mc[0x3]=_mm512_mullo_epi16(L1preds[k*6+3], mu[3]);
					mc[0x4]=_mm512_mullo_epi16(L1preds[k*6+4], mu[4]);
					mc[0x5]=_mm512_mullo_epi16(L1preds[k*6+5], mu[5]);
#else
					__mmask32 mask0=_mm512_movepi16_mask(mu[0]);//L1
					__mmask32 mask1=_mm512_movepi16_mask(mu[1]);
					__mmask32 mask2=_mm512_movepi16_mask(mu[2]);
					__mmask32 mask3=_mm512_movepi16_mask(mu[3]);
					__mmask32 mask4=_mm512_movepi16_mask(mu[4]);
					__mmask32 mask5=_mm512_movepi16_mask(mu[5]);
					__m512i zero=_mm512_setzero_si512();
					mc[0x0]=_mm512_mask_sub_epi16(L1preds[k*6+0], mask0, zero, L1preds[k*6+0]);
					mc[0x1]=_mm512_mask_sub_epi16(L1preds[k*6+1], mask1, zero, L1preds[k*6+1]);
					mc[0x2]=_mm512_mask_sub_epi16(L1preds[k*6+2], mask2, zero, L1preds[k*6+2]);
					mc[0x3]=_mm512_mask_sub_epi16(L1preds[k*6+3], mask3, zero, L1preds[k*6+3]);
					mc[0x4]=_mm512_mask_sub_epi16(L1preds[k*6+4], mask4, zero, L1preds[k*6+4]);
					mc[0x5]=_mm512_mask_sub_epi16(L1preds[k*6+5], mask5, zero, L1preds[k*6+5]);
#endif
					//16 -> 32	3 lo 3 hi registers
					mc[0x6]=_mm512_srai_epi32(mc[0x0], 16);
					mc[0x7]=_mm512_srai_epi32(mc[0x1], 16);
					mc[0x8]=_mm512_srai_epi32(mc[0x2], 16);
					mc[0x9]=_mm512_srai_epi32(mc[0x3], 16);
					mc[0xA]=_mm512_srai_epi32(mc[0x4], 16);
					mc[0xB]=_mm512_srai_epi32(mc[0x5], 16);

					mc[0x0]=_mm512_slli_epi32(mc[0x0], 16);
					mc[0x1]=_mm512_slli_epi32(mc[0x1], 16);
					mc[0x2]=_mm512_slli_epi32(mc[0x2], 16);
					mc[0x3]=_mm512_slli_epi32(mc[0x3], 16);
					mc[0x4]=_mm512_slli_epi32(mc[0x4], 16);
					mc[0x5]=_mm512_slli_epi32(mc[0x5], 16);

					mc[0x0]=_mm512_srai_epi32(mc[0x0], 16);
					mc[0x1]=_mm512_srai_epi32(mc[0x1], 16);
					mc[0x2]=_mm512_srai_epi32(mc[0x2], 16);
					mc[0x3]=_mm512_srai_epi32(mc[0x3], 16);
					mc[0x4]=_mm512_srai_epi32(mc[0x4], 16);
					mc[0x5]=_mm512_srai_epi32(mc[0x5], 16);
					mc[0x0]=_mm512_add_epi32(mc[0x0], _mm512_load_si512((__m512i*)L1weights+k*12+0x0));//update coeffs
					mc[0x1]=_mm512_add_epi32(mc[0x1], _mm512_load_si512((__m512i*)L1weights+k*12+0x1));
					mc[0x2]=_mm512_add_epi32(mc[0x2], _mm512_load_si512((__m512i*)L1weights+k*12+0x2));
					mc[0x3]=_mm512_add_epi32(mc[0x3], _mm512_load_si512((__m512i*)L1weights+k*12+0x3));
					mc[0x4]=_mm512_add_epi32(mc[0x4], _mm512_load_si512((__m512i*)L1weights+k*12+0x4));
					mc[0x5]=_mm512_add_epi32(mc[0x5], _mm512_load_si512((__m512i*)L1weights+k*12+0x5));
					mc[0x6]=_mm512_add_epi32(mc[0x6], _mm512_load_si512((__m512i*)L1weights+k*12+0x6));
					mc[0x7]=_mm512_add_epi32(mc[0x7], _mm512_load_si512((__m512i*)L1weights+k*12+0x7));
					mc[0x8]=_mm512_add_epi32(mc[0x8], _mm512_load_si512((__m512i*)L1weights+k*12+0x8));
					mc[0x9]=_mm512_add_epi32(mc[0x9], _mm512_load_si512((__m512i*)L1weights+k*12+0x9));
					mc[0xA]=_mm512_add_epi32(mc[0xA], _mm512_load_si512((__m512i*)L1weights+k*12+0xA));
					mc[0xB]=_mm512_add_epi32(mc[0xB], _mm512_load_si512((__m512i*)L1weights+k*12+0xB));
					_mm512_store_si512((__m512i*)L1weights+k*12+0x0, mc[0x0]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0x1, mc[0x1]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0x2, mc[0x2]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0x3, mc[0x3]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0x4, mc[0x4]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0x5, mc[0x5]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0x6, mc[0x6]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0x7, mc[0x7]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0x8, mc[0x8]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0x9, mc[0x9]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0xA, mc[0xA]);
					_mm512_store_si512((__m512i*)L1weights+k*12+0xB, mc[0xB]);
				}
			}
			//context update = (2*eW+(e<<3)+max(eNEE, eNEEE))>>2
			eNEEE[0]=_mm512_load_si512((__m512i*)rows[1]+0+(1+3*NROWS*NVAL)*NREG*NCH);
			eNEEE[1]=_mm512_load_si512((__m512i*)rows[1]+1+(1+3*NROWS*NVAL)*NREG*NCH);
			eNEEE[2]=_mm512_load_si512((__m512i*)rows[1]+2+(1+3*NROWS*NVAL)*NREG*NCH);
			eNEEE[3]=_mm512_load_si512((__m512i*)rows[1]+3+(1+3*NROWS*NVAL)*NREG*NCH);
			eNEEE[4]=_mm512_load_si512((__m512i*)rows[1]+4+(1+3*NROWS*NVAL)*NREG*NCH);
			eNEEE[5]=_mm512_load_si512((__m512i*)rows[1]+5+(1+3*NROWS*NVAL)*NREG*NCH);
			
			ecurr[0]=_mm512_slli_epi16(ecurr[0], GRBITS);
			ecurr[1]=_mm512_slli_epi16(ecurr[1], GRBITS);
			ecurr[2]=_mm512_slli_epi16(ecurr[2], GRBITS);
			ecurr[3]=_mm512_slli_epi16(ecurr[3], GRBITS);
			ecurr[4]=_mm512_slli_epi16(ecurr[4], GRBITS);
			ecurr[5]=_mm512_slli_epi16(ecurr[5], GRBITS);
			ecurr[0]=_mm512_avg_epu16(ecurr[0], _mm512_max_epi16(eNEE[0], eNEEE[0]));
			ecurr[1]=_mm512_avg_epu16(ecurr[1], _mm512_max_epi16(eNEE[1], eNEEE[1]));
			ecurr[2]=_mm512_avg_epu16(ecurr[2], _mm512_max_epi16(eNEE[2], eNEEE[2]));
			ecurr[3]=_mm512_avg_epu16(ecurr[3], _mm512_max_epi16(eNEE[3], eNEEE[3]));
			ecurr[4]=_mm512_avg_epu16(ecurr[4], _mm512_max_epi16(eNEE[4], eNEEE[4]));
			ecurr[5]=_mm512_avg_epu16(ecurr[5], _mm512_max_epi16(eNEE[5], eNEEE[5]));
			eW[0]=_mm512_avg_epu16(eW[0], ecurr[0]);
			eW[1]=_mm512_avg_epu16(eW[1], ecurr[1]);
			eW[2]=_mm512_avg_epu16(eW[2], ecurr[2]);
			eW[3]=_mm512_avg_epu16(eW[3], ecurr[3]);
			eW[4]=_mm512_avg_epu16(eW[4], ecurr[4]);
			eW[5]=_mm512_avg_epu16(eW[5], ecurr[5]);

			_mm512_store_si512((__m512i*)rows[0]+0+(1+0*NROWS*NVAL)*NREG*NCH, eW[0]);//store current contexts
			_mm512_store_si512((__m512i*)rows[0]+1+(1+0*NROWS*NVAL)*NREG*NCH, eW[1]);
			_mm512_store_si512((__m512i*)rows[0]+2+(1+0*NROWS*NVAL)*NREG*NCH, eW[2]);
			_mm512_store_si512((__m512i*)rows[0]+3+(1+0*NROWS*NVAL)*NREG*NCH, eW[3]);
			_mm512_store_si512((__m512i*)rows[0]+4+(1+0*NROWS*NVAL)*NREG*NCH, eW[4]);
			_mm512_store_si512((__m512i*)rows[0]+5+(1+0*NROWS*NVAL)*NREG*NCH, eW[5]);
			eNEE[0]=eNEEE[0];
			eNEE[1]=eNEEE[1];
			eNEE[2]=eNEEE[2];
			eNEE[3]=eNEEE[3];
			eNEE[4]=eNEEE[4];
			eNEE[5]=eNEEE[5];
			NW[0]=N[0];
			NW[1]=N[1];
			NW[2]=N[2];
			NW[3]=N[3];
			NW[4]=N[4];
			NW[5]=N[5];
			rows[0]+=NCODERS*NVAL*NCH*NROWS;
			rows[1]+=NCODERS*NVAL*NCH*NROWS;
			rows[2]+=NCODERS*NVAL*NCH*NROWS;
			rows[3]+=NCODERS*NVAL*NCH*NROWS;
			imptr+=NCH*NCODERS;
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
		_mm_free(L1state);
	if(fwd)//all rANS encoding is bwd-bwd
	{
		rANS_SIMD_SymInfo *syminfo=(rANS_SIMD_SymInfo*)CDF2syms;
		rANS_SIMD_SymInfo *rsyminfo=(rANS_SIMD_SymInfo*)CDF2syms+(ptrdiff_t)3*NCTX*256;
		int32_t *rhist=hists+(ptrdiff_t)3*NCTX*256;
		
		//encode remainder
		if(xremw||yremh)
		{
			for(int ky=0;ky<yremh;++ky)
				decorr1d(image+rowstride*(blockh*YCODERS+ky), iw, 3, bestrct, rhist);
			for(int kx=0;kx<xremw;++kx)
				decorr1d(image+qxbytes+3*kx, blockh*YCODERS, rowstride, bestrct, rhist);
		}

		//normalize/integrate hists
		for(int kc=0;kc<nctx;++kc)
			enc_hist2stats(hists+(ptrdiff_t)256*kc, syminfo+(ptrdiff_t)256*kc, &bypassmask, kc, 0);
			
		if(xremw||yremh)
		{
			uint32_t state=1<<(RANS_STATE_BITS-RANS_RENORM_BITS);
			for(int kx=xremw-1;kx>=0;--kx)
				encode1d(image+qxbytes+3*kx, blockh*YCODERS, rowstride, &state, &streamptr, image, rsyminfo);
			for(int ky=yremh-1;ky>=0;--ky)
				encode1d(image+rowstride*(blockh*YCODERS+ky), iw, 3, &state, &streamptr, image, rsyminfo);
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
		mstate[3]=mstate[2]=mstate[1]=mstate[0]=_mm512_set1_epi32(1<<(RANS_STATE_BITS-RANS_RENORM_BITS));
		uint16_t *ctxptr2=(uint16_t*)(interleaved+(isize<<1)-sizeof(int16_t[NCODERS]));
		for(int ky=blockh-1;ky>=0;--ky)
		{
#ifdef ESTIMATE_SIZE
			int kc=2;
#endif
			for(int kx=3*blockw-1;kx>=0;--kx)//blockw = iw/XCODERS
			{
				__m512i mmax[4], minvf[4], mcdf[4], mnegf_sh[4];
				{
					__m512i s0, s1, s2, s3;
					__m512i t0, t1, t2, t3;
#define SHUFFLE_PS(LO, HI, IMM8_HHLL) _mm512_castps_si512(_mm512_shuffle_ps(_mm512_castsi512_ps(LO), _mm512_castsi512_ps(HI), IMM8_HHLL))

					s0=_mm512_castsi128_si512(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0x0*2+0])));
					s1=_mm512_castsi128_si512(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0x1*2+0])));
					s2=_mm512_castsi128_si512(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0x2*2+0])));
					s3=_mm512_castsi128_si512(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0x3*2+0])));
					s0=_mm512_inserti64x2(s0,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0x4*2+0])), 1);
					s1=_mm512_inserti64x2(s1,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0x5*2+0])), 1);
					s2=_mm512_inserti64x2(s2,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0x6*2+0])), 1);
					s3=_mm512_inserti64x2(s3,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0x7*2+0])), 1);
					s0=_mm512_inserti64x2(s0,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0x8*2+0])), 2);
					s1=_mm512_inserti64x2(s1,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0x9*2+0])), 2);
					s2=_mm512_inserti64x2(s2,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0xA*2+0])), 2);
					s3=_mm512_inserti64x2(s3,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0xB*2+0])), 2);
					s0=_mm512_inserti64x2(s0,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0xC*2+0])), 3);
					s1=_mm512_inserti64x2(s1,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0xD*2+0])), 3);
					s2=_mm512_inserti64x2(s2,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0xE*2+0])), 3);
					s3=_mm512_inserti64x2(s3,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0xF*2+0])), 3);
					t0=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(1, 0, 1, 0));//_MM_TRANSPOSE4_PS
					t2=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(3, 2, 3, 2));
					t1=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(1, 0, 1, 0));
					t3=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(3, 2, 3, 2));
					mmax	[0]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(2, 0, 2, 0));
					minvf	[0]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(3, 1, 3, 1));
					mcdf	[0]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(2, 0, 2, 0));
					mnegf_sh[0]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(3, 1, 3, 1));
					
					s0=_mm512_castsi128_si512(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0x0*2+0])));
					s1=_mm512_castsi128_si512(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0x1*2+0])));
					s2=_mm512_castsi128_si512(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0x2*2+0])));
					s3=_mm512_castsi128_si512(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0x3*2+0])));
					s0=_mm512_inserti64x2(s0,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0x4*2+0])), 1);
					s1=_mm512_inserti64x2(s1,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0x5*2+0])), 1);
					s2=_mm512_inserti64x2(s2,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0x6*2+0])), 1);
					s3=_mm512_inserti64x2(s3,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0x7*2+0])), 1);
					s0=_mm512_inserti64x2(s0,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0x8*2+0])), 2);
					s1=_mm512_inserti64x2(s1,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0x9*2+0])), 2);
					s2=_mm512_inserti64x2(s2,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0xA*2+0])), 2);
					s3=_mm512_inserti64x2(s3,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0xB*2+0])), 2);
					s0=_mm512_inserti64x2(s0,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0xC*2+0])), 3);
					s1=_mm512_inserti64x2(s1,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0xD*2+0])), 3);
					s2=_mm512_inserti64x2(s2,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0xE*2+0])), 3);
					s3=_mm512_inserti64x2(s3,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0xF*2+0])), 3);
					t0=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(1, 0, 1, 0));//_MM_TRANSPOSE4_PS
					t2=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(3, 2, 3, 2));
					t1=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(1, 0, 1, 0));
					t3=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(3, 2, 3, 2));
					mmax	[1]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(2, 0, 2, 0));
					minvf	[1]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(3, 1, 3, 1));
					mcdf	[1]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(2, 0, 2, 0));
					mnegf_sh[1]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(3, 1, 3, 1));
					
					s0=_mm512_castsi128_si512(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0x0*2+1])));
					s1=_mm512_castsi128_si512(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0x1*2+1])));
					s2=_mm512_castsi128_si512(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0x2*2+1])));
					s3=_mm512_castsi128_si512(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0x3*2+1])));
					s0=_mm512_inserti64x2(s0,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0x4*2+1])), 1);
					s1=_mm512_inserti64x2(s1,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0x5*2+1])), 1);
					s2=_mm512_inserti64x2(s2,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0x6*2+1])), 1);
					s3=_mm512_inserti64x2(s3,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0x7*2+1])), 1);
					s0=_mm512_inserti64x2(s0,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0x8*2+1])), 2);
					s1=_mm512_inserti64x2(s1,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0x9*2+1])), 2);
					s2=_mm512_inserti64x2(s2,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0xA*2+1])), 2);
					s3=_mm512_inserti64x2(s3,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0xB*2+1])), 2);
					s0=_mm512_inserti64x2(s0,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0xC*2+1])), 3);
					s1=_mm512_inserti64x2(s1,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0xD*2+1])), 3);
					s2=_mm512_inserti64x2(s2,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0xE*2+1])), 3);
					s3=_mm512_inserti64x2(s3,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*32+0xF*2+1])), 3);
					t0=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(1, 0, 1, 0));//_MM_TRANSPOSE4_PS
					t2=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(3, 2, 3, 2));
					t1=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(1, 0, 1, 0));
					t3=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(3, 2, 3, 2));
					mmax	[2]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(2, 0, 2, 0));
					minvf	[2]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(3, 1, 3, 1));
					mcdf	[2]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(2, 0, 2, 0));
					mnegf_sh[2]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(3, 1, 3, 1));
					
					s0=_mm512_castsi128_si512(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0x0*2+1])));
					s1=_mm512_castsi128_si512(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0x1*2+1])));
					s2=_mm512_castsi128_si512(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0x2*2+1])));
					s3=_mm512_castsi128_si512(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0x3*2+1])));
					s0=_mm512_inserti64x2(s0,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0x4*2+1])), 1);
					s1=_mm512_inserti64x2(s1,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0x5*2+1])), 1);
					s2=_mm512_inserti64x2(s2,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0x6*2+1])), 1);
					s3=_mm512_inserti64x2(s3,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0x7*2+1])), 1);
					s0=_mm512_inserti64x2(s0,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0x8*2+1])), 2);
					s1=_mm512_inserti64x2(s1,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0x9*2+1])), 2);
					s2=_mm512_inserti64x2(s2,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0xA*2+1])), 2);
					s3=_mm512_inserti64x2(s3,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0xB*2+1])), 2);
					s0=_mm512_inserti64x2(s0,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0xC*2+1])), 3);
					s1=_mm512_inserti64x2(s1,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0xD*2+1])), 3);
					s2=_mm512_inserti64x2(s2,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0xE*2+1])), 3);
					s3=_mm512_inserti64x2(s3,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*32+0xF*2+1])), 3);
					t0=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(1, 0, 1, 0));//_MM_TRANSPOSE4_PS
					t2=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(3, 2, 3, 2));
					t1=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(1, 0, 1, 0));
					t3=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(3, 2, 3, 2));
					mmax	[3]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(2, 0, 2, 0));
					minvf	[3]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(3, 1, 3, 1));
					mcdf	[3]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(2, 0, 2, 0));
					mnegf_sh[3]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(3, 1, 3, 1));
					ctxptr2-=NCODERS;
				}
#ifdef ESTIMATE_SIZE
				{
					ALIGN(64) int anegf[NCODERS]={0};
					memcpy(anegf, mnegf_sh, sizeof(anegf));
					const double norm=1./(1<<PROBBITS);
					for(int k=0;k<NCODERS;++k)
					{
						int freq=(1<<PROBBITS)-(anegf[k]&0xFFFF);
						if((uint32_t)(freq-1)>=(uint32_t)((1<<PROBBITS)-1))
							CRASH("freq = %d", freq);
						esize[kc*NCODERS+k]-=log2(freq*norm)*0.125;
					}
				}
				--kc;
				if(kc<0)
					kc=2;
#endif
#ifdef ANS_VAL
				ansval_push(mstate, sizeof(int), NCODERS);
#endif

				//enc renorm		if(state > (freq<<(31-12))-1){write16(state); state>>=16;}
				__mmask16 mask0=_mm512_cmpgt_epi32_mask(mstate[0], mmax[0]);
				__mmask16 mask1=_mm512_cmpgt_epi32_mask(mstate[1], mmax[1]);
				__mmask16 mask2=_mm512_cmpgt_epi32_mask(mstate[2], mmax[2]);
				__mmask16 mask3=_mm512_cmpgt_epi32_mask(mstate[3], mmax[3]);
				int step0=_mm_popcnt_u32(mask0);
				int step1=_mm_popcnt_u32(mask1);
				int step2=_mm_popcnt_u32(mask2);
				int step3=_mm_popcnt_u32(mask3);
				__m256i renorm0=_mm512_cvtepi32_epi16(mstate[0]);
				__m256i renorm1=_mm512_cvtepi32_epi16(mstate[1]);
				__m256i renorm2=_mm512_cvtepi32_epi16(mstate[2]);
				__m256i renorm3=_mm512_cvtepi32_epi16(mstate[3]);
				renorm0=_mm256_maskz_compress_epi16(mask0, renorm0);
				renorm1=_mm256_maskz_compress_epi16(mask1, renorm1);
				renorm2=_mm256_maskz_compress_epi16(mask2, renorm2);
				renorm3=_mm256_maskz_compress_epi16(mask3, renorm3);
				streamptr-=step3*(int)sizeof(int16_t); _mm256_mask_storeu_epi16(streamptr, (uint16_t)((1<<step3)-1), renorm3);
				streamptr-=step2*(int)sizeof(int16_t); _mm256_mask_storeu_epi16(streamptr, (uint16_t)((1<<step2)-1), renorm2);
				streamptr-=step1*(int)sizeof(int16_t); _mm256_mask_storeu_epi16(streamptr, (uint16_t)((1<<step1)-1), renorm1);
				streamptr-=step0*(int)sizeof(int16_t); _mm256_mask_storeu_epi16(streamptr, (uint16_t)((1<<step0)-1), renorm0);
				mstate[0]=_mm512_mask_srli_epi32(mstate[0], mask0, mstate[0], 16);
				mstate[1]=_mm512_mask_srli_epi32(mstate[1], mask1, mstate[1], 16);
				mstate[2]=_mm512_mask_srli_epi32(mstate[2], mask2, mstate[2], 16);
				mstate[3]=_mm512_mask_srli_epi32(mstate[3], mask3, mstate[3], 16);
#ifdef ANS_VAL
				ansval_push(mstate, sizeof(int), NCODERS);
#endif

				//enc update		state += (state*invf>>32>>sh)*negf+cdf		state = state/freq<<12|(cdf+state%freq)
				{
					__m512i lo0=_mm512_mul_epu32(mstate[0], minvf[0]);//q = mulhi32(state, invf)
					__m512i lo1=_mm512_mul_epu32(mstate[1], minvf[1]);
					__m512i lo2=_mm512_mul_epu32(mstate[2], minvf[2]);
					__m512i lo3=_mm512_mul_epu32(mstate[3], minvf[3]);
					__m512i hi0=_mm512_mul_epu32(_mm512_srli_epi64(mstate[0], 32), _mm512_srli_epi64(minvf[0], 32));
					__m512i hi1=_mm512_mul_epu32(_mm512_srli_epi64(mstate[1], 32), _mm512_srli_epi64(minvf[1], 32));
					__m512i hi2=_mm512_mul_epu32(_mm512_srli_epi64(mstate[2], 32), _mm512_srli_epi64(minvf[2], 32));
					__m512i hi3=_mm512_mul_epu32(_mm512_srli_epi64(mstate[3], 32), _mm512_srli_epi64(minvf[3], 32));
					minvf[0]=_mm512_mask_blend_epi32(0xAAAA, _mm512_srli_epi64(lo0, 32), hi0);
					minvf[1]=_mm512_mask_blend_epi32(0xAAAA, _mm512_srli_epi64(lo1, 32), hi1);
					minvf[2]=_mm512_mask_blend_epi32(0xAAAA, _mm512_srli_epi64(lo2, 32), hi2);
					minvf[3]=_mm512_mask_blend_epi32(0xAAAA, _mm512_srli_epi64(lo3, 32), hi3);
				}
				{
					__m512i sh0=_mm512_srli_epi32(mnegf_sh[0], 16);
					__m512i sh1=_mm512_srli_epi32(mnegf_sh[1], 16);
					__m512i sh2=_mm512_srli_epi32(mnegf_sh[2], 16);
					__m512i sh3=_mm512_srli_epi32(mnegf_sh[3], 16);
					minvf[0]=_mm512_srlv_epi32(minvf[0], sh0);
					minvf[1]=_mm512_srlv_epi32(minvf[1], sh1);
					minvf[2]=_mm512_srlv_epi32(minvf[2], sh2);
					minvf[3]=_mm512_srlv_epi32(minvf[3], sh3);
				}
				mstate[0]=_mm512_add_epi32(mstate[0], mcdf[0]);
				mstate[1]=_mm512_add_epi32(mstate[1], mcdf[1]);
				mstate[2]=_mm512_add_epi32(mstate[2], mcdf[2]);
				mstate[3]=_mm512_add_epi32(mstate[3], mcdf[3]);
				{
					__m512i negf0=_mm512_maskz_mov_epi16(0x55555555, mnegf_sh[0]);
					__m512i negf1=_mm512_maskz_mov_epi16(0x55555555, mnegf_sh[1]);
					__m512i negf2=_mm512_maskz_mov_epi16(0x55555555, mnegf_sh[2]);
					__m512i negf3=_mm512_maskz_mov_epi16(0x55555555, mnegf_sh[3]);
					minvf[0]=_mm512_mullo_epi32(minvf[0], negf0);
					minvf[1]=_mm512_mullo_epi32(minvf[1], negf1);
					minvf[2]=_mm512_mullo_epi32(minvf[2], negf2);
					minvf[3]=_mm512_mullo_epi32(minvf[3], negf3);
				}
#ifdef ANS_VAL
				{
					__m512i one=_mm512_set1_epi32(1);
					mmax[0]=_mm512_add_epi32(mmax[0], one);
					mmax[1]=_mm512_add_epi32(mmax[1], one);
					mmax[2]=_mm512_add_epi32(mmax[2], one);
					mmax[3]=_mm512_add_epi32(mmax[3], one);
					mmax[0]=_mm512_srli_epi32(mmax[0], RANS_STATE_BITS-PROBBITS);
					mmax[1]=_mm512_srli_epi32(mmax[1], RANS_STATE_BITS-PROBBITS);
					mmax[2]=_mm512_srli_epi32(mmax[2], RANS_STATE_BITS-PROBBITS);
					mmax[3]=_mm512_srli_epi32(mmax[3], RANS_STATE_BITS-PROBBITS);
					mmax[2]=_mm512_slli_epi32(mmax[2], 16);
					mmax[3]=_mm512_slli_epi32(mmax[3], 16);
					mmax[0]=_mm512_or_si512(mmax[0], mmax[2]);
					mmax[1]=_mm512_or_si512(mmax[1], mmax[3]);
					ansval_push(mmax, sizeof(int16_t), NCODERS);
				}
#endif
				mstate[0]=_mm512_add_epi32(mstate[0], minvf[0]);
				mstate[1]=_mm512_add_epi32(mstate[1], minvf[1]);
				mstate[2]=_mm512_add_epi32(mstate[2], minvf[2]);
				mstate[3]=_mm512_add_epi32(mstate[3], minvf[3]);
			}
			//printf("%8d/%8d\r", blockh-1-ky, blockh);
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
			FILE *fdst=fopen(dstfn, "wb");
			if(!fdst)
			{
				CRASH("Cannot open \"%s\" for writing", dstfn);
				return 1;
			}
			ptrdiff_t csize2=0;
			csize2+=fwrite("V1", 1, 2, fdst);
			csize2+=fwrite(&iw, 1, 3, fdst);
			csize2+=fwrite(&ih, 1, 3, fdst);
			{
				int flags=bestrct<<2|effort;
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
			for(int k=0;k<3*NCODERS;++k)
			{
				etotal+=esize[k];
				printf("E %12.2lf\n", esize[k]);//
			}
			printf("Total estimate  %12.2lf bytes\n", etotal);
#endif
#ifdef LOUD
			printf("L1C AVX512 WH %d*%d  RCT %2d %s  effort %d  dist %3d  \"%s\"\n"
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
			(void)csize2;
			prof_checkpoint(csize2, "fwrite");
		}
		free(hists);
#ifdef SAVE_RESIDUALS
		{
			uint8_t *result=(uint8_t*)malloc(usize);
			if(!result)
			{
				CRASH("Alloc error");
				return 1;
			}
			memset(result, -128, usize);
			interleave_blocks_inv(residuals, iw, ih, result);
			{
				const char fn[]="20250605_0804pm.ppm";
				FILE *fdst2=fopen(fn, "wb");
				if(!fdst2)
				{
					CRASH("Cannot open \"%s\" for writing", fn);
					return 1;
				}
				fprintf(fdst2, "P6\n%d %d\n255\n", iw, ih);
				fwrite(result, 1, usize, fdst2);
				fclose(fdst2);
			}
			free(residuals);
			free(result);
		}
#endif
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
	_mm_free(pixels);
	_mm_free(CDF2syms);
	_mm_free(interleaved);
	free(image);

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
	return 0;
}