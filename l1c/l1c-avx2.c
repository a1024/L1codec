#ifdef _MSC_VER
#	ifndef _CRT_SECURE_NO_WARNINGS
#		define _CRT_SECURE_NO_WARNINGS
#	endif
#elif defined __linux__ && !defined _GNU_SOURCE
#	define _GNU_SOURCE
#	include<stddef.h>//ptrdiff_t
#endif
#include<math.h>
#include<immintrin.h>
#include<sys/stat.h>


	#define PROFILE_TIME		//should be on

#ifdef _MSC_VER
	#define PROFILE_SIZE
	#define LOUD			//size & time

//	#define ESTIMATE_SIZE		//DEBUG		checks for zero frequency, visualizes context usage

//	#define SAVE_RESIDUALS
//	#define PRINT_L1_BOUNDS
//	#define TEST_INTERLEAVE
#endif

//	#define ANALYSIS_GRAD
	#define ENABLE_RCT_EXTENSION
	#define INTERLEAVESIMD		//2.5x faster interleave
//	#define EMULATE_GATHER		//gather is a little faster


enum
{
	XCODERS=4,
	YCODERS=4,
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
};

#define COMMON_rANS
#include"common.h"
AWM_INLINE void gather32(int *dst, const int *src, const int *offsets)
{
#ifdef EMULATE_GATHER
	volatile int *ptr=dst;
	ptr[0]=src[offsets[0]];
	ptr[1]=src[offsets[1]];
	ptr[2]=src[offsets[2]];
	ptr[3]=src[offsets[3]];
	ptr[4]=src[offsets[4]];
	ptr[5]=src[offsets[5]];
	ptr[6]=src[offsets[6]];
	ptr[7]=src[offsets[7]];
#else
	_mm256_store_si256((__m256i*)dst, _mm256_i32gather_epi32(src, _mm256_load_si256((__m256i*)offsets), sizeof(int)));
#endif
}

AWM_INLINE void dec_yuv(
	__m256i *mstate,
	const __m256i *ctx0,
	const uint32_t *CDF2syms,
	const int *ans_permute,
	uint8_t **pstreamptr,
	const uint8_t *streamend,
	__m256i *syms
)
{
	const uint8_t *streamptr=*pstreamptr;
	__m256i decctx[2];
	{
		decctx[1]=_mm256_cvtepi16_epi32(_mm256_extracti128_si256(*ctx0, 1));
		decctx[0]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(*ctx0));
	}
	decctx[0]=_mm256_slli_epi32(decctx[0], PROBBITS);
	decctx[1]=_mm256_slli_epi32(decctx[1], PROBBITS);
#ifdef ANS_VAL
	ansval_check(mstate, sizeof(int), NCODERS);
#endif
	{
		__m256i mprobmask=_mm256_set1_epi32((1<<PROBBITS)-1);
		__m256i rem0=_mm256_and_si256(mstate[0], mprobmask);
		__m256i rem1=_mm256_and_si256(mstate[1], mprobmask);
		decctx[0]=_mm256_or_si256(decctx[0], rem0);
		decctx[1]=_mm256_or_si256(decctx[1], rem1);
	}
#ifdef ANS_VAL
	ALIGN(32) int debugctx[NCODERS];
	memcpy(debugctx, decctx, sizeof(int[NCODERS]));
#endif
	gather32((int*)(decctx+0), (const int*)CDF2syms, (int*)(decctx+0));//James Bonfield's CDF2sym: {freq<<20 | bias<<8 | sym}
	gather32((int*)(decctx+1), (const int*)CDF2syms, (int*)(decctx+1));
	//decctx[0]=_mm256_i32gather_epi32(statsptr, decctx[0], sizeof(int));
	//decctx[1]=_mm256_i32gather_epi32(statsptr, decctx[1], sizeof(int));

	//update		state = (state>>12)*freq+(rem-cdf)	rem-cdf is prebaked
	{
		__m256i mfreq0=_mm256_srli_epi32(decctx[0], PROBBITS+8);//1 <= freq <= 0xF01
		__m256i mfreq1=_mm256_srli_epi32(decctx[1], PROBBITS+8);
#ifdef ANS_VAL
		__m256i mdebugfreq[1];
		mdebugfreq[0]=_mm256_packus_epi32(mfreq0, mfreq1);
		mdebugfreq[0]=_mm256_permute4x64_epi64(mdebugfreq[0], _MM_SHUFFLE(3, 1, 2, 0));
		ALIGN(32) uint16_t freqs[NCODERS];
		memcpy(freqs, mdebugfreq, sizeof(freqs));
		ansval_check(freqs, sizeof(int16_t), NCODERS);
#endif
		mstate[0]=_mm256_srli_epi32(mstate[0], PROBBITS);
		mstate[1]=_mm256_srli_epi32(mstate[1], PROBBITS);
		mstate[0]=_mm256_mullo_epi32(mstate[0], mfreq0);//10 cycles
		mstate[1]=_mm256_mullo_epi32(mstate[1], mfreq1);
	}
	{
		__m256i mbias0=_mm256_slli_epi32(decctx[0], PROBBITS);
		__m256i mbias1=_mm256_slli_epi32(decctx[1], PROBBITS);
		mbias0=_mm256_srli_epi32(mbias0, 32-PROBBITS);
		mbias1=_mm256_srli_epi32(mbias1, 32-PROBBITS);
		mstate[0]=_mm256_add_epi32(mstate[0], mbias0);
		mstate[1]=_mm256_add_epi32(mstate[1], mbias1);
	}
	__m256i symmask=_mm256_set1_epi32(255);
	decctx[0]=_mm256_and_si256(decctx[0], symmask);
	decctx[1]=_mm256_and_si256(decctx[1], symmask);
	decctx[0]=_mm256_packus_epi16(decctx[0], decctx[1]);
	syms[0]=_mm256_permute4x64_epi64(decctx[0], _MM_SHUFFLE(3, 1, 2, 0));
#ifdef ANS_VAL
	ansval_check(mstate, sizeof(int), NCODERS);
#endif
	//renorm
	{
		__m256i smin=_mm256_set1_epi32(1<<(RANS_STATE_BITS-RANS_RENORM_BITS));
#ifdef _DEBUG
		if(streamptr>streamend)
			CRASH("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
		__m256i lo0=_mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i*)streamptr));
		__m256i cond0=_mm256_cmpgt_epi32(smin, mstate[0]);//FIXME this is signed comparison
		__m256i cond1=_mm256_cmpgt_epi32(smin, mstate[1]);
		int mask0=_mm256_movemask_ps(_mm256_castsi256_ps(cond0));
		int mask1=_mm256_movemask_ps(_mm256_castsi256_ps(cond1));
		__m256i idx0=_mm256_load_si256((const __m256i*)ans_permute+mask0);
		__m256i idx1=_mm256_load_si256((const __m256i*)ans_permute+mask1);
		mask0=_mm_popcnt_u32(mask0);
		mask1=_mm_popcnt_u32(mask1);
		__m256i renorm0=_mm256_slli_epi32(mstate[0], 16);
		__m256i renorm1=_mm256_slli_epi32(mstate[1], 16);
		streamptr+=mask0*sizeof(int16_t);
		lo0=_mm256_permutevar8x32_epi32(lo0, idx0);
#ifdef _DEBUG
		if(streamptr>streamend)
			CRASH("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
		__m256i lo1=_mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i*)streamptr));
		streamptr+=mask1*sizeof(int16_t);
		lo1=_mm256_permutevar8x32_epi32(lo1, idx1);
		renorm0=_mm256_or_si256(renorm0, lo0);
		renorm1=_mm256_or_si256(renorm1, lo1);

		mstate[0]=_mm256_blendv_epi8(mstate[0], renorm0, cond0);
		mstate[1]=_mm256_blendv_epi8(mstate[1], renorm1, cond1);
	}
	*pstreamptr=(uint8_t*)(size_t)streamptr;
}
AWM_INLINE void transpose16(__m128i *data)
{
#if 1
	__m128i a[16], b[16];
	a[0x0]=_mm_load_si128((__m128i*)data+0x0);
	a[0x1]=_mm_load_si128((__m128i*)data+0x1);
	a[0x2]=_mm_load_si128((__m128i*)data+0x2);
	a[0x3]=_mm_load_si128((__m128i*)data+0x3);
	a[0x4]=_mm_load_si128((__m128i*)data+0x4);
	a[0x5]=_mm_load_si128((__m128i*)data+0x5);
	a[0x6]=_mm_load_si128((__m128i*)data+0x6);
	a[0x7]=_mm_load_si128((__m128i*)data+0x7);
	a[0x8]=_mm_load_si128((__m128i*)data+0x8);
	a[0x9]=_mm_load_si128((__m128i*)data+0x9);
	a[0xA]=_mm_load_si128((__m128i*)data+0xA);
	a[0xB]=_mm_load_si128((__m128i*)data+0xB);
	a[0xC]=_mm_load_si128((__m128i*)data+0xC);
	a[0xD]=_mm_load_si128((__m128i*)data+0xD);
	a[0xE]=_mm_load_si128((__m128i*)data+0xE);
	a[0xF]=_mm_load_si128((__m128i*)data+0xF);

	b[0x0]=_mm_unpacklo_epi8(a[0x0], a[0x1]);
	b[0x1]=_mm_unpackhi_epi8(a[0x0], a[0x1]);
	b[0x2]=_mm_unpacklo_epi8(a[0x2], a[0x3]);
	b[0x3]=_mm_unpackhi_epi8(a[0x2], a[0x3]);
	b[0x4]=_mm_unpacklo_epi8(a[0x4], a[0x5]);
	b[0x5]=_mm_unpackhi_epi8(a[0x4], a[0x5]);
	b[0x6]=_mm_unpacklo_epi8(a[0x6], a[0x7]);
	b[0x7]=_mm_unpackhi_epi8(a[0x6], a[0x7]);
	b[0x8]=_mm_unpacklo_epi8(a[0x8], a[0x9]);
	b[0x9]=_mm_unpackhi_epi8(a[0x8], a[0x9]);
	b[0xA]=_mm_unpacklo_epi8(a[0xA], a[0xB]);
	b[0xB]=_mm_unpackhi_epi8(a[0xA], a[0xB]);
	b[0xC]=_mm_unpacklo_epi8(a[0xC], a[0xD]);
	b[0xD]=_mm_unpackhi_epi8(a[0xC], a[0xD]);
	b[0xE]=_mm_unpacklo_epi8(a[0xE], a[0xF]);
	b[0xF]=_mm_unpackhi_epi8(a[0xE], a[0xF]);

	a[0x0]=_mm_unpacklo_epi16(b[0x0], b[0x2]);
	a[0x1]=_mm_unpackhi_epi16(b[0x0], b[0x2]);
	a[0x2]=_mm_unpacklo_epi16(b[0x1], b[0x3]);
	a[0x3]=_mm_unpackhi_epi16(b[0x1], b[0x3]);
	a[0x4]=_mm_unpacklo_epi16(b[0x4], b[0x6]);
	a[0x5]=_mm_unpackhi_epi16(b[0x4], b[0x6]);
	a[0x6]=_mm_unpacklo_epi16(b[0x5], b[0x7]);
	a[0x7]=_mm_unpackhi_epi16(b[0x5], b[0x7]);
	a[0x8]=_mm_unpacklo_epi16(b[0x8], b[0xA]);
	a[0x9]=_mm_unpackhi_epi16(b[0x8], b[0xA]);
	a[0xA]=_mm_unpacklo_epi16(b[0x9], b[0xB]);
	a[0xB]=_mm_unpackhi_epi16(b[0x9], b[0xB]);
	a[0xC]=_mm_unpacklo_epi16(b[0xC], b[0xE]);
	a[0xD]=_mm_unpackhi_epi16(b[0xC], b[0xE]);
	a[0xE]=_mm_unpacklo_epi16(b[0xD], b[0xF]);
	a[0xF]=_mm_unpackhi_epi16(b[0xD], b[0xF]);

	b[0x0]=_mm_unpacklo_epi32(a[0x0], a[0x4]);
	b[0x1]=_mm_unpackhi_epi32(a[0x0], a[0x4]);
	b[0x2]=_mm_unpacklo_epi32(a[0x1], a[0x5]);
	b[0x3]=_mm_unpackhi_epi32(a[0x1], a[0x5]);
	b[0x4]=_mm_unpacklo_epi32(a[0x2], a[0x6]);
	b[0x5]=_mm_unpackhi_epi32(a[0x2], a[0x6]);
	b[0x6]=_mm_unpacklo_epi32(a[0x3], a[0x7]);
	b[0x7]=_mm_unpackhi_epi32(a[0x3], a[0x7]);
	b[0x8]=_mm_unpacklo_epi32(a[0x8], a[0xC]);
	b[0x9]=_mm_unpackhi_epi32(a[0x8], a[0xC]);
	b[0xA]=_mm_unpacklo_epi32(a[0x9], a[0xD]);
	b[0xB]=_mm_unpackhi_epi32(a[0x9], a[0xD]);
	b[0xC]=_mm_unpacklo_epi32(a[0xA], a[0xE]);
	b[0xD]=_mm_unpackhi_epi32(a[0xA], a[0xE]);
	b[0xE]=_mm_unpacklo_epi32(a[0xB], a[0xF]);
	b[0xF]=_mm_unpackhi_epi32(a[0xB], a[0xF]);

	a[0x0]=_mm_unpacklo_epi64(b[0x0], b[0x8]);
	a[0x1]=_mm_unpackhi_epi64(b[0x0], b[0x8]);
	a[0x2]=_mm_unpacklo_epi64(b[0x1], b[0x9]);
	a[0x3]=_mm_unpackhi_epi64(b[0x1], b[0x9]);
	a[0x4]=_mm_unpacklo_epi64(b[0x2], b[0xA]);
	a[0x5]=_mm_unpackhi_epi64(b[0x2], b[0xA]);
	a[0x6]=_mm_unpacklo_epi64(b[0x3], b[0xB]);
	a[0x7]=_mm_unpackhi_epi64(b[0x3], b[0xB]);
	a[0x8]=_mm_unpacklo_epi64(b[0x4], b[0xC]);
	a[0x9]=_mm_unpackhi_epi64(b[0x4], b[0xC]);
	a[0xA]=_mm_unpacklo_epi64(b[0x5], b[0xD]);
	a[0xB]=_mm_unpackhi_epi64(b[0x5], b[0xD]);
	a[0xC]=_mm_unpacklo_epi64(b[0x6], b[0xE]);
	a[0xD]=_mm_unpackhi_epi64(b[0x6], b[0xE]);
	a[0xE]=_mm_unpacklo_epi64(b[0x7], b[0xF]);
	a[0xF]=_mm_unpackhi_epi64(b[0x7], b[0xF]);

	_mm_store_si128((__m128i*)data+0x0, a[0x0]);
	_mm_store_si128((__m128i*)data+0x1, a[0x1]);
	_mm_store_si128((__m128i*)data+0x2, a[0x2]);
	_mm_store_si128((__m128i*)data+0x3, a[0x3]);
	_mm_store_si128((__m128i*)data+0x4, a[0x4]);
	_mm_store_si128((__m128i*)data+0x5, a[0x5]);
	_mm_store_si128((__m128i*)data+0x6, a[0x6]);
	_mm_store_si128((__m128i*)data+0x7, a[0x7]);
	_mm_store_si128((__m128i*)data+0x8, a[0x8]);
	_mm_store_si128((__m128i*)data+0x9, a[0x9]);
	_mm_store_si128((__m128i*)data+0xA, a[0xA]);
	_mm_store_si128((__m128i*)data+0xB, a[0xB]);
	_mm_store_si128((__m128i*)data+0xC, a[0xC]);
	_mm_store_si128((__m128i*)data+0xD, a[0xD]);
	_mm_store_si128((__m128i*)data+0xE, a[0xE]);
	_mm_store_si128((__m128i*)data+0xF, a[0xF]);
#endif

#if 0
	//https://pzemtsov.github.io/2014/10/01/how-to-transpose-a-16x16-matrix.html
	__m128i t0, t1, t2, t3;
#define SHUFFLE(LO, HI, IMM8_HHLL) _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(LO), _mm_castsi128_ps(HI), IMM8_HHLL))
#define TRANSPOSE4(S0, S1, S2, S3,  T0, T1, T2, T3,  D0, D1, D2, D3)\
	do\
	{\
		T0=SHUFFLE(S0, S1, _MM_SHUFFLE(1, 0, 1, 0));\
		T2=SHUFFLE(S0, S1, _MM_SHUFFLE(3, 2, 3, 2));\
		T1=SHUFFLE(S2, S3, _MM_SHUFFLE(1, 0, 1, 0));\
		T3=SHUFFLE(S2, S3, _MM_SHUFFLE(3, 2, 3, 2));\
		D0=SHUFFLE(T0, T1, _MM_SHUFFLE(2, 0, 2, 0));\
		D1=SHUFFLE(T0, T1, _MM_SHUFFLE(3, 1, 3, 1));\
		D2=SHUFFLE(T2, T3, _MM_SHUFFLE(2, 0, 2, 0));\
		D3=SHUFFLE(T2, T3, _MM_SHUFFLE(3, 1, 3, 1));\
	}while(0)
	TRANSPOSE4(data[0x0], data[0x1], data[0x2], data[0x3],  t0, t1, t2, t3,  data[0x0], data[0x1], data[0x2], data[0x3]);
	TRANSPOSE4(data[0x4], data[0x5], data[0x6], data[0x7],  t0, t1, t2, t3,  data[0x4], data[0x5], data[0x6], data[0x7]);
	TRANSPOSE4(data[0x8], data[0x9], data[0xA], data[0xB],  t0, t1, t2, t3,  data[0x8], data[0x9], data[0xA], data[0xB]);
	TRANSPOSE4(data[0xC], data[0xD], data[0xE], data[0xF],  t0, t1, t2, t3,  data[0xC], data[0xD], data[0xE], data[0xF]);
	t0=_mm_set_epi8(
		15, 11,  7,  3,
		14, 10,  6,  2,
		13,  9,  5,  1,
		12,  8,  4,  0
	);
	/*
	transpose the 4*4 registers themselves while shuffling
	exchange:
	1 <-> 4
	2 <-> 8
	3 <-> C
	6 <-> 9
	7 <-> D
	B <-> E
	*/
	//diagonals
	data[0x0]=_mm_shuffle_epi8(data[0x0], t0);
	data[0x5]=_mm_shuffle_epi8(data[0x5], t0);
	data[0xA]=_mm_shuffle_epi8(data[0xA], t0);
	data[0xF]=_mm_shuffle_epi8(data[0xF], t0);
	//nondiagonals
	t1=data[0x4]; data[0x4]=_mm_shuffle_epi8(data[0x1], t0); data[0x1]=_mm_shuffle_epi8(t1, t0);
	t1=data[0x8]; data[0x8]=_mm_shuffle_epi8(data[0x2], t0); data[0x2]=_mm_shuffle_epi8(t1, t0);
	t1=data[0xC]; data[0xC]=_mm_shuffle_epi8(data[0x3], t0); data[0x3]=_mm_shuffle_epi8(t1, t0);
	t1=data[0x9]; data[0x9]=_mm_shuffle_epi8(data[0x6], t0); data[0x6]=_mm_shuffle_epi8(t1, t0);
	t1=data[0xD]; data[0xD]=_mm_shuffle_epi8(data[0x7], t0); data[0x7]=_mm_shuffle_epi8(t1, t0);
	t1=data[0xE]; data[0xE]=_mm_shuffle_epi8(data[0xB], t0); data[0xB]=_mm_shuffle_epi8(t1, t0);
	TRANSPOSE4(data[0x0], data[0x1], data[0x2], data[0x3],  t0, t1, t2, t3,  data[0x0], data[0x1], data[0x2], data[0x3]);
	TRANSPOSE4(data[0x4], data[0x5], data[0x6], data[0x7],  t0, t1, t2, t3,  data[0x4], data[0x5], data[0x6], data[0x7]);
	TRANSPOSE4(data[0x8], data[0x9], data[0xA], data[0xB],  t0, t1, t2, t3,  data[0x8], data[0x9], data[0xA], data[0xB]);
	TRANSPOSE4(data[0xC], data[0xD], data[0xE], data[0xF],  t0, t1, t2, t3,  data[0xC], data[0xD], data[0xE], data[0xF]);
#undef  SHUFFLE
#undef  TRANSPOSE4
#endif
}
static void interleave_blocks_fwd(const uint8_t *original, int iw, int ih, uint8_t *interleaved)
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
	int SIMDxcount=blockxbytes&~((int)sizeof(__m128i[NCODERS])-1);
	__m256i slowinc=_mm256_set1_epi64x(sizeof(__m128i));
#endif
	uint8_t *fastptr=interleaved;
	ALIGN(32) const uint8_t *slowptrs[NCODERS]={0}, *slowptrs0[NCODERS]={0};
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
		for(;kx<SIMDxcount;kx+=(int)sizeof(__m128i[NCODERS]))
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
			__m128i block[16];
			block[0x0]=_mm_loadu_si128((__m128i*)slowptrs[0x00]);
			block[0x1]=_mm_loadu_si128((__m128i*)slowptrs[0x01]);
			block[0x2]=_mm_loadu_si128((__m128i*)slowptrs[0x02]);
			block[0x3]=_mm_loadu_si128((__m128i*)slowptrs[0x03]);
			block[0x4]=_mm_loadu_si128((__m128i*)slowptrs[0x04]);
			block[0x5]=_mm_loadu_si128((__m128i*)slowptrs[0x05]);
			block[0x6]=_mm_loadu_si128((__m128i*)slowptrs[0x06]);
			block[0x7]=_mm_loadu_si128((__m128i*)slowptrs[0x07]);
			block[0x8]=_mm_loadu_si128((__m128i*)slowptrs[0x08]);
			block[0x9]=_mm_loadu_si128((__m128i*)slowptrs[0x09]);
			block[0xA]=_mm_loadu_si128((__m128i*)slowptrs[0x0A]);
			block[0xB]=_mm_loadu_si128((__m128i*)slowptrs[0x0B]);
			block[0xC]=_mm_loadu_si128((__m128i*)slowptrs[0x0C]);
			block[0xD]=_mm_loadu_si128((__m128i*)slowptrs[0x0D]);
			block[0xE]=_mm_loadu_si128((__m128i*)slowptrs[0x0E]);
			block[0xF]=_mm_loadu_si128((__m128i*)slowptrs[0x0F]);
			transpose16(block);
			_mm_store_si128((__m128i*)fastptr+0x0, block[0x0]);
			_mm_store_si128((__m128i*)fastptr+0x1, block[0x1]);
			_mm_store_si128((__m128i*)fastptr+0x2, block[0x2]);
			_mm_store_si128((__m128i*)fastptr+0x3, block[0x3]);
			_mm_store_si128((__m128i*)fastptr+0x4, block[0x4]);
			_mm_store_si128((__m128i*)fastptr+0x5, block[0x5]);
			_mm_store_si128((__m128i*)fastptr+0x6, block[0x6]);
			_mm_store_si128((__m128i*)fastptr+0x7, block[0x7]);
			_mm_store_si128((__m128i*)fastptr+0x8, block[0x8]);
			_mm_store_si128((__m128i*)fastptr+0x9, block[0x9]);
			_mm_store_si128((__m128i*)fastptr+0xA, block[0xA]);
			_mm_store_si128((__m128i*)fastptr+0xB, block[0xB]);
			_mm_store_si128((__m128i*)fastptr+0xC, block[0xC]);
			_mm_store_si128((__m128i*)fastptr+0xD, block[0xD]);
			_mm_store_si128((__m128i*)fastptr+0xE, block[0xE]);
			_mm_store_si128((__m128i*)fastptr+0xF, block[0xF]);

			fastptr+=sizeof(__m128i[NCODERS]);
			{
				__m256i mp[4];
				mp[0]=_mm256_add_epi64(slowinc, _mm256_load_si256((__m256i*)slowptrs+0));
				mp[1]=_mm256_add_epi64(slowinc, _mm256_load_si256((__m256i*)slowptrs+1));
				mp[2]=_mm256_add_epi64(slowinc, _mm256_load_si256((__m256i*)slowptrs+2));
				mp[3]=_mm256_add_epi64(slowinc, _mm256_load_si256((__m256i*)slowptrs+3));
				_mm256_store_si256((__m256i*)slowptrs+0, mp[0]);
				_mm256_store_si256((__m256i*)slowptrs+1, mp[1]);
				_mm256_store_si256((__m256i*)slowptrs+2, mp[2]);
				_mm256_store_si256((__m256i*)slowptrs+3, mp[3]);
			}
		}
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
			fastptr+=NCODERS;
//#if defined __GNUC__ && !defined PROFILER
//#pragma GCC unroll 32
//#endif
//			for(int k=0;k<NCODERS;++k)
//				*fastptr++=*slowptrs[k]++;
		}
#endif
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
	int SIMDxcount=blockxbytes&~((int)sizeof(__m128i[NCODERS])-1);
	__m256i slowinc=_mm256_set1_epi64x(sizeof(__m128i));
#endif
	const uint8_t *fastptr=interleaved;
	ALIGN(32) uint8_t *slowptrs[NCODERS]={0}, *slowptrs0[NCODERS]={0};
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
		for(;kx<SIMDxcount;kx+=(int)sizeof(__m128i[NCODERS]))
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
			__m128i block[16];
			block[0x0]=_mm_load_si128((__m128i*)fastptr+0x0);
			block[0x1]=_mm_load_si128((__m128i*)fastptr+0x1);
			block[0x2]=_mm_load_si128((__m128i*)fastptr+0x2);
			block[0x3]=_mm_load_si128((__m128i*)fastptr+0x3);
			block[0x4]=_mm_load_si128((__m128i*)fastptr+0x4);
			block[0x5]=_mm_load_si128((__m128i*)fastptr+0x5);
			block[0x6]=_mm_load_si128((__m128i*)fastptr+0x6);
			block[0x7]=_mm_load_si128((__m128i*)fastptr+0x7);
			block[0x8]=_mm_load_si128((__m128i*)fastptr+0x8);
			block[0x9]=_mm_load_si128((__m128i*)fastptr+0x9);
			block[0xA]=_mm_load_si128((__m128i*)fastptr+0xA);
			block[0xB]=_mm_load_si128((__m128i*)fastptr+0xB);
			block[0xC]=_mm_load_si128((__m128i*)fastptr+0xC);
			block[0xD]=_mm_load_si128((__m128i*)fastptr+0xD);
			block[0xE]=_mm_load_si128((__m128i*)fastptr+0xE);
			block[0xF]=_mm_load_si128((__m128i*)fastptr+0xF);
			transpose16(block);
			_mm_storeu_si128((__m128i*)slowptrs[0x00], block[0x0]);
			_mm_storeu_si128((__m128i*)slowptrs[0x01], block[0x1]);
			_mm_storeu_si128((__m128i*)slowptrs[0x02], block[0x2]);
			_mm_storeu_si128((__m128i*)slowptrs[0x03], block[0x3]);
			_mm_storeu_si128((__m128i*)slowptrs[0x04], block[0x4]);
			_mm_storeu_si128((__m128i*)slowptrs[0x05], block[0x5]);
			_mm_storeu_si128((__m128i*)slowptrs[0x06], block[0x6]);
			_mm_storeu_si128((__m128i*)slowptrs[0x07], block[0x7]);
			_mm_storeu_si128((__m128i*)slowptrs[0x08], block[0x8]);
			_mm_storeu_si128((__m128i*)slowptrs[0x09], block[0x9]);
			_mm_storeu_si128((__m128i*)slowptrs[0x0A], block[0xA]);
			_mm_storeu_si128((__m128i*)slowptrs[0x0B], block[0xB]);
			_mm_storeu_si128((__m128i*)slowptrs[0x0C], block[0xC]);
			_mm_storeu_si128((__m128i*)slowptrs[0x0D], block[0xD]);
			_mm_storeu_si128((__m128i*)slowptrs[0x0E], block[0xE]);
			_mm_storeu_si128((__m128i*)slowptrs[0x0F], block[0xF]);

			fastptr+=sizeof(__m128i[NCODERS]);
			{
				__m256i mp[4];
				mp[0]=_mm256_add_epi64(slowinc, _mm256_load_si256((__m256i*)slowptrs+0));
				mp[1]=_mm256_add_epi64(slowinc, _mm256_load_si256((__m256i*)slowptrs+1));
				mp[2]=_mm256_add_epi64(slowinc, _mm256_load_si256((__m256i*)slowptrs+2));
				mp[3]=_mm256_add_epi64(slowinc, _mm256_load_si256((__m256i*)slowptrs+3));
				_mm256_store_si256((__m256i*)slowptrs+0, mp[0]);
				_mm256_store_si256((__m256i*)slowptrs+1, mp[1]);
				_mm256_store_si256((__m256i*)slowptrs+2, mp[2]);
				_mm256_store_si256((__m256i*)slowptrs+3, mp[3]);
			}
		}
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
			fastptr+=NCODERS;
//#if defined __GNUC__ && !defined PROFILER
//#pragma GCC unroll 32
//#endif
//			for(int k=0;k<NCODERS;++k)
//				*slowptrs[k]++=*fastptr++;
		}
#endif
		for(int k=0;k<NCODERS;++k)
			slowptrs0[k]+=rowstride;
	}
}
int codec_l1_avx2(int argc, char **argv)
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
	double esize[3*NCODERS]={0};
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
		image=(uint8_t*)malloc(cap+sizeof(__m256i));
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
			streamptr=streamstart=image+cap-(csize-cheadersize)-sizeof(__m256i);
			streamend=image+cap-sizeof(__m256i);
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
	uint8_t *interleaved=(uint8_t*)_mm_malloc(interleavedsize, sizeof(__m256i));
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
	uint32_t *CDF2syms=(uint32_t*)_mm_malloc(CDF2syms_size, sizeof(__m256i));

	int ans_permute_size=sizeof(__m256i[256]);
	int *ans_permute=(int*)_mm_malloc(ans_permute_size, sizeof(__m256i));

	psize=(blockw+2*XPAD)*(int)sizeof(int16_t[NCH*NROWS*NVAL*NCODERS]);//int16_t[blockw+2*XPAD][NCH*NROWS*NVAL*NCODERS]
	pixels=(int16_t*)_mm_malloc(psize, sizeof(__m256i));
	if((fwd&&!hists)||!CDF2syms||!ans_permute||!pixels)
	{
		CRASH("Alloc error");
		return 1;
	}
	memset(ans_permute, 0, ans_permute_size);//_mm256_permutevar8x32_epi32 can't clear elements like _mm256_shuffle_epi8 does
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
			ALIGN(32) int64_t counters[OCH_COUNT]={0};
			__m256i mcounters[OCH_COUNT];//64-bit
			__m128i half8=_mm_set1_epi8(-128);
			__m256i wordmask=_mm256_set1_epi64x(0xFFFF);
			memset(mcounters, 0, sizeof(mcounters));
			imptr=interleaved+isize;
			for(int ky=ystart;ky<blockh;ky+=ANALYSIS_YSTRIDE)
			{
				__m256i prev[OCH_COUNT];//16-bit
				memset(prev, 0, sizeof(prev));
				for(int kx=0;kx<blockw-1;kx+=ANALYSIS_XSTRIDE)
				{
					__m256i r=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_load_si128((__m128i*)imptr+0), half8));
					__m256i g=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_load_si128((__m128i*)imptr+1), half8));
					__m256i b=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_load_si128((__m128i*)imptr+2), half8));
#ifdef ANALYSIS_GRAD
					__m256i rN=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_load_si128((__m128i*)(imptr-ixbytes)+0), half8));
					__m256i gN=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_load_si128((__m128i*)(imptr-ixbytes)+1), half8));
					__m256i bN=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_load_si128((__m128i*)(imptr-ixbytes)+2), half8));
					r=_mm256_sub_epi16(r, rN);
					g=_mm256_sub_epi16(g, gN);
					b=_mm256_sub_epi16(b, bN);
#endif
					imptr+=3*NCODERS*ANALYSIS_XSTRIDE;
					r=_mm256_slli_epi16(r, 2);
					g=_mm256_slli_epi16(g, 2);
					b=_mm256_slli_epi16(b, 2);
					__m256i rg=_mm256_sub_epi16(r, g);
					__m256i gb=_mm256_sub_epi16(g, b);
					__m256i br=_mm256_sub_epi16(b, r);
#define UPDATE(IDXA, A0, IDXB, B0, IDXC, C0)\
	do\
	{\
		__m256i sa=A0;\
		__m256i sb=B0;\
		__m256i sc=C0;\
		__m256i ta=_mm256_sub_epi16(sa, prev[IDXA]);\
		__m256i tb=_mm256_sub_epi16(sb, prev[IDXB]);\
		__m256i tc=_mm256_sub_epi16(sc, prev[IDXC]);\
		prev[IDXA]=sa;\
		prev[IDXB]=sb;\
		prev[IDXC]=sc;\
		ta=_mm256_abs_epi16(ta);\
		tb=_mm256_abs_epi16(tb);\
		tc=_mm256_abs_epi16(tc);\
		ta=_mm256_add_epi16(ta, _mm256_srli_epi64(ta, 32));\
		tb=_mm256_add_epi16(tb, _mm256_srli_epi64(tb, 32));\
		tc=_mm256_add_epi16(tc, _mm256_srli_epi64(tc, 32));\
		ta=_mm256_add_epi16(ta, _mm256_srli_epi64(ta, 16));\
		tb=_mm256_add_epi16(tb, _mm256_srli_epi64(tb, 16));\
		tc=_mm256_add_epi16(tc, _mm256_srli_epi64(tc, 16));\
		mcounters[IDXA]=_mm256_add_epi64(mcounters[IDXA], _mm256_and_si256(ta, wordmask));\
		mcounters[IDXB]=_mm256_add_epi64(mcounters[IDXB], _mm256_and_si256(tb, wordmask));\
		mcounters[IDXC]=_mm256_add_epi64(mcounters[IDXC], _mm256_and_si256(tc, wordmask));\
	}while(0)
					UPDATE(OCH_YX00, r, OCH_Y0X0, g, OCH_Y00X, b);
					UPDATE(OCH_CX40, rg, OCH_C0X4, gb, OCH_C40X, br);
#ifdef ENABLE_RCT_EXTENSION
					UPDATE(
						OCH_CX31, _mm256_add_epi16(rg, _mm256_srai_epi16(gb, 2)),//r-(3*g+b)/4 = r-g-(b-g)/4
						OCH_C3X1, _mm256_add_epi16(rg, _mm256_srai_epi16(br, 2)),//g-(3*r+b)/4 = g-r-(b-r)/4
						OCH_C31X, _mm256_add_epi16(br, _mm256_srai_epi16(rg, 2)) //b-(3*r+g)/4 = b-r-(g-r)/4
					);
					UPDATE(
						OCH_CX13, _mm256_add_epi16(br, _mm256_srai_epi16(gb, 2)),//r-(g+3*b)/4 = r-b-(g-b)/4
						OCH_C1X3, _mm256_add_epi16(gb, _mm256_srai_epi16(br, 2)),//g-(r+3*b)/4 = g-b-(r-b)/4
						OCH_C13X, _mm256_add_epi16(gb, _mm256_srai_epi16(rg, 2)) //b-(r+3*g)/4 = b-g-(r-g)/4
					);
					UPDATE(
						OCH_CX22,_mm256_srai_epi16(_mm256_sub_epi16(rg, br), 1),//r-(g+b)/2 = (r-g + r-b)/2
						OCH_C2X2,_mm256_srai_epi16(_mm256_sub_epi16(gb, rg), 1),//g-(r+b)/2 = (g-r + g-b)/2
						OCH_C22X,_mm256_srai_epi16(_mm256_sub_epi16(br, gb), 1) //b-(r+g)/2 = (b-r + b-g)/2
					);
#endif
				}
				imptr+=ixbytes*(ANALYSIS_YSTRIDE-1);
			}
			for(int k=0;k<OCH_COUNT;++k)
			{
				ALIGN(32) int64_t temp[4]={0};
				_mm256_store_si256((__m256i*)temp, mcounters[k]);
				counters[k]=temp[0]+temp[1]+temp[2]+temp[3];
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
				printf("%-14s %12lld + %12lld + %12lld = %12lld%s\n"
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

//			bestrct=0;//
//#ifdef __GNUC__
//#error remove above
//#endif
			//printf("%2d ", bestrct);
			prof_checkpoint(usize, "analysis");
		}

		//generate encode permutations		eg: mask = MSB 0b11000101 LSB  ->  LO {x, x, x, x, 0, 2, 6, 7} HI
		for(int km=0;km<256;++km)
		{
			int *curr=ans_permute+((ptrdiff_t)km<<3);
			int kb2=7;
			for(int kb=7;kb>=0;--kb)
			{
				int bit=km>>kb&1;
				if(bit)
					curr[kb2--]=kb;
			}
		}
		prof_checkpoint(ans_permute_size, "gen permutation");
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

		//generate decode permutations		eg: mask = MSB 0b11000101 LSB  ->  LO {0, x, 1, x, x, x, 2, 3} HI
		for(int km=0;km<256;++km)
		{
			int *curr=ans_permute+((ptrdiff_t)km<<3);
			int idx=0;
			for(int kb=0;kb<8;++kb)
			{
				int bit=km>>kb&1;
				if(bit)
					curr[kb]=idx++;
			}
		}
		prof_checkpoint(ans_permute_size, "gen permutation");
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
		L1statesize=(int)sizeof(int[2*NCODERS*3*(L1_NPREDS3+1)]);//{preds, coeffs} * (NPREDS+{bias}) * 3 channels * NCODERS
		L1state=(int*)_mm_malloc(L1statesize, sizeof(__m256i));
		if(!L1state)
		{
			CRASH("Alloc error");
			return 1;
		}
		memset(L1state, 0, L1statesize);
	}
#ifdef PRINT_L1_BOUNDS
	int cmin=0, cmax=0;
	int bmin=0, bmax=0;
#endif
	const uint8_t *combination=rct_combinations[bestrct];
	int
		yidx=combination[II_PERM_Y]*NCODERS,
		uidx=combination[II_PERM_U]*NCODERS,
		vidx=combination[II_PERM_V]*NCODERS;
	__m256i uhelpmask=_mm256_set1_epi16(-(combination[II_COEFF_U_SUB_Y]!=0));
	__m256i vc0=_mm256_set1_epi16(combination[II_COEFF_V_SUB_Y]);
	__m256i vc1=_mm256_set1_epi16(combination[II_COEFF_V_SUB_U]);
	//int paddedwidth=blockw+2*XPAD;
	memset(pixels, 0, psize);
	__m256i mctxmax=_mm256_set1_epi16(NCTX-1);
	__m256i mctxuoffset=_mm256_set1_epi16(NCTX);
	__m256i mctxvoffset=_mm256_set1_epi16(NCTX*2);
	__m256i amin=_mm256_set1_epi16(-128);
	__m256i amax=_mm256_set1_epi16(127);
	__m128i half8=_mm_set1_epi8(-128);
	__m256i bytemask=_mm256_set1_epi16(255);
	__m256i wordmask=_mm256_set1_epi32(0xFFFF);
	__m256i myuv[3];
	__m256i dist_rcp=_mm256_set1_epi16(0x7FFF), mdist=_mm256_set1_epi16(1);
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
		dist_rcp=_mm256_set1_epi16(((1<<16)+dist-1)/dist);//x/dist  ->  {x*=inv; x=(x>>16)+((uint32_t)x>>31);}
		mdist=_mm256_set1_epi16(dist);
	}
	memset(myuv, 0, sizeof(myuv));
	uint8_t *ctxptr=interleaved;
	imptr=interleaved+(fwd?isize:0);
	__m256i mstate[2];
	__m256i *L1preds=effort?(__m256i*)L1state:0;
	int *L1weights=effort?(int*)(L1state+1*(ptrdiff_t)NCODERS*3*(L1_NPREDS3+1)):0;
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
		ALIGN(32) int16_t *rows[]=
		{
			pixels+(XPAD*NROWS*NCH+(ky-0LL+NROWS)%NROWS)*NVAL*NCODERS,
			pixels+(XPAD*NROWS*NCH+(ky-1LL+NROWS)%NROWS)*NVAL*NCODERS,
			pixels+(XPAD*NROWS*NCH+(ky-2LL+NROWS)%NROWS)*NVAL*NCODERS,
			pixels+(XPAD*NROWS*NCH+(ky-3LL+NROWS)%NROWS)*NVAL*NCODERS,
			//pixels+(paddedwidth*((ky-0LL+NROWS)%NROWS)+(ptrdiff_t)XPAD)*NCH*NVAL*NCODERS,
			//pixels+(paddedwidth*((ky-1LL+NROWS)%NROWS)+(ptrdiff_t)XPAD)*NCH*NVAL*NCODERS,
			//pixels+(paddedwidth*((ky-2LL+NROWS)%NROWS)+(ptrdiff_t)XPAD)*NCH*NVAL*NCODERS,
			//pixels+(paddedwidth*((ky-3LL+NROWS)%NROWS)+(ptrdiff_t)XPAD)*NCH*NVAL*NCODERS,
		};
		ALIGN(32) uint16_t syms[3*NCODERS]={0};
		__m256i NW[3], N[3], W[3];
		__m256i eW[3], ecurr[3], eNEE[3], eNEEE[3];
		memset(NW, 0, sizeof(NW));
		memset(N, 0, sizeof(N));
		memset(W, 0, sizeof(W));
		memset(eW, 0, sizeof(eW));
		memset(ecurr, 0, sizeof(ecurr));
		memset(eNEE, 0, sizeof(eNEE));
		memset(eNEEE, 0, sizeof(eNEEE));
		//int16_t[blockw+2*XPAD][NCH*NROWS*NVAL*NCODERS]	(__m256i*)rows[1]+E+(C+X*NCH)*NROWS*NVAL
		eNEE[0]=_mm256_load_si256((__m256i*)rows[1]+1+(0+2*NCH)*NROWS*NVAL);
		eNEE[1]=_mm256_load_si256((__m256i*)rows[1]+1+(1+2*NCH)*NROWS*NVAL);
		eNEE[2]=_mm256_load_si256((__m256i*)rows[1]+1+(2+2*NCH)*NROWS*NVAL);
		for(int kx=0;kx<blockw;++kx)
		{
			N[0]=_mm256_load_si256((__m256i*)rows[1]+0+(0+0*NCH)*NROWS*NVAL);//y neighbors
			N[1]=_mm256_load_si256((__m256i*)rows[1]+0+(1+0*NCH)*NROWS*NVAL);//u
			N[2]=_mm256_load_si256((__m256i*)rows[1]+0+(2+0*NCH)*NROWS*NVAL);//v
			__m256i
				predY, ctxY,
				predU, ctxU,
				predV, ctxV;
			__m256i predYUV0[3];
			{
				//context = FLOOR_LOG2(eW*eW+1)
				__m256i one=_mm256_set1_epi32(1);
				__m256i cy0=_mm256_and_si256(eW[0], wordmask), cy1=_mm256_srli_epi32(eW[0], 16);
				__m256i cu0=_mm256_and_si256(eW[1], wordmask), cu1=_mm256_srli_epi32(eW[1], 16);
				__m256i cv0=_mm256_and_si256(eW[2], wordmask), cv1=_mm256_srli_epi32(eW[2], 16);
				cy0=_mm256_mullo_epi32(cy0, cy0);
				cy1=_mm256_mullo_epi32(cy1, cy1);
				cu0=_mm256_mullo_epi32(cu0, cu0);
				cu1=_mm256_mullo_epi32(cu1, cu1);
				cv0=_mm256_mullo_epi32(cv0, cv0);
				cv1=_mm256_mullo_epi32(cv1, cv1);
				cy0=_mm256_add_epi32(cy0, one);
				cy1=_mm256_add_epi32(cy1, one);
				cu0=_mm256_add_epi32(cu0, one);
				cu1=_mm256_add_epi32(cu1, one);
				cv0=_mm256_add_epi32(cv0, one);
				cv1=_mm256_add_epi32(cv1, one);
				//FLOOR_LOG2_32x8(X) = _mm256_sub_epi32(_mm256_srli_epi32(_mm256_castps_si256(_mm256_cvtepi32_ps(X)), 23), _mm256_set1_epi32(127))
				cy0=_mm256_castps_si256(_mm256_cvtepi32_ps(cy0));
				cy1=_mm256_castps_si256(_mm256_cvtepi32_ps(cy1));
				cu0=_mm256_castps_si256(_mm256_cvtepi32_ps(cu0));
				cu1=_mm256_castps_si256(_mm256_cvtepi32_ps(cu1));
				cv0=_mm256_castps_si256(_mm256_cvtepi32_ps(cv0));
				cv1=_mm256_castps_si256(_mm256_cvtepi32_ps(cv1));
				cy0=_mm256_srli_epi32(cy0, 23);
				cy1=_mm256_srli_epi32(cy1, 23);
				cu0=_mm256_srli_epi32(cu0, 23);
				cu1=_mm256_srli_epi32(cu1, 23);
				cv0=_mm256_srli_epi32(cv0, 23);
				cv1=_mm256_srli_epi32(cv1, 23);
				__m256i expbias=_mm256_set1_epi32(127);
				cy0=_mm256_sub_epi32(cy0, expbias);
				cy1=_mm256_sub_epi32(cy1, expbias);
				cu0=_mm256_sub_epi32(cu0, expbias);
				cu1=_mm256_sub_epi32(cu1, expbias);
				cv0=_mm256_sub_epi32(cv0, expbias);
				cv1=_mm256_sub_epi32(cv1, expbias);
				cy1=_mm256_slli_epi32(cy1, 16);
				cu1=_mm256_slli_epi32(cu1, 16);
				cv1=_mm256_slli_epi32(cv1, 16);
				ctxY=_mm256_or_si256(cy0, cy1);
				ctxU=_mm256_or_si256(cu0, cu1);
				ctxV=_mm256_or_si256(cv0, cv1);
				ctxY=_mm256_min_epi16(ctxY, mctxmax);
				ctxU=_mm256_min_epi16(ctxU, mctxmax);
				ctxV=_mm256_min_epi16(ctxV, mctxmax);
			}
			{
				const int borderW=3;
				const int borderN=3;
				const int borderE=3;
				int cond_cg=
					(uint32_t)(kx-borderW)>=(uint32_t)(blockw-(borderW+borderE))||
					(uint32_t)(ky-borderN)>=(uint32_t)(blockh-borderN);
				__m256i mcg[3];
				__m256i ymin=_mm256_min_epi16(N[0], W[0]);
				__m256i ymax=_mm256_max_epi16(N[0], W[0]);
				__m256i umin=_mm256_min_epi16(N[1], W[1]);
				__m256i umax=_mm256_max_epi16(N[1], W[1]);
				__m256i vmin=_mm256_min_epi16(N[2], W[2]);
				__m256i vmax=_mm256_max_epi16(N[2], W[2]);
				predY=_mm256_add_epi16(N[0], W[0]);//N+W-NW
				predU=_mm256_add_epi16(N[1], W[1]);
				predV=_mm256_add_epi16(N[2], W[2]);
				predY=_mm256_sub_epi16(predY, NW[0]);
				predU=_mm256_sub_epi16(predU, NW[1]);
				predV=_mm256_sub_epi16(predV, NW[2]);
				mcg[0]=predY;
				mcg[1]=predU;
				mcg[2]=predV;

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
					L1preds[0*3+0]=predY;
					L1preds[0*3+1]=predU;
					L1preds[0*3+2]=predV;

					//2*N-NN
					L1preds[1*3+0]=_mm256_sub_epi16(_mm256_add_epi16(N[0], N[0]), _mm256_load_si256((__m256i*)rows[2]+0+(0+0*NCH)*NROWS*NVAL));
					L1preds[1*3+1]=_mm256_sub_epi16(_mm256_add_epi16(N[1], N[1]), _mm256_load_si256((__m256i*)rows[2]+0+(1+0*NCH)*NROWS*NVAL));
					L1preds[1*3+2]=_mm256_sub_epi16(_mm256_add_epi16(N[2], N[2]), _mm256_load_si256((__m256i*)rows[2]+0+(2+0*NCH)*NROWS*NVAL));

				//	//N
				//	L1preds[1*3+0]=N[0];
				//	L1preds[1*3+1]=N[1];
				//	L1preds[1*3+2]=N[2];

					//W
					L1preds[2*3+0]=W[0];
					L1preds[2*3+1]=W[1];
					L1preds[2*3+2]=W[2];

					//NE
					L1preds[3*3+0]=_mm256_load_si256((__m256i*)rows[1]+0+(0+1*NCH)*NROWS*NVAL);
					L1preds[3*3+1]=_mm256_load_si256((__m256i*)rows[1]+0+(1+1*NCH)*NROWS*NVAL);
					L1preds[3*3+2]=_mm256_load_si256((__m256i*)rows[1]+0+(2+1*NCH)*NROWS*NVAL);
					

					//mix
					__m256i mp[6], t[6];
					mp[0]=_mm256_setzero_si256();
					mp[1]=_mm256_setzero_si256();
					mp[2]=_mm256_setzero_si256();
					mp[3]=_mm256_setzero_si256();
					mp[4]=_mm256_setzero_si256();
					mp[5]=_mm256_setzero_si256();
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 4
#endif
					for(int k=0;k<L1_NPREDS1;++k)
					{
						//16 -> 32		3 lo 3 hi registers
						t[0]=_mm256_slli_epi32(L1preds[k*3+0], 16);
						t[1]=_mm256_slli_epi32(L1preds[k*3+1], 16);
						t[2]=_mm256_slli_epi32(L1preds[k*3+2], 16);
						t[3]=_mm256_srai_epi32(L1preds[k*3+0], 16);
						t[4]=_mm256_srai_epi32(L1preds[k*3+1], 16);
						t[5]=_mm256_srai_epi32(L1preds[k*3+2], 16);
						t[0]=_mm256_srai_epi32(t[0], 16);
						t[1]=_mm256_srai_epi32(t[1], 16);
						t[2]=_mm256_srai_epi32(t[2], 16);
						t[0]=_mm256_mullo_epi32(t[0], _mm256_load_si256((__m256i*)L1weights+k*6+0));
						t[1]=_mm256_mullo_epi32(t[1], _mm256_load_si256((__m256i*)L1weights+k*6+1));
						t[2]=_mm256_mullo_epi32(t[2], _mm256_load_si256((__m256i*)L1weights+k*6+2));
						t[3]=_mm256_mullo_epi32(t[3], _mm256_load_si256((__m256i*)L1weights+k*6+3));
						t[4]=_mm256_mullo_epi32(t[4], _mm256_load_si256((__m256i*)L1weights+k*6+4));
						t[5]=_mm256_mullo_epi32(t[5], _mm256_load_si256((__m256i*)L1weights+k*6+5));
						mp[0]=_mm256_add_epi32(mp[0], t[0]);
						mp[1]=_mm256_add_epi32(mp[1], t[1]);
						mp[2]=_mm256_add_epi32(mp[2], t[2]);
						mp[3]=_mm256_add_epi32(mp[3], t[3]);
						mp[4]=_mm256_add_epi32(mp[4], t[4]);
						mp[5]=_mm256_add_epi32(mp[5], t[5]);
					}
					//__m256i rcon=_mm256_set1_epi32((1<<L1_SH1)-1);
					//mp[0]=_mm256_add_epi32(mp[0], _mm256_and_si256(rcon, _mm256_srli_epi32(mp[0], 31)));//rounding to zero
					//mp[1]=_mm256_add_epi32(mp[1], _mm256_and_si256(rcon, _mm256_srli_epi32(mp[1], 31)));
					//mp[2]=_mm256_add_epi32(mp[2], _mm256_and_si256(rcon, _mm256_srli_epi32(mp[2], 31)));
					//mp[3]=_mm256_add_epi32(mp[3], _mm256_and_si256(rcon, _mm256_srli_epi32(mp[3], 31)));
					//mp[4]=_mm256_add_epi32(mp[4], _mm256_and_si256(rcon, _mm256_srli_epi32(mp[4], 31)));
					//mp[5]=_mm256_add_epi32(mp[5], _mm256_and_si256(rcon, _mm256_srli_epi32(mp[5], 31)));

					__m256i rcon=_mm256_set1_epi32(1<<L1_SH1>>1);
					mp[0]=_mm256_add_epi32(mp[0], rcon);//rounding to nearest
					mp[1]=_mm256_add_epi32(mp[1], rcon);
					mp[2]=_mm256_add_epi32(mp[2], rcon);
					mp[3]=_mm256_add_epi32(mp[3], rcon);
					mp[4]=_mm256_add_epi32(mp[4], rcon);
					mp[5]=_mm256_add_epi32(mp[5], rcon);

					mp[0]=_mm256_srai_epi32(mp[0], L1_SH1);
					mp[1]=_mm256_srai_epi32(mp[1], L1_SH1);
					mp[2]=_mm256_srai_epi32(mp[2], L1_SH1);
					mp[3]=_mm256_slli_epi32(mp[3], 16-L1_SH1);
					mp[4]=_mm256_slli_epi32(mp[4], 16-L1_SH1);
					mp[5]=_mm256_slli_epi32(mp[5], 16-L1_SH1);
					//32 -> 16
					predY=_mm256_blend_epi16(mp[0], mp[3], 0xAA);
					predU=_mm256_blend_epi16(mp[1], mp[4], 0xAA);
					predV=_mm256_blend_epi16(mp[2], mp[5], 0xAA);


					//loosen pred range
					if(!cond_cg)
					{
						t[0]=_mm256_load_si256((__m256i*)rows[1]+0+(0+1*NCH)*NROWS*NVAL);//NE
						t[1]=_mm256_load_si256((__m256i*)rows[1]+0+(1+1*NCH)*NROWS*NVAL);
						t[2]=_mm256_load_si256((__m256i*)rows[1]+0+(2+1*NCH)*NROWS*NVAL);
						ymin=_mm256_min_epi16(ymin, t[0]);
						ymax=_mm256_max_epi16(ymax, t[0]);
						umin=_mm256_min_epi16(umin, t[1]);
						umax=_mm256_max_epi16(umax, t[1]);
						vmin=_mm256_min_epi16(vmin, t[2]);
						vmax=_mm256_max_epi16(vmax, t[2]);
						t[0]=_mm256_load_si256((__m256i*)rows[1]+0+(0+3*NCH)*NROWS*NVAL);//NEEE
						t[1]=_mm256_load_si256((__m256i*)rows[1]+0+(1+3*NCH)*NROWS*NVAL);
						t[2]=_mm256_load_si256((__m256i*)rows[1]+0+(2+3*NCH)*NROWS*NVAL);
						ymin=_mm256_min_epi16(ymin, t[0]);
						ymax=_mm256_max_epi16(ymax, t[0]);
						umin=_mm256_min_epi16(umin, t[1]);
						umax=_mm256_max_epi16(umax, t[1]);
						vmin=_mm256_min_epi16(vmin, t[2]);
						vmax=_mm256_max_epi16(vmax, t[2]);
					}
				}
				else if(effort==2)
				{
					__m256i cache[3];
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
					L1preds[0*3+0]=N[0];
					L1preds[0*3+1]=N[1];
					L1preds[0*3+2]=N[2];

					//W
					L1preds[1*3+0]=W[0];
					L1preds[1*3+1]=W[1];
					L1preds[1*3+2]=W[2];

					//3*(N-NN)+NNN
					cache[0]=_mm256_sub_epi16(N[0], _mm256_load_si256((__m256i*)rows[2]+0+(0+0*NCH)*NROWS*NVAL));//N-NN
					cache[1]=_mm256_sub_epi16(N[1], _mm256_load_si256((__m256i*)rows[2]+0+(1+0*NCH)*NROWS*NVAL));
					cache[2]=_mm256_sub_epi16(N[2], _mm256_load_si256((__m256i*)rows[2]+0+(2+0*NCH)*NROWS*NVAL));
					cache[0]=_mm256_add_epi16(cache[0], _mm256_slli_epi16(cache[0], 1));//*3
					cache[1]=_mm256_add_epi16(cache[1], _mm256_slli_epi16(cache[1], 1));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_slli_epi16(cache[2], 1));
					L1preds[2*3+0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[3]+0+(0+0*NCH)*NROWS*NVAL));//+NNN
					L1preds[2*3+1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[3]+0+(1+0*NCH)*NROWS*NVAL));
					L1preds[2*3+2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[3]+0+(2+0*NCH)*NROWS*NVAL));

					//3*(W-WW)+WWW
					cache[0]=_mm256_sub_epi16(W[0], _mm256_load_si256((__m256i*)rows[0]+0+(0-2*NCH)*NROWS*NVAL));//W-WW
					cache[1]=_mm256_sub_epi16(W[1], _mm256_load_si256((__m256i*)rows[0]+0+(1-2*NCH)*NROWS*NVAL));
					cache[2]=_mm256_sub_epi16(W[2], _mm256_load_si256((__m256i*)rows[0]+0+(2-2*NCH)*NROWS*NVAL));
					cache[0]=_mm256_add_epi16(cache[0], _mm256_slli_epi16(cache[0], 1));//*3
					cache[1]=_mm256_add_epi16(cache[1], _mm256_slli_epi16(cache[1], 1));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_slli_epi16(cache[2], 1));
					L1preds[3*3+0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[0]+0+(0-3*NCH)*NROWS*NVAL));//+WWW
					L1preds[3*3+1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[0]+0+(1-3*NCH)*NROWS*NVAL));
					L1preds[3*3+2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[0]+0+(2-3*NCH)*NROWS*NVAL));

					//W+NE-N
					cache[0]=_mm256_sub_epi16(W[0], N[0]);
					cache[1]=_mm256_sub_epi16(W[1], N[1]);
					cache[2]=_mm256_sub_epi16(W[2], N[2]);
					L1preds[4*3+0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[1]+0+(0+1*NCH)*NROWS*NVAL));//+NE
					L1preds[4*3+1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[1]+0+(1+1*NCH)*NROWS*NVAL));
					L1preds[4*3+2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[1]+0+(2+1*NCH)*NROWS*NVAL));

					//N+W-NW
					L1preds[5*3+0]=predY;
					L1preds[5*3+1]=predU;
					L1preds[5*3+2]=predV;

					//N+NE-NNE
					cache[0]=_mm256_add_epi16(N[0], _mm256_load_si256((__m256i*)rows[1]+0+(0+1*NCH)*NROWS*NVAL));//N+NE
					cache[1]=_mm256_add_epi16(N[1], _mm256_load_si256((__m256i*)rows[1]+0+(1+1*NCH)*NROWS*NVAL));
					cache[2]=_mm256_add_epi16(N[2], _mm256_load_si256((__m256i*)rows[1]+0+(2+1*NCH)*NROWS*NVAL));
					L1preds[6*3+0]=_mm256_sub_epi16(cache[0], _mm256_load_si256((__m256i*)rows[2]+0+(0+1*NCH)*NROWS*NVAL));//NNE
					L1preds[6*3+1]=_mm256_sub_epi16(cache[1], _mm256_load_si256((__m256i*)rows[2]+0+(1+1*NCH)*NROWS*NVAL));
					L1preds[6*3+2]=_mm256_sub_epi16(cache[2], _mm256_load_si256((__m256i*)rows[2]+0+(2+1*NCH)*NROWS*NVAL));

					//(WWWW+WWW+NNN+NEE+NEEE+NEEEE-(NW+N))>>2
					cache[0]=_mm256_load_si256((__m256i*)rows[0]+0+(0-4*NCH)*NROWS*NVAL);//WWWW
					cache[1]=_mm256_load_si256((__m256i*)rows[0]+0+(1-4*NCH)*NROWS*NVAL);
					cache[2]=_mm256_load_si256((__m256i*)rows[0]+0+(2-4*NCH)*NROWS*NVAL);
					cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[0]+0+(0-3*NCH)*NROWS*NVAL));//+WWW
					cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[0]+0+(1-3*NCH)*NROWS*NVAL));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[0]+0+(2-3*NCH)*NROWS*NVAL));
					cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[3]+0+(0+0*NCH)*NROWS*NVAL));//+NNN
					cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[3]+0+(1+0*NCH)*NROWS*NVAL));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[3]+0+(2+0*NCH)*NROWS*NVAL));
					cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[1]+0+(0+2*NCH)*NROWS*NVAL));//+NEE
					cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[1]+0+(1+2*NCH)*NROWS*NVAL));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[1]+0+(2+2*NCH)*NROWS*NVAL));
					cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[1]+0+(0+3*NCH)*NROWS*NVAL));//+NEEE
					cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[1]+0+(1+3*NCH)*NROWS*NVAL));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[1]+0+(2+3*NCH)*NROWS*NVAL));
					cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[1]+0+(0+4*NCH)*NROWS*NVAL));//+NEEEE
					cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[1]+0+(1+4*NCH)*NROWS*NVAL));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[1]+0+(2+4*NCH)*NROWS*NVAL));
					cache[0]=_mm256_sub_epi16(cache[0], _mm256_add_epi16(N[0], NW[0]));
					cache[1]=_mm256_sub_epi16(cache[1], _mm256_add_epi16(N[1], NW[1]));
					cache[2]=_mm256_sub_epi16(cache[2], _mm256_add_epi16(N[2], NW[2]));
					L1preds[7*3+0]=_mm256_srai_epi16(cache[0], 2);
					L1preds[7*3+1]=_mm256_srai_epi16(cache[1], 2);
					L1preds[7*3+2]=_mm256_srai_epi16(cache[2], 2);


					//mix
					__m256i mp[6], t[6];
					mp[0]=_mm256_setzero_si256();
					mp[1]=_mm256_setzero_si256();
					mp[2]=_mm256_setzero_si256();
					mp[3]=_mm256_setzero_si256();
					mp[4]=_mm256_setzero_si256();
					mp[5]=_mm256_setzero_si256();
#ifdef PRINT_L1_BOUNDS
					for(int k2=0;k2<6*8;++k2)
					{
						int bias=L1coeffs[L1_NPREDS3*6*8+k2];
						if(bmin>bias)bmin=bias;
						if(bmax<bias)bmax=bias;
					}
#endif
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 8
#endif
					for(int k=0;k<L1_NPREDS2;++k)
					{
						//16 -> 32		3 lo 3 hi registers
						t[0]=_mm256_slli_epi32(L1preds[k*3+0], 16);
						t[1]=_mm256_slli_epi32(L1preds[k*3+1], 16);
						t[2]=_mm256_slli_epi32(L1preds[k*3+2], 16);
						t[3]=_mm256_srai_epi32(L1preds[k*3+0], 16);
						t[4]=_mm256_srai_epi32(L1preds[k*3+1], 16);
						t[5]=_mm256_srai_epi32(L1preds[k*3+2], 16);
						t[0]=_mm256_srai_epi32(t[0], 16);
						t[1]=_mm256_srai_epi32(t[1], 16);
						t[2]=_mm256_srai_epi32(t[2], 16);
						t[0]=_mm256_mullo_epi32(t[0], _mm256_load_si256((__m256i*)L1weights+k*6+0));
						t[1]=_mm256_mullo_epi32(t[1], _mm256_load_si256((__m256i*)L1weights+k*6+1));
						t[2]=_mm256_mullo_epi32(t[2], _mm256_load_si256((__m256i*)L1weights+k*6+2));
						t[3]=_mm256_mullo_epi32(t[3], _mm256_load_si256((__m256i*)L1weights+k*6+3));
						t[4]=_mm256_mullo_epi32(t[4], _mm256_load_si256((__m256i*)L1weights+k*6+4));
						t[5]=_mm256_mullo_epi32(t[5], _mm256_load_si256((__m256i*)L1weights+k*6+5));
						mp[0]=_mm256_add_epi32(mp[0], t[0]);
						mp[1]=_mm256_add_epi32(mp[1], t[1]);
						mp[2]=_mm256_add_epi32(mp[2], t[2]);
						mp[3]=_mm256_add_epi32(mp[3], t[3]);
						mp[4]=_mm256_add_epi32(mp[4], t[4]);
						mp[5]=_mm256_add_epi32(mp[5], t[5]);
					}
					__m256i rcon=_mm256_set1_epi32(1<<L1_SH2>>1);
					mp[0]=_mm256_add_epi32(mp[0], rcon);//rounding to nearest
					mp[1]=_mm256_add_epi32(mp[1], rcon);
					mp[2]=_mm256_add_epi32(mp[2], rcon);
					mp[3]=_mm256_add_epi32(mp[3], rcon);
					mp[4]=_mm256_add_epi32(mp[4], rcon);
					mp[5]=_mm256_add_epi32(mp[5], rcon);

					mp[0]=_mm256_srai_epi32(mp[0], L1_SH2);
					mp[1]=_mm256_srai_epi32(mp[1], L1_SH2);
					mp[2]=_mm256_srai_epi32(mp[2], L1_SH2);
					mp[3]=_mm256_srai_epi32(mp[3], L1_SH2-16);
					mp[4]=_mm256_srai_epi32(mp[4], L1_SH2-16);
					mp[5]=_mm256_srai_epi32(mp[5], L1_SH2-16);
					//32 -> 16
					predY=_mm256_blend_epi16(mp[0], mp[3], 0xAA);
					predU=_mm256_blend_epi16(mp[1], mp[4], 0xAA);
					predV=_mm256_blend_epi16(mp[2], mp[5], 0xAA);


					//loosen pred range
					if(!cond_cg)
					{
						cache[0]=_mm256_load_si256((__m256i*)rows[1]+0+(0+1*NCH)*NROWS*NVAL);//NE
						cache[1]=_mm256_load_si256((__m256i*)rows[1]+0+(1+1*NCH)*NROWS*NVAL);
						cache[2]=_mm256_load_si256((__m256i*)rows[1]+0+(2+1*NCH)*NROWS*NVAL);
						ymin=_mm256_min_epi16(ymin, cache[0]);
						ymax=_mm256_max_epi16(ymax, cache[0]);
						umin=_mm256_min_epi16(umin, cache[1]);
						umax=_mm256_max_epi16(umax, cache[1]);
						vmin=_mm256_min_epi16(vmin, cache[2]);
						vmax=_mm256_max_epi16(vmax, cache[2]);
						cache[0]=_mm256_load_si256((__m256i*)rows[1]+0+(0+3*NCH)*NROWS*NVAL);//NEEE
						cache[1]=_mm256_load_si256((__m256i*)rows[1]+0+(1+3*NCH)*NROWS*NVAL);
						cache[2]=_mm256_load_si256((__m256i*)rows[1]+0+(2+3*NCH)*NROWS*NVAL);
						ymin=_mm256_min_epi16(ymin, cache[0]);
						ymax=_mm256_max_epi16(ymax, cache[0]);
						umin=_mm256_min_epi16(umin, cache[1]);
						umax=_mm256_max_epi16(umax, cache[1]);
						vmin=_mm256_min_epi16(vmin, cache[2]);
						vmax=_mm256_max_epi16(vmax, cache[2]);
					}
				}
				else if(effort==3)
				{
					__m256i cache[3];
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
					L1preds[0*3+0]=predY;
					L1preds[0*3+1]=predU;
					L1preds[0*3+2]=predV;

					//N
					L1preds[1*3+0]=N[0];
					L1preds[1*3+1]=N[1];
					L1preds[1*3+2]=N[2];

					//W
					L1preds[2*3+0]=W[0];
					L1preds[2*3+1]=W[1];
					L1preds[2*3+2]=W[2];

					//W+NE-N
					cache[0]=_mm256_add_epi16(W[0], _mm256_load_si256((__m256i*)rows[1]+0+(0+1*NCH)*NROWS*NVAL));
					cache[1]=_mm256_add_epi16(W[1], _mm256_load_si256((__m256i*)rows[1]+0+(1+1*NCH)*NROWS*NVAL));
					cache[2]=_mm256_add_epi16(W[2], _mm256_load_si256((__m256i*)rows[1]+0+(2+1*NCH)*NROWS*NVAL));
					L1preds[3*3+0]=_mm256_sub_epi16(cache[0], N[0]);
					L1preds[3*3+1]=_mm256_sub_epi16(cache[1], N[1]);
					L1preds[3*3+2]=_mm256_sub_epi16(cache[2], N[2]);

					//3*(N-NN)+NNN
					cache[0]=_mm256_sub_epi16(N[0], _mm256_load_si256((__m256i*)rows[2]+0+(0+0*NCH)*NROWS*NVAL));//N-NN
					cache[1]=_mm256_sub_epi16(N[1], _mm256_load_si256((__m256i*)rows[2]+0+(1+0*NCH)*NROWS*NVAL));
					cache[2]=_mm256_sub_epi16(N[2], _mm256_load_si256((__m256i*)rows[2]+0+(2+0*NCH)*NROWS*NVAL));
					cache[0]=_mm256_add_epi16(cache[0], _mm256_slli_epi16(cache[0], 1));//*3
					cache[1]=_mm256_add_epi16(cache[1], _mm256_slli_epi16(cache[1], 1));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_slli_epi16(cache[2], 1));
					L1preds[4*3+0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[3]+0+(0+0*NCH)*NROWS*NVAL));//+NNN
					L1preds[4*3+1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[3]+0+(1+0*NCH)*NROWS*NVAL));
					L1preds[4*3+2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[3]+0+(2+0*NCH)*NROWS*NVAL));

					//3*(W-WW)+WWW
					cache[0]=_mm256_sub_epi16(W[0], _mm256_load_si256((__m256i*)rows[0]+0+(0-2*NCH)*NROWS*NVAL));//W-WW
					cache[1]=_mm256_sub_epi16(W[1], _mm256_load_si256((__m256i*)rows[0]+0+(1-2*NCH)*NROWS*NVAL));
					cache[2]=_mm256_sub_epi16(W[2], _mm256_load_si256((__m256i*)rows[0]+0+(2-2*NCH)*NROWS*NVAL));
					cache[0]=_mm256_add_epi16(cache[0], _mm256_slli_epi16(cache[0], 1));//*3
					cache[1]=_mm256_add_epi16(cache[1], _mm256_slli_epi16(cache[1], 1));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_slli_epi16(cache[2], 1));
					L1preds[5*3+0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[0]+0+(0-3*NCH)*NROWS*NVAL));//+WWW
					L1preds[5*3+1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[0]+0+(1-3*NCH)*NROWS*NVAL));
					L1preds[5*3+2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[0]+0+(2-3*NCH)*NROWS*NVAL));

					//N+NE-NNE
					cache[0]=_mm256_add_epi16(N[0], _mm256_load_si256((__m256i*)rows[1]+0+(0+1*NCH)*NROWS*NVAL));//N+NE
					cache[1]=_mm256_add_epi16(N[1], _mm256_load_si256((__m256i*)rows[1]+0+(1+1*NCH)*NROWS*NVAL));
					cache[2]=_mm256_add_epi16(N[2], _mm256_load_si256((__m256i*)rows[1]+0+(2+1*NCH)*NROWS*NVAL));
					L1preds[6*3+0]=_mm256_sub_epi16(cache[0], _mm256_load_si256((__m256i*)rows[2]+0+(0+1*NCH)*NROWS*NVAL));//NNE
					L1preds[6*3+1]=_mm256_sub_epi16(cache[1], _mm256_load_si256((__m256i*)rows[2]+0+(1+1*NCH)*NROWS*NVAL));
					L1preds[6*3+2]=_mm256_sub_epi16(cache[2], _mm256_load_si256((__m256i*)rows[2]+0+(2+1*NCH)*NROWS*NVAL));

					//NEE
					L1preds[7*3+0]=_mm256_load_si256((__m256i*)rows[1]+0+(0+2*NCH)*NROWS*NVAL);
					L1preds[7*3+1]=_mm256_load_si256((__m256i*)rows[1]+0+(1+2*NCH)*NROWS*NVAL);
					L1preds[7*3+2]=_mm256_load_si256((__m256i*)rows[1]+0+(2+2*NCH)*NROWS*NVAL);

					//NN
					L1preds[8*3+0]=_mm256_load_si256((__m256i*)rows[2]+0+(0+0*NCH)*NROWS*NVAL);
					L1preds[8*3+1]=_mm256_load_si256((__m256i*)rows[2]+0+(1+0*NCH)*NROWS*NVAL);
					L1preds[8*3+2]=_mm256_load_si256((__m256i*)rows[2]+0+(2+0*NCH)*NROWS*NVAL);

					//WW
					L1preds[9*3+0]=_mm256_load_si256((__m256i*)rows[0]+0+(0-2*NCH)*NROWS*NVAL);
					L1preds[9*3+1]=_mm256_load_si256((__m256i*)rows[0]+0+(1-2*NCH)*NROWS*NVAL);
					L1preds[9*3+2]=_mm256_load_si256((__m256i*)rows[0]+0+(2-2*NCH)*NROWS*NVAL);

					//2*N-NN
					L1preds[10*3+0]=_mm256_sub_epi16(_mm256_slli_epi16(N[0], 1), _mm256_load_si256((__m256i*)rows[2]+0+(0+0*NCH)*NROWS*NVAL));
					L1preds[10*3+1]=_mm256_sub_epi16(_mm256_slli_epi16(N[1], 1), _mm256_load_si256((__m256i*)rows[2]+0+(1+0*NCH)*NROWS*NVAL));
					L1preds[10*3+2]=_mm256_sub_epi16(_mm256_slli_epi16(N[2], 1), _mm256_load_si256((__m256i*)rows[2]+0+(2+0*NCH)*NROWS*NVAL));

					//2*W-WW
					L1preds[11*3+0]=_mm256_sub_epi16(_mm256_slli_epi16(W[0], 1), _mm256_load_si256((__m256i*)rows[0]+0+(0-2*NCH)*NROWS*NVAL));
					L1preds[11*3+1]=_mm256_sub_epi16(_mm256_slli_epi16(W[1], 1), _mm256_load_si256((__m256i*)rows[0]+0+(1-2*NCH)*NROWS*NVAL));
					L1preds[11*3+2]=_mm256_sub_epi16(_mm256_slli_epi16(W[2], 1), _mm256_load_si256((__m256i*)rows[0]+0+(2-2*NCH)*NROWS*NVAL));

					//NEEE
					L1preds[12*3+0]=_mm256_load_si256((__m256i*)rows[1]+0+(0+3*NCH)*NROWS*NVAL);
					L1preds[12*3+1]=_mm256_load_si256((__m256i*)rows[1]+0+(1+3*NCH)*NROWS*NVAL);
					L1preds[12*3+2]=_mm256_load_si256((__m256i*)rows[1]+0+(2+3*NCH)*NROWS*NVAL);

					//NEEEE
					L1preds[13*3+0]=_mm256_load_si256((__m256i*)rows[1]+0+(0+4*NCH)*NROWS*NVAL);
					L1preds[13*3+1]=_mm256_load_si256((__m256i*)rows[1]+0+(1+4*NCH)*NROWS*NVAL);
					L1preds[13*3+2]=_mm256_load_si256((__m256i*)rows[1]+0+(2+4*NCH)*NROWS*NVAL);

					//NNWW
					L1preds[14*3+0]=_mm256_load_si256((__m256i*)rows[2]+0+(0-2*NCH)*NROWS*NVAL);
					L1preds[14*3+1]=_mm256_load_si256((__m256i*)rows[2]+0+(1-2*NCH)*NROWS*NVAL);
					L1preds[14*3+2]=_mm256_load_si256((__m256i*)rows[2]+0+(2-2*NCH)*NROWS*NVAL);

					//NNEE
					L1preds[15*3+0]=_mm256_load_si256((__m256i*)rows[2]+0+(0+2*NCH)*NROWS*NVAL);
					L1preds[15*3+1]=_mm256_load_si256((__m256i*)rows[2]+0+(1+2*NCH)*NROWS*NVAL);
					L1preds[15*3+2]=_mm256_load_si256((__m256i*)rows[2]+0+(2+2*NCH)*NROWS*NVAL);

					//N+NW-NNW
					cache[0]=_mm256_add_epi16(N[0], NW[0]);//N+NW
					cache[1]=_mm256_add_epi16(N[1], NW[1]);
					cache[2]=_mm256_add_epi16(N[2], NW[2]);
					L1preds[16*3+0]=_mm256_sub_epi16(cache[0], _mm256_load_si256((__m256i*)rows[2]+0+(0-1*NCH)*NROWS*NVAL));//NNW
					L1preds[16*3+1]=_mm256_sub_epi16(cache[1], _mm256_load_si256((__m256i*)rows[2]+0+(1-1*NCH)*NROWS*NVAL));
					L1preds[16*3+2]=_mm256_sub_epi16(cache[2], _mm256_load_si256((__m256i*)rows[2]+0+(2-1*NCH)*NROWS*NVAL));

					//W+NW-NWW
					cache[0]=_mm256_add_epi16(W[0], NW[0]);//W+NW
					cache[1]=_mm256_add_epi16(W[1], NW[1]);
					cache[2]=_mm256_add_epi16(W[2], NW[2]);
					L1preds[17*3+0]=_mm256_sub_epi16(cache[0], _mm256_load_si256((__m256i*)rows[1]+0+(0-2*NCH)*NROWS*NVAL));//NWW
					L1preds[17*3+1]=_mm256_sub_epi16(cache[1], _mm256_load_si256((__m256i*)rows[1]+0+(1-2*NCH)*NROWS*NVAL));
					L1preds[17*3+2]=_mm256_sub_epi16(cache[2], _mm256_load_si256((__m256i*)rows[1]+0+(2-2*NCH)*NROWS*NVAL));

					//(WWWW+NEEEE)>>1
					cache[0]=_mm256_add_epi16(_mm256_load_si256((__m256i*)rows[0]+0+(0-4*NCH)*NROWS*NVAL), _mm256_load_si256((__m256i*)rows[1]+0+(0+4*NCH)*NROWS*NVAL));//WWWW+NEEEE
					cache[1]=_mm256_add_epi16(_mm256_load_si256((__m256i*)rows[0]+0+(1-4*NCH)*NROWS*NVAL), _mm256_load_si256((__m256i*)rows[1]+0+(1+4*NCH)*NROWS*NVAL));
					cache[2]=_mm256_add_epi16(_mm256_load_si256((__m256i*)rows[0]+0+(2-4*NCH)*NROWS*NVAL), _mm256_load_si256((__m256i*)rows[1]+0+(2+4*NCH)*NROWS*NVAL));
					L1preds[18*3+0]=_mm256_srai_epi16(cache[0], 1);
					L1preds[18*3+1]=_mm256_srai_epi16(cache[1], 1);
					L1preds[18*3+2]=_mm256_srai_epi16(cache[2], 1);

					//(WWW+NNN+NEEE-NW)>>1
					cache[0]=_mm256_add_epi16(_mm256_load_si256((__m256i*)rows[0]+0+(0-3*NCH)*NROWS*NVAL), _mm256_load_si256((__m256i*)rows[3]+0+(0+0*NCH)*NROWS*NVAL));//WWW+NNN
					cache[1]=_mm256_add_epi16(_mm256_load_si256((__m256i*)rows[0]+0+(1-3*NCH)*NROWS*NVAL), _mm256_load_si256((__m256i*)rows[3]+0+(1+0*NCH)*NROWS*NVAL));
					cache[2]=_mm256_add_epi16(_mm256_load_si256((__m256i*)rows[0]+0+(2-3*NCH)*NROWS*NVAL), _mm256_load_si256((__m256i*)rows[3]+0+(2+0*NCH)*NROWS*NVAL));
					cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[1]+0+(0+3*NCH)*NROWS*NVAL));//+NEEE
					cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[1]+0+(1+3*NCH)*NROWS*NVAL));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[1]+0+(2+3*NCH)*NROWS*NVAL));
					cache[0]=_mm256_sub_epi16(cache[0], NW[0]);//-NW
					cache[1]=_mm256_sub_epi16(cache[1], NW[1]);
					cache[2]=_mm256_sub_epi16(cache[2], NW[2]);
					L1preds[19*3+0]=_mm256_srai_epi16(cache[0], 1);
					L1preds[19*3+1]=_mm256_srai_epi16(cache[1], 1);
					L1preds[19*3+2]=_mm256_srai_epi16(cache[2], 1);


					//mix
					__m256i mp[6], t[6];
					mp[0]=_mm256_setzero_si256();
					mp[1]=_mm256_setzero_si256();
					mp[2]=_mm256_setzero_si256();
					mp[3]=_mm256_setzero_si256();
					mp[4]=_mm256_setzero_si256();
					mp[5]=_mm256_setzero_si256();
#ifdef PRINT_L1_BOUNDS
					for(int k2=0;k2<6*8;++k2)
					{
						int bias=L1coeffs[L1_NPREDS3*6*8+k2];
						if(bmin>bias)bmin=bias;
						if(bmax<bias)bmax=bias;
					}
#endif
//#if defined __GNUC__ && !defined PROFILER
//#pragma GCC unroll 11
//#endif
					for(int k=0;k<L1_NPREDS3;++k)
					{
						//16 -> 32		3 lo 3 hi registers
						t[0]=_mm256_slli_epi32(L1preds[k*3+0], 16);
						t[1]=_mm256_slli_epi32(L1preds[k*3+1], 16);
						t[2]=_mm256_slli_epi32(L1preds[k*3+2], 16);
						t[3]=_mm256_srai_epi32(L1preds[k*3+0], 16);
						t[4]=_mm256_srai_epi32(L1preds[k*3+1], 16);
						t[5]=_mm256_srai_epi32(L1preds[k*3+2], 16);
						t[0]=_mm256_srai_epi32(t[0], 16);
						t[1]=_mm256_srai_epi32(t[1], 16);
						t[2]=_mm256_srai_epi32(t[2], 16);
						t[0]=_mm256_mullo_epi32(t[0], _mm256_load_si256((__m256i*)L1weights+k*6+0));
						t[1]=_mm256_mullo_epi32(t[1], _mm256_load_si256((__m256i*)L1weights+k*6+1));
						t[2]=_mm256_mullo_epi32(t[2], _mm256_load_si256((__m256i*)L1weights+k*6+2));
						t[3]=_mm256_mullo_epi32(t[3], _mm256_load_si256((__m256i*)L1weights+k*6+3));
						t[4]=_mm256_mullo_epi32(t[4], _mm256_load_si256((__m256i*)L1weights+k*6+4));
						t[5]=_mm256_mullo_epi32(t[5], _mm256_load_si256((__m256i*)L1weights+k*6+5));
						mp[0]=_mm256_add_epi32(mp[0], t[0]);
						mp[1]=_mm256_add_epi32(mp[1], t[1]);
						mp[2]=_mm256_add_epi32(mp[2], t[2]);
						mp[3]=_mm256_add_epi32(mp[3], t[3]);
						mp[4]=_mm256_add_epi32(mp[4], t[4]);
						mp[5]=_mm256_add_epi32(mp[5], t[5]);
					}
					__m256i rcon=_mm256_set1_epi32(1<<L1_SH3>>1);
					mp[0]=_mm256_add_epi32(mp[0], rcon);//rounding to nearest
					mp[1]=_mm256_add_epi32(mp[1], rcon);
					mp[2]=_mm256_add_epi32(mp[2], rcon);
					mp[3]=_mm256_add_epi32(mp[3], rcon);
					mp[4]=_mm256_add_epi32(mp[4], rcon);
					mp[5]=_mm256_add_epi32(mp[5], rcon);

					mp[0]=_mm256_srai_epi32(mp[0], L1_SH3);
					mp[1]=_mm256_srai_epi32(mp[1], L1_SH3);
					mp[2]=_mm256_srai_epi32(mp[2], L1_SH3);
					mp[3]=_mm256_srai_epi32(mp[3], L1_SH3-16);
					mp[4]=_mm256_srai_epi32(mp[4], L1_SH3-16);
					mp[5]=_mm256_srai_epi32(mp[5], L1_SH3-16);
					//32 -> 16
					predY=_mm256_blend_epi16(mp[0], mp[3], 0xAA);
					predU=_mm256_blend_epi16(mp[1], mp[4], 0xAA);
					predV=_mm256_blend_epi16(mp[2], mp[5], 0xAA);


					//loosen pred range
					if(!cond_cg)
					{
						cache[0]=_mm256_load_si256((__m256i*)rows[1]+0+(0+1*NCH)*NROWS*NVAL);//NE
						cache[1]=_mm256_load_si256((__m256i*)rows[1]+0+(1+1*NCH)*NROWS*NVAL);
						cache[2]=_mm256_load_si256((__m256i*)rows[1]+0+(2+1*NCH)*NROWS*NVAL);
						ymin=_mm256_min_epi16(ymin, cache[0]);
						ymax=_mm256_max_epi16(ymax, cache[0]);
						umin=_mm256_min_epi16(umin, cache[1]);
						umax=_mm256_max_epi16(umax, cache[1]);
						vmin=_mm256_min_epi16(vmin, cache[2]);
						vmax=_mm256_max_epi16(vmax, cache[2]);
						cache[0]=_mm256_load_si256((__m256i*)rows[1]+0+(0+3*NCH)*NROWS*NVAL);//NEEE	crisp DIV2K -0.18%, noisy GDCC +0.13%
						cache[1]=_mm256_load_si256((__m256i*)rows[1]+0+(1+3*NCH)*NROWS*NVAL);
						cache[2]=_mm256_load_si256((__m256i*)rows[1]+0+(2+3*NCH)*NROWS*NVAL);
						ymin=_mm256_min_epi16(ymin, cache[0]);
						ymax=_mm256_max_epi16(ymax, cache[0]);
						umin=_mm256_min_epi16(umin, cache[1]);
						umax=_mm256_max_epi16(umax, cache[1]);
						vmin=_mm256_min_epi16(vmin, cache[2]);
						vmax=_mm256_max_epi16(vmax, cache[2]);
					}
				}
				predYUV0[0]=predY;
				predYUV0[1]=predU;
				predYUV0[2]=predV;
				if(cond_cg)
				{
					predY=mcg[0];
					predU=mcg[1];
					predV=mcg[2];
				}

				predY=_mm256_max_epi16(predY, ymin);
				predU=_mm256_max_epi16(predU, umin);
				predV=_mm256_max_epi16(predV, vmin);
				predY=_mm256_min_epi16(predY, ymax);
				predU=_mm256_min_epi16(predU, umax);
				predV=_mm256_min_epi16(predV, vmax);
			}
			__m256i msyms, moffset;
			if(fwd)
			{
				__m256i ctxblendmask=_mm256_set1_epi16(255);
				myuv[0]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+yidx)), half8));//load yuv
				myuv[1]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+uidx)), half8));
				myuv[2]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+vidx)), half8));
				
				//decorrelate Y
				msyms=_mm256_sub_epi16(myuv[0], predY);//sub pred
				if(dist>1)
				{
#ifdef ENABLE_GUIDE
					__m256i yuv0=myuv[0];
#endif
					__m256i t0=_mm256_srai_epi16(msyms, 15);
					msyms=_mm256_mulhi_epi16(msyms, dist_rcp);
					msyms=_mm256_sub_epi16(msyms, t0);
					myuv[0]=_mm256_mullo_epi16(msyms, mdist);
					myuv[0]=_mm256_add_epi16(myuv[0], predY);
					myuv[0]=_mm256_max_epi16(myuv[0], amin);
					myuv[0]=_mm256_min_epi16(myuv[0], amax);
#ifdef ENABLE_GUIDE
					__m128i modified=_mm_packs_epi16(_mm256_extracti128_si256(myuv[0], 0), _mm256_extracti128_si256(myuv[0], 1));
					modified=_mm_xor_si128(modified, _mm_set1_epi8(-128));
					_mm_storeu_si128((__m128i*)(g_image+(imptr-interleaved-isize)+yidx), modified);

					__m256i mdiff=_mm256_sub_epi16(yuv0, myuv[0]);
					mdiff=_mm256_abs_epi16(mdiff);
					mdiff=_mm256_mullo_epi16(mdiff, mdiff);
					ALIGN(32) uint16_t diff[16];
					_mm256_store_si256((__m256i*)diff, mdiff);
					for(int k=0;k<16;++k)
						g_sqe[0]+=diff[k];
#endif
				}
				W[0]=myuv[0];
#ifdef SAVE_RESIDUALS
				{
					ptrdiff_t idx=imptr-interleaved-isize+(ptrdiff_t)0*NCODERS;
					int16_t syms2[16];
					_mm256_storeu_si256((__m256i*)syms2, msyms);
					for(int k=0;k<16;++k)
						residuals[idx+k]=(uint8_t)(syms2[k]+128);
				}
#endif
				ecurr[0]=_mm256_xor_si256(_mm256_slli_epi16(msyms, 1), _mm256_srai_epi16(msyms, 15));//ecurr = pack_sign(yuv-pred)
				msyms=_mm256_sub_epi16(msyms, amin);
				ctxY=_mm256_slli_epi16(ctxY, 8);
				ctxY=_mm256_blendv_epi8(ctxY, msyms, ctxblendmask);
				_mm256_store_si256((__m256i*)syms+0, ctxY);
				_mm256_store_si256((__m256i*)ctxptr+0, ctxY);//store Y  ctx|residuals
				
				//decorrelate U
				moffset=_mm256_and_si256(myuv[0], uhelpmask);
				predU=_mm256_add_epi16(predU, moffset);
				predU=_mm256_max_epi16(predU, amin);
				predU=_mm256_min_epi16(predU, amax);

				msyms=_mm256_sub_epi16(myuv[1], predU);
				if(dist>1)
				{
#ifdef ENABLE_GUIDE
					__m256i yuv0=myuv[1];
#endif
					__m256i t0=_mm256_srai_epi16(msyms, 15);
					msyms=_mm256_mulhi_epi16(msyms, dist_rcp);
					msyms=_mm256_sub_epi16(msyms, t0);
					myuv[1]=_mm256_mullo_epi16(msyms, mdist);
					myuv[1]=_mm256_add_epi16(myuv[1], predU);
					myuv[1]=_mm256_max_epi16(myuv[1], amin);
					myuv[1]=_mm256_min_epi16(myuv[1], amax);
#ifdef ENABLE_GUIDE
					__m128i modified=_mm_packs_epi16(_mm256_extracti128_si256(myuv[1], 0), _mm256_extracti128_si256(myuv[1], 1));
					modified=_mm_xor_si128(modified, _mm_set1_epi8(-128));
					_mm_storeu_si128((__m128i*)(g_image+(imptr-interleaved-isize)+uidx), modified);

					__m256i mdiff=_mm256_sub_epi16(yuv0, myuv[1]);
					mdiff=_mm256_abs_epi16(mdiff);
					mdiff=_mm256_mullo_epi16(mdiff, mdiff);
					ALIGN(32) uint16_t diff[16];
					_mm256_store_si256((__m256i*)diff, mdiff);
					for(int k=0;k<16;++k)
						g_sqe[1]+=diff[k];
#endif
				}
#ifdef SAVE_RESIDUALS
				{
					ptrdiff_t idx=imptr-interleaved-isize+(ptrdiff_t)1*NCODERS;
					int16_t syms2[16];
					_mm256_storeu_si256((__m256i*)syms2, msyms);
					for(int k=0;k<16;++k)
						residuals[idx+k]=(uint8_t)(syms2[k]+128);
				}
#endif
				ecurr[1]=_mm256_xor_si256(_mm256_slli_epi16(msyms, 1), _mm256_srai_epi16(msyms, 15));
				msyms=_mm256_sub_epi16(msyms, amin);
				ctxU=_mm256_add_epi16(ctxU, mctxuoffset);
				ctxU=_mm256_slli_epi16(ctxU, 8);
				ctxU=_mm256_blendv_epi8(ctxU, msyms, ctxblendmask);
				_mm256_store_si256((__m256i*)syms+1, ctxU);
				_mm256_store_si256((__m256i*)ctxptr+1, ctxU);//store U  ctx|residuals
				W[1]=_mm256_sub_epi16(myuv[1], moffset);
				
				//decorrelate V
				moffset=_mm256_mullo_epi16(vc0, myuv[0]);
				moffset=_mm256_add_epi16(moffset, _mm256_mullo_epi16(vc1, myuv[1]));
				moffset=_mm256_srai_epi16(moffset, 2);
				predV=_mm256_add_epi16(predV, moffset);
				predV=_mm256_max_epi16(predV, amin);
				predV=_mm256_min_epi16(predV, amax);

				msyms=_mm256_sub_epi16(myuv[2], predV);
				if(dist>1)
				{
#ifdef ENABLE_GUIDE
					__m256i yuv0=myuv[2];
#endif
					__m256i t0=_mm256_srai_epi16(msyms, 15);
					msyms=_mm256_mulhi_epi16(msyms, dist_rcp);
					msyms=_mm256_sub_epi16(msyms, t0);
					myuv[2]=_mm256_mullo_epi16(msyms, mdist);
					myuv[2]=_mm256_add_epi16(myuv[2], predV);
					myuv[2]=_mm256_max_epi16(myuv[2], amin);
					myuv[2]=_mm256_min_epi16(myuv[2], amax);
#ifdef ENABLE_GUIDE
					__m128i modified=_mm_packs_epi16(_mm256_extracti128_si256(myuv[2], 0), _mm256_extracti128_si256(myuv[2], 1));
					modified=_mm_xor_si128(modified, _mm_set1_epi8(-128));
					_mm_storeu_si128((__m128i*)(g_image+(imptr-interleaved-isize)+vidx), modified);

					__m256i mdiff=_mm256_sub_epi16(yuv0, myuv[2]);
					mdiff=_mm256_abs_epi16(mdiff);
					mdiff=_mm256_mullo_epi16(mdiff, mdiff);
					ALIGN(32) uint16_t diff[16];
					_mm256_store_si256((__m256i*)diff, mdiff);
					for(int k=0;k<16;++k)
						g_sqe[2]+=diff[k];
#endif
				}
#ifdef SAVE_RESIDUALS
				{
					ptrdiff_t idx=imptr-interleaved-isize+(ptrdiff_t)2*NCODERS;
					int16_t syms2[16];
					_mm256_storeu_si256((__m256i*)syms2, msyms);
					for(int k=0;k<16;++k)
						residuals[idx+k]=(uint8_t)(syms2[k]+128);
				}
#endif
				ecurr[2]=_mm256_xor_si256(_mm256_slli_epi16(msyms, 1), _mm256_srai_epi16(msyms, 15));
				msyms=_mm256_sub_epi16(msyms, amin);
				ctxV=_mm256_add_epi16(ctxV, mctxvoffset);
				ctxV=_mm256_slli_epi16(ctxV, 8);
				ctxV=_mm256_blendv_epi8(ctxV, msyms, ctxblendmask);
				_mm256_store_si256((__m256i*)syms+2, ctxV);
				_mm256_store_si256((__m256i*)ctxptr+2, ctxV);//store V  ctx|residuals		ctxptr+NCODERS*(C*2+R)
				W[2]=_mm256_sub_epi16(myuv[2], moffset);
#if 1
				int *pa, *pb, *pc, va, vb, vc;
				pa=hists+syms[0*16+0x0]; pb=hists+syms[1*16+0x0]; pc=hists+syms[2*16+0x0]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
				pa=hists+syms[0*16+0x1]; pb=hists+syms[1*16+0x1]; pc=hists+syms[2*16+0x1]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
				pa=hists+syms[0*16+0x2]; pb=hists+syms[1*16+0x2]; pc=hists+syms[2*16+0x2]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
				pa=hists+syms[0*16+0x3]; pb=hists+syms[1*16+0x3]; pc=hists+syms[2*16+0x3]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
				pa=hists+syms[0*16+0x4]; pb=hists+syms[1*16+0x4]; pc=hists+syms[2*16+0x4]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
				pa=hists+syms[0*16+0x5]; pb=hists+syms[1*16+0x5]; pc=hists+syms[2*16+0x5]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
				pa=hists+syms[0*16+0x6]; pb=hists+syms[1*16+0x6]; pc=hists+syms[2*16+0x6]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
				pa=hists+syms[0*16+0x7]; pb=hists+syms[1*16+0x7]; pc=hists+syms[2*16+0x7]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
				pa=hists+syms[0*16+0x8]; pb=hists+syms[1*16+0x8]; pc=hists+syms[2*16+0x8]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
				pa=hists+syms[0*16+0x9]; pb=hists+syms[1*16+0x9]; pc=hists+syms[2*16+0x9]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
				pa=hists+syms[0*16+0xA]; pb=hists+syms[1*16+0xA]; pc=hists+syms[2*16+0xA]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
				pa=hists+syms[0*16+0xB]; pb=hists+syms[1*16+0xB]; pc=hists+syms[2*16+0xB]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
				pa=hists+syms[0*16+0xC]; pb=hists+syms[1*16+0xC]; pc=hists+syms[2*16+0xC]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
				pa=hists+syms[0*16+0xD]; pb=hists+syms[1*16+0xD]; pc=hists+syms[2*16+0xD]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
				pa=hists+syms[0*16+0xE]; pb=hists+syms[1*16+0xE]; pc=hists+syms[2*16+0xE]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
				pa=hists+syms[0*16+0xF]; pb=hists+syms[1*16+0xF]; pc=hists+syms[2*16+0xF]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
#endif
#if 0
				++hists[syms[0*16+0x0]];
				++hists[syms[1*16+0x0]];
				++hists[syms[2*16+0x0]];
				++hists[syms[0*16+0x1]];
				++hists[syms[1*16+0x1]];
				++hists[syms[2*16+0x1]];
				++hists[syms[0*16+0x2]];
				++hists[syms[1*16+0x2]];
				++hists[syms[2*16+0x2]];
				++hists[syms[0*16+0x3]];
				++hists[syms[1*16+0x3]];
				++hists[syms[2*16+0x3]];
				++hists[syms[0*16+0x4]];
				++hists[syms[1*16+0x4]];
				++hists[syms[2*16+0x4]];
				++hists[syms[0*16+0x5]];
				++hists[syms[1*16+0x5]];
				++hists[syms[2*16+0x5]];
				++hists[syms[0*16+0x6]];
				++hists[syms[1*16+0x6]];
				++hists[syms[2*16+0x6]];
				++hists[syms[0*16+0x7]];
				++hists[syms[1*16+0x7]];
				++hists[syms[2*16+0x7]];
				++hists[syms[0*16+0x8]];
				++hists[syms[1*16+0x8]];
				++hists[syms[2*16+0x8]];
				++hists[syms[0*16+0x9]];
				++hists[syms[1*16+0x9]];
				++hists[syms[2*16+0x9]];
				++hists[syms[0*16+0xA]];
				++hists[syms[1*16+0xA]];
				++hists[syms[2*16+0xA]];
				++hists[syms[0*16+0xB]];
				++hists[syms[1*16+0xB]];
				++hists[syms[2*16+0xB]];
				++hists[syms[0*16+0xC]];
				++hists[syms[1*16+0xC]];
				++hists[syms[2*16+0xC]];
				++hists[syms[0*16+0xD]];
				++hists[syms[1*16+0xD]];
				++hists[syms[2*16+0xD]];
				++hists[syms[0*16+0xE]];
				++hists[syms[1*16+0xE]];
				++hists[syms[2*16+0xE]];
				++hists[syms[0*16+0xF]];
				++hists[syms[1*16+0xF]];
				++hists[syms[2*16+0xF]];
#endif
				ctxptr+=sizeof(int16_t[3][NCODERS]);
			}
			else
			{
				//decode main
				__m128i msyms8;
				
				//yuv = (char)(sym+pred-128)	= (uint8_t)(sym+pred)-128
				dec_yuv(mstate, &ctxY, CDF2syms+((ptrdiff_t)NCTX*0<<PROBBITS), ans_permute, &streamptr, streamend, myuv+0);//residuals from [0 ~ 255]
				dec_yuv(mstate, &ctxU, CDF2syms+((ptrdiff_t)NCTX*1<<PROBBITS), ans_permute, &streamptr, streamend, myuv+1);
				dec_yuv(mstate, &ctxV, CDF2syms+((ptrdiff_t)NCTX*2<<PROBBITS), ans_permute, &streamptr, streamend, myuv+2);

//#ifdef ANS_VAL
//				ALIGN(32) uint8_t debugvals[6*NCODERS];
//				msyms8=_mm256_packus_epi16(ctxY0, ctxY1);
//				msyms8=_mm256_permute4x64_epi64(msyms8, _MM_SHUFFLE(3, 1, 2, 0));
//				_mm256_store_si256((__m256i*)debugvals+0, msyms8);
//				msyms8=_mm256_packus_epi16(myuv[0], myuv[1]);
//				msyms8=_mm256_permute4x64_epi64(msyms8, _MM_SHUFFLE(3, 1, 2, 0));
//				_mm256_store_si256((__m256i*)debugvals+1, msyms8);
//				ansval_check(debugvals+0*NCODERS, 1, NCODERS);//contexts
//				ansval_check(debugvals+1*NCODERS, 1, NCODERS);//residuals
//#endif
				
				//reconstruct Y
				if(dist>1)
				{
					msyms=_mm256_add_epi16(myuv[0], amin);
					myuv[0]=_mm256_mullo_epi16(msyms, mdist);
					myuv[0]=_mm256_add_epi16(myuv[0], predY);
					myuv[0]=_mm256_max_epi16(myuv[0], amin);
					myuv[0]=_mm256_min_epi16(myuv[0], amax);
				}
				else
				{
					myuv[0]=_mm256_add_epi16(myuv[0], predY);
					myuv[0]=_mm256_and_si256(myuv[0], bytemask);
					myuv[0]=_mm256_add_epi16(myuv[0], amin);
					msyms=_mm256_sub_epi16(myuv[0], predY);//sub pred
				}
				ecurr[0]=_mm256_xor_si256(_mm256_slli_epi16(msyms, 1), _mm256_srai_epi16(msyms, 15));
				msyms8=_mm_packs_epi16(_mm256_extracti128_si256(myuv[0], 0), _mm256_extracti128_si256(myuv[0], 1));
				msyms8=_mm_xor_si128(msyms8, _mm_set1_epi8(-128));
				_mm_store_si128((__m128i*)(imptr+yidx), msyms8);//store Y bytes
#ifdef ENABLE_GUIDE
				if(memcmp(imptr+yidx, g_image+(imptr-interleaved)+yidx, NCODERS))
				{
					printf("original  decoded  original-decoded  XYC0 %d %d %d\n", kx, ky, yidx);
					for(int k=0;k<NCODERS;++k)
						printf("0x%02X  0x%02X  %4d\n",
							g_image[imptr-interleaved+yidx+k],
							imptr[yidx+k],
							g_image[imptr-interleaved+yidx+k]-imptr[yidx+k]
						);
					CRASH("guide error XYC0 %d %d %d/%d", kx, ky, yidx, NCODERS);
				}
#endif
//#ifdef ANS_VAL
//				ALIGN(32) uint8_t actx[3*NCODERS];
//				msyms8=_mm256_packus_epi16(ctxY0, ctxY1);
//				msyms8=_mm256_permute4x64_epi64(msyms8, _MM_SHUFFLE(3, 1, 2, 0));
//				_mm256_storeu_si256((__m256i*)actx+0, msyms8);
//				ansval_check(actx+0*NCODERS, 1, NCODERS);
//				ansval_check(imptr+kx+yidx, 1, NCODERS);
//#endif
				W[0]=myuv[0];


				//reconstruct U
				moffset=_mm256_and_si256(myuv[0], uhelpmask);
				predU=_mm256_add_epi16(predU, moffset);
				predU=_mm256_max_epi16(predU, amin);
				predU=_mm256_min_epi16(predU, amax);
//#ifdef ANS_VAL
//				msyms8=_mm256_packus_epi16(ctxU0, ctxU1);
//				msyms8=_mm256_permute4x64_epi64(msyms8, _MM_SHUFFLE(3, 1, 2, 0));
//				_mm256_store_si256((__m256i*)debugvals+2, msyms8);
//				msyms8=_mm256_packus_epi16(myuv[2], myuv[3]);
//				msyms8=_mm256_permute4x64_epi64(msyms8, _MM_SHUFFLE(3, 1, 2, 0));
//				_mm256_store_si256((__m256i*)debugvals+3, msyms8);
//				ansval_check(debugvals+2*NCODERS, 1, NCODERS);//contexts
//				ansval_check(debugvals+3*NCODERS, 1, NCODERS);//residuals
//#endif
				
				if(dist>1)
				{
					msyms=_mm256_add_epi16(myuv[1], amin);
					myuv[1]=_mm256_mullo_epi16(msyms, mdist);
					myuv[1]=_mm256_add_epi16(myuv[1], predU);
					myuv[1]=_mm256_max_epi16(myuv[1], amin);
					myuv[1]=_mm256_min_epi16(myuv[1], amax);
				}
				else
				{
					myuv[1]=_mm256_add_epi16(myuv[1], predU);
					myuv[1]=_mm256_and_si256(myuv[1], bytemask);
					myuv[1]=_mm256_add_epi16(myuv[1], amin);
					msyms=_mm256_sub_epi16(myuv[1], predU);//sub pred
				}
				ecurr[1]=_mm256_xor_si256(_mm256_slli_epi16(msyms, 1), _mm256_srai_epi16(msyms, 15));
				msyms8=_mm_packs_epi16(_mm256_extracti128_si256(myuv[1], 0), _mm256_extracti128_si256(myuv[1], 1));
				msyms8=_mm_xor_si128(msyms8, _mm_set1_epi8(-128));
				_mm_store_si128((__m128i*)(imptr+uidx), msyms8);//store U bytes
#ifdef ENABLE_GUIDE
				if(memcmp(imptr+uidx, g_image+(imptr-interleaved)+uidx, NCODERS))
				{
					printf("original  decoded  original-decoded  XYC1 %d %d %d\n", kx, ky, uidx);
					for(int k=0;k<NCODERS;++k)
						printf("0x%02X  0x%02X  %4d\n",
							g_image[imptr-interleaved+uidx+k],
							imptr[uidx+k],
							g_image[imptr-interleaved+uidx+k]-imptr[uidx+k]
						);
					CRASH("guide error XYC1 %d %d %d/%d", kx, ky, uidx, NCODERS);
				}
#endif
//#ifdef ANS_VAL
//				msyms8=_mm256_packus_epi16(ctxU0, ctxU1);
//				msyms8=_mm256_permute4x64_epi64(msyms8, _MM_SHUFFLE(3, 1, 2, 0));
//				_mm256_storeu_si256((__m256i*)actx+1, msyms8);
//				ansval_check(actx+1*NCODERS, 1, NCODERS);
//				ansval_check(imptr+kx+uidx, 1, NCODERS);
//#endif
				W[1]=_mm256_sub_epi16(myuv[1], moffset);//subtract Uoffset from U
				

				//reconstruct V
				moffset=_mm256_mullo_epi16(vc0, myuv[0]);
				moffset=_mm256_add_epi16(moffset, _mm256_mullo_epi16(vc1, myuv[1]));
				moffset=_mm256_srai_epi16(moffset, 2);
				predV=_mm256_add_epi16(predV, moffset);
				predV=_mm256_max_epi16(predV, amin);
				predV=_mm256_min_epi16(predV, amax);
				
				if(dist>1)
				{
					msyms=_mm256_add_epi16(myuv[2], amin);
					myuv[2]=_mm256_mullo_epi16(msyms, mdist);
					myuv[2]=_mm256_add_epi16(myuv[2], predV);
					myuv[2]=_mm256_max_epi16(myuv[2], amin);
					myuv[2]=_mm256_min_epi16(myuv[2], amax);
				}
				else
				{
					myuv[2]=_mm256_add_epi16(myuv[2], predV);
					myuv[2]=_mm256_and_si256(myuv[2], bytemask);
					myuv[2]=_mm256_add_epi16(myuv[2], amin);
					msyms=_mm256_sub_epi16(myuv[2], predV);
				}
				ecurr[2]=_mm256_xor_si256(_mm256_slli_epi16(msyms, 1), _mm256_srai_epi16(msyms, 15));
				msyms8=_mm_packs_epi16(_mm256_extracti128_si256(myuv[2], 0), _mm256_extracti128_si256(myuv[2], 1));
				msyms8=_mm_xor_si128(msyms8, _mm_set1_epi8(-128));
				_mm_store_si128((__m128i*)(imptr+vidx), msyms8);//store V bytes
#ifdef ENABLE_GUIDE
				if(memcmp(imptr+vidx, g_image+(imptr-interleaved)+vidx, NCODERS))
				{
					printf("original  decoded  original-decoded  XYC2 %d %d %d\n", kx, ky, vidx);
					for(int k=0;k<NCODERS;++k)
						printf("0x%02X  0x%02X  %4d\n",
							g_image[imptr-interleaved+vidx+k],
							imptr[vidx+k],
							g_image[imptr-interleaved+vidx+k]-imptr[vidx+k]
						);
					CRASH("guide error XYC2 %d %d %d/%d", kx, ky, vidx, NCODERS);
				}
#endif
				W[2]=_mm256_sub_epi16(myuv[2], moffset);//subtract Voffset from V
			}
			_mm256_store_si256((__m256i*)rows[0]+0+(0+0*NCH)*NROWS*NVAL, W[0]);//store neighbors
			_mm256_store_si256((__m256i*)rows[0]+0+(1+0*NCH)*NROWS*NVAL, W[1]);
			_mm256_store_si256((__m256i*)rows[0]+0+(2+0*NCH)*NROWS*NVAL, W[2]);
			if(effort==1)//update
			{
				__m256i mu[3];

				mu[0]=_mm256_sub_epi16(W[0], predYUV0[0]);
				mu[1]=_mm256_sub_epi16(W[1], predYUV0[1]);
				mu[2]=_mm256_sub_epi16(W[2], predYUV0[2]);
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 4
#endif
				for(int k=0;k<L1_NPREDS1;++k)//update
				{
					__m256i mc[6];
					mc[0]=_mm256_sign_epi16(L1preds[k*3+0], mu[0]);//L1
					mc[1]=_mm256_sign_epi16(L1preds[k*3+1], mu[1]);
					mc[2]=_mm256_sign_epi16(L1preds[k*3+2], mu[2]);
					//mc[0]=_mm256_mullo_epi16(L1preds[k*3+0], mu[0]);//L2
					//mc[1]=_mm256_mullo_epi16(L1preds[k*3+1], mu[1]);
					//mc[2]=_mm256_mullo_epi16(L1preds[k*3+2], mu[2]);

					//16 -> 32	3 lo 3 hi registers
					mc[3]=_mm256_srai_epi32(mc[0], 16);
					mc[4]=_mm256_srai_epi32(mc[1], 16);
					mc[5]=_mm256_srai_epi32(mc[2], 16);
					mc[0]=_mm256_slli_epi32(mc[0], 16);
					mc[1]=_mm256_slli_epi32(mc[1], 16);
					mc[2]=_mm256_slli_epi32(mc[2], 16);
					mc[0]=_mm256_srai_epi32(mc[0], 16);
					mc[1]=_mm256_srai_epi32(mc[1], 16);
					mc[2]=_mm256_srai_epi32(mc[2], 16);
					mc[0]=_mm256_add_epi32(mc[0], _mm256_load_si256((__m256i*)L1weights+k*6+0));//update coeffs
					mc[1]=_mm256_add_epi32(mc[1], _mm256_load_si256((__m256i*)L1weights+k*6+1));
					mc[2]=_mm256_add_epi32(mc[2], _mm256_load_si256((__m256i*)L1weights+k*6+2));
					mc[3]=_mm256_add_epi32(mc[3], _mm256_load_si256((__m256i*)L1weights+k*6+3));
					mc[4]=_mm256_add_epi32(mc[4], _mm256_load_si256((__m256i*)L1weights+k*6+4));
					mc[5]=_mm256_add_epi32(mc[5], _mm256_load_si256((__m256i*)L1weights+k*6+5));
					_mm256_store_si256((__m256i*)L1weights+k*6+0, mc[0]);
					_mm256_store_si256((__m256i*)L1weights+k*6+1, mc[1]);
					_mm256_store_si256((__m256i*)L1weights+k*6+2, mc[2]);
					_mm256_store_si256((__m256i*)L1weights+k*6+3, mc[3]);
					_mm256_store_si256((__m256i*)L1weights+k*6+4, mc[4]);
					_mm256_store_si256((__m256i*)L1weights+k*6+5, mc[5]);
				}
			}
			else if(effort==2)//update
			{
				__m256i mu[3];

				mu[0]=_mm256_sub_epi16(W[0], predYUV0[0]);
				mu[1]=_mm256_sub_epi16(W[1], predYUV0[1]);
				mu[2]=_mm256_sub_epi16(W[2], predYUV0[2]);
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 8
#endif
				for(int k=0;k<L1_NPREDS2;++k)//update
				{
					__m256i mc[6];
					mc[0]=_mm256_sign_epi16(L1preds[k*3+0], mu[0]);
					mc[1]=_mm256_sign_epi16(L1preds[k*3+1], mu[1]);
					mc[2]=_mm256_sign_epi16(L1preds[k*3+2], mu[2]);
					//16 -> 32	3 lo 3 hi registers
					mc[3]=_mm256_srai_epi32(mc[0], 16);
					mc[4]=_mm256_srai_epi32(mc[1], 16);
					mc[5]=_mm256_srai_epi32(mc[2], 16);
					mc[0]=_mm256_slli_epi32(mc[0], 16);
					mc[1]=_mm256_slli_epi32(mc[1], 16);
					mc[2]=_mm256_slli_epi32(mc[2], 16);
					mc[0]=_mm256_srai_epi32(mc[0], 16);
					mc[1]=_mm256_srai_epi32(mc[1], 16);
					mc[2]=_mm256_srai_epi32(mc[2], 16);
					mc[0]=_mm256_add_epi32(mc[0], _mm256_load_si256((__m256i*)L1weights+k*6+0));//update coeffs
					mc[1]=_mm256_add_epi32(mc[1], _mm256_load_si256((__m256i*)L1weights+k*6+1));
					mc[2]=_mm256_add_epi32(mc[2], _mm256_load_si256((__m256i*)L1weights+k*6+2));
					mc[3]=_mm256_add_epi32(mc[3], _mm256_load_si256((__m256i*)L1weights+k*6+3));
					mc[4]=_mm256_add_epi32(mc[4], _mm256_load_si256((__m256i*)L1weights+k*6+4));
					mc[5]=_mm256_add_epi32(mc[5], _mm256_load_si256((__m256i*)L1weights+k*6+5));
					_mm256_store_si256((__m256i*)L1weights+k*6+0, mc[0]);
					_mm256_store_si256((__m256i*)L1weights+k*6+1, mc[1]);
					_mm256_store_si256((__m256i*)L1weights+k*6+2, mc[2]);
					_mm256_store_si256((__m256i*)L1weights+k*6+3, mc[3]);
					_mm256_store_si256((__m256i*)L1weights+k*6+4, mc[4]);
					_mm256_store_si256((__m256i*)L1weights+k*6+5, mc[5]);
				}
			}
			else if(effort==3)//update
			{
				__m256i mu[3];

				mu[0]=_mm256_sub_epi16(W[0], predYUV0[0]);
				mu[1]=_mm256_sub_epi16(W[1], predYUV0[1]);
				mu[2]=_mm256_sub_epi16(W[2], predYUV0[2]);
//#if defined __GNUC__ && !defined PROFILER
//#pragma GCC unroll 11
//#endif
				for(int k=0;k<L1_NPREDS3;++k)//update
				{
					__m256i mc[6];
					mc[0]=_mm256_sign_epi16(L1preds[k*3+0], mu[0]);
					mc[1]=_mm256_sign_epi16(L1preds[k*3+1], mu[1]);
					mc[2]=_mm256_sign_epi16(L1preds[k*3+2], mu[2]);
					//16 -> 32	3 lo 3 hi registers
					mc[3]=_mm256_srai_epi32(mc[0], 16);
					mc[4]=_mm256_srai_epi32(mc[1], 16);
					mc[5]=_mm256_srai_epi32(mc[2], 16);
					mc[0]=_mm256_slli_epi32(mc[0], 16);
					mc[1]=_mm256_slli_epi32(mc[1], 16);
					mc[2]=_mm256_slli_epi32(mc[2], 16);
					mc[0]=_mm256_srai_epi32(mc[0], 16);
					mc[1]=_mm256_srai_epi32(mc[1], 16);
					mc[2]=_mm256_srai_epi32(mc[2], 16);
					mc[0]=_mm256_add_epi32(mc[0], _mm256_load_si256((__m256i*)L1weights+k*6+0));//update coeffs
					mc[1]=_mm256_add_epi32(mc[1], _mm256_load_si256((__m256i*)L1weights+k*6+1));
					mc[2]=_mm256_add_epi32(mc[2], _mm256_load_si256((__m256i*)L1weights+k*6+2));
					mc[3]=_mm256_add_epi32(mc[3], _mm256_load_si256((__m256i*)L1weights+k*6+3));
					mc[4]=_mm256_add_epi32(mc[4], _mm256_load_si256((__m256i*)L1weights+k*6+4));
					mc[5]=_mm256_add_epi32(mc[5], _mm256_load_si256((__m256i*)L1weights+k*6+5));
					_mm256_store_si256((__m256i*)L1weights+k*6+0, mc[0]);
					_mm256_store_si256((__m256i*)L1weights+k*6+1, mc[1]);
					_mm256_store_si256((__m256i*)L1weights+k*6+2, mc[2]);
					_mm256_store_si256((__m256i*)L1weights+k*6+3, mc[3]);
					_mm256_store_si256((__m256i*)L1weights+k*6+4, mc[4]);
					_mm256_store_si256((__m256i*)L1weights+k*6+5, mc[5]);
				}
			}
			//context update = (2*eW+(e<<3)+max(eNEE, eNEEE))>>2
			eNEEE[0]=_mm256_load_si256((__m256i*)rows[1]+1+(0+3*NCH)*NROWS*NVAL);
			eNEEE[1]=_mm256_load_si256((__m256i*)rows[1]+1+(1+3*NCH)*NROWS*NVAL);
			eNEEE[2]=_mm256_load_si256((__m256i*)rows[1]+1+(2+3*NCH)*NROWS*NVAL);
			
			ecurr[0]=_mm256_slli_epi16(ecurr[0], GRBITS);
			ecurr[1]=_mm256_slli_epi16(ecurr[1], GRBITS);
			ecurr[2]=_mm256_slli_epi16(ecurr[2], GRBITS);
			ecurr[0]=_mm256_avg_epu16(ecurr[0], _mm256_max_epi16(eNEE[0], eNEEE[0]));
			ecurr[1]=_mm256_avg_epu16(ecurr[1], _mm256_max_epi16(eNEE[1], eNEEE[1]));
			ecurr[2]=_mm256_avg_epu16(ecurr[2], _mm256_max_epi16(eNEE[2], eNEEE[2]));
			eW[0]=_mm256_avg_epu16(eW[0], ecurr[0]);
			eW[1]=_mm256_avg_epu16(eW[1], ecurr[1]);
			eW[2]=_mm256_avg_epu16(eW[2], ecurr[2]);

			_mm256_store_si256((__m256i*)rows[0]+1+(0+0*NCH)*NROWS*NVAL, eW[0]);//store current contexts
			_mm256_store_si256((__m256i*)rows[0]+1+(1+0*NCH)*NROWS*NVAL, eW[1]);
			_mm256_store_si256((__m256i*)rows[0]+1+(2+0*NCH)*NROWS*NVAL, eW[2]);
			eNEE[0]=eNEEE[0];
			eNEE[1]=eNEEE[1];
			eNEE[2]=eNEEE[2];
			NW[0]=N[0];
			NW[1]=N[1];
			NW[2]=N[2];
			rows[0]+=NCH*NROWS*NVAL*NCODERS;
			rows[1]+=NCH*NROWS*NVAL*NCODERS;
			rows[2]+=NCH*NROWS*NVAL*NCODERS;
			rows[3]+=NCH*NROWS*NVAL*NCODERS;
			imptr+=3*NCODERS;
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

#ifdef PRINT_L1_BOUNDS
	printf("Coeff %8d ~ %8d\n", cmin, cmax);
	printf("Bias  %8d ~ %8d\n", bmin, bmax);
#endif
	if(effort)
		_mm_free(L1state);
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
			enc_hist2stats(hists+(ptrdiff_t)256*kc, syminfo+(ptrdiff_t)256*kc, &bypassmask, kc, 0);
			
		if(xremw||yremh)//encode remainder
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
		mstate[1]=mstate[0]=_mm256_set1_epi32(1<<(RANS_STATE_BITS-RANS_RENORM_BITS));
		uint16_t *ctxptr2=(uint16_t*)(interleaved+(isize<<1));
		for(int ky=blockh-1;ky>=0;--ky)
		{
#ifdef ESTIMATE_SIZE
			int kc=2;
#endif
			for(int kx=3*blockw-1;kx>=0;--kx)//blockw = iw/XCODERS
			{
				__m256i mmax[2], minvf[2], mcdf[2], mnegf_sh[2];
				{
					__m256i s0, s1, s2, s3;
					__m256i t0, t1, t2, t3;
#define SHUFFLE_PS(LO, HI, IMM8_HHLL) _mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(LO), _mm256_castsi256_ps(HI), IMM8_HHLL))

					ctxptr2-=NCODERS;

					s0=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*8+0])));
					s1=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*8+1])));
					s2=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*8+2])));
					s3=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*8+3])));
					s0=_mm256_inserti128_si256(s0,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*8+4])), 1);
					s1=_mm256_inserti128_si256(s1,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*8+5])), 1);
					s2=_mm256_inserti128_si256(s2,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*8+6])), 1);
					s3=_mm256_inserti128_si256(s3,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*8+7])), 1);
					t0=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(1, 0, 1, 0));//_MM_TRANSPOSE4_PS
					t2=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(3, 2, 3, 2));
					t1=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(1, 0, 1, 0));
					t3=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(3, 2, 3, 2));
					mmax	[0]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(2, 0, 2, 0));
					minvf	[0]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(3, 1, 3, 1));
					mcdf	[0]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(2, 0, 2, 0));
					mnegf_sh[0]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(3, 1, 3, 1));

					s0=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*8+0])));
					s1=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*8+1])));
					s2=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*8+2])));
					s3=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*8+3])));
					s0=_mm256_inserti128_si256(s0,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*8+4])), 1);
					s1=_mm256_inserti128_si256(s1,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*8+5])), 1);
					s2=_mm256_inserti128_si256(s2,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*8+6])), 1);
					s3=_mm256_inserti128_si256(s3,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*8+7])), 1);
					t0=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(1, 0, 1, 0));
					t2=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(3, 2, 3, 2));
					t1=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(1, 0, 1, 0));
					t3=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(3, 2, 3, 2));
					mmax	[1]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(2, 0, 2, 0));
					minvf	[1]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(3, 1, 3, 1));
					mcdf	[1]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(2, 0, 2, 0));
					mnegf_sh[1]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(3, 1, 3, 1));
				}
#ifdef ESTIMATE_SIZE
				{
					ALIGN(32) int anegf[NCODERS]={0};
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
				
				if(!ky&&!kx)//
					printf("");

				//enc renorm		if(state>(freq<<(31-12))-1){write16(state); state>>=16;}
				__m256i cond0=_mm256_cmpgt_epi32(mstate[0], mmax[0]);
				__m256i cond1=_mm256_cmpgt_epi32(mstate[1], mmax[1]);
				int mask0=_mm256_movemask_ps(_mm256_castsi256_ps(cond0));
				int mask1=_mm256_movemask_ps(_mm256_castsi256_ps(cond1));
				__m256i idx0=_mm256_load_si256((__m256i*)ans_permute+mask0);
				__m256i idx1=_mm256_load_si256((__m256i*)ans_permute+mask1);
				__m256i emit0=_mm256_permutevar8x32_epi32(_mm256_and_si256(mstate[0], cond0), idx0);
				__m256i emit1=_mm256_permutevar8x32_epi32(_mm256_and_si256(mstate[1], cond1), idx1);
				emit0=_mm256_and_si256(emit0, _mm256_set1_epi32(0xFFFF));
				emit1=_mm256_and_si256(emit1, _mm256_set1_epi32(0xFFFF));
				emit0=_mm256_packus_epi32(emit0, emit1);
				emit0=_mm256_permute4x64_epi64(emit0, _MM_SHUFFLE(3, 1, 2, 0));
				__m128i e1=_mm256_extractf128_si256(emit0, 1);
				__m128i e0=_mm256_castsi256_si128(emit0);
				mask1=_mm_popcnt_u32(mask1);
				mask0=_mm_popcnt_u32(mask0);
#ifdef _DEBUG
				if(streamptr-(2*((ptrdiff_t)mask0+mask1)+sizeof(__m128i))<=image)
					CRASH("OOB ptr %016zX <= %016zX", streamptr, image);
#endif
				_mm_storeu_si128((__m128i*)streamptr-1, e1); streamptr-=mask1*sizeof(int16_t);
				_mm_storeu_si128((__m128i*)streamptr-1, e0); streamptr-=mask0*sizeof(int16_t);
				{
					__m256i state0=_mm256_srli_epi32(mstate[0], 16);
					__m256i state1=_mm256_srli_epi32(mstate[1], 16);
					mstate[0]=_mm256_blendv_epi8(mstate[0], state0, cond0);
					mstate[1]=_mm256_blendv_epi8(mstate[1], state1, cond1);
				}
#ifdef ANS_VAL
				ansval_push(mstate, sizeof(int), NCODERS);
#endif
				//enc update		state += (state*invf>>32>>sh)*negf+cdf		state = state/freq<<12|(cdf+state%freq)
				{
					__m256i lo0=_mm256_mul_epu32(mstate[0], minvf[0]);//q = mulhi32(state, invf)
					__m256i lo1=_mm256_mul_epu32(mstate[1], minvf[1]);
					__m256i hi0=_mm256_mul_epu32(_mm256_srli_epi64(mstate[0], 32), _mm256_srli_epi64(minvf[0], 32));
					__m256i hi1=_mm256_mul_epu32(_mm256_srli_epi64(mstate[1], 32), _mm256_srli_epi64(minvf[1], 32));
					minvf[0]=_mm256_blend_epi32(_mm256_srli_epi64(lo0, 32), hi0, 0xAA);
					minvf[1]=_mm256_blend_epi32(_mm256_srli_epi64(lo1, 32), hi1, 0xAA);
				}
				{
					__m256i sh0=_mm256_srli_epi32(mnegf_sh[0], 16);
					__m256i sh1=_mm256_srli_epi32(mnegf_sh[1], 16);
					minvf[0]=_mm256_srlv_epi32(minvf[0], sh0);
					minvf[1]=_mm256_srlv_epi32(minvf[1], sh1);
				}
				mstate[0]=_mm256_add_epi32(mstate[0], mcdf[0]);
				mstate[1]=_mm256_add_epi32(mstate[1], mcdf[1]);
				{
					__m256i lomask=_mm256_set1_epi32(0xFFFF);
					__m256i negf0=_mm256_and_si256(mnegf_sh[0], lomask);
					__m256i negf1=_mm256_and_si256(mnegf_sh[1], lomask);
					minvf[0]=_mm256_mullo_epi32(minvf[0], negf0);
					minvf[1]=_mm256_mullo_epi32(minvf[1], negf1);
#ifdef ANS_VAL
					__m256i mM=_mm256_set1_epi32(1<<PROBBITS);
					negf0=_mm256_sub_epi16(mM, negf0);
					negf1=_mm256_sub_epi16(mM, negf1);
					negf0=_mm256_packus_epi32(negf0, negf1);
					negf0=_mm256_permute4x64_epi64(negf0, _MM_SHUFFLE(3, 1, 2, 0));
					ALIGN(32) uint16_t freqs[NCODERS];
					_mm256_store_si256((__m256i*)freqs, negf0);
					ansval_push(freqs, sizeof(int16_t), NCODERS);
#endif
				}
				mstate[0]=_mm256_add_epi32(mstate[0], minvf[0]);
				mstate[1]=_mm256_add_epi32(mstate[1], minvf[1]);
#ifdef ANS_VAL
				ansval_push(mstate, sizeof(int), NCODERS);
#endif
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
			for(int k=0;k<3*NCODERS;++k)
			{
				etotal+=esize[k];
				printf("E %12.2lf\n", esize[k]);//
			}
			printf("Total estimate  %12.2lf bytes\n", etotal);
#endif
#ifdef LOUD
			printf("L1C AVX2 WH %d*%d  RCT %2d %s  effort %d  dist %3d  \"%s\"\n"
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
				const char fn[]="20250603_0628PM.PPM";
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
	_mm_free(ans_permute);
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
