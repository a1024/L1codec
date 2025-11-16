#if defined _MSC_VER && !defined _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<emmintrin.h>
#ifdef _MSC_VER
#include<intrin.h>
#else
#include<x86intrin.h>
#endif
static const char file[]=__FILE__;


	#define PROFILE_TIME		//should be on

#ifdef _DEBUG
	#define PROFILE_SIZE
	#define LOUD			//size & time

//	#define ESTIMATE_SIZE		//DEBUG		checks for zero frequency, visualizes context usage
	#define ENABLE_GUIDE		//DEBUG		checks interleaved pixels
//	#define ANS_VAL			//DEBUG

//	#define PRINT_SHIFTBOUNDS
//	#define TEST_INTERLEAVE
#endif

	#define ENABLE_RCT_EXTENSION
	#define INTERLEAVESIMD		//2.5x faster interleave


#define ANALYSIS_XSTRIDE 2
#define ANALYSIS_YSTRIDE 2

#define DEFAULT_EFFORT_LEVEL 2
#define L1_NPREDS1 4
#define L1_NPREDS2 8
#define L1_NPREDS3 11
#define L1_SH1 15
#define L1_SH2 17
#define L1_SH3 17

//3*17+3=54 contexts
#define GRBITS 3
#define NCTX 18		//18*3+3 = 57 total

#define XCODERS 4	//xrem 1~3 cols
#define YCODERS 4	//yrem 1~3 rows

#define NCODERS 16

#define PROBBITS 12	//12 bit max	James Bonfield's CDF2sym: {freq<<20 | bias<<8 | sym}

#define RANS_STATE_BITS 31
#define RANS_RENORM_BITS 16

#include"entropy.h"
#ifdef ENABLE_GUIDE
static int g_iw=0, g_ih=0;
static unsigned char *g_image=0;
static double g_sqe[3]={0};
static void guide_save(const unsigned char *image, int iw, int ih)
{
	int size=3*iw*ih;
	g_iw=iw;
	g_ih=ih;
	g_image=(unsigned char*)malloc(size);
	if(!g_image)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memcpy(g_image, image, size);
}
static void guide_check(const unsigned char *image, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx);
	if(memcmp(image+idx, g_image+idx, 3))
	{
		LOG_ERROR("");
		printf("");
	}
}
#else
#define guide_save(...)
#define guide_check(...)
#endif
#ifdef PROFILE_SIZE
#include<stdarg.h>
static void profile_size(const unsigned char *dstbwdptr, const char *msg, ...)
{
	static ptrdiff_t size=0;
	static const unsigned char *prev=0;
	if(prev)
	{
		ptrdiff_t diff=prev-dstbwdptr;
		size+=diff;
		printf("%10d (%+10d) bytes", (int)size, (int)diff);
		if(msg)
		{
			va_list args;
			va_start(args, msg);
			printf("  ");
			vprintf(msg, args);
			va_end(args);
		}
		printf("\n");
	}
	prev=dstbwdptr;
}
#else
#define profile_size(...)
#endif
#ifdef PROFILE_TIME
typedef struct _SpeedProfilerInfo
{
	double t;
	ptrdiff_t size;
	const char *msg;
} SpeedProfilerInfo;
#define PROF_CAP 128
static double prof_timestamp=0;
static SpeedProfilerInfo prof_data[PROF_CAP]={0};
static int prof_count=0;
static void prof_checkpoint(ptrdiff_t size, const char *msg)
{
	double t2=time_sec();
	if(prof_timestamp)
	{
		SpeedProfilerInfo *info=prof_data+prof_count++;
		if(prof_count>PROF_CAP)
		{
			LOG_ERROR("Profiler OOB");
			return;
		}
		info->t=t2-prof_timestamp;
		info->size=size;
		info->msg=msg;
		//double delta=t2-t;
		//printf("%16.12lf sec", delta);
		//if(size)
		//	printf(" %12.6lf MB/s %10d bytes", size/(delta*1024*1024), (int)size);
		//if(msg)
		//	printf(" %s", msg);
		//printf("\n");
	}
	prof_timestamp=t2;
}
static void prof_print(ptrdiff_t usize)
{
	static char buf[2048]={0};
	double timesum=0, tmax=0;
	int k;
	int prev=0;
	double csum=0;
	int colors[128]={0};
	const int scale=5;

	for(k=0;k<prof_count;++k)
	{
		double t=prof_data[k].t;
		timesum+=t;
		if(tmax<t)
			tmax=t;
	}
	srand((unsigned)__rdtsc());
	colorgen(colors, prof_count, 64, 300, 100);
	//colorgen0(colors, prof_count, 0xC0C0C0);
	printf("1 char = %d ms\n", scale);
	printf("|");
	for(k=0;k<prof_count;++k)
	{
		int curr, space, len;

		SpeedProfilerInfo *info=prof_data+k;
		csum+=info->t;
		curr=(int)(csum*1000/scale);//fixed scale
		space=curr-prev;
		len=0;
		if(info->msg)
			len=(int)strlen(info->msg);
		if(space>2047)//printf("%*s", HUGE, ""); CRASHES
			space=2047;
		if(info->msg&&space>=len)
		{
			int labelstart=(space-len)>>1;
			int labelend=labelstart+len;

			memset(buf, '-', labelstart);
			buf[labelstart]=0;
			colorprintf(colors[k], colors[k], buf);
			colorprintf(COLORPRINTF_TXT_DEFAULT, colors[k], "%s", info->msg);
			memset(buf, '-', (ptrdiff_t)space-labelend);
			buf[space-labelend]=0;
			colorprintf(colors[k], colors[k], buf);
		//	colorprintf(COLORPRINTF_TXT_DEFAULT, colors[k], "%*s%s%*s", labelstart, "", info->msg, space-labelend, "");
		}
		else
		{
			memset(buf, '-', space);
			buf[space]=0;
			colorprintf(colors[k], colors[k], buf);
		}
		//	colorprintf(COLORPRINTF_TXT_DEFAULT, colors[k], "%*s", space, "");
		printf("|");
		prev=curr;
	}
	printf("\n");
	for(k=0;k<prof_count;++k)
	{
		SpeedProfilerInfo *info=prof_data+k;
		//int nstars=(int)(info->t/tmax*64+0.5);
	//	colorprintf(COLORPRINTF_TXT_DEFAULT, colors[k], "  %c", k+'A');
		printf("%16.7lf ms %8.4lf%% ", info->t*1000, 100.*info->t/timesum);
	//	printf("  %c: %16.12lf sec %8.4lf%% ", k+'A', info->t, info->t/timesum);
		if(info->size)
			printf(" %16.6lf MB/s %10d bytes ", info->size/(info->t*1024*1024), (int)info->size);
		if(info->msg)
			colorprintf(COLORPRINTF_TXT_DEFAULT, colors[k], "%-20s", info->msg);
		else// if(nstars)
			colorprintf(COLORPRINTF_TXT_DEFAULT, colors[k], "%-20s", "");
		//for(int k2=0;k2<nstars;++k2)
		//	printf("*");
		printf("\n");
	}
	printf("\n");
	printf("%16.7lf ms %12.6lf MB/s Total\n", timesum*1000, usize/(timesum*1024*1024));
	printf("\n");
	prof_count=0;
	prof_timestamp=0;
}
#else
#define prof_checkpoint(...)
#define prof_print(...)
#endif

#ifndef ENABLE_RCT_EXTENSION
#define OCHLIST\
	OCH(Y400) OCH(Y040) OCH(Y004)\
	OCH(CX40) OCH(C0X4) OCH(C40X)
#endif
#ifdef ENABLE_RCT_EXTENSION
#define OCHLIST\
	OCH(Y400) OCH(Y040) OCH(Y004)\
	OCH(CX40) OCH(C0X4) OCH(C40X)\
	OCH(CX31) OCH(C3X1) OCH(C31X)\
	OCH(CX13) OCH(C1X3) OCH(C13X)\
	OCH(CX22) OCH(C2X2) OCH(C22X)
#if 0
#define OCHLIST\
	OCH(Y400) OCH(Y040) OCH(Y004)\
	OCH(Y310) OCH(Y031) OCH(Y103)\
	OCH(Y301) OCH(Y130) OCH(Y013)\
	OCH(Y211) OCH(Y121) OCH(Y112)\
	OCH(CX40) OCH(C0X4) OCH(C40X)\
	OCH(CX31) OCH(C3X1) OCH(C31X)\
	OCH(CX13) OCH(C1X3) OCH(C13X)\
	OCH(CX22) OCH(C2X2) OCH(C22X)
#endif
#endif
typedef enum _OCHIndex
{
#define OCH(X) OCH_##X,
	OCHLIST
#undef  OCH
	OCH_COUNT,

	OCH_R=OCH_Y400,
	OCH_G=OCH_Y040,
	OCH_B=OCH_Y004,
	OCH_C4X0=OCH_CX40,
	OCH_C04X=OCH_C0X4,
	OCH_CX04=OCH_C40X,
	OCH_BG=OCH_C04X,
	OCH_BR=OCH_C40X,
	OCH_RG=OCH_CX40,
	OCH_RB=OCH_CX04,
	OCH_GB=OCH_C0X4,
	OCH_GR=OCH_C4X0,
#ifdef ENABLE_RCT_EXTENSION
	OCH_R1=OCH_CX13,
	OCH_G1=OCH_C3X1,
	OCH_B1=OCH_C13X,
	OCH_R2=OCH_CX22,
	OCH_G2=OCH_C2X2,
	OCH_B2=OCH_C22X,
	OCH_R3=OCH_CX31,
	OCH_G3=OCH_C1X3,
	OCH_B3=OCH_C31X,
#endif
} OCHIndex;
static const char *och_names[]=
{
#define OCH(X) #X,
	OCHLIST
#undef  OCH
};
typedef enum _RCTInfoIdx
{
	II_OCH_Y,
	II_OCH_U,
	II_OCH_V,

	II_PERM_Y,
	II_PERM_U,
	II_PERM_V,

	II_COEFF_U_SUB_Y,
	II_COEFF_V_SUB_Y,
	II_COEFF_V_SUB_U,

//	II_COEFF_Y_SUB_U,
//	II_COEFF_Y_SUB_V,
//	II_COEFF_U_SUB_V_NBLI,
//	II_COEFF_V_SUB_U_NBLI,

	II_COUNT,
} RCTInfoIdx;
//YUV = RCT * RGB	watch out for permutation in last row
//luma: averaging	chroma: subtraction
//example: _400_40X_3X1 == [1 0 0; -1 0 1; -3/4 1 -1/4]
#ifndef ENABLE_RCT_EXTENSION
#define RCTLIST\
	RCT(_400_0X0_00X,	OCH_R,		OCH_G,		OCH_B,		0, 1, 2,	0,  0, 0)\
	RCT(_400_0X0_04X,	OCH_R,		OCH_G,		OCH_BG,		0, 1, 2,	0,  0, 4)\
	RCT(_400_0X0_40X,	OCH_R,		OCH_G,		OCH_BR,		0, 1, 2,	0,  4, 0)\
	RCT(_040_00X_X40,	OCH_G,		OCH_B,		OCH_RG,		1, 2, 0,	0,  4, 0)\
	RCT(_040_00X_X04,	OCH_G,		OCH_B,		OCH_RB,		1, 2, 0,	0,  0, 4)\
	RCT(_004_X00_4X0,	OCH_B,		OCH_R,		OCH_GR,		2, 0, 1,	0,  0, 4)\
	RCT(_004_X00_0X4,	OCH_B,		OCH_R,		OCH_GB,		2, 0, 1,	0,  4, 0)\
	RCT(_040_04X_X40,	OCH_G,		OCH_BG,		OCH_RG,		1, 2, 0,	4,  4, 0)\
	RCT(_040_04X_X04,	OCH_G,		OCH_BG,		OCH_RB,		1, 2, 0,	4,  0, 4)\
	RCT(_040_X40_40X,	OCH_G,		OCH_RG,		OCH_BR,		1, 0, 2,	4,  0, 4)\
	RCT(_004_X04_0X4,	OCH_B,		OCH_RB,		OCH_GB,		2, 0, 1,	4,  4, 0)\
	RCT(_004_X04_4X0,	OCH_B,		OCH_RB,		OCH_GR,		2, 0, 1,	4,  0, 4)\
	RCT(_004_0X4_X40,	OCH_B,		OCH_GB,		OCH_RG,		2, 1, 0,	4,  0, 4)\
	RCT(_400_4X0_40X,	OCH_R,		OCH_GR,		OCH_BR,		0, 1, 2,	4,  4, 0)\
	RCT(_400_4X0_04X,	OCH_R,		OCH_GR,		OCH_BG,		0, 1, 2,	4,  0, 4)\
	RCT(_400_40X_0X4,	OCH_R,		OCH_BR,		OCH_GB,		0, 2, 1,	4,  0, 4)
#endif
#ifdef ENABLE_RCT_EXTENSION
#define RCTLIST\
	RCT(_400_0X0_00X,	OCH_R,		OCH_G,		OCH_B,		0, 1, 2,	0,  0, 0)\
	RCT(_400_0X0_04X,	OCH_R,		OCH_G,		OCH_BG,		0, 1, 2,	0,  0, 4)\
	RCT(_400_0X0_40X,	OCH_R,		OCH_G,		OCH_BR,		0, 1, 2,	0,  4, 0)\
	RCT(_040_00X_X40,	OCH_G,		OCH_B,		OCH_RG,		1, 2, 0,	0,  4, 0)\
	RCT(_040_00X_X04,	OCH_G,		OCH_B,		OCH_RB,		1, 2, 0,	0,  0, 4)\
	RCT(_004_X00_4X0,	OCH_B,		OCH_R,		OCH_GR,		2, 0, 1,	0,  0, 4)\
	RCT(_004_X00_0X4,	OCH_B,		OCH_R,		OCH_GB,		2, 0, 1,	0,  4, 0)\
	RCT(_040_04X_X40,	OCH_G,		OCH_BG,		OCH_RG,		1, 2, 0,	4,  4, 0)\
	RCT(_040_04X_X04,	OCH_G,		OCH_BG,		OCH_RB,		1, 2, 0,	4,  0, 4)\
	RCT(_040_X40_40X,	OCH_G,		OCH_RG,		OCH_BR,		1, 0, 2,	4,  0, 4)\
	RCT(_004_X04_0X4,	OCH_B,		OCH_RB,		OCH_GB,		2, 0, 1,	4,  4, 0)\
	RCT(_004_X04_4X0,	OCH_B,		OCH_RB,		OCH_GR,		2, 0, 1,	4,  0, 4)\
	RCT(_004_0X4_X40,	OCH_B,		OCH_GB,		OCH_RG,		2, 1, 0,	4,  0, 4)\
	RCT(_400_4X0_40X,	OCH_R,		OCH_GR,		OCH_BR,		0, 1, 2,	4,  4, 0)\
	RCT(_400_4X0_04X,	OCH_R,		OCH_GR,		OCH_BG,		0, 1, 2,	4,  0, 4)\
	RCT(_400_40X_0X4,	OCH_R,		OCH_BR,		OCH_GB,		0, 2, 1,	4,  0, 4)\
	RCT(_400_0X0_13X,	OCH_R,		OCH_G,		OCH_B1,		0, 1, 2,	0,  1, 3)\
	RCT(_400_4X0_13X,	OCH_R,		OCH_GR,		OCH_B1,		0, 1, 2,	4,  1, 3)\
	RCT(_400_00X_3X1,	OCH_R,		OCH_B,		OCH_G1,		0, 2, 1,	0,  3, 1)\
	RCT(_400_40X_3X1,	OCH_R,		OCH_BR,		OCH_G1,		0, 2, 1,	4,  3, 1)\
	RCT(_040_00X_X13,	OCH_G,		OCH_B,		OCH_R1,		1, 2, 0,	0,  1, 3)\
	RCT(_040_04X_X13,	OCH_G,		OCH_BG,		OCH_R1,		1, 2, 0,	4,  1, 3)\
	RCT(_040_X40_13X,	OCH_G,		OCH_RG,		OCH_B1,		1, 0, 2,	4,  3, 1)\
	RCT(_004_X04_3X1,	OCH_B,		OCH_RB,		OCH_G1,		2, 0, 1,	4,  1, 3)\
	RCT(_004_04X_X13,	OCH_B,		OCH_GB,		OCH_R1,		2, 1, 0,	4,  3, 1)\
	RCT(_400_0X0_22X,	OCH_R,		OCH_G,		OCH_B2,		0, 1, 2,	0,  2, 2)\
	RCT(_400_4X0_22X,	OCH_R,		OCH_GR,		OCH_B2,		0, 1, 2,	4,  2, 2)\
	RCT(_400_00X_2X2,	OCH_R,		OCH_B,		OCH_G2,		0, 2, 1,	0,  2, 2)\
	RCT(_400_40X_2X2,	OCH_R,		OCH_BR,		OCH_G2,		0, 2, 1,	4,  2, 2)\
	RCT(_040_00X_X22,	OCH_G,		OCH_B,		OCH_R2,		1, 2, 0,	0,  2, 2)\
	RCT(_040_04X_X22,	OCH_G,		OCH_BG,		OCH_R2,		1, 2, 0,	4,  2, 2)\
	RCT(_040_X40_22X,	OCH_G,		OCH_RG,		OCH_B2,		1, 0, 2,	4,  2, 2)\
	RCT(_004_X04_2X2,	OCH_B,		OCH_RB,		OCH_G2,		2, 0, 1,	4,  2, 2)\
	RCT(_004_0X4_X22,	OCH_B,		OCH_GB,		OCH_R2,		2, 1, 0,	4,  2, 2)\
	RCT(_400_0X0_31X,	OCH_R,		OCH_G,		OCH_B3,		0, 1, 2,	0,  3, 1)\
	RCT(_400_4X0_31X,	OCH_R,		OCH_GR,		OCH_B3,		0, 1, 2,	4,  3, 1)\
	RCT(_400_00X_1X3,	OCH_R,		OCH_B,		OCH_G3,		0, 2, 1,	0,  1, 3)\
	RCT(_400_40X_1X3,	OCH_R,		OCH_BR,		OCH_G3,		0, 2, 1,	4,  1, 3)\
	RCT(_040_00X_X31,	OCH_G,		OCH_B,		OCH_R3,		1, 2, 0,	0,  3, 1)\
	RCT(_040_04X_X31,	OCH_G,		OCH_BG,		OCH_R3,		1, 2, 0,	4,  3, 1)\
	RCT(_040_X40_31X,	OCH_G,		OCH_RG,		OCH_B3,		1, 0, 2,	4,  1, 3)\
	RCT(_004_X04_1X3,	OCH_B,		OCH_RB,		OCH_G3,		2, 0, 1,	4,  3, 1)\
	RCT(_004_0X4_X31,	OCH_B,		OCH_GB,		OCH_R3,		2, 1, 0,	4,  1, 3)
#if 0
	RCT(_211_4X0_40X,	OCH_Y211,	OCH_C4X0,	OCH_C40X,	0, 1, 2,	4,  4, 0,	1, 1, 0, 0)\
	RCT(_211_4X0_31X,	OCH_Y211,	OCH_C4X0,	OCH_C31X,	0, 1, 2,	4,  4, 0,	1, 1, 0, 1)\
	RCT(_211_3X1_40X,	OCH_Y211,	OCH_C3X1,	OCH_C40X,	0, 1, 2,	4,  4, 0,	1, 1, 1, 0)\
	RCT(_310_4X0_40X,	OCH_Y310,	OCH_C4X0,	OCH_C40X,	0, 1, 2,	4,  4, 0,	1, 0, 0, 0)\
	RCT(_310_4X0_31X,	OCH_Y310,	OCH_C4X0,	OCH_C31X,	0, 1, 2,	4,  4, 0,	1, 0, 0, 1)\
	RCT(_310_3X1_40X,	OCH_Y310,	OCH_C3X1,	OCH_C40X,	0, 1, 2,	4,  4, 0,	1, 0, 1, 0)\
	RCT(_301_4X0_40X,	OCH_Y301,	OCH_C4X0,	OCH_C40X,	0, 1, 2,	4,  4, 0,	0, 1, 0, 0)\
	RCT(_301_4X0_31X,	OCH_Y301,	OCH_C4X0,	OCH_C31X,	0, 1, 2,	4,  4, 0,	0, 1, 0, 1)\
	RCT(_301_3X1_40X,	OCH_Y301,	OCH_C3X1,	OCH_C40X,	0, 1, 2,	4,  4, 0,	0, 1, 1, 0)\
	RCT(_121_04X_X40,	OCH_Y121,	OCH_C04X,	OCH_CX40,	1, 2, 0,	4,  4, 0,	1, 1, 0, 0)\
	RCT(_121_04X_X31,	OCH_Y121,	OCH_C04X,	OCH_CX31,	1, 2, 0,	4,  4, 0,	1, 1, 0, 1)\
	RCT(_121_13X_X40,	OCH_Y121,	OCH_C13X,	OCH_CX40,	1, 2, 0,	4,  4, 0,	1, 1, 1, 0)\
	RCT(_031_04X_X40,	OCH_Y031,	OCH_C04X,	OCH_CX40,	1, 2, 0,	4,  4, 0,	0, 1, 0, 0)\
	RCT(_031_04X_X31,	OCH_Y031,	OCH_C04X,	OCH_CX31,	1, 2, 0,	4,  4, 0,	0, 1, 0, 1)\
	RCT(_031_13X_X40,	OCH_Y031,	OCH_C13X,	OCH_CX40,	1, 2, 0,	4,  4, 0,	0, 1, 1, 0)\
	RCT(_130_40X_X40,	OCH_Y130,	OCH_C04X,	OCH_CX40,	1, 2, 0,	4,  4, 0,	1, 0, 0, 0)\
	RCT(_130_40X_X31,	OCH_Y130,	OCH_C04X,	OCH_CX31,	1, 2, 0,	4,  4, 0,	1, 0, 0, 1)\
	RCT(_130_31X_X40,	OCH_Y130,	OCH_C13X,	OCH_CX40,	1, 2, 0,	4,  4, 0,	1, 0, 1, 0)\
	RCT(_112_X04_0X4,	OCH_Y112,	OCH_CX04,	OCH_C0X4,	2, 0, 1,	4,  4, 0,	1, 1, 0, 0)\
	RCT(_112_X04_1X3,	OCH_Y112,	OCH_CX04,	OCH_C1X3,	2, 0, 1,	4,  4, 0,	1, 1, 0, 1)\
	RCT(_112_X13_0X4,	OCH_Y112,	OCH_CX13,	OCH_C0X4,	2, 0, 1,	4,  4, 0,	1, 1, 1, 0)\
	RCT(_013_X04_0X4,	OCH_Y013,	OCH_CX04,	OCH_C0X4,	2, 0, 1,	4,  4, 0,	0, 1, 0, 0)\
	RCT(_013_X04_1X3,	OCH_Y013,	OCH_CX04,	OCH_C1X3,	2, 0, 1,	4,  4, 0,	0, 1, 0, 1)\
	RCT(_013_X13_0X4,	OCH_Y013,	OCH_CX13,	OCH_C0X4,	2, 0, 1,	4,  4, 0,	0, 1, 1, 0)\
	RCT(_103_X40_0X4,	OCH_Y103,	OCH_CX04,	OCH_C0X4,	2, 0, 1,	4,  4, 0,	1, 0, 0, 0)\
	RCT(_103_X40_1X3,	OCH_Y103,	OCH_CX04,	OCH_C1X3,	2, 0, 1,	4,  4, 0,	1, 0, 0, 1)\
	RCT(_103_X31_0X4,	OCH_Y103,	OCH_CX13,	OCH_C0X4,	2, 0, 1,	4,  4, 0,	1, 0, 1, 0)
#endif
#endif
typedef enum _RCTIndex
{
#define RCT(LABEL, ...) RCT_##LABEL,
	RCTLIST
#undef  RCT
	RCT_COUNT,
} RCTIndex;
static const unsigned char rct_combinations[RCT_COUNT][II_COUNT]=
{
#define RCT(LABEL, ...) {__VA_ARGS__},
	RCTLIST
#undef  RCT
};
static const char *rct_names[RCT_COUNT]=
{
#define RCT(LABEL, ...) #LABEL,
	RCTLIST
#undef  RCT
};


AWM_INLINE void gather32(int *dst, const int *src, const int *offsets)
{
	dst[0]=src[offsets[0]];
	dst[1]=src[offsets[1]];
	dst[2]=src[offsets[2]];
	dst[3]=src[offsets[3]];
}

#ifdef PRINT_SHIFTBOUNDS
static int minsh=0x7FFFFFFF, maxsh=0;
#endif

//https://github.com/rygorous/ryg_rans
//https://github.com/samtools/htscodecs
typedef struct _rANS_SIMD_SymInfo	//16 bytes/level	4KB/ctx = 1<<12 bytes
{
	unsigned smax, invf, cdf;
	unsigned short negf, sh;
} rANS_SIMD_SymInfo;
static void enc_hist2stats(int *hist, rANS_SIMD_SymInfo *syminfo, unsigned long long *ctxmask, int ctxidx)
{
	int sum=0, count=0, ks, rare;
	for(ks=0;ks<256;++ks)
	{
		int freq=hist[ks];
		sum+=freq;
		count+=freq!=0;
	}
	rare=sum<12*256/8;
	*ctxmask|=(unsigned long long)rare<<ctxidx;
#ifdef ESTIMATE_SIZE
	int count0=count, sum0=sum;
#endif
	if(rare)
	{
		for(ks=0;ks<256;++ks)//bypass
			hist[ks]=1;
		sum=256;
		count=256;
	}
	else if(count==1)//disallow degenerate distribution
	{
		for(ks=0;ks<256;++ks)
		{
			int freq=hist[ks];
			if(freq==(1<<PROBBITS))
			{
				--freq;
				if(!ks)
					++hist[ks+1];
				else
					++hist[ks-1];
				break;
			}
		}
		count=2;
	}
	{
		int sum2, ks2;

		for(ks=0, ks2=0, sum2=0;ks<256;++ks)//absent symbols get zero freqs
		{
			int freq=hist[ks];
			hist[ks]=(int)(sum2*((1ULL<<PROBBITS)-count)/sum)+ks2;
			ks2+=freq!=0;
			sum2+=freq;
		}
		//for(ks=0, sum2=0;ks<256;++ks)//never allows zero freqs	INEFFICIENT
		//{
		//	int freq=hist[ks];
		//	hist[ks]=(int)(sum2*((1ULL<<PROBBITS)-256)/sum)+ks;
		//	sum2+=freq;
		//}
	}
#ifdef ESTIMATE_SIZE
	double e=sum0;
	if(count==count0)
	{
		double norm=1./0x1000;
		e=0;
		for(int ks=0;ks<256;++ks)//estimate
		{
			int freq=(ks<256-1?hist[ks+1]:1<<PROBBITS)-hist[ks];
			if(freq)
			{
				double p=freq*norm;
				e-=p*log2(p);
			}
			if(e!=e)
				LOG_ERROR("");
		}
		e*=sum/8.;
	}
	if(ctxidx&&!(ctxidx%NCTX))
		printf("\n");
	printf("%c  ctx %3d  %12.2lf / %9d bytes%8.2lf%%  %3d %s",
		ctxidx<3*NCTX?"YUV"[ctxidx/NCTX]:"yuv"[ctxidx-3*NCTX],
		ctxidx, e, sum0, 100.*e/sum0, count0, count==count0?"levels":"bypass"
	);
	if(count==count0&&count<256)
	{
		printf(" %3d", count);
		int fmax=0;
		for(int ks=0;ks<256;++ks)
		{
			int freq=(ks<256-1?hist[ks+1]:1<<PROBBITS)-hist[ks];
			if(fmax<freq)
				fmax=freq;
		}
		for(int ks=0;ks<256;++ks)
		{
			int freq=(ks<256-1?hist[ks+1]:1<<PROBBITS)-hist[ks];
			if(!(ks&15))
				printf(" ");

			int shade=48+freq*(255-48)/fmax;
			colorprintf(shade<<16|shade<<8|shade, freq?0x808080:COLORPRINTF_BK_DEFAULT, "%c", "0123456789ABCDEF"[ks&15]);
			//int shade=freq*255/fmax;
			//colorprintf(freq?0xFFFFFF:0x808080, shade<<16|0<<8|shade, "%c", "0123456789ABCDEF"[ks&15]);

			//printf("%c", freq?"0123456789ABCDEF"[ks&15]:'-');
		}
#if 0
		int printmissing=count>128, printcount=printmissing?256-count:count;
		if(printmissing)
			printf(" MISSING %3d: ", printcount);
		else
			printf("            : ");
		//printf(" %3d %-7s: ", printcount, printmissing?"MISSING":"       ");
		for(int ks=0, printed=0;ks<256;++ks)
		{
			int ks2=((ks>>1^-(ks&1))+128)&255;
			int freq=(ks2<256-1?hist[ks2+1]:1<<PROBBITS)-hist[ks2];
			if(printmissing!=(freq!=0))
			{
				printf(" %02X", ks2);
				//++printed;
				//if(printed&1)
				//	printf("%02x", ks2);
				//else
				//	printf("%02X", ks2);
				//if(printed>=90)
				//{
				//	printf("...%+4d more", printcount-printed);
				//	break;
				//}
			}
		}
#endif
	}
	printf("\n");
#if 0
	if(ctxidx==26)
	{
		const int amplitude=512;
		printf("Context %d: (1 star = %d steps)\n", ctxidx, 4096/amplitude);
		for(int ks=0;ks<256;++ks)
		{
			int freq=(ks<256-1?hist[ks+1]:1<<PROBBITS)-hist[ks], nstars=freq*amplitude>>PROBBITS;
			printf("%3d %4d ", ks, freq);
			for(int k=0;k<nstars;++k)
				printf("*");
			printf("\n");
		}
	}
#endif
#endif
	{
		int next=1<<PROBBITS;
		for(ks=255;ks>=0;--ks)
		{
			rANS_SIMD_SymInfo *info=syminfo+ks;
			int curr=hist[ks];
			int freq=next-curr;
			next=curr;
			hist[ks]=freq;
			info->smax=(freq<<(RANS_STATE_BITS-PROBBITS))-1;//rescale freq to match the rANS state, and decrement to use _mm_cmpgt_epi32 instead of '>='
			info->cdf=curr;
			info->negf=(1<<PROBBITS)-freq;
			//encoding:  state  =  q<<16|(cdf+r)
			//div-free:  state  =  q*M+cdf+state-q*freq  =  state+q*(M-freq)+cdf  =  state+(state*invf>>sh)*(M-freq)+cdf
			//sh = FLOOR_LOG2(freq)+32
			//invf = ceil(2^sh/freq)		state is 31 bits
			if(freq<2)
			{
				//freq=1
				//ideally  q = state*inv(1)>>sh(1) = state*2^32>>32
				//here  q' = state*(2^32-1)>>32 = floor(state-state/2^32) = state-1  if  1 <= x < 2^32
				//enc  state = (state/1)*M+cdf+state%1  =  state+q*(M-1)+cdf
				//but  q' = state-1
				//so  state = state+(state-1+1)*(M-1)+cdf  =  state+q'*(M-1)+(cdf+M-1)
				info->sh=0;
				info->invf=0xFFFFFFFF;
				info->cdf+=(1<<PROBBITS)-1;
			}
			else
			{
				unsigned long long inv;

				info->sh=FLOOR_LOG2(freq);//eg: x/2 = x*0x80000000>>32>>0
				inv=((0x100000000ULL<<info->sh)+freq-1)/freq;
				info->invf=(unsigned)inv;
				if(inv>0xFFFFFFFF)
				{
					--info->sh;
					info->invf=(unsigned)(inv>>1);
				}
			}
#ifdef PRINT_SHIFTBOUNDS
			if(minsh>info->sh)
				minsh=info->sh;
			if(maxsh<info->sh)
				maxsh=info->sh;
#endif
			info->sh=1<<(PROBBITS-1-info->sh);
		}
	}
}
static void enc_packhist(BitPackerLIFO *ec, const int *hist, unsigned long long ctxmask, int ctxidx)//histogram must be normalized to PROBBITS, with spike at 128
{
	unsigned short CDF[257];
	int ks;

	if(ctxmask>>ctxidx&1)
		return;
	{
		int sum=0;
		for(ks=0;ks<256;++ks)//integrage to zigzag CDF to be packed backwards
		{
			int sym=((ks>>1^-(ks&1))+128)&255;
			int freq=hist[sym];
			CDF[ks]=sum;//separate buffer for faster access in 2nd loop
			sum+=freq;
		}
		CDF[256]=1<<PROBBITS;
	}
	{
		int cdfW=CDF[0];
		int CDFlevels=1<<PROBBITS;
		int startsym=0;
		for(ks=1;ks<=256;++ks)//push GR.k
		{
			int next=CDF[ks], freq=next-cdfW;
			int nbypass=FLOOR_LOG2(CDFlevels);
			if(ks>1)
				nbypass-=7;
			if(nbypass<0)
				nbypass=0;
			CDF[ks]=nbypass<<PROBBITS|freq;
			cdfW=next;
			CDFlevels-=freq;
			startsym=ks;
			if(!CDFlevels)
				break;
		}
		for(ks=startsym;ks>0;--ks)//encode GR
		{
			int freq, nbypass, nzeros, bypass;

			freq=CDF[ks];
			nbypass=freq>>PROBBITS;
			freq&=(1<<PROBBITS)-1;
			nzeros=freq>>nbypass, bypass=freq&((1<<nbypass)-1);
			if(nbypass)
				bitpacker_enc(ec, nbypass, bypass);
			bitpacker_enc(ec, 1, 1);
			while(nzeros)
			{
				bitpacker_enc(ec, 1, 0);
				--nzeros;
			}
#ifdef ANS_VAL
			ansval_push(&ks, sizeof(ks), 1);
#endif
		}
	}
}
static void dec_unpackhist(BitPackerLIFO *ec, unsigned *CDF2sym, unsigned long long ctxmask, int ctxidx)
{
	unsigned short hist[257];
	int ks;

	if(ctxmask>>ctxidx&1)//rare context
	{
		for(ks=0;ks<256;++ks)//bypass
			hist[ks]=(1<<PROBBITS)/256;
	}
	else
	{
		unsigned short CDF[257]={0};
		int CDFlevels=1<<PROBBITS;
		CDF[0]=0;
		for(ks=0;ks<256;++ks)//decode GR
		{
			int freq, nbypass, ks2, bit;

			freq=-1;//stop bit doesn't count
			nbypass=FLOOR_LOG2(CDFlevels);
			ks2=ks+1;
			if(ks2>1)
				nbypass-=7;
			if(nbypass<0)
				nbypass=0;
#ifdef ANS_VAL
			ansval_check(&ks2, sizeof(ks2), 1);
#endif
			bit=0;
			do
			{
				bit=bitpacker_dec(ec, 1);
				++freq;
			}while(!bit);
			if(nbypass)
				freq=freq<<nbypass|bitpacker_dec(ec, nbypass);

			CDF[ks]=freq;
			CDFlevels-=freq;
			if(CDFlevels<=0)
			{
#ifdef _DEBUG
				if(CDFlevels<0)
					LOG_ERROR("CDF unpack error");
#endif
				break;
			}
		}
		if(CDFlevels)
			LOG_ERROR("CDF unpack error");
		for(ks=0;ks<256;++ks)//undo zigzag
		{
			int sym=((ks>>1^-(ks&1))+128)&255;
			hist[sym]=CDF[ks];
		}
	}
	{
		int sum=0;
		for(ks=0;ks<256;++ks)//integrate
		{
			int freq=hist[ks];
			hist[ks]=sum;
			sum+=freq;
		}
	}
	hist[256]=1<<PROBBITS;
	for(ks=0;ks<256;++ks)//CDF2sym contains {freq, (state&0xFFF)-cdf, sym}
	{
		int cdf, next, freq, val, ks2;

		cdf=hist[ks];
		next=hist[ks+1];
		freq=next-cdf;
		val=(freq<<PROBBITS|0)<<8|ks;
		for(ks2=cdf;ks2<next;++ks2, val+=1<<8)
			CDF2sym[ks2]=val;
	}
}

ALIGN(16) static const char ans_permute_enc[]=//pack		eg: mask = MSB 0b1010 LSB  ->  LO {x, x, 1, 3} HI
{
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,//0000
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  1,//0001
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  4,  5,//0010
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  1,  4,  5,//0011
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  8,  9,//0100
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  1,  8,  9,//0101
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  4,  5,  8,  9,//0110
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  1,  4,  5,  8,  9,//0111
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 12, 13,//1000
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  1, 12, 13,//1001
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  4,  5, 12, 13,//1010
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  1,  4,  5, 12, 13,//1011
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  8,  9, 12, 13,//1100
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  1,  8,  9, 12, 13,//1101
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  4,  5,  8,  9, 12, 13,//1110
	-1, -1, -1, -1, -1, -1, -1, -1,  0,  1,  4,  5,  8,  9, 12, 13,//1111
};
ALIGN(16) static const char ans_permute_dec[]=//unpack		eg: mask = MSB 0b1010 LSB  ->  LO {x, 0, x, 1} HI
{
	-1, -1, -1, -1,		-1, -1, -1, -1,		-1, -1, -1, -1,		-1, -1, -1, -1,//0000
	 0,  1, -1, -1,		-1, -1, -1, -1,		-1, -1, -1, -1,		-1, -1, -1, -1,//0001
	-1, -1, -1, -1,		 0,  1, -1, -1,		-1, -1, -1, -1,		-1, -1, -1, -1,//0010
	 0,  1, -1, -1,		 2,  3, -1, -1,		-1, -1, -1, -1,		-1, -1, -1, -1,//0011
	-1, -1, -1, -1,		-1, -1, -1, -1,		 0,  1, -1, -1,		-1, -1, -1, -1,//0100
	 0,  1, -1, -1,		-1, -1, -1, -1,		 2,  3, -1, -1,		-1, -1, -1, -1,//0101
	-1, -1, -1, -1,		 0,  1, -1, -1,		 2,  3, -1, -1,		-1, -1, -1, -1,//0110
	 0,  1, -1, -1,		 2,  3, -1, -1,		 4,  5, -1, -1,		-1, -1, -1, -1,//0111
	-1, -1, -1, -1,		-1, -1, -1, -1,		-1, -1, -1, -1,		 0,  1, -1, -1,//1000
	 0,  1, -1, -1,		-1, -1, -1, -1,		-1, -1, -1, -1,		 2,  3, -1, -1,//1001
	-1, -1, -1, -1,		 0,  1, -1, -1,		-1, -1, -1, -1,		 2,  3, -1, -1,//1010
	 0,  1, -1, -1,		 2,  3, -1, -1,		-1, -1, -1, -1,		 4,  5, -1, -1,//1011
	-1, -1, -1, -1,		-1, -1, -1, -1,		 0,  1, -1, -1,		 2,  3, -1, -1,//1100
	 0,  1, -1, -1,		-1, -1, -1, -1,		 2,  3, -1, -1,		 4,  5, -1, -1,//1101
	-1, -1, -1, -1,		 0,  1, -1, -1,		 2,  3, -1, -1,		 4,  5, -1, -1,//1110
	 0,  1, -1, -1,		 2,  3, -1, -1,		 4,  5, -1, -1,		 6,  7, -1, -1,//1111
};
static const unsigned char popcnt_table[]=
{
	0*2,//0000
	1*2,//0001
	1*2,//0010
	2*2,//0011
	1*2,//0100
	2*2,//0101
	2*2,//0110
	3*2,//0111
	1*2,//1000
	2*2,//1001
	2*2,//1010
	3*2,//1011
	2*2,//1100
	3*2,//1101
	3*2,//1110
	4*2,//1111
};
AWM_INLINE void dec_yuv(
	__m128i *mstate,
	const int kc,
	const __m128i *ctx0,
	const int *CDF2syms,
	unsigned char **pstreamptr,
	const unsigned char *streamend,
	__m128i *syms
)
{
	const unsigned char *streamptr=*pstreamptr;
	__m128i decctx[4];
	decctx[1]=_mm_cvtepi16_epi32(_mm_shuffle_epi32(ctx0[0], _MM_SHUFFLE(1, 0, 3, 2)));
	decctx[0]=_mm_cvtepi16_epi32(ctx0[0]);
	decctx[3]=_mm_cvtepi16_epi32(_mm_shuffle_epi32(ctx0[1], _MM_SHUFFLE(1, 0, 3, 2)));
	decctx[2]=_mm_cvtepi16_epi32(ctx0[1]);
	decctx[0]=_mm_slli_epi32(decctx[0], PROBBITS);
	decctx[1]=_mm_slli_epi32(decctx[1], PROBBITS);
	decctx[2]=_mm_slli_epi32(decctx[2], PROBBITS);
	decctx[3]=_mm_slli_epi32(decctx[3], PROBBITS);
#ifdef ANS_VAL
	ansval_check(mstate, sizeof(int), NCODERS);
#endif
	{
		__m128i mprobmask=_mm_set1_epi32((1<<PROBBITS)-1);
		__m128i rem0=_mm_and_si128(mstate[0], mprobmask);
		__m128i rem1=_mm_and_si128(mstate[1], mprobmask);
		__m128i rem2=_mm_and_si128(mstate[2], mprobmask);
		__m128i rem3=_mm_and_si128(mstate[3], mprobmask);
		decctx[0]=_mm_or_si128(decctx[0], rem0);
		decctx[1]=_mm_or_si128(decctx[1], rem1);
		decctx[2]=_mm_or_si128(decctx[2], rem2);
		decctx[3]=_mm_or_si128(decctx[3], rem3);
	}
#ifdef ANS_VAL
	ALIGN(32) int debugctx[NCODERS];
	memcpy(debugctx, decctx, sizeof(int[NCODERS]));
#endif
	{
		const int *statsptr=CDF2syms+((ptrdiff_t)NCTX*kc<<PROBBITS);
		gather32((int*)(decctx+0), statsptr, (int*)(decctx+0));
		gather32((int*)(decctx+1), statsptr, (int*)(decctx+1));
		gather32((int*)(decctx+2), statsptr, (int*)(decctx+2));
		gather32((int*)(decctx+3), statsptr, (int*)(decctx+3));
	}

	//update		state = (state>>12)*freq+(rem-cdf)	rem-cdf is prebaked
	{
		__m128i mfreq0=_mm_srli_epi32(decctx[0], PROBBITS+8);//1 <= freq <= 0xF01
		__m128i mfreq1=_mm_srli_epi32(decctx[1], PROBBITS+8);
		__m128i mfreq2=_mm_srli_epi32(decctx[2], PROBBITS+8);
		__m128i mfreq3=_mm_srli_epi32(decctx[3], PROBBITS+8);
#ifdef ANS_VAL
		{
			__m128i mdebugfreq[2];
			ALIGN(32) unsigned short freqs[NCODERS];
			mdebugfreq[0]=_mm_packus_epi32(mfreq0, mfreq1);
			mdebugfreq[1]=_mm_packus_epi32(mfreq2, mfreq3);
			memcpy(freqs, mdebugfreq, sizeof(freqs));
			ansval_check(freqs, sizeof(short), NCODERS);
		}
#endif
		mstate[0]=_mm_srli_epi32(mstate[0], PROBBITS);
		mstate[1]=_mm_srli_epi32(mstate[1], PROBBITS);
		mstate[2]=_mm_srli_epi32(mstate[2], PROBBITS);
		mstate[3]=_mm_srli_epi32(mstate[3], PROBBITS);
		mstate[0]=_mm_mullo_epi32(mstate[0], mfreq0);//10 cycles
		mstate[1]=_mm_mullo_epi32(mstate[1], mfreq1);
		mstate[2]=_mm_mullo_epi32(mstate[2], mfreq2);
		mstate[3]=_mm_mullo_epi32(mstate[3], mfreq3);
	}
	{
		__m128i mbias0=_mm_slli_epi32(decctx[0], PROBBITS);
		__m128i mbias1=_mm_slli_epi32(decctx[1], PROBBITS);
		__m128i mbias2=_mm_slli_epi32(decctx[2], PROBBITS);
		__m128i mbias3=_mm_slli_epi32(decctx[3], PROBBITS);
		mbias0=_mm_srli_epi32(mbias0, 32-PROBBITS);
		mbias1=_mm_srli_epi32(mbias1, 32-PROBBITS);
		mbias2=_mm_srli_epi32(mbias2, 32-PROBBITS);
		mbias3=_mm_srli_epi32(mbias3, 32-PROBBITS);
		mstate[0]=_mm_add_epi32(mstate[0], mbias0);
		mstate[1]=_mm_add_epi32(mstate[1], mbias1);
		mstate[2]=_mm_add_epi32(mstate[2], mbias2);
		mstate[3]=_mm_add_epi32(mstate[3], mbias3);
	}
	{
		__m128i symmask=_mm_set1_epi32(255);
		decctx[0]=_mm_and_si128(decctx[0], symmask);
		decctx[1]=_mm_and_si128(decctx[1], symmask);
		decctx[2]=_mm_and_si128(decctx[2], symmask);
		decctx[3]=_mm_and_si128(decctx[3], symmask);
	}
	syms[0]=_mm_packus_epi16(decctx[0], decctx[1]);
	syms[1]=_mm_packus_epi16(decctx[2], decctx[3]);
#ifdef ANS_VAL
	ansval_check(mstate, sizeof(int), NCODERS);
#endif
	//renorm	if(state<(1<<(31-16)))state=state<<16|read16();
	{
		__m128i
			lo0, cond0, idx0, renorm0,
			lo1, cond1, idx1, renorm1,
			lo2, cond2, idx2, renorm2,
			lo3, cond3, idx3, renorm3;
		int mask0, mask1, mask2, mask3;
		{
			__m128i smin=_mm_set1_epi32(1<<(RANS_STATE_BITS-RANS_RENORM_BITS));
			cond0=_mm_cmpgt_epi32(smin, mstate[0]);//signed comparison!
			cond1=_mm_cmpgt_epi32(smin, mstate[1]);
			cond2=_mm_cmpgt_epi32(smin, mstate[2]);
			cond3=_mm_cmpgt_epi32(smin, mstate[3]);
		}
		mask0=_mm_movemask_ps(_mm_castsi128_ps(cond0));
		mask1=_mm_movemask_ps(_mm_castsi128_ps(cond1));
		mask2=_mm_movemask_ps(_mm_castsi128_ps(cond2));
		mask3=_mm_movemask_ps(_mm_castsi128_ps(cond3));
		idx0=_mm_load_si128((__m128i*)ans_permute_dec+mask0);
		idx1=_mm_load_si128((__m128i*)ans_permute_dec+mask1);
		idx2=_mm_load_si128((__m128i*)ans_permute_dec+mask2);
		idx3=_mm_load_si128((__m128i*)ans_permute_dec+mask3);
		mask0=popcnt_table[mask0];
		mask1=popcnt_table[mask1];
		mask2=popcnt_table[mask2];
		mask3=popcnt_table[mask3];
#ifdef _DEBUG
		if(streamptr+mask0+mask1+mask2+mask3>streamend)
			LOG_ERROR("OOB ptr %016zX >= %016zX", (ptrdiff_t)streamptr+mask0+mask1+mask2+mask3, (ptrdiff_t)streamend);
#endif
		lo0=_mm_loadu_si128((__m128i*)streamptr); streamptr+=mask0;
		lo1=_mm_loadu_si128((__m128i*)streamptr); streamptr+=mask1;
		lo2=_mm_loadu_si128((__m128i*)streamptr); streamptr+=mask2;
		lo3=_mm_loadu_si128((__m128i*)streamptr); streamptr+=mask3;
		lo0=_mm_shuffle_epi8(lo0, idx0);
		lo1=_mm_shuffle_epi8(lo1, idx1);
		lo2=_mm_shuffle_epi8(lo2, idx2);
		lo3=_mm_shuffle_epi8(lo3, idx3);
		renorm0=_mm_slli_epi32(mstate[0], 16);
		renorm1=_mm_slli_epi32(mstate[1], 16);
		renorm2=_mm_slli_epi32(mstate[2], 16);
		renorm3=_mm_slli_epi32(mstate[3], 16);
		renorm0=_mm_or_si128(renorm0, lo0);
		renorm1=_mm_or_si128(renorm1, lo1);
		renorm2=_mm_or_si128(renorm2, lo2);
		renorm3=_mm_or_si128(renorm3, lo3);
		mstate[0]=_mm_blendv_epi8(mstate[0], renorm0, cond0);
		mstate[1]=_mm_blendv_epi8(mstate[1], renorm1, cond1);
		mstate[2]=_mm_blendv_epi8(mstate[2], renorm2, cond2);
		mstate[3]=_mm_blendv_epi8(mstate[3], renorm3, cond3);
	}
	*pstreamptr=(unsigned char*)(size_t)streamptr;
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
static void interleave_blocks_fwd(const unsigned char *original, int iw, int ih, unsigned char *interleaved)
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
	__m128i slowinc=sizeof(void*)==8?_mm_set_epi32(0, sizeof(__m128i), 0, sizeof(__m128i)):_mm_set1_epi32(sizeof(__m128i));
#endif
	unsigned char *fastptr=interleaved;
	ALIGN(32) const unsigned char *slowptrs[NCODERS]={0}, *slowptrs0[NCODERS]={0};
	int kx, ky;

	for(ky=0;ky<YCODERS;++ky)//spread slow pointers
	{
		for(kx=0;kx<XCODERS;++kx)
			slowptrs0[XCODERS*ky+kx]=original+3*(iw*ixyblockh*ky+ixyblockw*kx);
	}
	for(ky=0;ky<ixyblockh;++ky)//interleave
	{
		kx=0;
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

			slowptrs[NCODERS]		fastptr  (aligned because NCODERS == sizeof(__m128i[2]))
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
			if(sizeof(void*)==8)
			{
				__m128i mp[4];
				mp[0]=_mm_add_epi64(slowinc, _mm_load_si128((__m128i*)slowptrs+0));
				mp[1]=_mm_add_epi64(slowinc, _mm_load_si128((__m128i*)slowptrs+1));
				mp[2]=_mm_add_epi64(slowinc, _mm_load_si128((__m128i*)slowptrs+2));
				mp[3]=_mm_add_epi64(slowinc, _mm_load_si128((__m128i*)slowptrs+3));
				_mm_store_si128((__m128i*)slowptrs+0, mp[0]);
				_mm_store_si128((__m128i*)slowptrs+1, mp[1]);
				_mm_store_si128((__m128i*)slowptrs+2, mp[2]);
				_mm_store_si128((__m128i*)slowptrs+3, mp[3]);
				mp[0]=_mm_add_epi64(slowinc, _mm_load_si128((__m128i*)slowptrs+0+4));
				mp[1]=_mm_add_epi64(slowinc, _mm_load_si128((__m128i*)slowptrs+1+4));
				mp[2]=_mm_add_epi64(slowinc, _mm_load_si128((__m128i*)slowptrs+2+4));
				mp[3]=_mm_add_epi64(slowinc, _mm_load_si128((__m128i*)slowptrs+3+4));
				_mm_store_si128((__m128i*)slowptrs+0+4, mp[0]);
				_mm_store_si128((__m128i*)slowptrs+1+4, mp[1]);
				_mm_store_si128((__m128i*)slowptrs+2+4, mp[2]);
				_mm_store_si128((__m128i*)slowptrs+3+4, mp[3]);
			}
			else
			{
				__m128i mp[4];
				mp[0]=_mm_add_epi32(slowinc, _mm_load_si128((__m128i*)slowptrs+0));
				mp[1]=_mm_add_epi32(slowinc, _mm_load_si128((__m128i*)slowptrs+1));
				mp[2]=_mm_add_epi32(slowinc, _mm_load_si128((__m128i*)slowptrs+2));
				mp[3]=_mm_add_epi32(slowinc, _mm_load_si128((__m128i*)slowptrs+3));
				_mm_store_si128((__m128i*)slowptrs+0, mp[0]);
				_mm_store_si128((__m128i*)slowptrs+1, mp[1]);
				_mm_store_si128((__m128i*)slowptrs+2, mp[2]);
				_mm_store_si128((__m128i*)slowptrs+3, mp[3]);
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
		{
			int k;
			for(k=0;k<NCODERS;++k)
				slowptrs0[k]+=rowstride;
		}
	}
}
static void interleave_blocks_inv(const unsigned char *interleaved, int iw, int ih, unsigned char *original)
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
	__m128i slowinc=sizeof(void*)==8?_mm_set_epi32(0, sizeof(__m128i), 0, sizeof(__m128i)):_mm_set1_epi32(sizeof(__m128i));
#endif
	const unsigned char *fastptr=interleaved;
	ALIGN(32) unsigned char *slowptrs[NCODERS]={0}, *slowptrs0[NCODERS]={0};
	int kx, ky;

	for(ky=0;ky<YCODERS;++ky)//spread slow pointers
	{
		for(kx=0;kx<XCODERS;++kx)
			slowptrs0[XCODERS*ky+kx]=original+3*(iw*ixyblockh*ky+ixyblockw*kx);
	}
	for(ky=0;ky<ixyblockh;++ky)//interleave
	{
		kx=0;
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

			slowptrs[NCODERS]		fastptr  (aligned because NCODERS == sizeof(__m128i[2]))
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
			if(sizeof(void*)==8)
			{
				__m128i mp[4];
				mp[0]=_mm_add_epi64(slowinc, _mm_load_si128((__m128i*)slowptrs+0));
				mp[1]=_mm_add_epi64(slowinc, _mm_load_si128((__m128i*)slowptrs+1));
				mp[2]=_mm_add_epi64(slowinc, _mm_load_si128((__m128i*)slowptrs+2));
				mp[3]=_mm_add_epi64(slowinc, _mm_load_si128((__m128i*)slowptrs+3));
				_mm_store_si128((__m128i*)slowptrs+0, mp[0]);
				_mm_store_si128((__m128i*)slowptrs+1, mp[1]);
				_mm_store_si128((__m128i*)slowptrs+2, mp[2]);
				_mm_store_si128((__m128i*)slowptrs+3, mp[3]);
				mp[0]=_mm_add_epi64(slowinc, _mm_load_si128((__m128i*)slowptrs+0+4));
				mp[1]=_mm_add_epi64(slowinc, _mm_load_si128((__m128i*)slowptrs+1+4));
				mp[2]=_mm_add_epi64(slowinc, _mm_load_si128((__m128i*)slowptrs+2+4));
				mp[3]=_mm_add_epi64(slowinc, _mm_load_si128((__m128i*)slowptrs+3+4));
				_mm_store_si128((__m128i*)slowptrs+0+4, mp[0]);
				_mm_store_si128((__m128i*)slowptrs+1+4, mp[1]);
				_mm_store_si128((__m128i*)slowptrs+2+4, mp[2]);
				_mm_store_si128((__m128i*)slowptrs+3+4, mp[3]);
			}
			else
			{
				__m128i mp[4];
				mp[0]=_mm_add_epi32(slowinc, _mm_load_si128((__m128i*)slowptrs+0));
				mp[1]=_mm_add_epi32(slowinc, _mm_load_si128((__m128i*)slowptrs+1));
				mp[2]=_mm_add_epi32(slowinc, _mm_load_si128((__m128i*)slowptrs+2));
				mp[3]=_mm_add_epi32(slowinc, _mm_load_si128((__m128i*)slowptrs+3));
				_mm_store_si128((__m128i*)slowptrs+0, mp[0]);
				_mm_store_si128((__m128i*)slowptrs+1, mp[1]);
				_mm_store_si128((__m128i*)slowptrs+2, mp[2]);
				_mm_store_si128((__m128i*)slowptrs+3, mp[3]);
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
		{
			int k;
			for(k=0;k<NCODERS;++k)
				slowptrs0[k]+=rowstride;
		}
	}
}
static void save_ppm(const char *fn, const unsigned char *image, int iw, int ih)
{
	FILE *fdst=fopen(fn, "wb");
	if(!fdst)
	{
		LOG_ERROR("Cannot open \"%s\" for writing", fn);
		return;
	}
	fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
	fwrite(image, 1, (ptrdiff_t)3*iw*ih, fdst);
	fclose(fdst);
}
static void decorr1d(unsigned char *data, int count, int bytestride, int bestrct, int *rhist)
{
	const unsigned char *combination=rct_combinations[bestrct];
	int yidx=combination[II_PERM_Y];
	int uidx=combination[II_PERM_U];
	int vidx=combination[II_PERM_V];
	int ufromy=-(combination[II_COEFF_U_SUB_Y]!=0);
	int vc0=combination[II_COEFF_V_SUB_Y];
	int vc1=combination[II_COEFF_V_SUB_U];

	unsigned char *ptr=data;
	int prevy=0, prevu=0, prevv=0, offset=0;
	int k, vpred;
	for(k=0;k<count;++k)
	{
		int y=ptr[yidx]-128;
		int u=ptr[uidx]-128;
		int v=ptr[vidx]-128;
		int sym;
		ptr[0]=sym=(unsigned char)(y-prevy+128);
		++rhist[256*0+sym];
		prevy=y;

		offset=y&ufromy;
		prevu+=offset;
		CLAMP2(prevu, -128, 127);
		ptr[1]=sym=(unsigned char)(u-prevu+128);
		++rhist[256*1+sym];
		prevu=u-offset;

		offset=vc0*y+vc1*u;
		vpred=(prevv+offset)>>2;
		CLAMP2(vpred, -128, 127);
		ptr[2]=sym=(unsigned char)(v-vpred+128);
		++rhist[256*2+sym];
		prevv=4*v-offset;
		ptr+=bytestride;
	}
}
static void encode1d(unsigned char *data, int count, int bytestride, unsigned *pstate, unsigned char **pstreamptr, const unsigned char *streamend, const rANS_SIMD_SymInfo *rsyminfo)
{
	unsigned char *streamptr=*pstreamptr;
	unsigned state=*pstate;
	unsigned char *ptr=data+(count-(ptrdiff_t)1)*bytestride;
	const rANS_SIMD_SymInfo *info=0;
	int k;
	for(k=0;k<count;++k)
	{
		info=rsyminfo+ptr[2]+256*2;
		if(state>info->smax)
		{
			streamptr-=2;
#ifdef _DEBUG
			if(streamptr<=streamend)//"streamend" is buffer start
				LOG_ERROR("OOB ptr %016zX <= %016zX", streamptr, streamend);
#endif
			*(unsigned short*)streamptr=(unsigned short)state;
			state>>=RANS_RENORM_BITS;
		}
		//state += ((state*invf>>32)*(1<<(11-sh))>>11)*negf+cdf
		state+=(((unsigned long long)state*info->invf>>32)*info->sh>>(PROBBITS-1))*info->negf+info->cdf;
		//state+=((unsigned long long)state*info->invf>>32>>info->sh)*info->negf+info->cdf;
#ifdef ANS_VAL
		ansval_push(&state, sizeof(state), 1);
#endif

		info=rsyminfo+ptr[1]+256*1;
		if(state>info->smax)
		{
			streamptr-=2;
#ifdef _DEBUG
			if(streamptr<=streamend)
				LOG_ERROR("OOB ptr %016zX <= %016zX", streamptr, streamend);
#endif
			*(unsigned short*)streamptr=(unsigned short)state;
			state>>=RANS_RENORM_BITS;
		}
		state+=(((unsigned long long)state*info->invf>>32)*info->sh>>(PROBBITS-1))*info->negf+info->cdf;
#ifdef ANS_VAL
		ansval_push(&state, sizeof(state), 1);
#endif

		info=rsyminfo+ptr[0]+256*0;
		if(state>info->smax)
		{
			streamptr-=2;
#ifdef _DEBUG
			if(streamptr<=streamend)
				LOG_ERROR("OOB ptr %016zX <= %016zX", streamptr, streamend);
#endif
			*(unsigned short*)streamptr=(unsigned short)state;
			state>>=RANS_RENORM_BITS;
		}
		state+=(((unsigned long long)state*info->invf>>32)*info->sh>>(PROBBITS-1))*info->negf+info->cdf;
#ifdef ANS_VAL
		ansval_push(&state, sizeof(state), 1);
#endif
		ptr-=bytestride;
	}
	*pstreamptr=streamptr;
	*pstate=state;
}
static void decode1d(unsigned char *data, int count, int bytestride, int bestrct, unsigned *pstate, const unsigned char **pstreamptr, const unsigned char *streamend, unsigned *rCDF2syms)
{
	const unsigned char *combination=rct_combinations[bestrct];
	int yidx=combination[II_PERM_Y];
	int uidx=combination[II_PERM_U];
	int vidx=combination[II_PERM_V];
	int ufromy=-(combination[II_COEFF_U_SUB_Y]!=0);
	int vc0=combination[II_COEFF_V_SUB_Y];
	int vc1=combination[II_COEFF_V_SUB_U];

	const unsigned char *streamptr=*pstreamptr;
	unsigned state=*pstate;
	unsigned char *ptr=data;
	int prevy=0, prevu=0, prevv=0, offset=0;
	int y=0, u=0, v=0;
	int k, vpred;
	for(k=0;k<count;++k)
	{
		unsigned info;

		//yuv = (char)(error+N-128)
		info=rCDF2syms[0<<PROBBITS|(state&((1<<PROBBITS)-1))];
		y=(char)(info+prevy-128);
		prevy=y;
#ifdef ANS_VAL
		ansval_check(&state, sizeof(state), 1);
#endif
		state=(state>>PROBBITS)*(info>>(PROBBITS+8))+(info<<PROBBITS>>(32-PROBBITS));
		if(state<(1<<(RANS_STATE_BITS-RANS_RENORM_BITS)))
		{
#ifdef _DEBUG
			if(streamptr>streamend)
				LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
			state=state<<16|*(unsigned short*)streamptr;
			streamptr+=2;
		}

		offset=y&ufromy;
		prevu+=offset;
		CLAMP2(prevu, -128, 127);
		info=rCDF2syms[1<<PROBBITS|(state&((1<<PROBBITS)-1))];
		u=(char)(info+prevu-128);
		prevu=u-offset;
#ifdef ANS_VAL
		ansval_check(&state, sizeof(state), 1);
#endif
		state=(state>>PROBBITS)*(info>>(PROBBITS+8))+(info<<PROBBITS>>(32-PROBBITS));
		if(state<(1<<(RANS_STATE_BITS-RANS_RENORM_BITS)))
		{
#ifdef _DEBUG
			if(streamptr>streamend)
				LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
			state=state<<16|*(unsigned short*)streamptr;
			streamptr+=2;
		}

		offset=vc0*y+vc1*u;
		vpred=(prevv+offset)>>2;
		CLAMP2(vpred, -128, 127);
		info=rCDF2syms[2<<PROBBITS|(state&((1<<PROBBITS)-1))];
		v=(char)(info+vpred-128);
		prevv=4*v-offset;
#ifdef ANS_VAL
		ansval_check(&state, sizeof(state), 1);
#endif
		state=(state>>PROBBITS)*(info>>(PROBBITS+8))+(info<<PROBBITS>>(32-PROBBITS));
		if(state<(1<<(RANS_STATE_BITS-RANS_RENORM_BITS)))
		{
#ifdef _DEBUG
			if(streamptr>streamend)
				LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
			state=state<<16|*(unsigned short*)streamptr;
			streamptr+=2;
		}
		ptr[yidx]=y+128;
		ptr[uidx]=u+128;
		ptr[vidx]=v+128;
		ptr+=bytestride;
	}
	*pstreamptr=streamptr;
	*pstate=state;
}
int l1_codec(int argc, char **argv)
{
	const char *srcfn, *dstfn;
	int effort, dist;
#ifdef ESTIMATE_SIZE
	double esize[3*NCODERS]={0};
#endif
#ifdef LOUD
	double t=time_sec();
#endif
	int fwd=0, iw=0, ih=0, rowstride=0;
	ptrdiff_t usize=0, cap=0;
	unsigned char *image=0, *imptr=0, *streamptr=0, *streamstart=0, *streamend=0;
	int psize=0;
	short *pixels=0;
	ptrdiff_t cheadersize=0, csize=0;
	int blockw, blockh;
	int qxbytes;//iw/XCODERS*XCODERS*3
	int ixcount, ixbytes;//ix = interleaved circular buffer width		iw/XCODERS*NCODERS
	int xremw, yremh;
	int xrembytes;
	ptrdiff_t isize;
	ptrdiff_t interleavedsize;//fwd ? interleave residuals & context : pack residuals
	unsigned char *interleaved=0;
	int bestrct=0, npreds=0, sh=0, profile=0;
	unsigned long long ctxmask=0;//3*NCTX+3 = 54 flags	0: rare context (bypass)  1: emit stats
	const int hsize=(int)sizeof(int[3*NCTX<<8]);//3 channels
	int *hists=0;//fwd-only
	const int rhsize=(int)sizeof(int[3*256]);
	int *rhist=0;
	int CDF2syms_size=0;
	unsigned *CDF2syms=0;
	int rCDF2syms_size;
	unsigned *rCDF2syms=0;
	int L1statesize=0;
	int *L1state=0;
	const unsigned char *combination=0;
	int yidx, uidx, vidx;
	__m128i uhelpmask, vc0, vc1;
	int paddedwidth;
	__m128i mctxmax=_mm_set1_epi16(NCTX-1);
	__m128i mctxuoffset=_mm_set1_epi16(NCTX);
	__m128i mctxvoffset=_mm_set1_epi16(NCTX*2);
	__m128i amin=_mm_set1_epi16(-128);
	__m128i amax=_mm_set1_epi16(127);
	__m128i half8=_mm_set1_epi8(-128);
	__m128i bytemask=_mm_set1_epi16(255);
	__m128i wordmask=_mm_set1_epi32(0xFFFF);
	__m128i myuv[6];//y y u u v v
	__m128i dist_rcp=_mm_set1_epi16(0x7FFF), mdist=_mm_set1_epi16(1);
	unsigned char *ctxptr=0;
	__m128i mstate[4];
	__m128i *L1preds=0;
	int *L1weights=0;
	int kx, ky;

	if(argc!=3&&argc!=4&&argc!=5)
	{
		printf(
			"Usage: \"%s\"  input  output  [Effort]  [Dist]    Encode/decode.\n"
			"  Effort  =  0 CG / 1~3 L1 | 4 Profiler.\n"
			"  Dist    =  lossy distortion. 4 <= Dist <= 16.\n"
			"Built on %s %s\n"
			, argv[0]
			, __DATE__, __TIME__
		);
		return 1;
	}
	srcfn=argv[1];
	dstfn=argv[2];
	effort=argc<4?DEFAULT_EFFORT_LEVEL:atoi(argv[3]);
	profile=effort>>2;
	effort&=3;
	dist=argc<5?1:atoi(argv[4]);
	if(dist>1)
		CLAMP2(dist, 4, 16);
	prof_checkpoint(0, 0);
	if(!srcfn||!dstfn)
	{
		LOG_ERROR("Codec requires both source and destination filenames");
		return 1;
	}
	{
		int c=0;
		FILE *fsrc=fopen(srcfn, "rb");
		if(!fsrc)
		{
			LOG_ERROR("Cannot open \"%s\"", srcfn);
			return 1;
		}
		fread(&c, 1, 2, fsrc);
		fwd=c==('P'|'6'<<8);
		if(!fwd&&c!=('3'|'2'<<8))
		{
			LOG_ERROR("Unsupported file \"%s\"", srcfn);
			return 1;
		}
		if(fwd)
		{
			int nread, vmax;
#ifdef LOUD
			print_timestamp("%Y-%m-%d_%H%M%S\n");
#endif
			c=fgetc(fsrc);
			if(c!='\n')
			{
				LOG_ERROR("Invalid PPM file");
				return 1;
			}
			nread=fscanf(fsrc, "%d %d", &iw, &ih);
			if(nread!=2)
			{
				LOG_ERROR("Unsupported PPM file");
				return 1;
			}
			vmax=0;
			nread=fscanf(fsrc, "%d", &vmax);
			if(nread!=1||vmax!=255)
			{
				LOG_ERROR("Unsupported PPM file");
				return 1;
			}
			c=fgetc(fsrc);
			if(c!='\n')
			{
				LOG_ERROR("Invalid PPM file");
				return 1;
			}
		}
		else
		{
			fread(&iw, 1, 4, fsrc);
			fread(&ih, 1, 4, fsrc);
			cheadersize=ftell(fsrc);
		}
		if(iw<1||ih<1)
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
		rowstride=3*iw;
		usize=(ptrdiff_t)3*iw*ih;
		cap=(ptrdiff_t)4*iw*ih;
		image=(unsigned char*)malloc(cap+sizeof(__m128i));
		if(!image)
		{
			LOG_ERROR("Alloc error");
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
			csize=get_filesize(srcfn);
			streamptr=streamstart=image+cap-(csize-cheadersize)-sizeof(__m128i);
			streamend=image+cap-sizeof(__m128i);
			fread(streamstart, 1, csize-cheadersize, fsrc);//read stream
		}
		fclose(fsrc);
	}
	prof_checkpoint(fwd?usize:csize, "fread");
	blockw=iw/XCODERS;
	blockh=ih/YCODERS;
	qxbytes=blockw*XCODERS*3;//iw/XCODERS*XCODERS*3
	ixcount=blockw*NCODERS;
	ixbytes=3*ixcount;//ix = interleaved circular buffer width		iw/XCODERS*NCODERS
	xremw=iw-blockw*XCODERS;
	yremh=ih-blockh*YCODERS;
	xrembytes=3*xremw;
	isize=(ptrdiff_t)ixbytes*blockh;
	interleavedsize=isize<<fwd;//fwd ? interleave residuals & context : pack residuals
	interleaved=(unsigned char*)_mm_malloc(interleavedsize, sizeof(__m128i));
	if(!interleaved)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	(void)xrembytes;
	bestrct=0, npreds=0, sh=0;
	ctxmask=0;//3*NCTX+3 = 54 flags	0: rare context (bypass)  1: emit stats
	hists=fwd?(int*)malloc(hsize):0;//fwd-only
	rhist=fwd?(int*)malloc(rhsize):0;

	CDF2syms_size=(int)sizeof(int[3*NCTX<<PROBBITS]);
	if(fwd)//DIV-free rANS encoder reuses these as SIMD symbol info
		CDF2syms_size=(int)sizeof(rANS_SIMD_SymInfo[3*NCTX<<8]);
	CDF2syms=(unsigned*)_mm_malloc(CDF2syms_size, sizeof(__m128i));

	rCDF2syms_size=(int)sizeof(int[3<<PROBBITS]);
	if(fwd)
		rCDF2syms_size=(int)sizeof(rANS_SIMD_SymInfo[3<<8]);
	rCDF2syms=(unsigned*)_mm_malloc(rCDF2syms_size, sizeof(__m128i));

	psize=(int)sizeof(short[4*6*NCODERS])*(blockw+16);//4 padded rows  *  {Y*NCODERS, U*NCODERS, V*NCODERS,  eY*NCODERS, eU*NCODERS, eV*NCODERS} = 2*3*32 = 192 channels  ~48*iw bytes
	pixels=(short*)_mm_malloc(psize, sizeof(__m128i));//~188 KB for 4K/12MP
	if((fwd&&(!hists||!rhist))||!CDF2syms||!rCDF2syms||!pixels)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(pixels, 0, psize);
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
			LOG_ERROR("ERROR");
		printf("SUCCESS\n");
		exit(0);
#endif
		interleave_blocks_fwd(image, iw, ih, interleaved+isize);//reuse memory: read 8-bit pixel, write 16-bit context<<8|residual
		guide_save(interleaved+isize, ixcount, blockh);
		prof_checkpoint(usize, "interleave");
		if(dist<=1)
		{
			int kx, ky;
			ALIGN(32) long long counters[OCH_COUNT]={0};
			__m128i mcounters[OCH_COUNT];//64-bit
			__m128i half8=_mm_set1_epi8(-128);
			__m128i wordmask=_mm_set_epi32(0, 0xFFFF, 0, 0xFFFF);
			memset(mcounters, 0, sizeof(mcounters));
			imptr=interleaved+isize;
			for(ky=0;ky<blockh;ky+=ANALYSIS_YSTRIDE)//analysis
			{
				__m128i prev[OCH_COUNT][2];//16-bit
				memset(prev, 0, sizeof(prev));
				for(kx=0;kx<ixbytes-3*NCODERS;kx+=3*NCODERS*ANALYSIS_XSTRIDE)
				{
					__m128i rg0, gb0, br0;
					__m128i rg1, gb1, br1;
					__m128i r0=_mm_add_epi8(_mm_load_si128((__m128i*)imptr+0), half8);
					__m128i g0=_mm_add_epi8(_mm_load_si128((__m128i*)imptr+1), half8);
					__m128i b0=_mm_add_epi8(_mm_load_si128((__m128i*)imptr+2), half8);
					__m128i r1=_mm_cvtepi8_epi16(_mm_shuffle_epi32(r0, _MM_SHUFFLE(1, 0, 3, 2)));
					__m128i g1=_mm_cvtepi8_epi16(_mm_shuffle_epi32(g0, _MM_SHUFFLE(1, 0, 3, 2)));
					__m128i b1=_mm_cvtepi8_epi16(_mm_shuffle_epi32(b0, _MM_SHUFFLE(1, 0, 3, 2)));
					r0=_mm_cvtepi8_epi16(r0);
					g0=_mm_cvtepi8_epi16(g0);
					b0=_mm_cvtepi8_epi16(b0);
					imptr+=3*NCODERS*ANALYSIS_XSTRIDE;
					r0=_mm_slli_epi16(r0, 2);
					g0=_mm_slli_epi16(g0, 2);
					b0=_mm_slli_epi16(b0, 2);
					r1=_mm_slli_epi16(r1, 2);
					g1=_mm_slli_epi16(g1, 2);
					b1=_mm_slli_epi16(b1, 2);
					rg0=_mm_sub_epi16(r0, g0);
					gb0=_mm_sub_epi16(g0, b0);
					br0=_mm_sub_epi16(b0, r0);
					rg1=_mm_sub_epi16(r1, g1);
					gb1=_mm_sub_epi16(g1, b1);
					br1=_mm_sub_epi16(b1, r1);
#ifdef ENABLE_RCT_EXTENSION
#endif
#define UPDATE(IDXA, IDXB, IDXC, A0, B0, C0, LANE)\
	do\
	{\
		__m128i t0=A0, t1=B0, t2=C0;\
		__m128i ta=_mm_sub_epi16(t0, prev[IDXA][LANE]);\
		__m128i tb=_mm_sub_epi16(t1, prev[IDXB][LANE]);\
		__m128i tc=_mm_sub_epi16(t2, prev[IDXC][LANE]);\
		prev[IDXA][LANE]=t0;\
		prev[IDXB][LANE]=t1;\
		prev[IDXC][LANE]=t2;\
		ta=_mm_abs_epi16(ta);\
		tb=_mm_abs_epi16(tb);\
		tc=_mm_abs_epi16(tc);\
		ta=_mm_add_epi16(ta, _mm_srli_epi64(ta, 32));\
		tb=_mm_add_epi16(tb, _mm_srli_epi64(tb, 32));\
		tc=_mm_add_epi16(tc, _mm_srli_epi64(tc, 32));\
		ta=_mm_add_epi16(ta, _mm_srli_epi64(ta, 16));\
		tb=_mm_add_epi16(tb, _mm_srli_epi64(tb, 16));\
		tc=_mm_add_epi16(tc, _mm_srli_epi64(tc, 16));\
		mcounters[IDXA]=_mm_add_epi64(mcounters[IDXA], _mm_and_si128(ta, wordmask));\
		mcounters[IDXB]=_mm_add_epi64(mcounters[IDXB], _mm_and_si128(tb, wordmask));\
		mcounters[IDXC]=_mm_add_epi64(mcounters[IDXC], _mm_and_si128(tc, wordmask));\
	}while(0)
					UPDATE(OCH_Y400, OCH_Y040, OCH_Y004, r0, g0, b0, 0);
					UPDATE(OCH_Y400, OCH_Y040, OCH_Y004, r1, g1, b1, 1);
					UPDATE(OCH_CX40, OCH_C0X4, OCH_C40X, rg0, gb0, br0, 0);
					UPDATE(OCH_CX40, OCH_C0X4, OCH_C40X, rg1, gb1, br1, 1);
#ifdef ENABLE_RCT_EXTENSION
					UPDATE(OCH_CX31, OCH_C3X1, OCH_C31X,
						_mm_add_epi16(rg0, _mm_srai_epi16(gb0, 2)),//r-(3*g+b)/4 = r-g-(b-g)/4
						_mm_add_epi16(rg0, _mm_srai_epi16(br0, 2)),//g-(3*r+b)/4 = g-r-(b-r)/4
						_mm_add_epi16(br0, _mm_srai_epi16(rg0, 2)),//b-(3*r+g)/4 = b-r-(g-r)/4
						0
					);
					UPDATE(OCH_CX31, OCH_C3X1, OCH_C31X,
						_mm_add_epi16(rg1, _mm_srai_epi16(gb1, 2)),
						_mm_add_epi16(rg1, _mm_srai_epi16(br1, 2)),
						_mm_add_epi16(br1, _mm_srai_epi16(rg1, 2)),
						1
					);
					UPDATE(OCH_CX13, OCH_C1X3, OCH_C13X,
						_mm_add_epi16(br0, _mm_srai_epi16(gb0, 2)),//r-(g+3*b)/4 = r-b-(g-b)/4
						_mm_add_epi16(gb0, _mm_srai_epi16(br0, 2)),//g-(r+3*b)/4 = g-b-(r-b)/4
						_mm_add_epi16(gb0, _mm_srai_epi16(rg0, 2)),//b-(r+3*g)/4 = b-g-(r-g)/4
						0
					);
					UPDATE(OCH_CX13, OCH_C1X3, OCH_C13X,
						_mm_add_epi16(br1, _mm_srai_epi16(gb1, 2)),
						_mm_add_epi16(gb1, _mm_srai_epi16(br1, 2)),
						_mm_add_epi16(gb1, _mm_srai_epi16(rg1, 2)),
						1
					);
					UPDATE(OCH_CX22, OCH_C2X2, OCH_C22X,
						_mm_srai_epi16(_mm_sub_epi16(rg0, br0), 1),//r-(g+b)/2 = (r-g + r-b)/2
						_mm_srai_epi16(_mm_sub_epi16(gb0, rg0), 1),//g-(r+b)/2 = (g-r + g-b)/2
						_mm_srai_epi16(_mm_sub_epi16(br0, gb0), 1),//b-(r+g)/2 = (b-r + b-g)/2
						0
					);
					UPDATE(OCH_CX22, OCH_C2X2, OCH_C22X,
						_mm_srai_epi16(_mm_sub_epi16(rg1, br1), 1),
						_mm_srai_epi16(_mm_sub_epi16(gb1, rg1), 1),
						_mm_srai_epi16(_mm_sub_epi16(br1, gb1), 1),
						1
					);
#endif
				}
				imptr+=ixbytes*(ANALYSIS_YSTRIDE-1);
			}
			{
				int k;
				for(k=0;k<OCH_COUNT;++k)
				{
					ALIGN(32) int64_t temp[2]={0};
					_mm_store_si128((__m128i*)temp, mcounters[k]);
					counters[k]=temp[0]+temp[1];
				}
			}
			{
				long long minerr=0;
				int kt;
				for(kt=0;kt<RCT_COUNT;++kt)
				{
					const unsigned char *rct=rct_combinations[kt];
					long long currerr=
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
			}

//			bestrct=0;//
//#ifdef __GNUC__
//#error remove above
//#endif
			//printf("%2d ", bestrct);
			prof_checkpoint(usize, "analysis");
		}
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
#ifndef LOUD
		if(profile)
#endif
			printf("%s  NPREDS=%d  %d bytes\n", rct_names[bestrct], npreds, (int)usize);
	}
	else
	{
		//decode flags, stats
		int flags=*streamptr++;
		effort=flags&3;
		flags>>=2;
		bestrct=flags%RCT_COUNT;
		dist=*streamptr++;
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

		ctxmask=*(unsigned long long*)streamptr;
		streamptr+=8;
		{
			BitPackerLIFO ec;
			int kc;

			bitpacker_dec_init(&ec, streamptr, streamend);
			for(kc=0;kc<3*NCTX;++kc)
				dec_unpackhist(&ec, CDF2syms+((ptrdiff_t)kc<<PROBBITS), ctxmask, kc);
			if(xremw||yremh)
			{
				dec_unpackhist(&ec, rCDF2syms+((ptrdiff_t)0<<PROBBITS), ctxmask, 3*NCTX+0);
				dec_unpackhist(&ec, rCDF2syms+((ptrdiff_t)1<<PROBBITS), ctxmask, 3*NCTX+1);
				dec_unpackhist(&ec, rCDF2syms+((ptrdiff_t)2<<PROBBITS), ctxmask, 3*NCTX+2);
			}
			streamptr=(unsigned char*)(size_t)ec.srcfwdptr;
		}
		prof_checkpoint((ptrdiff_t)CDF2syms_size+rCDF2syms_size, "unpack histograms");
	}
	if(effort)
	{
		L1statesize=(int)sizeof(int[2*NCODERS*3*(L1_NPREDS3+1)]);//{preds, coeffs} * (NPREDS+{bias}) * 3 channels * NCODERS
		L1state=(int*)_mm_malloc(L1statesize, sizeof(__m128i));
		if(!L1state)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memset(L1state, 0, L1statesize);
	}
	combination=rct_combinations[bestrct];
	yidx=combination[II_PERM_Y]*NCODERS;
	uidx=combination[II_PERM_U]*NCODERS;
	vidx=combination[II_PERM_V]*NCODERS;
	uhelpmask=_mm_set1_epi16(-(combination[II_COEFF_U_SUB_Y]!=0));
	vc0=_mm_set1_epi16(combination[II_COEFF_V_SUB_Y]);
	vc1=_mm_set1_epi16(combination[II_COEFF_V_SUB_U]);
	paddedwidth=blockw+16;
	if(dist>1)
	{
		dist_rcp=_mm_set1_epi16(((1<<16)+dist-1)/dist);//x/dist  ->  {x*=inv; x=(x>>16)+((unsigned)x>>31);}
		mdist=_mm_set1_epi16(dist);
	}
	memset(myuv, 0, sizeof(myuv));
	ctxptr=interleaved;
	imptr=interleaved+(fwd?isize:0);
	L1preds=effort?(__m128i*)L1state:0;
	L1weights=effort?(int*)(L1state+1*(ptrdiff_t)NCODERS*3*(L1_NPREDS3+1)):0;
	if(effort)
		FILLMEM(L1weights, (1<<sh)/npreds, npreds*sizeof(int[6*8]), sizeof(int));
	//	FILLMEM(L1weights, (1<<sh)/npreds, (npreds+1)*sizeof(int[6*8]), sizeof(int));//bias
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
			LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
		memcpy(mstate, streamptr, sizeof(mstate));
		streamptr+=sizeof(mstate);
	}
	for(ky=0;ky<blockh;++ky)//main coding loop
	{
		ALIGN(32) short *rows[]=
		{
			pixels+(paddedwidth*((ky-0LL)&3)+8LL)*6*NCODERS,
			pixels+(paddedwidth*((ky-1LL)&3)+8LL)*6*NCODERS,
			pixels+(paddedwidth*((ky-2LL)&3)+8LL)*6*NCODERS,
			pixels+(paddedwidth*((ky-3LL)&3)+8LL)*6*NCODERS,
		};
		ALIGN(32) unsigned short syms[3*NCODERS]={0};
		__m128i NW[6], N[6], W[6];
		__m128i eW[6], ecurr[6], eNEE[6], eNEEE[6];

		memset(NW, 0, sizeof(NW));
		memset(N, 0, sizeof(N));
		memset(W, 0, sizeof(W));
		memset(eW, 0, sizeof(eW));
		memset(ecurr, 0, sizeof(ecurr));
		memset(eNEE, 0, sizeof(eNEE));
		memset(eNEEE, 0, sizeof(eNEEE));
		//                    (__m128i*)rows[-Y]+E+C+X*12
		eNEE[0]=_mm_load_si128((__m128i*)rows[1]+6+0+2*12);
		eNEE[1]=_mm_load_si128((__m128i*)rows[1]+6+1+2*12);
		eNEE[2]=_mm_load_si128((__m128i*)rows[1]+6+2+2*12);
		eNEE[3]=_mm_load_si128((__m128i*)rows[1]+6+3+2*12);
		eNEE[4]=_mm_load_si128((__m128i*)rows[1]+6+4+2*12);
		eNEE[5]=_mm_load_si128((__m128i*)rows[1]+6+5+2*12);
		for(kx=0;kx<ixbytes;kx+=3*NCODERS)
		{
			__m128i
				predY[2], ctxY[2],
				predU[2], ctxU[2],
				predV[2], ctxV[2];
			__m128i predYUV0[6];
			__m128i msyms[2], moffset[2];

			//                    (__m128i*)rows[-Y]+E+C+X*12
			N[0]=_mm_load_si128((__m128i*)rows[1]+0+0+0*12);//y neighbors
			N[1]=_mm_load_si128((__m128i*)rows[1]+0+1+0*12);
			N[2]=_mm_load_si128((__m128i*)rows[1]+0+2+0*12);//u
			N[3]=_mm_load_si128((__m128i*)rows[1]+0+3+0*12);
			N[4]=_mm_load_si128((__m128i*)rows[1]+0+4+0*12);//v
			N[5]=_mm_load_si128((__m128i*)rows[1]+0+5+0*12);
			{
				//context = FLOOR_LOG2(eW*eW+1)
				__m128i one=_mm_set1_epi32(1);
				__m128i cy0=_mm_and_si128(eW[0], wordmask), cy1=_mm_srli_epi32(eW[0], 16);
				__m128i cy2=_mm_and_si128(eW[1], wordmask), cy3=_mm_srli_epi32(eW[1], 16);
				__m128i cu0=_mm_and_si128(eW[2], wordmask), cu1=_mm_srli_epi32(eW[2], 16);
				__m128i cu2=_mm_and_si128(eW[3], wordmask), cu3=_mm_srli_epi32(eW[3], 16);
				__m128i cv0=_mm_and_si128(eW[4], wordmask), cv1=_mm_srli_epi32(eW[4], 16);
				__m128i cv2=_mm_and_si128(eW[5], wordmask), cv3=_mm_srli_epi32(eW[5], 16);
				cy0=_mm_mullo_epi32(cy0, cy0);
				cy1=_mm_mullo_epi32(cy1, cy1);
				cy2=_mm_mullo_epi32(cy2, cy2);
				cy3=_mm_mullo_epi32(cy3, cy3);
				cu0=_mm_mullo_epi32(cu0, cu0);
				cu1=_mm_mullo_epi32(cu1, cu1);
				cu2=_mm_mullo_epi32(cu2, cu2);
				cu3=_mm_mullo_epi32(cu3, cu3);
				cv0=_mm_mullo_epi32(cv0, cv0);
				cv1=_mm_mullo_epi32(cv1, cv1);
				cv2=_mm_mullo_epi32(cv2, cv2);
				cv3=_mm_mullo_epi32(cv3, cv3);
				cy0=_mm_add_epi32(cy0, one);
				cy1=_mm_add_epi32(cy1, one);
				cy2=_mm_add_epi32(cy2, one);
				cy3=_mm_add_epi32(cy3, one);
				cu0=_mm_add_epi32(cu0, one);
				cu1=_mm_add_epi32(cu1, one);
				cu2=_mm_add_epi32(cu2, one);
				cu3=_mm_add_epi32(cu3, one);
				cv0=_mm_add_epi32(cv0, one);
				cv1=_mm_add_epi32(cv1, one);
				cv2=_mm_add_epi32(cv2, one);
				cv3=_mm_add_epi32(cv3, one);
				//FLOOR_LOG2_32x8(X) = _mm_sub_epi32(_mm_srli_epi32(_mm_castps_si128(_mm_cvtepi32_ps(X)), 23), _mm_set1_epi32(127))
				cy0=_mm_castps_si128(_mm_cvtepi32_ps(cy0));
				cy1=_mm_castps_si128(_mm_cvtepi32_ps(cy1));
				cy2=_mm_castps_si128(_mm_cvtepi32_ps(cy2));
				cy3=_mm_castps_si128(_mm_cvtepi32_ps(cy3));
				cu0=_mm_castps_si128(_mm_cvtepi32_ps(cu0));
				cu1=_mm_castps_si128(_mm_cvtepi32_ps(cu1));
				cu2=_mm_castps_si128(_mm_cvtepi32_ps(cu2));
				cu3=_mm_castps_si128(_mm_cvtepi32_ps(cu3));
				cv0=_mm_castps_si128(_mm_cvtepi32_ps(cv0));
				cv1=_mm_castps_si128(_mm_cvtepi32_ps(cv1));
				cv2=_mm_castps_si128(_mm_cvtepi32_ps(cv2));
				cv3=_mm_castps_si128(_mm_cvtepi32_ps(cv3));
				cy0=_mm_srli_epi32(cy0, 23);
				cy1=_mm_srli_epi32(cy1, 23);
				cy2=_mm_srli_epi32(cy2, 23);
				cy3=_mm_srli_epi32(cy3, 23);
				cu0=_mm_srli_epi32(cu0, 23);
				cu1=_mm_srli_epi32(cu1, 23);
				cu2=_mm_srli_epi32(cu2, 23);
				cu3=_mm_srli_epi32(cu3, 23);
				cv0=_mm_srli_epi32(cv0, 23);
				cv1=_mm_srli_epi32(cv1, 23);
				cv2=_mm_srli_epi32(cv2, 23);
				cv3=_mm_srli_epi32(cv3, 23);
				{
					__m128i expbias=_mm_set1_epi32(127);
					cy0=_mm_sub_epi32(cy0, expbias);
					cy1=_mm_sub_epi32(cy1, expbias);
					cy2=_mm_sub_epi32(cy2, expbias);
					cy3=_mm_sub_epi32(cy3, expbias);
					cu0=_mm_sub_epi32(cu0, expbias);
					cu1=_mm_sub_epi32(cu1, expbias);
					cu2=_mm_sub_epi32(cu2, expbias);
					cu3=_mm_sub_epi32(cu3, expbias);
					cv0=_mm_sub_epi32(cv0, expbias);
					cv1=_mm_sub_epi32(cv1, expbias);
					cv2=_mm_sub_epi32(cv2, expbias);
					cv3=_mm_sub_epi32(cv3, expbias);
				}
				cy1=_mm_slli_epi32(cy1, 16);
				cu1=_mm_slli_epi32(cu1, 16);
				cv1=_mm_slli_epi32(cv1, 16);
				cy3=_mm_slli_epi32(cy3, 16);
				cu3=_mm_slli_epi32(cu3, 16);
				cv3=_mm_slli_epi32(cv3, 16);
				ctxY[0]=_mm_or_si128(cy0, cy1);
				ctxU[0]=_mm_or_si128(cu0, cu1);
				ctxV[0]=_mm_or_si128(cv0, cv1);
				ctxY[1]=_mm_or_si128(cy2, cy3);
				ctxU[1]=_mm_or_si128(cu2, cu3);
				ctxV[1]=_mm_or_si128(cv2, cv3);
				ctxY[0]=_mm_min_epi16(ctxY[0], mctxmax);
				ctxY[1]=_mm_min_epi16(ctxY[1], mctxmax);
				ctxU[0]=_mm_min_epi16(ctxU[0], mctxmax);
				ctxU[1]=_mm_min_epi16(ctxU[1], mctxmax);
				ctxV[0]=_mm_min_epi16(ctxV[0], mctxmax);
				ctxV[1]=_mm_min_epi16(ctxV[1], mctxmax);
			}
			{
				const int borderW=3;
				const int borderN=3;
				const int borderE=3;
				int cond_cg=(unsigned)(kx-3*NCODERS*borderW)>=(unsigned)(ixbytes-3*NCODERS*(borderW+borderE))
					||(unsigned)(ky-borderN)>=(unsigned)(blockh-borderN);
				__m128i
					ymin[2], ymax[2],
					umin[2], umax[2],
					vmin[2], vmax[2];
				__m128i mcg[6];

				ymin[0]=_mm_min_epi16(N[0*2+0], W[0*2+0]);	ymin[1]=_mm_min_epi16(N[0*2+1], W[0*2+1]);
				ymax[0]=_mm_max_epi16(N[0*2+0], W[0*2+0]);	ymax[1]=_mm_max_epi16(N[0*2+1], W[0*2+1]);
				umin[0]=_mm_min_epi16(N[1*2+0], W[1*2+0]);	umin[1]=_mm_min_epi16(N[1*2+1], W[1*2+1]);
				umax[0]=_mm_max_epi16(N[1*2+0], W[1*2+0]);	umax[1]=_mm_max_epi16(N[1*2+1], W[1*2+1]);
				vmin[0]=_mm_min_epi16(N[2*2+0], W[2*2+0]);	vmin[1]=_mm_min_epi16(N[2*2+1], W[2*2+1]);
				vmax[0]=_mm_max_epi16(N[2*2+0], W[2*2+0]);	vmax[1]=_mm_max_epi16(N[2*2+1], W[2*2+1]);
				predY[0]=_mm_add_epi16(N[0*2+0], W[0*2+0]);	predY[1]=_mm_add_epi16(N[0*2+1], W[0*2+1]);//N+W-NW
				predU[0]=_mm_add_epi16(N[1*2+0], W[1*2+0]);	predU[1]=_mm_add_epi16(N[1*2+1], W[1*2+1]);
				predV[0]=_mm_add_epi16(N[2*2+0], W[2*2+0]);	predV[1]=_mm_add_epi16(N[2*2+1], W[2*2+1]);
				predY[0]=_mm_sub_epi16(predY[0], NW[0*2+0]);	predY[1]=_mm_sub_epi16(predY[1], NW[0*2+1]);
				predU[0]=_mm_sub_epi16(predU[0], NW[1*2+0]);	predU[1]=_mm_sub_epi16(predU[1], NW[1*2+1]);
				predV[0]=_mm_sub_epi16(predV[0], NW[2*2+0]);	predV[1]=_mm_sub_epi16(predV[1], NW[2*2+1]);
				mcg[0]=predY[0];
				mcg[1]=predY[1];
				mcg[2]=predU[0];
				mcg[3]=predU[1];
				mcg[4]=predV[0];
				mcg[5]=predV[1];

				if(effort==1)//predict
				{
					__m128i mp[12], t[12];
					int kp;
					/*
					effort 1
					0	N+W-NW
					1	N
					2	NE
					3	W
					*/

					//N+W-NW
					L1preds[0*6+0]=predY[0];
					L1preds[0*6+1]=predY[1];
					L1preds[0*6+2]=predU[0];
					L1preds[0*6+3]=predU[1];
					L1preds[0*6+4]=predV[0];
					L1preds[0*6+5]=predV[1];

					//N
					L1preds[1*6+0]=N[0];
					L1preds[1*6+1]=N[1];
					L1preds[1*6+2]=N[2];
					L1preds[1*6+3]=N[3];
					L1preds[1*6+4]=N[4];
					L1preds[1*6+5]=N[5];

					//NE
					L1preds[2*6+0]=_mm_load_si128((__m128i*)rows[1]+0+0+1*12);
					L1preds[2*6+1]=_mm_load_si128((__m128i*)rows[1]+0+1+1*12);
					L1preds[2*6+2]=_mm_load_si128((__m128i*)rows[1]+0+2+1*12);
					L1preds[2*6+3]=_mm_load_si128((__m128i*)rows[1]+0+3+1*12);
					L1preds[2*6+4]=_mm_load_si128((__m128i*)rows[1]+0+4+1*12);
					L1preds[2*6+5]=_mm_load_si128((__m128i*)rows[1]+0+5+1*12);

					//W
					L1preds[3*6+0]=W[0];
					L1preds[3*6+1]=W[1];
					L1preds[3*6+2]=W[2];
					L1preds[3*6+3]=W[3];
					L1preds[3*6+4]=W[4];
					L1preds[3*6+5]=W[5];


					//mix
					mp[0x0]=_mm_setzero_si128();
					mp[0x1]=_mm_setzero_si128();
					mp[0x2]=_mm_setzero_si128();
					mp[0x3]=_mm_setzero_si128();
					mp[0x4]=_mm_setzero_si128();
					mp[0x5]=_mm_setzero_si128();
					mp[0x6]=_mm_setzero_si128();
					mp[0x7]=_mm_setzero_si128();
					mp[0x8]=_mm_setzero_si128();
					mp[0x9]=_mm_setzero_si128();
					mp[0xA]=_mm_setzero_si128();
					mp[0xB]=_mm_setzero_si128();
					for(kp=0;kp<L1_NPREDS1;++kp)
					{
						//signed 16 -> 32	3 lo 3 hi registers
						t[0x0]=_mm_slli_epi32(L1preds[kp*6+0], 16);//y lo
						t[0x1]=_mm_slli_epi32(L1preds[kp*6+1], 16);
						t[0x2]=_mm_slli_epi32(L1preds[kp*6+2], 16);//u
						t[0x3]=_mm_slli_epi32(L1preds[kp*6+3], 16);
						t[0x4]=_mm_slli_epi32(L1preds[kp*6+4], 16);//v
						t[0x5]=_mm_slli_epi32(L1preds[kp*6+5], 16);

						t[0x6]=_mm_srai_epi32(L1preds[kp*6+0], 16);//y hi
						t[0x7]=_mm_srai_epi32(L1preds[kp*6+1], 16);
						t[0x8]=_mm_srai_epi32(L1preds[kp*6+2], 16);//u
						t[0x9]=_mm_srai_epi32(L1preds[kp*6+3], 16);
						t[0xA]=_mm_srai_epi32(L1preds[kp*6+4], 16);//v
						t[0xB]=_mm_srai_epi32(L1preds[kp*6+5], 16);
						t[0x0]=_mm_srai_epi32(t[0x0], 16);
						t[0x1]=_mm_srai_epi32(t[0x1], 16);
						t[0x2]=_mm_srai_epi32(t[0x2], 16);
						t[0x3]=_mm_srai_epi32(t[0x3], 16);
						t[0x4]=_mm_srai_epi32(t[0x4], 16);
						t[0x5]=_mm_srai_epi32(t[0x5], 16);
						t[0x0]=_mm_mullo_epi32(t[0x0], _mm_load_si128((__m128i*)L1weights+kp*12+0x0));
						t[0x1]=_mm_mullo_epi32(t[0x1], _mm_load_si128((__m128i*)L1weights+kp*12+0x1));
						t[0x2]=_mm_mullo_epi32(t[0x2], _mm_load_si128((__m128i*)L1weights+kp*12+0x2));
						t[0x3]=_mm_mullo_epi32(t[0x3], _mm_load_si128((__m128i*)L1weights+kp*12+0x3));
						t[0x4]=_mm_mullo_epi32(t[0x4], _mm_load_si128((__m128i*)L1weights+kp*12+0x4));
						t[0x5]=_mm_mullo_epi32(t[0x5], _mm_load_si128((__m128i*)L1weights+kp*12+0x5));
						t[0x6]=_mm_mullo_epi32(t[0x6], _mm_load_si128((__m128i*)L1weights+kp*12+0x6));
						t[0x7]=_mm_mullo_epi32(t[0x7], _mm_load_si128((__m128i*)L1weights+kp*12+0x7));
						t[0x8]=_mm_mullo_epi32(t[0x8], _mm_load_si128((__m128i*)L1weights+kp*12+0x8));
						t[0x9]=_mm_mullo_epi32(t[0x9], _mm_load_si128((__m128i*)L1weights+kp*12+0x9));
						t[0xA]=_mm_mullo_epi32(t[0xA], _mm_load_si128((__m128i*)L1weights+kp*12+0xA));
						t[0xB]=_mm_mullo_epi32(t[0xB], _mm_load_si128((__m128i*)L1weights+kp*12+0xB));
						mp[0x0]=_mm_add_epi32(mp[0x0], t[0x0]);
						mp[0x1]=_mm_add_epi32(mp[0x1], t[0x1]);
						mp[0x2]=_mm_add_epi32(mp[0x2], t[0x2]);
						mp[0x3]=_mm_add_epi32(mp[0x3], t[0x3]);
						mp[0x4]=_mm_add_epi32(mp[0x4], t[0x4]);
						mp[0x5]=_mm_add_epi32(mp[0x5], t[0x5]);
						mp[0x6]=_mm_add_epi32(mp[0x6], t[0x6]);
						mp[0x7]=_mm_add_epi32(mp[0x7], t[0x7]);
						mp[0x8]=_mm_add_epi32(mp[0x8], t[0x8]);
						mp[0x9]=_mm_add_epi32(mp[0x9], t[0x9]);
						mp[0xA]=_mm_add_epi32(mp[0xA], t[0xA]);
						mp[0xB]=_mm_add_epi32(mp[0xB], t[0xB]);
					}
					{
						__m128i rcon=_mm_set1_epi32(1<<L1_SH1>>1);
						mp[0x0]=_mm_add_epi32(mp[0x0], rcon);//rounding to nearest
						mp[0x1]=_mm_add_epi32(mp[0x1], rcon);
						mp[0x2]=_mm_add_epi32(mp[0x2], rcon);
						mp[0x3]=_mm_add_epi32(mp[0x3], rcon);
						mp[0x4]=_mm_add_epi32(mp[0x4], rcon);
						mp[0x5]=_mm_add_epi32(mp[0x5], rcon);
						mp[0x6]=_mm_add_epi32(mp[0x6], rcon);
						mp[0x7]=_mm_add_epi32(mp[0x7], rcon);
						mp[0x8]=_mm_add_epi32(mp[0x8], rcon);
						mp[0x9]=_mm_add_epi32(mp[0x9], rcon);
						mp[0xA]=_mm_add_epi32(mp[0xA], rcon);
						mp[0xB]=_mm_add_epi32(mp[0xB], rcon);
					}
					mp[0x0]=_mm_srai_epi32(mp[0x0], L1_SH1);
					mp[0x1]=_mm_srai_epi32(mp[0x1], L1_SH1);
					mp[0x2]=_mm_srai_epi32(mp[0x2], L1_SH1);
					mp[0x3]=_mm_srai_epi32(mp[0x3], L1_SH1);
					mp[0x4]=_mm_srai_epi32(mp[0x4], L1_SH1);
					mp[0x5]=_mm_srai_epi32(mp[0x5], L1_SH1);
					mp[0x6]=_mm_slli_epi32(mp[0x6], 16-L1_SH1);//y hi
					mp[0x7]=_mm_slli_epi32(mp[0x7], 16-L1_SH1);
					mp[0x8]=_mm_slli_epi32(mp[0x8], 16-L1_SH1);//u hi
					mp[0x9]=_mm_slli_epi32(mp[0x9], 16-L1_SH1);
					mp[0xA]=_mm_slli_epi32(mp[0xA], 16-L1_SH1);//v hi
					mp[0xB]=_mm_slli_epi32(mp[0xB], 16-L1_SH1);
					//32 -> 16
					predY[0]=_mm_blend_epi16(mp[0x0], mp[0x6], 0xAA);
					predY[1]=_mm_blend_epi16(mp[0x1], mp[0x7], 0xAA);
					predU[0]=_mm_blend_epi16(mp[0x2], mp[0x8], 0xAA);
					predU[1]=_mm_blend_epi16(mp[0x3], mp[0x9], 0xAA);
					predV[0]=_mm_blend_epi16(mp[0x4], mp[0xA], 0xAA);
					predV[1]=_mm_blend_epi16(mp[0x5], mp[0xB], 0xAA);


					//loosen pred range
					if(!cond_cg)
					{
						t[0]=_mm_load_si128((__m128i*)rows[1]+0+0+1*12);//NE
						t[1]=_mm_load_si128((__m128i*)rows[1]+0+1+1*12);
						t[2]=_mm_load_si128((__m128i*)rows[1]+0+2+1*12);
						t[3]=_mm_load_si128((__m128i*)rows[1]+0+3+1*12);
						t[4]=_mm_load_si128((__m128i*)rows[1]+0+4+1*12);
						t[5]=_mm_load_si128((__m128i*)rows[1]+0+5+1*12);
						ymin[0]=_mm_min_epi16(ymin[0], t[0*2+0]); ymin[1]=_mm_min_epi16(ymin[1], t[0*2+1]);
						ymax[0]=_mm_max_epi16(ymax[0], t[0*2+0]); ymax[1]=_mm_max_epi16(ymax[1], t[0*2+1]);
						umin[0]=_mm_min_epi16(umin[0], t[1*2+0]); umin[1]=_mm_min_epi16(umin[1], t[1*2+1]);
						umax[0]=_mm_max_epi16(umax[0], t[1*2+0]); umax[1]=_mm_max_epi16(umax[1], t[1*2+1]);
						vmin[0]=_mm_min_epi16(vmin[0], t[2*2+0]); vmin[1]=_mm_min_epi16(vmin[1], t[2*2+1]);
						vmax[0]=_mm_max_epi16(vmax[0], t[2*2+0]); vmax[1]=_mm_max_epi16(vmax[1], t[2*2+1]);
						t[0]=_mm_load_si128((__m128i*)rows[1]+0+0+3*12);//NEEE
						t[1]=_mm_load_si128((__m128i*)rows[1]+0+1+3*12);
						t[2]=_mm_load_si128((__m128i*)rows[1]+0+2+3*12);
						t[3]=_mm_load_si128((__m128i*)rows[1]+0+3+3*12);
						t[4]=_mm_load_si128((__m128i*)rows[1]+0+4+3*12);
						t[5]=_mm_load_si128((__m128i*)rows[1]+0+5+3*12);
						ymin[0]=_mm_min_epi16(ymin[0], t[0*2+0]); ymin[1]=_mm_min_epi16(ymin[1], t[0*2+1]);
						ymax[0]=_mm_max_epi16(ymax[0], t[0*2+0]); ymax[1]=_mm_max_epi16(ymax[1], t[0*2+1]);
						umin[0]=_mm_min_epi16(umin[0], t[1*2+0]); umin[1]=_mm_min_epi16(umin[1], t[1*2+1]);
						umax[0]=_mm_max_epi16(umax[0], t[1*2+0]); umax[1]=_mm_max_epi16(umax[1], t[1*2+1]);
						vmin[0]=_mm_min_epi16(vmin[0], t[2*2+0]); vmin[1]=_mm_min_epi16(vmin[1], t[2*2+1]);
						vmax[0]=_mm_max_epi16(vmax[0], t[2*2+0]); vmax[1]=_mm_max_epi16(vmax[1], t[2*2+1]);
					}
				}
				else if(effort==2)
				{
					__m128i cache[6], mp[12], t[12];
					int kp;
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
					cache[0]=_mm_sub_epi16(N[0], _mm_load_si128((__m128i*)rows[2]+0+0+0*12));//N-NN
					cache[1]=_mm_sub_epi16(N[1], _mm_load_si128((__m128i*)rows[2]+0+1+0*12));
					cache[2]=_mm_sub_epi16(N[2], _mm_load_si128((__m128i*)rows[2]+0+2+0*12));
					cache[3]=_mm_sub_epi16(N[3], _mm_load_si128((__m128i*)rows[2]+0+3+0*12));
					cache[4]=_mm_sub_epi16(N[4], _mm_load_si128((__m128i*)rows[2]+0+4+0*12));
					cache[5]=_mm_sub_epi16(N[5], _mm_load_si128((__m128i*)rows[2]+0+5+0*12));
					cache[0]=_mm_add_epi16(cache[0], _mm_slli_epi16(cache[0], 1));//*3
					cache[1]=_mm_add_epi16(cache[1], _mm_slli_epi16(cache[1], 1));
					cache[2]=_mm_add_epi16(cache[2], _mm_slli_epi16(cache[2], 1));
					cache[3]=_mm_add_epi16(cache[3], _mm_slli_epi16(cache[3], 1));
					cache[4]=_mm_add_epi16(cache[4], _mm_slli_epi16(cache[4], 1));
					cache[5]=_mm_add_epi16(cache[5], _mm_slli_epi16(cache[5], 1));
					L1preds[2*6+0]=_mm_add_epi16(cache[0], _mm_load_si128((__m128i*)rows[3]+0+0+0*12));//+NNN
					L1preds[2*6+1]=_mm_add_epi16(cache[1], _mm_load_si128((__m128i*)rows[3]+0+1+0*12));
					L1preds[2*6+2]=_mm_add_epi16(cache[2], _mm_load_si128((__m128i*)rows[3]+0+2+0*12));
					L1preds[2*6+3]=_mm_add_epi16(cache[3], _mm_load_si128((__m128i*)rows[3]+0+3+0*12));
					L1preds[2*6+4]=_mm_add_epi16(cache[4], _mm_load_si128((__m128i*)rows[3]+0+4+0*12));
					L1preds[2*6+5]=_mm_add_epi16(cache[5], _mm_load_si128((__m128i*)rows[3]+0+5+0*12));

					//3*(W-WW)+WWW
					cache[0]=_mm_sub_epi16(W[0], _mm_load_si128((__m128i*)rows[0]+0+0-2*12));//W-WW
					cache[1]=_mm_sub_epi16(W[1], _mm_load_si128((__m128i*)rows[0]+0+1-2*12));
					cache[2]=_mm_sub_epi16(W[2], _mm_load_si128((__m128i*)rows[0]+0+2-2*12));
					cache[3]=_mm_sub_epi16(W[3], _mm_load_si128((__m128i*)rows[0]+0+3-2*12));
					cache[4]=_mm_sub_epi16(W[4], _mm_load_si128((__m128i*)rows[0]+0+4-2*12));
					cache[5]=_mm_sub_epi16(W[5], _mm_load_si128((__m128i*)rows[0]+0+5-2*12));
					cache[0]=_mm_add_epi16(cache[0], _mm_slli_epi16(cache[0], 1));//*3
					cache[1]=_mm_add_epi16(cache[1], _mm_slli_epi16(cache[1], 1));
					cache[2]=_mm_add_epi16(cache[2], _mm_slli_epi16(cache[2], 1));
					cache[3]=_mm_add_epi16(cache[3], _mm_slli_epi16(cache[3], 1));
					cache[4]=_mm_add_epi16(cache[4], _mm_slli_epi16(cache[4], 1));
					cache[5]=_mm_add_epi16(cache[5], _mm_slli_epi16(cache[5], 1));
					L1preds[3*6+0]=_mm_add_epi16(cache[0], _mm_load_si128((__m128i*)rows[0]+0+0-3*12));//+WWW
					L1preds[3*6+1]=_mm_add_epi16(cache[1], _mm_load_si128((__m128i*)rows[0]+0+1-3*12));
					L1preds[3*6+2]=_mm_add_epi16(cache[2], _mm_load_si128((__m128i*)rows[0]+0+2-3*12));
					L1preds[3*6+3]=_mm_add_epi16(cache[3], _mm_load_si128((__m128i*)rows[0]+0+3-3*12));//+WWW
					L1preds[3*6+4]=_mm_add_epi16(cache[4], _mm_load_si128((__m128i*)rows[0]+0+4-3*12));
					L1preds[3*6+5]=_mm_add_epi16(cache[5], _mm_load_si128((__m128i*)rows[0]+0+5-3*12));

					//W+NE-N
					cache[0]=_mm_sub_epi16(W[0], N[0]);
					cache[1]=_mm_sub_epi16(W[1], N[1]);
					cache[2]=_mm_sub_epi16(W[2], N[2]);
					cache[3]=_mm_sub_epi16(W[3], N[3]);
					cache[4]=_mm_sub_epi16(W[4], N[4]);
					cache[5]=_mm_sub_epi16(W[5], N[5]);
					L1preds[4*6+0]=_mm_add_epi16(cache[0], _mm_load_si128((__m128i*)rows[1]+0+0+1*12));//+NE
					L1preds[4*6+1]=_mm_add_epi16(cache[1], _mm_load_si128((__m128i*)rows[1]+0+1+1*12));
					L1preds[4*6+2]=_mm_add_epi16(cache[2], _mm_load_si128((__m128i*)rows[1]+0+2+1*12));
					L1preds[4*6+3]=_mm_add_epi16(cache[3], _mm_load_si128((__m128i*)rows[1]+0+3+1*12));//+NE
					L1preds[4*6+4]=_mm_add_epi16(cache[4], _mm_load_si128((__m128i*)rows[1]+0+4+1*12));
					L1preds[4*6+5]=_mm_add_epi16(cache[5], _mm_load_si128((__m128i*)rows[1]+0+5+1*12));

					//N+W-NW
					L1preds[5*6+0]=predY[0];
					L1preds[5*6+1]=predY[1];
					L1preds[5*6+2]=predU[0];
					L1preds[5*6+3]=predU[1];
					L1preds[5*6+4]=predV[0];
					L1preds[5*6+5]=predV[1];

					//N+NE-NNE
					cache[0]=_mm_add_epi16(N[0], _mm_load_si128((__m128i*)rows[1]+0+0+1*12));//N+NE
					cache[1]=_mm_add_epi16(N[1], _mm_load_si128((__m128i*)rows[1]+0+1+1*12));
					cache[2]=_mm_add_epi16(N[2], _mm_load_si128((__m128i*)rows[1]+0+2+1*12));
					cache[3]=_mm_add_epi16(N[3], _mm_load_si128((__m128i*)rows[1]+0+3+1*12));
					cache[4]=_mm_add_epi16(N[4], _mm_load_si128((__m128i*)rows[1]+0+4+1*12));
					cache[5]=_mm_add_epi16(N[5], _mm_load_si128((__m128i*)rows[1]+0+5+1*12));
					L1preds[6*6+0]=_mm_sub_epi16(cache[0], _mm_load_si128((__m128i*)rows[2]+0+0+1*12));//NNE
					L1preds[6*6+1]=_mm_sub_epi16(cache[1], _mm_load_si128((__m128i*)rows[2]+0+1+1*12));
					L1preds[6*6+2]=_mm_sub_epi16(cache[2], _mm_load_si128((__m128i*)rows[2]+0+2+1*12));
					L1preds[6*6+3]=_mm_sub_epi16(cache[3], _mm_load_si128((__m128i*)rows[2]+0+3+1*12));
					L1preds[6*6+4]=_mm_sub_epi16(cache[4], _mm_load_si128((__m128i*)rows[2]+0+4+1*12));
					L1preds[6*6+5]=_mm_sub_epi16(cache[5], _mm_load_si128((__m128i*)rows[2]+0+5+1*12));

					//(WWWW+WWW+NNN+NEE+NEEE+NEEEE-(NW+N))>>2
					cache[0]=_mm_load_si128((__m128i*)rows[0]+0+0-4*12);//WWWW
					cache[1]=_mm_load_si128((__m128i*)rows[0]+0+1-4*12);
					cache[2]=_mm_load_si128((__m128i*)rows[0]+0+2-4*12);
					cache[3]=_mm_load_si128((__m128i*)rows[0]+0+3-4*12);
					cache[4]=_mm_load_si128((__m128i*)rows[0]+0+4-4*12);
					cache[5]=_mm_load_si128((__m128i*)rows[0]+0+5-4*12);
					cache[0]=_mm_add_epi16(cache[0], _mm_load_si128((__m128i*)rows[0]+0+0-3*12));//+WWW
					cache[1]=_mm_add_epi16(cache[1], _mm_load_si128((__m128i*)rows[0]+0+1-3*12));
					cache[2]=_mm_add_epi16(cache[2], _mm_load_si128((__m128i*)rows[0]+0+2-3*12));
					cache[3]=_mm_add_epi16(cache[3], _mm_load_si128((__m128i*)rows[0]+0+3-3*12));
					cache[4]=_mm_add_epi16(cache[4], _mm_load_si128((__m128i*)rows[0]+0+4-3*12));
					cache[5]=_mm_add_epi16(cache[5], _mm_load_si128((__m128i*)rows[0]+0+5-3*12));
					cache[0]=_mm_add_epi16(cache[0], _mm_load_si128((__m128i*)rows[3]+0+0+0*12));//+NNN
					cache[1]=_mm_add_epi16(cache[1], _mm_load_si128((__m128i*)rows[3]+0+1+0*12));
					cache[2]=_mm_add_epi16(cache[2], _mm_load_si128((__m128i*)rows[3]+0+2+0*12));
					cache[3]=_mm_add_epi16(cache[3], _mm_load_si128((__m128i*)rows[3]+0+3+0*12));
					cache[4]=_mm_add_epi16(cache[4], _mm_load_si128((__m128i*)rows[3]+0+4+0*12));
					cache[5]=_mm_add_epi16(cache[5], _mm_load_si128((__m128i*)rows[3]+0+5+0*12));
					cache[0]=_mm_add_epi16(cache[0], _mm_load_si128((__m128i*)rows[1]+0+0+2*12));//+NEE
					cache[1]=_mm_add_epi16(cache[1], _mm_load_si128((__m128i*)rows[1]+0+1+2*12));
					cache[2]=_mm_add_epi16(cache[2], _mm_load_si128((__m128i*)rows[1]+0+2+2*12));
					cache[3]=_mm_add_epi16(cache[3], _mm_load_si128((__m128i*)rows[1]+0+3+2*12));
					cache[4]=_mm_add_epi16(cache[4], _mm_load_si128((__m128i*)rows[1]+0+4+2*12));
					cache[5]=_mm_add_epi16(cache[5], _mm_load_si128((__m128i*)rows[1]+0+5+2*12));
					cache[0]=_mm_add_epi16(cache[0], _mm_load_si128((__m128i*)rows[1]+0+0+3*12));//+NEEE
					cache[1]=_mm_add_epi16(cache[1], _mm_load_si128((__m128i*)rows[1]+0+1+3*12));
					cache[2]=_mm_add_epi16(cache[2], _mm_load_si128((__m128i*)rows[1]+0+2+3*12));
					cache[3]=_mm_add_epi16(cache[3], _mm_load_si128((__m128i*)rows[1]+0+3+3*12));
					cache[4]=_mm_add_epi16(cache[4], _mm_load_si128((__m128i*)rows[1]+0+4+3*12));
					cache[5]=_mm_add_epi16(cache[5], _mm_load_si128((__m128i*)rows[1]+0+5+3*12));
					cache[0]=_mm_add_epi16(cache[0], _mm_load_si128((__m128i*)rows[1]+0+0+4*12));//+NEEEE
					cache[1]=_mm_add_epi16(cache[1], _mm_load_si128((__m128i*)rows[1]+0+1+4*12));
					cache[2]=_mm_add_epi16(cache[2], _mm_load_si128((__m128i*)rows[1]+0+2+4*12));
					cache[3]=_mm_add_epi16(cache[3], _mm_load_si128((__m128i*)rows[1]+0+3+4*12));
					cache[4]=_mm_add_epi16(cache[4], _mm_load_si128((__m128i*)rows[1]+0+4+4*12));
					cache[5]=_mm_add_epi16(cache[5], _mm_load_si128((__m128i*)rows[1]+0+5+4*12));
					cache[0]=_mm_sub_epi16(cache[0], _mm_add_epi16(N[0], NW[0]));
					cache[1]=_mm_sub_epi16(cache[1], _mm_add_epi16(N[1], NW[1]));
					cache[2]=_mm_sub_epi16(cache[2], _mm_add_epi16(N[2], NW[2]));
					cache[3]=_mm_sub_epi16(cache[3], _mm_add_epi16(N[3], NW[3]));
					cache[4]=_mm_sub_epi16(cache[4], _mm_add_epi16(N[4], NW[4]));
					cache[5]=_mm_sub_epi16(cache[5], _mm_add_epi16(N[5], NW[5]));
					L1preds[7*6+0]=_mm_srai_epi16(cache[0], 2);
					L1preds[7*6+1]=_mm_srai_epi16(cache[1], 2);
					L1preds[7*6+2]=_mm_srai_epi16(cache[2], 2);
					L1preds[7*6+3]=_mm_srai_epi16(cache[3], 2);
					L1preds[7*6+4]=_mm_srai_epi16(cache[4], 2);
					L1preds[7*6+5]=_mm_srai_epi16(cache[5], 2);
					

					//mix
					mp[0x0]=_mm_setzero_si128();
					mp[0x1]=_mm_setzero_si128();
					mp[0x2]=_mm_setzero_si128();
					mp[0x3]=_mm_setzero_si128();
					mp[0x4]=_mm_setzero_si128();
					mp[0x5]=_mm_setzero_si128();
					mp[0x6]=_mm_setzero_si128();
					mp[0x7]=_mm_setzero_si128();
					mp[0x8]=_mm_setzero_si128();
					mp[0x9]=_mm_setzero_si128();
					mp[0xA]=_mm_setzero_si128();
					mp[0xB]=_mm_setzero_si128();
					for(kp=0;kp<L1_NPREDS2;++kp)
					{
						//signed 16 -> 32	6 lo 6 hi registers
						t[0x0]=_mm_slli_epi32(L1preds[kp*6+0], 16);//y lo
						t[0x1]=_mm_slli_epi32(L1preds[kp*6+1], 16);
						t[0x2]=_mm_slli_epi32(L1preds[kp*6+2], 16);//u
						t[0x3]=_mm_slli_epi32(L1preds[kp*6+3], 16);
						t[0x4]=_mm_slli_epi32(L1preds[kp*6+4], 16);//v
						t[0x5]=_mm_slli_epi32(L1preds[kp*6+5], 16);

						t[0x6]=_mm_srai_epi32(L1preds[kp*6+0], 16);//y hi
						t[0x7]=_mm_srai_epi32(L1preds[kp*6+1], 16);
						t[0x8]=_mm_srai_epi32(L1preds[kp*6+2], 16);//u
						t[0x9]=_mm_srai_epi32(L1preds[kp*6+3], 16);
						t[0xA]=_mm_srai_epi32(L1preds[kp*6+4], 16);//v
						t[0xB]=_mm_srai_epi32(L1preds[kp*6+5], 16);
						t[0x0]=_mm_srai_epi32(t[0x0], 16);
						t[0x1]=_mm_srai_epi32(t[0x1], 16);
						t[0x2]=_mm_srai_epi32(t[0x2], 16);
						t[0x3]=_mm_srai_epi32(t[0x3], 16);
						t[0x4]=_mm_srai_epi32(t[0x4], 16);
						t[0x5]=_mm_srai_epi32(t[0x5], 16);
						t[0x0]=_mm_mullo_epi32(t[0x0], _mm_load_si128((__m128i*)L1weights+kp*12+0x0));
						t[0x1]=_mm_mullo_epi32(t[0x1], _mm_load_si128((__m128i*)L1weights+kp*12+0x1));
						t[0x2]=_mm_mullo_epi32(t[0x2], _mm_load_si128((__m128i*)L1weights+kp*12+0x2));
						t[0x3]=_mm_mullo_epi32(t[0x3], _mm_load_si128((__m128i*)L1weights+kp*12+0x3));
						t[0x4]=_mm_mullo_epi32(t[0x4], _mm_load_si128((__m128i*)L1weights+kp*12+0x4));
						t[0x5]=_mm_mullo_epi32(t[0x5], _mm_load_si128((__m128i*)L1weights+kp*12+0x5));
						t[0x6]=_mm_mullo_epi32(t[0x6], _mm_load_si128((__m128i*)L1weights+kp*12+0x6));
						t[0x7]=_mm_mullo_epi32(t[0x7], _mm_load_si128((__m128i*)L1weights+kp*12+0x7));
						t[0x8]=_mm_mullo_epi32(t[0x8], _mm_load_si128((__m128i*)L1weights+kp*12+0x8));
						t[0x9]=_mm_mullo_epi32(t[0x9], _mm_load_si128((__m128i*)L1weights+kp*12+0x9));
						t[0xA]=_mm_mullo_epi32(t[0xA], _mm_load_si128((__m128i*)L1weights+kp*12+0xA));
						t[0xB]=_mm_mullo_epi32(t[0xB], _mm_load_si128((__m128i*)L1weights+kp*12+0xB));
						mp[0x0]=_mm_add_epi32(mp[0x0], t[0x0]);
						mp[0x1]=_mm_add_epi32(mp[0x1], t[0x1]);
						mp[0x2]=_mm_add_epi32(mp[0x2], t[0x2]);
						mp[0x3]=_mm_add_epi32(mp[0x3], t[0x3]);
						mp[0x4]=_mm_add_epi32(mp[0x4], t[0x4]);
						mp[0x5]=_mm_add_epi32(mp[0x5], t[0x5]);
						mp[0x6]=_mm_add_epi32(mp[0x6], t[0x6]);
						mp[0x7]=_mm_add_epi32(mp[0x7], t[0x7]);
						mp[0x8]=_mm_add_epi32(mp[0x8], t[0x8]);
						mp[0x9]=_mm_add_epi32(mp[0x9], t[0x9]);
						mp[0xA]=_mm_add_epi32(mp[0xA], t[0xA]);
						mp[0xB]=_mm_add_epi32(mp[0xB], t[0xB]);
					}
					{
						__m128i rcon=_mm_set1_epi32(1<<L1_SH2>>1);
						mp[0x0]=_mm_add_epi32(mp[0x0], rcon);//rounding to nearest
						mp[0x1]=_mm_add_epi32(mp[0x1], rcon);
						mp[0x2]=_mm_add_epi32(mp[0x2], rcon);
						mp[0x3]=_mm_add_epi32(mp[0x3], rcon);
						mp[0x4]=_mm_add_epi32(mp[0x4], rcon);
						mp[0x5]=_mm_add_epi32(mp[0x5], rcon);
						mp[0x6]=_mm_add_epi32(mp[0x6], rcon);
						mp[0x7]=_mm_add_epi32(mp[0x7], rcon);
						mp[0x8]=_mm_add_epi32(mp[0x8], rcon);
						mp[0x9]=_mm_add_epi32(mp[0x9], rcon);
						mp[0xA]=_mm_add_epi32(mp[0xA], rcon);
						mp[0xB]=_mm_add_epi32(mp[0xB], rcon);
					}

					mp[0x0]=_mm_srai_epi32(mp[0x0], L1_SH2);
					mp[0x1]=_mm_srai_epi32(mp[0x1], L1_SH2);
					mp[0x2]=_mm_srai_epi32(mp[0x2], L1_SH2);
					mp[0x3]=_mm_srai_epi32(mp[0x3], L1_SH2);
					mp[0x4]=_mm_srai_epi32(mp[0x4], L1_SH2);
					mp[0x5]=_mm_srai_epi32(mp[0x5], L1_SH2);
					mp[0x6]=_mm_srai_epi32(mp[0x6], L1_SH2-16);//y hi
					mp[0x7]=_mm_srai_epi32(mp[0x7], L1_SH2-16);
					mp[0x8]=_mm_srai_epi32(mp[0x8], L1_SH2-16);//u hi
					mp[0x9]=_mm_srai_epi32(mp[0x9], L1_SH2-16);
					mp[0xA]=_mm_srai_epi32(mp[0xA], L1_SH2-16);//v hi
					mp[0xB]=_mm_srai_epi32(mp[0xB], L1_SH2-16);
					//32 -> 16
					predY[0]=_mm_blend_epi16(mp[0x0], mp[0x6], 0xAA);
					predY[1]=_mm_blend_epi16(mp[0x1], mp[0x7], 0xAA);
					predU[0]=_mm_blend_epi16(mp[0x2], mp[0x8], 0xAA);
					predU[1]=_mm_blend_epi16(mp[0x3], mp[0x9], 0xAA);
					predV[0]=_mm_blend_epi16(mp[0x4], mp[0xA], 0xAA);
					predV[1]=_mm_blend_epi16(mp[0x5], mp[0xB], 0xAA);


					//loosen pred range
					if(!cond_cg)
					{
						t[0]=_mm_load_si128((__m128i*)rows[1]+0+0+1*12);//NE
						t[1]=_mm_load_si128((__m128i*)rows[1]+0+1+1*12);
						t[2]=_mm_load_si128((__m128i*)rows[1]+0+2+1*12);
						t[3]=_mm_load_si128((__m128i*)rows[1]+0+3+1*12);
						t[4]=_mm_load_si128((__m128i*)rows[1]+0+4+1*12);
						t[5]=_mm_load_si128((__m128i*)rows[1]+0+5+1*12);
						ymin[0]=_mm_min_epi16(ymin[0], t[0*2+0]); ymin[1]=_mm_min_epi16(ymin[1], t[0*2+1]);
						ymax[0]=_mm_max_epi16(ymax[0], t[0*2+0]); ymax[1]=_mm_max_epi16(ymax[1], t[0*2+1]);
						umin[0]=_mm_min_epi16(umin[0], t[1*2+0]); umin[1]=_mm_min_epi16(umin[1], t[1*2+1]);
						umax[0]=_mm_max_epi16(umax[0], t[1*2+0]); umax[1]=_mm_max_epi16(umax[1], t[1*2+1]);
						vmin[0]=_mm_min_epi16(vmin[0], t[2*2+0]); vmin[1]=_mm_min_epi16(vmin[1], t[2*2+1]);
						vmax[0]=_mm_max_epi16(vmax[0], t[2*2+0]); vmax[1]=_mm_max_epi16(vmax[1], t[2*2+1]);
						t[0]=_mm_load_si128((__m128i*)rows[1]+0+0+3*12);//NEEE
						t[1]=_mm_load_si128((__m128i*)rows[1]+0+1+3*12);
						t[2]=_mm_load_si128((__m128i*)rows[1]+0+2+3*12);
						t[3]=_mm_load_si128((__m128i*)rows[1]+0+3+3*12);
						t[4]=_mm_load_si128((__m128i*)rows[1]+0+4+3*12);
						t[5]=_mm_load_si128((__m128i*)rows[1]+0+5+3*12);
						ymin[0]=_mm_min_epi16(ymin[0], t[0*2+0]); ymin[1]=_mm_min_epi16(ymin[1], t[0*2+1]);
						ymax[0]=_mm_max_epi16(ymax[0], t[0*2+0]); ymax[1]=_mm_max_epi16(ymax[1], t[0*2+1]);
						umin[0]=_mm_min_epi16(umin[0], t[1*2+0]); umin[1]=_mm_min_epi16(umin[1], t[1*2+1]);
						umax[0]=_mm_max_epi16(umax[0], t[1*2+0]); umax[1]=_mm_max_epi16(umax[1], t[1*2+1]);
						vmin[0]=_mm_min_epi16(vmin[0], t[2*2+0]); vmin[1]=_mm_min_epi16(vmin[1], t[2*2+1]);
						vmax[0]=_mm_max_epi16(vmax[0], t[2*2+0]); vmax[1]=_mm_max_epi16(vmax[1], t[2*2+1]);
					}
				}
				else if(effort==3)
				{
					__m128i cache[6], mp[12], t[12];
					int kp;
					/*
					effort 3
					0	N
					1	W
					2	NNN
					3	WWW
					4	NEEE
					5	3*(N-NN)+NNN
					6	3*(W-WW)+WWW
					7	W+NE-N
					8	N+W-NW
					9	N+NE-NNE
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

					//NNN
					L1preds[2*6+0]=_mm_load_si128((__m128i*)rows[3]+0+0+0*12);
					L1preds[2*6+1]=_mm_load_si128((__m128i*)rows[3]+0+1+0*12);
					L1preds[2*6+2]=_mm_load_si128((__m128i*)rows[3]+0+2+0*12);
					L1preds[2*6+3]=_mm_load_si128((__m128i*)rows[3]+0+3+0*12);
					L1preds[2*6+4]=_mm_load_si128((__m128i*)rows[3]+0+4+0*12);
					L1preds[2*6+5]=_mm_load_si128((__m128i*)rows[3]+0+5+0*12);

					//WWW
					L1preds[3*6+0]=_mm_load_si128((__m128i*)rows[0]+0+0-3*12);
					L1preds[3*6+1]=_mm_load_si128((__m128i*)rows[0]+0+1-3*12);
					L1preds[3*6+2]=_mm_load_si128((__m128i*)rows[0]+0+2-3*12);
					L1preds[3*6+3]=_mm_load_si128((__m128i*)rows[0]+0+3-3*12);
					L1preds[3*6+4]=_mm_load_si128((__m128i*)rows[0]+0+4-3*12);
					L1preds[3*6+5]=_mm_load_si128((__m128i*)rows[0]+0+5-3*12);

					//NEEE
					L1preds[4*6+0]=_mm_load_si128((__m128i*)rows[1]+0+0+3*12);
					L1preds[4*6+1]=_mm_load_si128((__m128i*)rows[1]+0+1+3*12);
					L1preds[4*6+2]=_mm_load_si128((__m128i*)rows[1]+0+2+3*12);
					L1preds[4*6+3]=_mm_load_si128((__m128i*)rows[1]+0+3+3*12);
					L1preds[4*6+4]=_mm_load_si128((__m128i*)rows[1]+0+4+3*12);
					L1preds[4*6+5]=_mm_load_si128((__m128i*)rows[1]+0+5+3*12);

					//3*(N-NN)+NNN
					cache[0]=_mm_sub_epi16(N[0], _mm_load_si128((__m128i*)rows[2]+0+0+0*12));//N-NN
					cache[1]=_mm_sub_epi16(N[1], _mm_load_si128((__m128i*)rows[2]+0+1+0*12));
					cache[2]=_mm_sub_epi16(N[2], _mm_load_si128((__m128i*)rows[2]+0+2+0*12));
					cache[3]=_mm_sub_epi16(N[3], _mm_load_si128((__m128i*)rows[2]+0+3+0*12));
					cache[4]=_mm_sub_epi16(N[4], _mm_load_si128((__m128i*)rows[2]+0+4+0*12));
					cache[5]=_mm_sub_epi16(N[5], _mm_load_si128((__m128i*)rows[2]+0+5+0*12));
					cache[0]=_mm_add_epi16(cache[0], _mm_slli_epi16(cache[0], 1));//*3
					cache[1]=_mm_add_epi16(cache[1], _mm_slli_epi16(cache[1], 1));
					cache[2]=_mm_add_epi16(cache[2], _mm_slli_epi16(cache[2], 1));
					cache[3]=_mm_add_epi16(cache[3], _mm_slli_epi16(cache[3], 1));
					cache[4]=_mm_add_epi16(cache[4], _mm_slli_epi16(cache[4], 1));
					cache[5]=_mm_add_epi16(cache[5], _mm_slli_epi16(cache[5], 1));
					L1preds[5*6+0]=_mm_add_epi16(cache[0], _mm_load_si128((__m128i*)rows[3]+0+0+0*12));//+NNN
					L1preds[5*6+1]=_mm_add_epi16(cache[1], _mm_load_si128((__m128i*)rows[3]+0+1+0*12));
					L1preds[5*6+2]=_mm_add_epi16(cache[2], _mm_load_si128((__m128i*)rows[3]+0+2+0*12));
					L1preds[5*6+3]=_mm_add_epi16(cache[3], _mm_load_si128((__m128i*)rows[3]+0+3+0*12));
					L1preds[5*6+4]=_mm_add_epi16(cache[4], _mm_load_si128((__m128i*)rows[3]+0+4+0*12));
					L1preds[5*6+5]=_mm_add_epi16(cache[5], _mm_load_si128((__m128i*)rows[3]+0+5+0*12));

					//3*(W-WW)+WWW
					cache[0]=_mm_sub_epi16(W[0], _mm_load_si128((__m128i*)rows[0]+0+0-2*12));//W-WW
					cache[1]=_mm_sub_epi16(W[1], _mm_load_si128((__m128i*)rows[0]+0+1-2*12));
					cache[2]=_mm_sub_epi16(W[2], _mm_load_si128((__m128i*)rows[0]+0+2-2*12));
					cache[3]=_mm_sub_epi16(W[3], _mm_load_si128((__m128i*)rows[0]+0+3-2*12));
					cache[4]=_mm_sub_epi16(W[4], _mm_load_si128((__m128i*)rows[0]+0+4-2*12));
					cache[5]=_mm_sub_epi16(W[5], _mm_load_si128((__m128i*)rows[0]+0+5-2*12));
					cache[0]=_mm_add_epi16(cache[0], _mm_slli_epi16(cache[0], 1));//*3
					cache[1]=_mm_add_epi16(cache[1], _mm_slli_epi16(cache[1], 1));
					cache[2]=_mm_add_epi16(cache[2], _mm_slli_epi16(cache[2], 1));
					cache[3]=_mm_add_epi16(cache[3], _mm_slli_epi16(cache[3], 1));
					cache[4]=_mm_add_epi16(cache[4], _mm_slli_epi16(cache[4], 1));
					cache[5]=_mm_add_epi16(cache[5], _mm_slli_epi16(cache[5], 1));
					L1preds[6*6+0]=_mm_add_epi16(cache[0], _mm_load_si128((__m128i*)rows[0]+0+0-3*12));//+WWW
					L1preds[6*6+1]=_mm_add_epi16(cache[1], _mm_load_si128((__m128i*)rows[0]+0+1-3*12));
					L1preds[6*6+2]=_mm_add_epi16(cache[2], _mm_load_si128((__m128i*)rows[0]+0+2-3*12));
					L1preds[6*6+3]=_mm_add_epi16(cache[3], _mm_load_si128((__m128i*)rows[0]+0+3-3*12));
					L1preds[6*6+4]=_mm_add_epi16(cache[4], _mm_load_si128((__m128i*)rows[0]+0+4-3*12));
					L1preds[6*6+5]=_mm_add_epi16(cache[5], _mm_load_si128((__m128i*)rows[0]+0+5-3*12));

					//W+NE-N
					cache[0]=_mm_sub_epi16(W[0], N[0]);
					cache[1]=_mm_sub_epi16(W[1], N[1]);
					cache[2]=_mm_sub_epi16(W[2], N[2]);
					cache[3]=_mm_sub_epi16(W[3], N[3]);
					cache[4]=_mm_sub_epi16(W[4], N[4]);
					cache[5]=_mm_sub_epi16(W[5], N[5]);
					L1preds[7*6+0]=_mm_add_epi16(cache[0], _mm_load_si128((__m128i*)rows[1]+0+0+1*12));//+NE
					L1preds[7*6+1]=_mm_add_epi16(cache[1], _mm_load_si128((__m128i*)rows[1]+0+1+1*12));
					L1preds[7*6+2]=_mm_add_epi16(cache[2], _mm_load_si128((__m128i*)rows[1]+0+2+1*12));
					L1preds[7*6+3]=_mm_add_epi16(cache[3], _mm_load_si128((__m128i*)rows[1]+0+3+1*12));
					L1preds[7*6+4]=_mm_add_epi16(cache[4], _mm_load_si128((__m128i*)rows[1]+0+4+1*12));
					L1preds[7*6+5]=_mm_add_epi16(cache[5], _mm_load_si128((__m128i*)rows[1]+0+5+1*12));

					//N+W-NW
					L1preds[8*6+0]=predY[0];
					L1preds[8*6+1]=predY[1];
					L1preds[8*6+2]=predU[0];
					L1preds[8*6+3]=predU[1];
					L1preds[8*6+4]=predV[0];
					L1preds[8*6+5]=predV[1];

					//N+NE-NNE
					cache[0]=_mm_add_epi16(N[0], _mm_load_si128((__m128i*)rows[1]+0+0+1*12));//N+NE
					cache[1]=_mm_add_epi16(N[1], _mm_load_si128((__m128i*)rows[1]+0+1+1*12));
					cache[2]=_mm_add_epi16(N[2], _mm_load_si128((__m128i*)rows[1]+0+2+1*12));
					cache[3]=_mm_add_epi16(N[3], _mm_load_si128((__m128i*)rows[1]+0+3+1*12));
					cache[4]=_mm_add_epi16(N[4], _mm_load_si128((__m128i*)rows[1]+0+4+1*12));
					cache[5]=_mm_add_epi16(N[5], _mm_load_si128((__m128i*)rows[1]+0+5+1*12));
					L1preds[9*6+0]=_mm_sub_epi16(cache[0], _mm_load_si128((__m128i*)rows[2]+0+0+1*12));//NNE
					L1preds[9*6+1]=_mm_sub_epi16(cache[1], _mm_load_si128((__m128i*)rows[2]+0+1+1*12));
					L1preds[9*6+2]=_mm_sub_epi16(cache[2], _mm_load_si128((__m128i*)rows[2]+0+2+1*12));
					L1preds[9*6+3]=_mm_sub_epi16(cache[3], _mm_load_si128((__m128i*)rows[2]+0+3+1*12));
					L1preds[9*6+4]=_mm_sub_epi16(cache[4], _mm_load_si128((__m128i*)rows[2]+0+4+1*12));
					L1preds[9*6+5]=_mm_sub_epi16(cache[5], _mm_load_si128((__m128i*)rows[2]+0+5+1*12));
					
					//(WWWW+WWW+NNN+NNEE+NEEE+NEEEE-(N+W))>>2
					cache[0]=_mm_load_si128((__m128i*)rows[0]+0+0-4*12);//WWWW
					cache[1]=_mm_load_si128((__m128i*)rows[0]+0+1-4*12);
					cache[2]=_mm_load_si128((__m128i*)rows[0]+0+2-4*12);
					cache[3]=_mm_load_si128((__m128i*)rows[0]+0+3-4*12);
					cache[4]=_mm_load_si128((__m128i*)rows[0]+0+4-4*12);
					cache[5]=_mm_load_si128((__m128i*)rows[0]+0+5-4*12);
					cache[0]=_mm_add_epi16(cache[0], _mm_load_si128((__m128i*)rows[0]+0+0-3*12));//+WWW
					cache[1]=_mm_add_epi16(cache[1], _mm_load_si128((__m128i*)rows[0]+0+1-3*12));
					cache[2]=_mm_add_epi16(cache[2], _mm_load_si128((__m128i*)rows[0]+0+2-3*12));
					cache[3]=_mm_add_epi16(cache[3], _mm_load_si128((__m128i*)rows[0]+0+3-3*12));
					cache[4]=_mm_add_epi16(cache[4], _mm_load_si128((__m128i*)rows[0]+0+4-3*12));
					cache[5]=_mm_add_epi16(cache[5], _mm_load_si128((__m128i*)rows[0]+0+5-3*12));
					cache[0]=_mm_add_epi16(cache[0], _mm_load_si128((__m128i*)rows[3]+0+0+0*12));//+NNN
					cache[1]=_mm_add_epi16(cache[1], _mm_load_si128((__m128i*)rows[3]+0+1+0*12));
					cache[2]=_mm_add_epi16(cache[2], _mm_load_si128((__m128i*)rows[3]+0+2+0*12));
					cache[3]=_mm_add_epi16(cache[3], _mm_load_si128((__m128i*)rows[3]+0+3+0*12));
					cache[4]=_mm_add_epi16(cache[4], _mm_load_si128((__m128i*)rows[3]+0+4+0*12));
					cache[5]=_mm_add_epi16(cache[5], _mm_load_si128((__m128i*)rows[3]+0+5+0*12));
					cache[0]=_mm_add_epi16(cache[0], _mm_load_si128((__m128i*)rows[2]+0+0+2*12));//+NNEE
					cache[1]=_mm_add_epi16(cache[1], _mm_load_si128((__m128i*)rows[2]+0+1+2*12));
					cache[2]=_mm_add_epi16(cache[2], _mm_load_si128((__m128i*)rows[2]+0+2+2*12));
					cache[3]=_mm_add_epi16(cache[3], _mm_load_si128((__m128i*)rows[2]+0+3+2*12));
					cache[4]=_mm_add_epi16(cache[4], _mm_load_si128((__m128i*)rows[2]+0+4+2*12));
					cache[5]=_mm_add_epi16(cache[5], _mm_load_si128((__m128i*)rows[2]+0+5+2*12));
					cache[0]=_mm_add_epi16(cache[0], _mm_load_si128((__m128i*)rows[1]+0+0+3*12));//+NEEE
					cache[1]=_mm_add_epi16(cache[1], _mm_load_si128((__m128i*)rows[1]+0+1+3*12));
					cache[2]=_mm_add_epi16(cache[2], _mm_load_si128((__m128i*)rows[1]+0+2+3*12));
					cache[3]=_mm_add_epi16(cache[3], _mm_load_si128((__m128i*)rows[1]+0+3+3*12));
					cache[4]=_mm_add_epi16(cache[4], _mm_load_si128((__m128i*)rows[1]+0+4+3*12));
					cache[5]=_mm_add_epi16(cache[5], _mm_load_si128((__m128i*)rows[1]+0+5+3*12));
					cache[0]=_mm_add_epi16(cache[0], _mm_load_si128((__m128i*)rows[1]+0+0+4*12));//+NEEEE
					cache[1]=_mm_add_epi16(cache[1], _mm_load_si128((__m128i*)rows[1]+0+1+4*12));
					cache[2]=_mm_add_epi16(cache[2], _mm_load_si128((__m128i*)rows[1]+0+2+4*12));
					cache[3]=_mm_add_epi16(cache[3], _mm_load_si128((__m128i*)rows[1]+0+3+4*12));
					cache[4]=_mm_add_epi16(cache[4], _mm_load_si128((__m128i*)rows[1]+0+4+4*12));
					cache[5]=_mm_add_epi16(cache[5], _mm_load_si128((__m128i*)rows[1]+0+5+4*12));
					cache[0]=_mm_sub_epi16(cache[0], _mm_add_epi16(N[0], W[0]));
					cache[1]=_mm_sub_epi16(cache[1], _mm_add_epi16(N[1], W[1]));
					cache[2]=_mm_sub_epi16(cache[2], _mm_add_epi16(N[2], W[2]));
					cache[3]=_mm_sub_epi16(cache[3], _mm_add_epi16(N[3], W[3]));
					cache[4]=_mm_sub_epi16(cache[4], _mm_add_epi16(N[4], W[4]));
					cache[5]=_mm_sub_epi16(cache[5], _mm_add_epi16(N[5], W[5]));
					L1preds[10*6+0]=_mm_srai_epi16(cache[0], 2);
					L1preds[10*6+1]=_mm_srai_epi16(cache[1], 2);
					L1preds[10*6+2]=_mm_srai_epi16(cache[2], 2);
					L1preds[10*6+3]=_mm_srai_epi16(cache[3], 2);
					L1preds[10*6+4]=_mm_srai_epi16(cache[4], 2);
					L1preds[10*6+5]=_mm_srai_epi16(cache[5], 2);


					//mix
					mp[0x0]=_mm_setzero_si128();
					mp[0x1]=_mm_setzero_si128();
					mp[0x2]=_mm_setzero_si128();
					mp[0x3]=_mm_setzero_si128();
					mp[0x4]=_mm_setzero_si128();
					mp[0x5]=_mm_setzero_si128();
					mp[0x6]=_mm_setzero_si128();
					mp[0x7]=_mm_setzero_si128();
					mp[0x8]=_mm_setzero_si128();
					mp[0x9]=_mm_setzero_si128();
					mp[0xA]=_mm_setzero_si128();
					mp[0xB]=_mm_setzero_si128();
					for(kp=0;kp<L1_NPREDS3;++kp)
					{
						//signed 16 -> 32	6 lo 6 hi registers
						t[0x0]=_mm_slli_epi32(L1preds[kp*6+0], 16);//y lo
						t[0x1]=_mm_slli_epi32(L1preds[kp*6+1], 16);
						t[0x2]=_mm_slli_epi32(L1preds[kp*6+2], 16);//u
						t[0x3]=_mm_slli_epi32(L1preds[kp*6+3], 16);
						t[0x4]=_mm_slli_epi32(L1preds[kp*6+4], 16);//v
						t[0x5]=_mm_slli_epi32(L1preds[kp*6+5], 16);

						t[0x6]=_mm_srai_epi32(L1preds[kp*6+0], 16);//y hi
						t[0x7]=_mm_srai_epi32(L1preds[kp*6+1], 16);
						t[0x8]=_mm_srai_epi32(L1preds[kp*6+2], 16);//u
						t[0x9]=_mm_srai_epi32(L1preds[kp*6+3], 16);
						t[0xA]=_mm_srai_epi32(L1preds[kp*6+4], 16);//v
						t[0xB]=_mm_srai_epi32(L1preds[kp*6+5], 16);
						t[0x0]=_mm_srai_epi32(t[0x0], 16);
						t[0x1]=_mm_srai_epi32(t[0x1], 16);
						t[0x2]=_mm_srai_epi32(t[0x2], 16);
						t[0x3]=_mm_srai_epi32(t[0x3], 16);
						t[0x4]=_mm_srai_epi32(t[0x4], 16);
						t[0x5]=_mm_srai_epi32(t[0x5], 16);
						t[0x0]=_mm_mullo_epi32(t[0x0], _mm_load_si128((__m128i*)L1weights+kp*12+0x0));
						t[0x1]=_mm_mullo_epi32(t[0x1], _mm_load_si128((__m128i*)L1weights+kp*12+0x1));
						t[0x2]=_mm_mullo_epi32(t[0x2], _mm_load_si128((__m128i*)L1weights+kp*12+0x2));
						t[0x3]=_mm_mullo_epi32(t[0x3], _mm_load_si128((__m128i*)L1weights+kp*12+0x3));
						t[0x4]=_mm_mullo_epi32(t[0x4], _mm_load_si128((__m128i*)L1weights+kp*12+0x4));
						t[0x5]=_mm_mullo_epi32(t[0x5], _mm_load_si128((__m128i*)L1weights+kp*12+0x5));
						t[0x6]=_mm_mullo_epi32(t[0x6], _mm_load_si128((__m128i*)L1weights+kp*12+0x6));
						t[0x7]=_mm_mullo_epi32(t[0x7], _mm_load_si128((__m128i*)L1weights+kp*12+0x7));
						t[0x8]=_mm_mullo_epi32(t[0x8], _mm_load_si128((__m128i*)L1weights+kp*12+0x8));
						t[0x9]=_mm_mullo_epi32(t[0x9], _mm_load_si128((__m128i*)L1weights+kp*12+0x9));
						t[0xA]=_mm_mullo_epi32(t[0xA], _mm_load_si128((__m128i*)L1weights+kp*12+0xA));
						t[0xB]=_mm_mullo_epi32(t[0xB], _mm_load_si128((__m128i*)L1weights+kp*12+0xB));
						mp[0x0]=_mm_add_epi32(mp[0x0], t[0x0]);
						mp[0x1]=_mm_add_epi32(mp[0x1], t[0x1]);
						mp[0x2]=_mm_add_epi32(mp[0x2], t[0x2]);
						mp[0x3]=_mm_add_epi32(mp[0x3], t[0x3]);
						mp[0x4]=_mm_add_epi32(mp[0x4], t[0x4]);
						mp[0x5]=_mm_add_epi32(mp[0x5], t[0x5]);
						mp[0x6]=_mm_add_epi32(mp[0x6], t[0x6]);
						mp[0x7]=_mm_add_epi32(mp[0x7], t[0x7]);
						mp[0x8]=_mm_add_epi32(mp[0x8], t[0x8]);
						mp[0x9]=_mm_add_epi32(mp[0x9], t[0x9]);
						mp[0xA]=_mm_add_epi32(mp[0xA], t[0xA]);
						mp[0xB]=_mm_add_epi32(mp[0xB], t[0xB]);
					}
					{
						__m128i rcon=_mm_set1_epi32(1<<L1_SH3>>1);
						mp[0x0]=_mm_add_epi32(mp[0x0], rcon);//rounding to nearest
						mp[0x1]=_mm_add_epi32(mp[0x1], rcon);
						mp[0x2]=_mm_add_epi32(mp[0x2], rcon);
						mp[0x3]=_mm_add_epi32(mp[0x3], rcon);
						mp[0x4]=_mm_add_epi32(mp[0x4], rcon);
						mp[0x5]=_mm_add_epi32(mp[0x5], rcon);
						mp[0x6]=_mm_add_epi32(mp[0x6], rcon);
						mp[0x7]=_mm_add_epi32(mp[0x7], rcon);
						mp[0x8]=_mm_add_epi32(mp[0x8], rcon);
						mp[0x9]=_mm_add_epi32(mp[0x9], rcon);
						mp[0xA]=_mm_add_epi32(mp[0xA], rcon);
						mp[0xB]=_mm_add_epi32(mp[0xB], rcon);
					}

					mp[0x0]=_mm_srai_epi32(mp[0x0], L1_SH3);
					mp[0x1]=_mm_srai_epi32(mp[0x1], L1_SH3);
					mp[0x2]=_mm_srai_epi32(mp[0x2], L1_SH3);
					mp[0x3]=_mm_srai_epi32(mp[0x3], L1_SH3);
					mp[0x4]=_mm_srai_epi32(mp[0x4], L1_SH3);
					mp[0x5]=_mm_srai_epi32(mp[0x5], L1_SH3);
					mp[0x6]=_mm_srai_epi32(mp[0x6], L1_SH3-16);//y hi
					mp[0x7]=_mm_srai_epi32(mp[0x7], L1_SH3-16);
					mp[0x8]=_mm_srai_epi32(mp[0x8], L1_SH3-16);//u hi
					mp[0x9]=_mm_srai_epi32(mp[0x9], L1_SH3-16);
					mp[0xA]=_mm_srai_epi32(mp[0xA], L1_SH3-16);//v hi
					mp[0xB]=_mm_srai_epi32(mp[0xB], L1_SH3-16);
					//32 -> 16
					predY[0]=_mm_blend_epi16(mp[0x0], mp[0x6], 0xAA);
					predY[1]=_mm_blend_epi16(mp[0x1], mp[0x7], 0xAA);
					predU[0]=_mm_blend_epi16(mp[0x2], mp[0x8], 0xAA);
					predU[1]=_mm_blend_epi16(mp[0x3], mp[0x9], 0xAA);
					predV[0]=_mm_blend_epi16(mp[0x4], mp[0xA], 0xAA);
					predV[1]=_mm_blend_epi16(mp[0x5], mp[0xB], 0xAA);


					//loosen pred range
					if(!cond_cg)
					{
						t[0]=_mm_load_si128((__m128i*)rows[1]+0+0+1*12);//NE
						t[1]=_mm_load_si128((__m128i*)rows[1]+0+1+1*12);
						t[2]=_mm_load_si128((__m128i*)rows[1]+0+2+1*12);
						t[3]=_mm_load_si128((__m128i*)rows[1]+0+3+1*12);
						t[4]=_mm_load_si128((__m128i*)rows[1]+0+4+1*12);
						t[5]=_mm_load_si128((__m128i*)rows[1]+0+5+1*12);
						ymin[0]=_mm_min_epi16(ymin[0], t[0*2+0]); ymin[1]=_mm_min_epi16(ymin[1], t[0*2+1]);
						ymax[0]=_mm_max_epi16(ymax[0], t[0*2+0]); ymax[1]=_mm_max_epi16(ymax[1], t[0*2+1]);
						umin[0]=_mm_min_epi16(umin[0], t[1*2+0]); umin[1]=_mm_min_epi16(umin[1], t[1*2+1]);
						umax[0]=_mm_max_epi16(umax[0], t[1*2+0]); umax[1]=_mm_max_epi16(umax[1], t[1*2+1]);
						vmin[0]=_mm_min_epi16(vmin[0], t[2*2+0]); vmin[1]=_mm_min_epi16(vmin[1], t[2*2+1]);
						vmax[0]=_mm_max_epi16(vmax[0], t[2*2+0]); vmax[1]=_mm_max_epi16(vmax[1], t[2*2+1]);
						t[0]=_mm_load_si128((__m128i*)rows[1]+0+0+3*12);//NEEE
						t[1]=_mm_load_si128((__m128i*)rows[1]+0+1+3*12);
						t[2]=_mm_load_si128((__m128i*)rows[1]+0+2+3*12);
						t[3]=_mm_load_si128((__m128i*)rows[1]+0+3+3*12);
						t[4]=_mm_load_si128((__m128i*)rows[1]+0+4+3*12);
						t[5]=_mm_load_si128((__m128i*)rows[1]+0+5+3*12);
						ymin[0]=_mm_min_epi16(ymin[0], t[0*2+0]); ymin[1]=_mm_min_epi16(ymin[1], t[0*2+1]);
						ymax[0]=_mm_max_epi16(ymax[0], t[0*2+0]); ymax[1]=_mm_max_epi16(ymax[1], t[0*2+1]);
						umin[0]=_mm_min_epi16(umin[0], t[1*2+0]); umin[1]=_mm_min_epi16(umin[1], t[1*2+1]);
						umax[0]=_mm_max_epi16(umax[0], t[1*2+0]); umax[1]=_mm_max_epi16(umax[1], t[1*2+1]);
						vmin[0]=_mm_min_epi16(vmin[0], t[2*2+0]); vmin[1]=_mm_min_epi16(vmin[1], t[2*2+1]);
						vmax[0]=_mm_max_epi16(vmax[0], t[2*2+0]); vmax[1]=_mm_max_epi16(vmax[1], t[2*2+1]);
					}
				}
				predYUV0[0]=predY[0];
				predYUV0[1]=predY[1];
				predYUV0[2]=predU[0];
				predYUV0[3]=predU[1];
				predYUV0[4]=predV[0];
				predYUV0[5]=predV[1];
				if(cond_cg)
				{
					predY[0]=mcg[0];
					predY[1]=mcg[1];
					predU[0]=mcg[2];
					predU[1]=mcg[3];
					predV[0]=mcg[4];
					predV[1]=mcg[5];
				}

				predY[0]=_mm_max_epi16(predY[0], ymin[0]); predY[1]=_mm_max_epi16(predY[1], ymin[1]);
				predU[0]=_mm_max_epi16(predU[0], umin[0]); predU[1]=_mm_max_epi16(predU[1], umin[1]);
				predV[0]=_mm_max_epi16(predV[0], vmin[0]); predV[1]=_mm_max_epi16(predV[1], vmin[1]);
				predY[0]=_mm_min_epi16(predY[0], ymax[0]); predY[1]=_mm_min_epi16(predY[1], ymax[1]);
				predU[0]=_mm_min_epi16(predU[0], umax[0]); predU[1]=_mm_min_epi16(predU[1], umax[1]);
				predV[0]=_mm_min_epi16(predV[0], vmax[0]); predV[1]=_mm_min_epi16(predV[1], vmax[1]);
			}
			if(dist>1)
			{
				__m128i val[6], tmp;
				if(fwd)
				{
					myuv[0]=_mm_add_epi8(_mm_load_si128((__m128i*)(imptr+0*NCODERS)), half8);//load rgb
					myuv[2]=_mm_add_epi8(_mm_load_si128((__m128i*)(imptr+1*NCODERS)), half8);
					myuv[4]=_mm_add_epi8(_mm_load_si128((__m128i*)(imptr+2*NCODERS)), half8);
					myuv[1]=_mm_cvtepi8_epi16(_mm_shuffle_epi32(myuv[0], _MM_SHUFFLE(1, 0, 3, 2)));
					myuv[3]=_mm_cvtepi8_epi16(_mm_shuffle_epi32(myuv[2], _MM_SHUFFLE(1, 0, 3, 2)));
					myuv[5]=_mm_cvtepi8_epi16(_mm_shuffle_epi32(myuv[4], _MM_SHUFFLE(1, 0, 3, 2)));
					myuv[0]=_mm_cvtepi8_epi16(myuv[0]);
					myuv[2]=_mm_cvtepi8_epi16(myuv[2]);
					myuv[4]=_mm_cvtepi8_epi16(myuv[4]);

					//val=(curr-(int)pred)/dist
					val[0]=_mm_sub_epi16(myuv[0], predY[0]);
					val[1]=_mm_sub_epi16(myuv[1], predV[1]);
					val[2]=_mm_sub_epi16(myuv[2], predU[0]);
					val[3]=_mm_sub_epi16(myuv[3], predU[1]);
					val[4]=_mm_sub_epi16(myuv[4], predV[0]);
					val[5]=_mm_sub_epi16(myuv[5], predV[1]);
					{
						__m128i cond[6];
						cond[0]=_mm_srai_epi16(val[0], 15);
						cond[1]=_mm_srai_epi16(val[1], 15);
						cond[2]=_mm_srai_epi16(val[2], 15);
						cond[3]=_mm_srai_epi16(val[3], 15);
						cond[4]=_mm_srai_epi16(val[4], 15);
						cond[5]=_mm_srai_epi16(val[5], 15);
						val[0]=_mm_mulhi_epi16(val[0], dist_rcp);
						val[1]=_mm_mulhi_epi16(val[1], dist_rcp);
						val[2]=_mm_mulhi_epi16(val[2], dist_rcp);
						val[3]=_mm_mulhi_epi16(val[3], dist_rcp);
						val[4]=_mm_mulhi_epi16(val[4], dist_rcp);
						val[5]=_mm_mulhi_epi16(val[5], dist_rcp);
						val[0]=_mm_sub_epi16(val[0], cond[0]);//(x*inv>>16)-(x>>31)
						val[1]=_mm_sub_epi16(val[1], cond[1]);
						val[2]=_mm_sub_epi16(val[2], cond[2]);
						val[3]=_mm_sub_epi16(val[3], cond[3]);
						val[4]=_mm_sub_epi16(val[4], cond[4]);
						val[5]=_mm_sub_epi16(val[5], cond[5]);
					}

					//curr=dist*val+pred
					myuv[0]=_mm_mullo_epi16(val[0], mdist);
					myuv[1]=_mm_mullo_epi16(val[1], mdist);
					myuv[2]=_mm_mullo_epi16(val[2], mdist);
					myuv[3]=_mm_mullo_epi16(val[3], mdist);
					myuv[4]=_mm_mullo_epi16(val[4], mdist);
					myuv[5]=_mm_mullo_epi16(val[5], mdist);

					myuv[0]=_mm_add_epi16(myuv[0], predY[0]);
					myuv[1]=_mm_add_epi16(myuv[1], predY[1]);
					myuv[2]=_mm_add_epi16(myuv[2], predU[0]);
					myuv[3]=_mm_add_epi16(myuv[3], predU[1]);
					myuv[4]=_mm_add_epi16(myuv[4], predV[0]);
					myuv[5]=_mm_add_epi16(myuv[5], predV[1]);

					tmp=_mm_set1_epi16(-128);
					myuv[0]=_mm_max_epi16(myuv[0], tmp);
					myuv[1]=_mm_max_epi16(myuv[1], tmp);
					myuv[2]=_mm_max_epi16(myuv[2], tmp);
					myuv[3]=_mm_max_epi16(myuv[3], tmp);
					myuv[4]=_mm_max_epi16(myuv[4], tmp);
					myuv[5]=_mm_max_epi16(myuv[5], tmp);
					tmp=_mm_set1_epi16(127);
					myuv[0]=_mm_min_epi16(myuv[0], tmp);
					myuv[1]=_mm_min_epi16(myuv[1], tmp);
					myuv[2]=_mm_min_epi16(myuv[2], tmp);
					myuv[3]=_mm_min_epi16(myuv[3], tmp);
					myuv[4]=_mm_min_epi16(myuv[4], tmp);
					myuv[5]=_mm_min_epi16(myuv[5], tmp);
					
					//post-RCT
					val[0*2+0]=_mm_sub_epi16(val[0*2+0], val[1*2+0]);
					val[0*2+1]=_mm_sub_epi16(val[0*2+1], val[1*2+1]);

					val[2*2+0]=_mm_sub_epi16(val[2*2+0], val[1*2+0]);
					val[2*2+1]=_mm_sub_epi16(val[2*2+1], val[1*2+1]);

					val[1*2+0]=_mm_add_epi16(val[1*2+0], _mm_srai_epi16(_mm_add_epi16(val[0*2+0], val[2*2+0]), 2));
					val[1*2+1]=_mm_add_epi16(val[1*2+1], _mm_srai_epi16(_mm_add_epi16(val[0*2+1], val[2*2+1]), 2));

					ecurr[0]=val[0];
					ecurr[1]=val[1];
					ecurr[2]=val[2];
					ecurr[3]=val[3];
					ecurr[4]=val[4];
					ecurr[5]=val[5];
#ifdef _DEBUG
					{
						int k;
						for(k=0;k<48;++k)
						{
							int v=((short*)val)[k];
							if((unsigned)(v+128)>=256)
								LOG_ERROR("");
						}
					}
#endif
					tmp=_mm_set1_epi16(128);
					val[0]=_mm_add_epi16(val[0], tmp);
					val[1]=_mm_add_epi16(val[1], tmp);
					val[2]=_mm_add_epi16(val[2], tmp);
					val[3]=_mm_add_epi16(val[3], tmp);
					val[4]=_mm_add_epi16(val[4], tmp);
					val[5]=_mm_add_epi16(val[5], tmp);
					tmp=_mm_set1_epi16(255);
					val[0]=_mm_and_si128(val[0], tmp);
					val[1]=_mm_and_si128(val[1], tmp);
					val[2]=_mm_and_si128(val[2], tmp);
					val[3]=_mm_and_si128(val[3], tmp);
					val[4]=_mm_and_si128(val[4], tmp);
					val[5]=_mm_and_si128(val[5], tmp);

					ctxU[0]=_mm_add_epi16(ctxU[0], mctxuoffset);
					ctxU[1]=_mm_add_epi16(ctxU[1], mctxuoffset);
					ctxV[0]=_mm_add_epi16(ctxV[0], mctxvoffset);
					ctxV[1]=_mm_add_epi16(ctxV[1], mctxvoffset);
					ctxY[0]=_mm_slli_epi16(ctxY[0], 8);
					ctxY[1]=_mm_slli_epi16(ctxY[1], 8);
					ctxU[0]=_mm_slli_epi16(ctxU[0], 8);
					ctxU[1]=_mm_slli_epi16(ctxU[1], 8);
					ctxV[0]=_mm_slli_epi16(ctxV[0], 8);
					ctxV[1]=_mm_slli_epi16(ctxV[1], 8);
					ctxY[0]=_mm_or_si128(ctxY[0], val[0]);
					ctxY[1]=_mm_or_si128(ctxY[1], val[1]);
					ctxU[0]=_mm_or_si128(ctxU[0], val[2]);
					ctxU[1]=_mm_or_si128(ctxU[1], val[3]);
					ctxV[0]=_mm_or_si128(ctxV[0], val[4]);
					ctxV[1]=_mm_or_si128(ctxV[1], val[5]);
					_mm_store_si128((__m128i*)syms+0, ctxY[0]);
					_mm_store_si128((__m128i*)syms+1, ctxY[1]);
					_mm_store_si128((__m128i*)syms+2, ctxU[0]);
					_mm_store_si128((__m128i*)syms+3, ctxU[1]);
					_mm_store_si128((__m128i*)syms+4, ctxV[0]);
					_mm_store_si128((__m128i*)syms+5, ctxV[1]);
					_mm_store_si128((__m128i*)ctxptr+0, ctxY[0]);
					_mm_store_si128((__m128i*)ctxptr+1, ctxY[1]);
					_mm_store_si128((__m128i*)ctxptr+2, ctxU[0]);
					_mm_store_si128((__m128i*)ctxptr+3, ctxU[1]);
					_mm_store_si128((__m128i*)ctxptr+4, ctxV[0]);
					_mm_store_si128((__m128i*)ctxptr+5, ctxV[1]);

#if 0
					for(int k=0;k<48;++k)
					{
						if((unsigned)syms[k]>=(unsigned)(18*3<<8))
							LOG_ERROR("");
					}
#endif
					{
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
					}
					ctxptr+=sizeof(short[3][NCODERS]);
				}
				else
				{
					__m128i msyms24[3];

					dec_yuv(mstate, 0, ctxY, (int*)CDF2syms, &streamptr, streamend, myuv+0*2);
					dec_yuv(mstate, 1, ctxU, (int*)CDF2syms, &streamptr, streamend, myuv+1*2);
					dec_yuv(mstate, 2, ctxV, (int*)CDF2syms, &streamptr, streamend, myuv+2*2);
					
					tmp=_mm_set1_epi16(128);
					val[0]=_mm_sub_epi16(myuv[0], tmp);
					val[1]=_mm_sub_epi16(myuv[1], tmp);
					val[2]=_mm_sub_epi16(myuv[2], tmp);
					val[3]=_mm_sub_epi16(myuv[3], tmp);
					val[4]=_mm_sub_epi16(myuv[4], tmp);
					val[5]=_mm_sub_epi16(myuv[5], tmp);
					ecurr[0]=val[0];
					ecurr[1]=val[1];
					ecurr[2]=val[2];
					ecurr[3]=val[3];
					ecurr[4]=val[4];
					ecurr[5]=val[5];
					
					//post-invRCT
					val[1*2+0]=_mm_sub_epi16(val[1*2+0], _mm_srai_epi16(_mm_add_epi16(val[0*2+0], val[2*2+0]), 2));
					val[1*2+1]=_mm_sub_epi16(val[1*2+1], _mm_srai_epi16(_mm_add_epi16(val[0*2+1], val[2*2+1]), 2));
					val[0*2+0]=_mm_add_epi16(val[0*2+0], val[1*2+0]);
					val[0*2+1]=_mm_add_epi16(val[0*2+1], val[1*2+1]);

					val[2*2+0]=_mm_add_epi16(val[2*2+0], val[1*2+0]);
					val[2*2+1]=_mm_add_epi16(val[2*2+1], val[1*2+1]);

					//val=dist*curr+pred
					val[0]=_mm_mullo_epi16(val[0], mdist);
					val[1]=_mm_mullo_epi16(val[1], mdist);
					val[2]=_mm_mullo_epi16(val[2], mdist);
					val[3]=_mm_mullo_epi16(val[3], mdist);
					val[4]=_mm_mullo_epi16(val[4], mdist);
					val[5]=_mm_mullo_epi16(val[5], mdist);
					val[0]=_mm_add_epi16(val[0], predY[0]);
					val[1]=_mm_add_epi16(val[1], predY[1]);
					val[2]=_mm_add_epi16(val[2], predU[0]);
					val[3]=_mm_add_epi16(val[3], predU[1]);
					val[4]=_mm_add_epi16(val[4], predV[0]);
					val[5]=_mm_add_epi16(val[5], predV[1]);

					//clamp
					tmp=_mm_set1_epi16(-128);
					val[0]=_mm_max_epi16(val[0], tmp);
					val[1]=_mm_max_epi16(val[1], tmp);
					val[2]=_mm_max_epi16(val[2], tmp);
					val[3]=_mm_max_epi16(val[3], tmp);
					val[4]=_mm_max_epi16(val[4], tmp);
					val[5]=_mm_max_epi16(val[5], tmp);
					tmp=_mm_set1_epi16(127);
					val[0]=_mm_min_epi16(val[0], tmp);
					val[1]=_mm_min_epi16(val[1], tmp);
					val[2]=_mm_min_epi16(val[2], tmp);
					val[3]=_mm_min_epi16(val[3], tmp);
					val[4]=_mm_min_epi16(val[4], tmp);
					val[5]=_mm_min_epi16(val[5], tmp);
					myuv[0]=val[0];
					myuv[1]=val[1];
					myuv[2]=val[2];
					myuv[3]=val[3];
					myuv[4]=val[4];
					myuv[5]=val[5];

					//16 -> 8
					msyms24[0]=_mm_packs_epi16(val[0], val[1]);
					msyms24[1]=_mm_packs_epi16(val[2], val[3]);
					msyms24[2]=_mm_packs_epi16(val[4], val[5]);
					tmp=_mm_set1_epi8(-128);
					msyms24[0]=_mm_xor_si128(msyms24[0], tmp);
					msyms24[1]=_mm_xor_si128(msyms24[1], tmp);
					msyms24[2]=_mm_xor_si128(msyms24[2], tmp);
					_mm_store_si128((__m128i*)(imptr+0*NCODERS), msyms24[0]);//store rgb
					_mm_store_si128((__m128i*)(imptr+1*NCODERS), msyms24[1]);
					_mm_store_si128((__m128i*)(imptr+2*NCODERS), msyms24[2]);
#ifdef ENABLE_GUIDE
					{
						int kc, k;

						for(kc=0;kc<3;++kc)
						{
							for(k=0;k<NCODERS;++k)
							{
								double diff=imptr[kc*NCODERS+k]-g_image[imptr-interleaved+kc*NCODERS+k];
								g_sqe[kc]+=diff*diff;
							}
						}
					}
#endif
				}
				_mm_store_si128((__m128i*)rows[0]+0+0+0*12, myuv[0]);
				_mm_store_si128((__m128i*)rows[0]+0+1+0*12, myuv[1]);
				_mm_store_si128((__m128i*)rows[0]+0+2+0*12, myuv[2]);
				_mm_store_si128((__m128i*)rows[0]+0+3+0*12, myuv[3]);
				_mm_store_si128((__m128i*)rows[0]+0+4+0*12, myuv[4]);
				_mm_store_si128((__m128i*)rows[0]+0+5+0*12, myuv[5]);

				W[0]=myuv[0];
				W[1]=myuv[1];
				W[2]=myuv[2];
				W[3]=myuv[3];
				W[4]=myuv[4];
				W[5]=myuv[5];
				ecurr[0]=_mm_xor_si128(_mm_slli_epi16(ecurr[0], 1), _mm_srai_epi16(ecurr[0], 15));
				ecurr[1]=_mm_xor_si128(_mm_slli_epi16(ecurr[1], 1), _mm_srai_epi16(ecurr[1], 15));
				ecurr[2]=_mm_xor_si128(_mm_slli_epi16(ecurr[2], 1), _mm_srai_epi16(ecurr[2], 15));
				ecurr[3]=_mm_xor_si128(_mm_slli_epi16(ecurr[3], 1), _mm_srai_epi16(ecurr[3], 15));
				ecurr[4]=_mm_xor_si128(_mm_slli_epi16(ecurr[4], 1), _mm_srai_epi16(ecurr[4], 15));
				ecurr[5]=_mm_xor_si128(_mm_slli_epi16(ecurr[5], 1), _mm_srai_epi16(ecurr[5], 15));
			}
			else if(fwd)
			{
				__m128i ctxblendmask=_mm_set1_epi16(255);
				myuv[0]=_mm_add_epi8(_mm_load_si128((__m128i*)(imptr+yidx)), half8);//load yuv
				myuv[2]=_mm_add_epi8(_mm_load_si128((__m128i*)(imptr+uidx)), half8);
				myuv[4]=_mm_add_epi8(_mm_load_si128((__m128i*)(imptr+vidx)), half8);
				myuv[1]=_mm_cvtepi8_epi16(_mm_shuffle_epi32(myuv[0], _MM_SHUFFLE(1, 0, 3, 2)));
				myuv[3]=_mm_cvtepi8_epi16(_mm_shuffle_epi32(myuv[2], _MM_SHUFFLE(1, 0, 3, 2)));
				myuv[5]=_mm_cvtepi8_epi16(_mm_shuffle_epi32(myuv[4], _MM_SHUFFLE(1, 0, 3, 2)));
				myuv[0]=_mm_cvtepi8_epi16(myuv[0]);
				myuv[2]=_mm_cvtepi8_epi16(myuv[2]);
				myuv[4]=_mm_cvtepi8_epi16(myuv[4]);

				//encode Y
				W[0]=myuv[0*2+0];
				W[1]=myuv[0*2+1];
				_mm_store_si128((__m128i*)rows[0]+0+0+0*12, myuv[0*2+0]);//store Y neighbors
				_mm_store_si128((__m128i*)rows[0]+0+1+0*12, myuv[0*2+1]);
				msyms[0]=_mm_sub_epi16(myuv[0*2+0], predY[0]);//sub pred
				msyms[1]=_mm_sub_epi16(myuv[0*2+1], predY[1]);
				ecurr[0]=_mm_xor_si128(_mm_slli_epi16(msyms[0], 1), _mm_srai_epi16(msyms[0], 15));//ecurr = pack_sign(yuv-pred)
				ecurr[1]=_mm_xor_si128(_mm_slli_epi16(msyms[1], 1), _mm_srai_epi16(msyms[1], 15));
				msyms[0]=_mm_sub_epi16(msyms[0], amin);
				msyms[1]=_mm_sub_epi16(msyms[1], amin);
				ctxY[0]=_mm_slli_epi16(ctxY[0], 8);
				ctxY[1]=_mm_slli_epi16(ctxY[1], 8);
				ctxY[0]=_mm_blendv_epi8(ctxY[0], msyms[0], ctxblendmask);
				ctxY[1]=_mm_blendv_epi8(ctxY[1], msyms[1], ctxblendmask);
				_mm_store_si128((__m128i*)syms+0, ctxY[0]);
				_mm_store_si128((__m128i*)syms+1, ctxY[1]);
				_mm_store_si128((__m128i*)ctxptr+0, ctxY[0]);//store Y  ctx|residuals
				_mm_store_si128((__m128i*)ctxptr+1, ctxY[1]);
				
				//encode U
				moffset[0]=_mm_and_si128(myuv[0*2+0], uhelpmask);
				moffset[1]=_mm_and_si128(myuv[0*2+1], uhelpmask);
				predU[0]=_mm_add_epi16(predU[0], moffset[0]);
				predU[1]=_mm_add_epi16(predU[1], moffset[1]);
				predU[0]=_mm_max_epi16(predU[0], amin);
				predU[1]=_mm_max_epi16(predU[1], amin);
				predU[0]=_mm_min_epi16(predU[0], amax);
				predU[1]=_mm_min_epi16(predU[1], amax);

				msyms[0]=_mm_sub_epi16(myuv[1*2+0], predU[0]);
				msyms[1]=_mm_sub_epi16(myuv[1*2+1], predU[1]);
				ecurr[1*2+0]=_mm_xor_si128(_mm_slli_epi16(msyms[0], 1), _mm_srai_epi16(msyms[0], 15));
				ecurr[1*2+1]=_mm_xor_si128(_mm_slli_epi16(msyms[1], 1), _mm_srai_epi16(msyms[1], 15));
				msyms[0]=_mm_sub_epi16(msyms[0], amin);
				msyms[1]=_mm_sub_epi16(msyms[1], amin);
				ctxU[0]=_mm_add_epi16(ctxU[0], mctxuoffset);
				ctxU[1]=_mm_add_epi16(ctxU[1], mctxuoffset);
				ctxU[0]=_mm_slli_epi16(ctxU[0], 8);
				ctxU[1]=_mm_slli_epi16(ctxU[1], 8);
				ctxU[0]=_mm_blendv_epi8(ctxU[0], msyms[0], ctxblendmask);
				ctxU[1]=_mm_blendv_epi8(ctxU[1], msyms[1], ctxblendmask);
				_mm_store_si128((__m128i*)syms+1*2+0, ctxU[0]);
				_mm_store_si128((__m128i*)syms+1*2+1, ctxU[1]);
				_mm_store_si128((__m128i*)ctxptr+1*2+0, ctxU[0]);//store U  ctx|residuals
				_mm_store_si128((__m128i*)ctxptr+1*2+1, ctxU[1]);
				W[1*2+0]=_mm_sub_epi16(myuv[1*2+0], moffset[0]);
				W[1*2+1]=_mm_sub_epi16(myuv[1*2+1], moffset[1]);
				_mm_store_si128((__m128i*)rows[0]+0+2+0*12, W[1*2+0]);//store U neighbors
				_mm_store_si128((__m128i*)rows[0]+0+3+0*12, W[1*2+1]);
				
				//encode V
				moffset[0]=_mm_mullo_epi16(vc0, myuv[0]);
				moffset[1]=_mm_mullo_epi16(vc0, myuv[1]);
				moffset[0]=_mm_add_epi16(moffset[0], _mm_mullo_epi16(vc1, myuv[2]));
				moffset[1]=_mm_add_epi16(moffset[1], _mm_mullo_epi16(vc1, myuv[3]));
				moffset[0]=_mm_srai_epi16(moffset[0], 2);
				moffset[1]=_mm_srai_epi16(moffset[1], 2);
				predV[0]=_mm_add_epi16(predV[0], moffset[0]);
				predV[1]=_mm_add_epi16(predV[1], moffset[1]);
				predV[0]=_mm_max_epi16(predV[0], amin);
				predV[1]=_mm_max_epi16(predV[1], amin);
				predV[0]=_mm_min_epi16(predV[0], amax);
				predV[1]=_mm_min_epi16(predV[1], amax);

				msyms[0]=_mm_sub_epi16(myuv[2*2+0], predV[0]);
				msyms[1]=_mm_sub_epi16(myuv[2*2+1], predV[1]);
				ecurr[2*2+0]=_mm_xor_si128(_mm_slli_epi16(msyms[0], 1), _mm_srai_epi16(msyms[0], 15));
				ecurr[2*2+1]=_mm_xor_si128(_mm_slli_epi16(msyms[1], 1), _mm_srai_epi16(msyms[1], 15));
				msyms[0]=_mm_sub_epi16(msyms[0], amin);
				msyms[1]=_mm_sub_epi16(msyms[1], amin);
				ctxV[0]=_mm_add_epi16(ctxV[0], mctxvoffset);
				ctxV[1]=_mm_add_epi16(ctxV[1], mctxvoffset);
				ctxV[0]=_mm_slli_epi16(ctxV[0], 8);
				ctxV[1]=_mm_slli_epi16(ctxV[1], 8);
				ctxV[0]=_mm_blendv_epi8(ctxV[0], msyms[0], ctxblendmask);
				ctxV[1]=_mm_blendv_epi8(ctxV[1], msyms[1], ctxblendmask);
				_mm_store_si128((__m128i*)syms+2*2+0, ctxV[0]);
				_mm_store_si128((__m128i*)syms+2*2+1, ctxV[1]);
				_mm_store_si128((__m128i*)ctxptr+2*2+0, ctxV[0]);//store V  ctx|residuals		ctxptr+NCODERS*(C*2+R)
				_mm_store_si128((__m128i*)ctxptr+2*2+1, ctxV[1]);
				W[2*2+0]=_mm_sub_epi16(myuv[2*2+0], moffset[0]);
				W[2*2+1]=_mm_sub_epi16(myuv[2*2+1], moffset[1]);
				_mm_store_si128((__m128i*)rows[0]+0+4+0*12, W[2*2+0]);//store V neighbors
				_mm_store_si128((__m128i*)rows[0]+0+5+0*12, W[2*2+1]);

				{
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
				}
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
				ctxptr+=sizeof(short[3][NCODERS]);
			}
			else
			{
				//decode main
				__m128i msyms8;
				
				//decode Y
				dec_yuv(mstate, 0, ctxY, (int*)CDF2syms, &streamptr, streamend, myuv+0*2);//residuals from [0 ~ 255]
				dec_yuv(mstate, 1, ctxU, (int*)CDF2syms, &streamptr, streamend, myuv+1*2);
				dec_yuv(mstate, 2, ctxV, (int*)CDF2syms, &streamptr, streamend, myuv+2*2);

				//yuv = (char)(sym+pred-128)	= (unsigned char)(sym+pred)-128
				myuv[0*2+0]=_mm_add_epi16(myuv[0*2+0], predY[0]);
				myuv[0*2+1]=_mm_add_epi16(myuv[0*2+1], predY[1]);
				myuv[0*2+0]=_mm_and_si128(myuv[0*2+0], bytemask);
				myuv[0*2+1]=_mm_and_si128(myuv[0*2+1], bytemask);
				myuv[0*2+0]=_mm_add_epi16(myuv[0*2+0], amin);
				myuv[0*2+1]=_mm_add_epi16(myuv[0*2+1], amin);
				msyms[0]=_mm_sub_epi16(myuv[0*2+0], predY[0]);//sub pred
				msyms[1]=_mm_sub_epi16(myuv[0*2+1], predY[1]);
				ecurr[0*2+0]=_mm_xor_si128(_mm_slli_epi16(msyms[0], 1), _mm_srai_epi16(msyms[0], 15));
				ecurr[0*2+1]=_mm_xor_si128(_mm_slli_epi16(msyms[1], 1), _mm_srai_epi16(msyms[1], 15));
				msyms8=_mm_packs_epi16(myuv[0*2+0], myuv[0*2+1]);
				msyms8=_mm_xor_si128(msyms8, _mm_set1_epi8(-128));
				_mm_store_si128((__m128i*)(imptr+yidx), msyms8);//store Y bytes
#ifdef ENABLE_GUIDE
				if(memcmp(imptr+yidx, g_image+(imptr-interleaved)+yidx, NCODERS))
				{
					int k;
					printf("original  decoded  original-decoded  XYC0 %d %d %d\n", kx, ky, yidx);
					for(k=0;k<NCODERS;++k)
						printf("0x%02X  0x%02X  %4d\n",
							g_image[imptr-interleaved+yidx+k],
							imptr[yidx+k],
							g_image[imptr-interleaved+yidx+k]-imptr[yidx+k]
						);
					LOG_ERROR("guide error XYC0 %d %d %d/%d", kx, ky, yidx, NCODERS);
				}
#endif
				W[0]=myuv[0];
				W[1]=myuv[1];
				_mm_store_si128((__m128i*)rows[0]+0+0+0*12, myuv[0]);//store Y neighbors
				_mm_store_si128((__m128i*)rows[0]+0+1+0*12, myuv[1]);


				//decode U
				moffset[0]=_mm_and_si128(myuv[0*2+0], uhelpmask);
				moffset[1]=_mm_and_si128(myuv[0*2+1], uhelpmask);
				predU[0]=_mm_add_epi16(predU[0], moffset[0]);
				predU[1]=_mm_add_epi16(predU[1], moffset[1]);
				predU[0]=_mm_max_epi16(predU[0], amin);
				predU[1]=_mm_max_epi16(predU[1], amin);
				predU[0]=_mm_min_epi16(predU[0], amax);
				predU[1]=_mm_min_epi16(predU[1], amax);
				myuv[1*2+0]=_mm_add_epi16(myuv[1*2+0], predU[0]);
				myuv[1*2+1]=_mm_add_epi16(myuv[1*2+1], predU[1]);
				myuv[1*2+0]=_mm_and_si128(myuv[1*2+0], bytemask);
				myuv[1*2+1]=_mm_and_si128(myuv[1*2+1], bytemask);
				myuv[1*2+0]=_mm_add_epi16(myuv[1*2+0], amin);
				myuv[1*2+1]=_mm_add_epi16(myuv[1*2+1], amin);
				msyms[0]=_mm_sub_epi16(myuv[1*2+0], predU[0]);//sub pred
				msyms[1]=_mm_sub_epi16(myuv[1*2+1], predU[1]);
				ecurr[1*2+0]=_mm_xor_si128(_mm_slli_epi16(msyms[0], 1), _mm_srai_epi16(msyms[0], 15));
				ecurr[1*2+1]=_mm_xor_si128(_mm_slli_epi16(msyms[1], 1), _mm_srai_epi16(msyms[1], 15));
				msyms8=_mm_packs_epi16(myuv[1*2+0], myuv[1*2+1]);
				msyms8=_mm_xor_si128(msyms8, _mm_set1_epi8(-128));
				_mm_store_si128((__m128i*)(imptr+uidx), msyms8);//store U bytes
#ifdef ENABLE_GUIDE
				if(memcmp(imptr+uidx, g_image+(imptr-interleaved)+uidx, NCODERS))
				{
					int k;
					printf("original  decoded  original-decoded  XYC1 %d %d %d\n", kx, ky, uidx);
					for(k=0;k<NCODERS;++k)
						printf("0x%02X  0x%02X  %4d\n",
							g_image[imptr-interleaved+uidx+k],
							imptr[uidx+k],
							g_image[imptr-interleaved+uidx+k]-imptr[uidx+k]
						);
					LOG_ERROR("guide error XYC1 %d %d %d/%d", kx, ky, uidx, NCODERS);
				}
#endif
				W[1*2+0]=_mm_sub_epi16(myuv[1*2+0], moffset[0]);//subtract Uoffset from U
				W[1*2+1]=_mm_sub_epi16(myuv[1*2+1], moffset[1]);
				_mm_store_si128((__m128i*)rows[0]+0+2+0*12, W[1*2+0]);//store U neighbors
				_mm_store_si128((__m128i*)rows[0]+0+3+0*12, W[1*2+1]);
				

				//decode V
				moffset[0]=_mm_mullo_epi16(vc0, myuv[0*2+0]);
				moffset[1]=_mm_mullo_epi16(vc0, myuv[0*2+1]);
				moffset[0]=_mm_add_epi16(moffset[0], _mm_mullo_epi16(vc1, myuv[1*2+0]));
				moffset[1]=_mm_add_epi16(moffset[1], _mm_mullo_epi16(vc1, myuv[1*2+1]));
				moffset[0]=_mm_srai_epi16(moffset[0], 2);
				moffset[1]=_mm_srai_epi16(moffset[1], 2);
				predV[0]=_mm_add_epi16(predV[0], moffset[0]);
				predV[1]=_mm_add_epi16(predV[1], moffset[1]);
				predV[0]=_mm_max_epi16(predV[0], amin);
				predV[1]=_mm_max_epi16(predV[1], amin);
				predV[0]=_mm_min_epi16(predV[0], amax);
				predV[1]=_mm_min_epi16(predV[1], amax);
				myuv[2*2+0]=_mm_add_epi16(myuv[2*2+0], predV[0]);
				myuv[2*2+1]=_mm_add_epi16(myuv[2*2+1], predV[1]);
				myuv[2*2+0]=_mm_and_si128(myuv[2*2+0], bytemask);
				myuv[2*2+1]=_mm_and_si128(myuv[2*2+1], bytemask);
				myuv[2*2+0]=_mm_add_epi16(myuv[2*2+0], amin);
				myuv[2*2+1]=_mm_add_epi16(myuv[2*2+1], amin);
				msyms[0]=_mm_sub_epi16(myuv[2*2+0], predV[0]);
				msyms[1]=_mm_sub_epi16(myuv[2*2+1], predV[1]);
				ecurr[2*2+0]=_mm_xor_si128(_mm_slli_epi16(msyms[0], 1), _mm_srai_epi16(msyms[0], 15));
				ecurr[2*2+1]=_mm_xor_si128(_mm_slli_epi16(msyms[1], 1), _mm_srai_epi16(msyms[1], 15));
				msyms8=_mm_packs_epi16(myuv[2*2+0], myuv[2*2+1]);
				msyms8=_mm_xor_si128(msyms8, _mm_set1_epi8(-128));
				_mm_store_si128((__m128i*)(imptr+vidx), msyms8);//store V bytes
#ifdef ENABLE_GUIDE
				if(memcmp(imptr+vidx, g_image+(imptr-interleaved)+vidx, NCODERS))
				{
					int k;
					printf("original  decoded  original-decoded  XYC2 %d %d %d\n", kx, ky, vidx);
					for(k=0;k<NCODERS;++k)
						printf("0x%02X  0x%02X  %4d\n",
							g_image[imptr-interleaved+vidx+k],
							imptr[vidx+k],
							g_image[imptr-interleaved+vidx+k]-imptr[vidx+k]
						);
					LOG_ERROR("guide error XYC2 %d %d %d/%d", kx, ky, vidx, NCODERS);
				}
#endif
				W[2*2+0]=_mm_sub_epi16(myuv[2*2+0], moffset[0]);//subtract Voffset from V
				W[2*2+1]=_mm_sub_epi16(myuv[2*2+1], moffset[1]);
				_mm_store_si128((__m128i*)rows[0]+0+4+0*12, W[2*2+0]);//store V neighbors
				_mm_store_si128((__m128i*)rows[0]+0+5+0*12, W[2*2+1]);
			}
			if(effort==1)//update
			{
				__m128i mu[6];
				int k;

				mu[0]=_mm_sub_epi16(W[0], predYUV0[0]);
				mu[1]=_mm_sub_epi16(W[1], predYUV0[1]);
				mu[2]=_mm_sub_epi16(W[2], predYUV0[2]);
				mu[3]=_mm_sub_epi16(W[3], predYUV0[3]);
				mu[4]=_mm_sub_epi16(W[4], predYUV0[4]);
				mu[5]=_mm_sub_epi16(W[5], predYUV0[5]);
				for(k=0;k<L1_NPREDS1;++k)//update
				{
					__m128i mc[12];
					mc[0x0]=_mm_sign_epi16(L1preds[k*6+0], mu[0]);//L1
					mc[0x1]=_mm_sign_epi16(L1preds[k*6+1], mu[1]);
					mc[0x2]=_mm_sign_epi16(L1preds[k*6+2], mu[2]);
					mc[0x3]=_mm_sign_epi16(L1preds[k*6+3], mu[3]);
					mc[0x4]=_mm_sign_epi16(L1preds[k*6+4], mu[4]);
					mc[0x5]=_mm_sign_epi16(L1preds[k*6+5], mu[5]);
					//mc[0x0]=_mm_mullo_epi16(L1preds[k*6+0], mu[0]);//L2
					//mc[0x1]=_mm_mullo_epi16(L1preds[k*6+1], mu[1]);
					//mc[0x2]=_mm_mullo_epi16(L1preds[k*6+2], mu[2]);
					//mc[0x3]=_mm_mullo_epi16(L1preds[k*6+3], mu[3]);
					//mc[0x4]=_mm_mullo_epi16(L1preds[k*6+4], mu[4]);
					//mc[0x5]=_mm_mullo_epi16(L1preds[k*6+5], mu[5]);

					//signed 16 -> 32	6 lo 6 hi registers
					mc[0x6]=_mm_srai_epi32(mc[0x0], 16);
					mc[0x7]=_mm_srai_epi32(mc[0x1], 16);
					mc[0x8]=_mm_srai_epi32(mc[0x2], 16);
					mc[0x9]=_mm_srai_epi32(mc[0x3], 16);
					mc[0xA]=_mm_srai_epi32(mc[0x4], 16);
					mc[0xB]=_mm_srai_epi32(mc[0x5], 16);

					mc[0x0]=_mm_slli_epi32(mc[0x0], 16);
					mc[0x1]=_mm_slli_epi32(mc[0x1], 16);
					mc[0x2]=_mm_slli_epi32(mc[0x2], 16);
					mc[0x3]=_mm_slli_epi32(mc[0x3], 16);
					mc[0x4]=_mm_slli_epi32(mc[0x4], 16);
					mc[0x5]=_mm_slli_epi32(mc[0x5], 16);

					mc[0x0]=_mm_srai_epi32(mc[0x0], 16);
					mc[0x1]=_mm_srai_epi32(mc[0x1], 16);
					mc[0x2]=_mm_srai_epi32(mc[0x2], 16);
					mc[0x3]=_mm_srai_epi32(mc[0x3], 16);
					mc[0x4]=_mm_srai_epi32(mc[0x4], 16);
					mc[0x5]=_mm_srai_epi32(mc[0x5], 16);
					mc[0x0]=_mm_add_epi32(mc[0x0], _mm_load_si128((__m128i*)L1weights+k*12+0x0));//update coeffs
					mc[0x1]=_mm_add_epi32(mc[0x1], _mm_load_si128((__m128i*)L1weights+k*12+0x1));
					mc[0x2]=_mm_add_epi32(mc[0x2], _mm_load_si128((__m128i*)L1weights+k*12+0x2));
					mc[0x3]=_mm_add_epi32(mc[0x3], _mm_load_si128((__m128i*)L1weights+k*12+0x3));
					mc[0x4]=_mm_add_epi32(mc[0x4], _mm_load_si128((__m128i*)L1weights+k*12+0x4));
					mc[0x5]=_mm_add_epi32(mc[0x5], _mm_load_si128((__m128i*)L1weights+k*12+0x5));
					mc[0x6]=_mm_add_epi32(mc[0x6], _mm_load_si128((__m128i*)L1weights+k*12+0x6));
					mc[0x7]=_mm_add_epi32(mc[0x7], _mm_load_si128((__m128i*)L1weights+k*12+0x7));
					mc[0x8]=_mm_add_epi32(mc[0x8], _mm_load_si128((__m128i*)L1weights+k*12+0x8));
					mc[0x9]=_mm_add_epi32(mc[0x9], _mm_load_si128((__m128i*)L1weights+k*12+0x9));
					mc[0xA]=_mm_add_epi32(mc[0xA], _mm_load_si128((__m128i*)L1weights+k*12+0xA));
					mc[0xB]=_mm_add_epi32(mc[0xB], _mm_load_si128((__m128i*)L1weights+k*12+0xB));
					_mm_store_si128((__m128i*)L1weights+k*12+0x0, mc[0x0]);
					_mm_store_si128((__m128i*)L1weights+k*12+0x1, mc[0x1]);
					_mm_store_si128((__m128i*)L1weights+k*12+0x2, mc[0x2]);
					_mm_store_si128((__m128i*)L1weights+k*12+0x3, mc[0x3]);
					_mm_store_si128((__m128i*)L1weights+k*12+0x4, mc[0x4]);
					_mm_store_si128((__m128i*)L1weights+k*12+0x5, mc[0x5]);
					_mm_store_si128((__m128i*)L1weights+k*12+0x6, mc[0x6]);
					_mm_store_si128((__m128i*)L1weights+k*12+0x7, mc[0x7]);
					_mm_store_si128((__m128i*)L1weights+k*12+0x8, mc[0x8]);
					_mm_store_si128((__m128i*)L1weights+k*12+0x9, mc[0x9]);
					_mm_store_si128((__m128i*)L1weights+k*12+0xA, mc[0xA]);
					_mm_store_si128((__m128i*)L1weights+k*12+0xB, mc[0xB]);
				}
			}
			else if(effort==2)//update
			{
				__m128i mu[6];
				int k;

				mu[0]=_mm_sub_epi16(W[0], predYUV0[0]);
				mu[1]=_mm_sub_epi16(W[1], predYUV0[1]);
				mu[2]=_mm_sub_epi16(W[2], predYUV0[2]);
				mu[3]=_mm_sub_epi16(W[3], predYUV0[3]);
				mu[4]=_mm_sub_epi16(W[4], predYUV0[4]);
				mu[5]=_mm_sub_epi16(W[5], predYUV0[5]);
				for(k=0;k<L1_NPREDS2;++k)//update
				{
					__m128i mc[12];
					mc[0x0]=_mm_sign_epi16(L1preds[k*6+0], mu[0]);//L1
					mc[0x1]=_mm_sign_epi16(L1preds[k*6+1], mu[1]);
					mc[0x2]=_mm_sign_epi16(L1preds[k*6+2], mu[2]);
					mc[0x3]=_mm_sign_epi16(L1preds[k*6+3], mu[3]);
					mc[0x4]=_mm_sign_epi16(L1preds[k*6+4], mu[4]);
					mc[0x5]=_mm_sign_epi16(L1preds[k*6+5], mu[5]);
					//mc[0x0]=_mm_mullo_epi16(L1preds[k*6+0], mu[0]);//L2
					//mc[0x1]=_mm_mullo_epi16(L1preds[k*6+1], mu[1]);
					//mc[0x2]=_mm_mullo_epi16(L1preds[k*6+2], mu[2]);
					//mc[0x3]=_mm_mullo_epi16(L1preds[k*6+3], mu[3]);
					//mc[0x4]=_mm_mullo_epi16(L1preds[k*6+4], mu[4]);
					//mc[0x5]=_mm_mullo_epi16(L1preds[k*6+5], mu[5]);

					//signed 16 -> 32	6 lo 6 hi registers
					mc[0x6]=_mm_srai_epi32(mc[0x0], 16);
					mc[0x7]=_mm_srai_epi32(mc[0x1], 16);
					mc[0x8]=_mm_srai_epi32(mc[0x2], 16);
					mc[0x9]=_mm_srai_epi32(mc[0x3], 16);
					mc[0xA]=_mm_srai_epi32(mc[0x4], 16);
					mc[0xB]=_mm_srai_epi32(mc[0x5], 16);

					mc[0x0]=_mm_slli_epi32(mc[0x0], 16);
					mc[0x1]=_mm_slli_epi32(mc[0x1], 16);
					mc[0x2]=_mm_slli_epi32(mc[0x2], 16);
					mc[0x3]=_mm_slli_epi32(mc[0x3], 16);
					mc[0x4]=_mm_slli_epi32(mc[0x4], 16);
					mc[0x5]=_mm_slli_epi32(mc[0x5], 16);

					mc[0x0]=_mm_srai_epi32(mc[0x0], 16);
					mc[0x1]=_mm_srai_epi32(mc[0x1], 16);
					mc[0x2]=_mm_srai_epi32(mc[0x2], 16);
					mc[0x3]=_mm_srai_epi32(mc[0x3], 16);
					mc[0x4]=_mm_srai_epi32(mc[0x4], 16);
					mc[0x5]=_mm_srai_epi32(mc[0x5], 16);
					mc[0x0]=_mm_add_epi32(mc[0x0], _mm_load_si128((__m128i*)L1weights+k*12+0x0));//update coeffs
					mc[0x1]=_mm_add_epi32(mc[0x1], _mm_load_si128((__m128i*)L1weights+k*12+0x1));
					mc[0x2]=_mm_add_epi32(mc[0x2], _mm_load_si128((__m128i*)L1weights+k*12+0x2));
					mc[0x3]=_mm_add_epi32(mc[0x3], _mm_load_si128((__m128i*)L1weights+k*12+0x3));
					mc[0x4]=_mm_add_epi32(mc[0x4], _mm_load_si128((__m128i*)L1weights+k*12+0x4));
					mc[0x5]=_mm_add_epi32(mc[0x5], _mm_load_si128((__m128i*)L1weights+k*12+0x5));
					mc[0x6]=_mm_add_epi32(mc[0x6], _mm_load_si128((__m128i*)L1weights+k*12+0x6));
					mc[0x7]=_mm_add_epi32(mc[0x7], _mm_load_si128((__m128i*)L1weights+k*12+0x7));
					mc[0x8]=_mm_add_epi32(mc[0x8], _mm_load_si128((__m128i*)L1weights+k*12+0x8));
					mc[0x9]=_mm_add_epi32(mc[0x9], _mm_load_si128((__m128i*)L1weights+k*12+0x9));
					mc[0xA]=_mm_add_epi32(mc[0xA], _mm_load_si128((__m128i*)L1weights+k*12+0xA));
					mc[0xB]=_mm_add_epi32(mc[0xB], _mm_load_si128((__m128i*)L1weights+k*12+0xB));
					_mm_store_si128((__m128i*)L1weights+k*12+0x0, mc[0x0]);
					_mm_store_si128((__m128i*)L1weights+k*12+0x1, mc[0x1]);
					_mm_store_si128((__m128i*)L1weights+k*12+0x2, mc[0x2]);
					_mm_store_si128((__m128i*)L1weights+k*12+0x3, mc[0x3]);
					_mm_store_si128((__m128i*)L1weights+k*12+0x4, mc[0x4]);
					_mm_store_si128((__m128i*)L1weights+k*12+0x5, mc[0x5]);
					_mm_store_si128((__m128i*)L1weights+k*12+0x6, mc[0x6]);
					_mm_store_si128((__m128i*)L1weights+k*12+0x7, mc[0x7]);
					_mm_store_si128((__m128i*)L1weights+k*12+0x8, mc[0x8]);
					_mm_store_si128((__m128i*)L1weights+k*12+0x9, mc[0x9]);
					_mm_store_si128((__m128i*)L1weights+k*12+0xA, mc[0xA]);
					_mm_store_si128((__m128i*)L1weights+k*12+0xB, mc[0xB]);
				}
			}
			else if(effort==3)//update
			{
				__m128i mu[6];
				int k;

				mu[0]=_mm_sub_epi16(W[0], predYUV0[0]);
				mu[1]=_mm_sub_epi16(W[1], predYUV0[1]);
				mu[2]=_mm_sub_epi16(W[2], predYUV0[2]);
				mu[3]=_mm_sub_epi16(W[3], predYUV0[3]);
				mu[4]=_mm_sub_epi16(W[4], predYUV0[4]);
				mu[5]=_mm_sub_epi16(W[5], predYUV0[5]);
				for(k=0;k<L1_NPREDS3;++k)//update
				{
					__m128i mc[12];
					mc[0x0]=_mm_sign_epi16(L1preds[k*6+0], mu[0]);//L1
					mc[0x1]=_mm_sign_epi16(L1preds[k*6+1], mu[1]);
					mc[0x2]=_mm_sign_epi16(L1preds[k*6+2], mu[2]);
					mc[0x3]=_mm_sign_epi16(L1preds[k*6+3], mu[3]);
					mc[0x4]=_mm_sign_epi16(L1preds[k*6+4], mu[4]);
					mc[0x5]=_mm_sign_epi16(L1preds[k*6+5], mu[5]);
					//mc[0x0]=_mm_mullo_epi16(L1preds[k*6+0], mu[0]);//L2
					//mc[0x1]=_mm_mullo_epi16(L1preds[k*6+1], mu[1]);
					//mc[0x2]=_mm_mullo_epi16(L1preds[k*6+2], mu[2]);
					//mc[0x3]=_mm_mullo_epi16(L1preds[k*6+3], mu[3]);
					//mc[0x4]=_mm_mullo_epi16(L1preds[k*6+4], mu[4]);
					//mc[0x5]=_mm_mullo_epi16(L1preds[k*6+5], mu[5]);

					//signed 16 -> 32	6 lo 6 hi registers
					mc[0x6]=_mm_srai_epi32(mc[0x0], 16);
					mc[0x7]=_mm_srai_epi32(mc[0x1], 16);
					mc[0x8]=_mm_srai_epi32(mc[0x2], 16);
					mc[0x9]=_mm_srai_epi32(mc[0x3], 16);
					mc[0xA]=_mm_srai_epi32(mc[0x4], 16);
					mc[0xB]=_mm_srai_epi32(mc[0x5], 16);

					mc[0x0]=_mm_slli_epi32(mc[0x0], 16);
					mc[0x1]=_mm_slli_epi32(mc[0x1], 16);
					mc[0x2]=_mm_slli_epi32(mc[0x2], 16);
					mc[0x3]=_mm_slli_epi32(mc[0x3], 16);
					mc[0x4]=_mm_slli_epi32(mc[0x4], 16);
					mc[0x5]=_mm_slli_epi32(mc[0x5], 16);

					mc[0x0]=_mm_srai_epi32(mc[0x0], 16);
					mc[0x1]=_mm_srai_epi32(mc[0x1], 16);
					mc[0x2]=_mm_srai_epi32(mc[0x2], 16);
					mc[0x3]=_mm_srai_epi32(mc[0x3], 16);
					mc[0x4]=_mm_srai_epi32(mc[0x4], 16);
					mc[0x5]=_mm_srai_epi32(mc[0x5], 16);
					mc[0x0]=_mm_add_epi32(mc[0x0], _mm_load_si128((__m128i*)L1weights+k*12+0x0));//update coeffs
					mc[0x1]=_mm_add_epi32(mc[0x1], _mm_load_si128((__m128i*)L1weights+k*12+0x1));
					mc[0x2]=_mm_add_epi32(mc[0x2], _mm_load_si128((__m128i*)L1weights+k*12+0x2));
					mc[0x3]=_mm_add_epi32(mc[0x3], _mm_load_si128((__m128i*)L1weights+k*12+0x3));
					mc[0x4]=_mm_add_epi32(mc[0x4], _mm_load_si128((__m128i*)L1weights+k*12+0x4));
					mc[0x5]=_mm_add_epi32(mc[0x5], _mm_load_si128((__m128i*)L1weights+k*12+0x5));
					mc[0x6]=_mm_add_epi32(mc[0x6], _mm_load_si128((__m128i*)L1weights+k*12+0x6));
					mc[0x7]=_mm_add_epi32(mc[0x7], _mm_load_si128((__m128i*)L1weights+k*12+0x7));
					mc[0x8]=_mm_add_epi32(mc[0x8], _mm_load_si128((__m128i*)L1weights+k*12+0x8));
					mc[0x9]=_mm_add_epi32(mc[0x9], _mm_load_si128((__m128i*)L1weights+k*12+0x9));
					mc[0xA]=_mm_add_epi32(mc[0xA], _mm_load_si128((__m128i*)L1weights+k*12+0xA));
					mc[0xB]=_mm_add_epi32(mc[0xB], _mm_load_si128((__m128i*)L1weights+k*12+0xB));
					_mm_store_si128((__m128i*)L1weights+k*12+0x0, mc[0x0]);
					_mm_store_si128((__m128i*)L1weights+k*12+0x1, mc[0x1]);
					_mm_store_si128((__m128i*)L1weights+k*12+0x2, mc[0x2]);
					_mm_store_si128((__m128i*)L1weights+k*12+0x3, mc[0x3]);
					_mm_store_si128((__m128i*)L1weights+k*12+0x4, mc[0x4]);
					_mm_store_si128((__m128i*)L1weights+k*12+0x5, mc[0x5]);
					_mm_store_si128((__m128i*)L1weights+k*12+0x6, mc[0x6]);
					_mm_store_si128((__m128i*)L1weights+k*12+0x7, mc[0x7]);
					_mm_store_si128((__m128i*)L1weights+k*12+0x8, mc[0x8]);
					_mm_store_si128((__m128i*)L1weights+k*12+0x9, mc[0x9]);
					_mm_store_si128((__m128i*)L1weights+k*12+0xA, mc[0xA]);
					_mm_store_si128((__m128i*)L1weights+k*12+0xB, mc[0xB]);
				}
			}
			//context update = (2*eW+(e<<3)+max(eNEE, eNEEE))>>2
			eNEEE[0]=_mm_load_si128((__m128i*)rows[1]+6+0+3*12);
			eNEEE[1]=_mm_load_si128((__m128i*)rows[1]+6+1+3*12);
			eNEEE[2]=_mm_load_si128((__m128i*)rows[1]+6+2+3*12);
			eNEEE[3]=_mm_load_si128((__m128i*)rows[1]+6+3+3*12);
			eNEEE[4]=_mm_load_si128((__m128i*)rows[1]+6+4+3*12);
			eNEEE[5]=_mm_load_si128((__m128i*)rows[1]+6+5+3*12);

			ecurr[0]=_mm_slli_epi16(ecurr[0], GRBITS);
			ecurr[1]=_mm_slli_epi16(ecurr[1], GRBITS);
			ecurr[2]=_mm_slli_epi16(ecurr[2], GRBITS);
			ecurr[3]=_mm_slli_epi16(ecurr[3], GRBITS);
			ecurr[4]=_mm_slli_epi16(ecurr[4], GRBITS);
			ecurr[5]=_mm_slli_epi16(ecurr[5], GRBITS);
			ecurr[0]=_mm_avg_epu16(ecurr[0], _mm_max_epi16(eNEE[0], eNEEE[0]));
			ecurr[1]=_mm_avg_epu16(ecurr[1], _mm_max_epi16(eNEE[1], eNEEE[1]));
			ecurr[2]=_mm_avg_epu16(ecurr[2], _mm_max_epi16(eNEE[2], eNEEE[2]));
			ecurr[3]=_mm_avg_epu16(ecurr[3], _mm_max_epi16(eNEE[3], eNEEE[3]));
			ecurr[4]=_mm_avg_epu16(ecurr[4], _mm_max_epi16(eNEE[4], eNEEE[4]));
			ecurr[5]=_mm_avg_epu16(ecurr[5], _mm_max_epi16(eNEE[5], eNEEE[5]));
			eW[0]=_mm_avg_epu16(eW[0], ecurr[0]);
			eW[1]=_mm_avg_epu16(eW[1], ecurr[1]);
			eW[2]=_mm_avg_epu16(eW[2], ecurr[2]);
			eW[3]=_mm_avg_epu16(eW[3], ecurr[3]);
			eW[4]=_mm_avg_epu16(eW[4], ecurr[4]);
			eW[5]=_mm_avg_epu16(eW[5], ecurr[5]);

			_mm_store_si128((__m128i*)rows[0]+6+0+0*12, eW[0]);//store current contexts
			_mm_store_si128((__m128i*)rows[0]+6+1+0*12, eW[1]);
			_mm_store_si128((__m128i*)rows[0]+6+2+0*12, eW[2]);
			_mm_store_si128((__m128i*)rows[0]+6+3+0*12, eW[3]);
			_mm_store_si128((__m128i*)rows[0]+6+4+0*12, eW[4]);
			_mm_store_si128((__m128i*)rows[0]+6+5+0*12, eW[5]);
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
			rows[0]+=6*NCODERS;
			rows[1]+=6*NCODERS;
			rows[2]+=6*NCODERS;
			rows[3]+=6*NCODERS;
			imptr+=3*NCODERS;
		}
	}
	prof_checkpoint(isize, "main");
#ifdef ENABLE_GUIDE
	if(dist>1&&!fwd)
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
		rANS_SIMD_SymInfo *rsyminfo=(rANS_SIMD_SymInfo*)rCDF2syms;
		int kx, ky;
		unsigned short *ctxptr2=0;

		//normalize/integrate hists
		{
			int kc;
			for(kc=0;kc<3*NCTX;++kc)
				enc_hist2stats(hists+(ptrdiff_t)256*kc, syminfo+(ptrdiff_t)256*kc, &ctxmask, kc);
		}
		
		//encode remainder
		if(xremw||yremh)
		{
			memset(rhist, 0, rhsize);
			for(ky=0;ky<yremh;++ky)
				decorr1d(image+rowstride*(blockh*YCODERS+ky), iw, 3, bestrct, rhist);
			for(kx=0;kx<xremw;++kx)
				decorr1d(image+qxbytes+3*kx, blockh*YCODERS, rowstride, bestrct, rhist);
			enc_hist2stats(rhist+(ptrdiff_t)256*0, rsyminfo+(ptrdiff_t)256*0, &ctxmask, 3*NCTX+0);
			enc_hist2stats(rhist+(ptrdiff_t)256*1, rsyminfo+(ptrdiff_t)256*1, &ctxmask, 3*NCTX+1);
			enc_hist2stats(rhist+(ptrdiff_t)256*2, rsyminfo+(ptrdiff_t)256*2, &ctxmask, 3*NCTX+2);
			
			{
				unsigned state=1<<(RANS_STATE_BITS-RANS_RENORM_BITS);
				for(kx=xremw-1;kx>=0;--kx)
					encode1d(image+qxbytes+3*kx, blockh*YCODERS, rowstride, &state, &streamptr, image, rsyminfo);
				for(ky=yremh-1;ky>=0;--ky)
					encode1d(image+rowstride*(blockh*YCODERS+ky), iw, 3, &state, &streamptr, image, rsyminfo);
				//flush
				streamptr-=4;
#ifdef _DEBUG
				if(streamptr<=image)
					LOG_ERROR("OOB ptr %016zX <= %016zX", streamptr, image);
#endif
				*(unsigned*)streamptr=state;
			}
			prof_checkpoint(usize-isize, "encode remainder");
			profile_size(streamptr, "/ %9d bytes remainders", (int)(usize-isize));
		}

		//encode main
		mstate[3]=mstate[2]=mstate[1]=mstate[0]=_mm_set1_epi32(1<<(RANS_STATE_BITS-RANS_RENORM_BITS));
		ctxptr2=(unsigned short*)(interleaved+(isize<<1)-sizeof(short[NCODERS]));
		for(ky=blockh-1;ky>=0;--ky)
		{
#ifdef ESTIMATE_SIZE
			int kc=2;
#endif
			for(kx=ixbytes-NCODERS;kx>=0;kx-=NCODERS)//ixbytes = iw/XCODERS*NCODERS*3
			{
				__m128i mmax[4], minvf[4], mcdf[4], mnegf_sh[4];
				{
					__m128i s0, s1, s2, s3;
					__m128i t0, t1, t2, t3;
#define SHUFFLE_PS(LO, HI, IMM8_HHLL) _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(LO), _mm_castsi128_ps(HI), IMM8_HHLL))

					s0=_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*4+0]));
					s1=_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*4+1]));
					s2=_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*4+2]));
					s3=_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*4+3]));
					t0=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(1, 0, 1, 0));//_MM_TRANSPOSE4_PS
					t2=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(3, 2, 3, 2));
					t1=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(1, 0, 1, 0));
					t3=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(3, 2, 3, 2));
					mmax	[0]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(2, 0, 2, 0));
					minvf	[0]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(3, 1, 3, 1));
					mcdf	[0]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(2, 0, 2, 0));
					mnegf_sh[0]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(3, 1, 3, 1));

					s0=_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*4+0]));
					s1=_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*4+1]));
					s2=_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*4+2]));
					s3=_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*4+3]));
					t0=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(1, 0, 1, 0));
					t2=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(3, 2, 3, 2));
					t1=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(1, 0, 1, 0));
					t3=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(3, 2, 3, 2));
					mmax	[1]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(2, 0, 2, 0));
					minvf	[1]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(3, 1, 3, 1));
					mcdf	[1]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(2, 0, 2, 0));
					mnegf_sh[1]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(3, 1, 3, 1));

					s0=_mm_load_si128((__m128i*)(syminfo+ctxptr2[2*4+0]));
					s1=_mm_load_si128((__m128i*)(syminfo+ctxptr2[2*4+1]));
					s2=_mm_load_si128((__m128i*)(syminfo+ctxptr2[2*4+2]));
					s3=_mm_load_si128((__m128i*)(syminfo+ctxptr2[2*4+3]));
					t0=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(1, 0, 1, 0));
					t2=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(3, 2, 3, 2));
					t1=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(1, 0, 1, 0));
					t3=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(3, 2, 3, 2));
					mmax	[2]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(2, 0, 2, 0));
					minvf	[2]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(3, 1, 3, 1));
					mcdf	[2]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(2, 0, 2, 0));
					mnegf_sh[2]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(3, 1, 3, 1));

					s0=_mm_load_si128((__m128i*)(syminfo+ctxptr2[3*4+0]));
					s1=_mm_load_si128((__m128i*)(syminfo+ctxptr2[3*4+1]));
					s2=_mm_load_si128((__m128i*)(syminfo+ctxptr2[3*4+2]));
					s3=_mm_load_si128((__m128i*)(syminfo+ctxptr2[3*4+3]));
					t0=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(1, 0, 1, 0));
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
					ALIGN(32) int freqs[NCODERS]={0};
					memcpy(freqs, mmax, sizeof(freqs));
					const double norm=1./(1<<PROBBITS);
					for(int k=0;k<NCODERS;++k)
					{
						int freq=(freqs[k]+1)>>(RANS_STATE_BITS-PROBBITS);
						if((unsigned)(freq-1)>=(unsigned)((1<<PROBBITS)-1))
							LOG_ERROR("freq = %d", freq);
						esize[kc*NCODERS+k]-=log2(freq*norm)*0.125;
					}
				}
				--kc;
				if(kc<0)
					kc=2;
#endif
				//enc renorm		if(state>(freq<<(31-12))-1){write16(state); state>>=16;}
				{
					__m128i cond0=_mm_cmpgt_epi32(mstate[0], mmax[0]);
					__m128i cond1=_mm_cmpgt_epi32(mstate[1], mmax[1]);
					__m128i cond2=_mm_cmpgt_epi32(mstate[2], mmax[2]);
					__m128i cond3=_mm_cmpgt_epi32(mstate[3], mmax[3]);
					int mask0=_mm_movemask_ps(_mm_castsi128_ps(cond0));
					int mask1=_mm_movemask_ps(_mm_castsi128_ps(cond1));
					int mask2=_mm_movemask_ps(_mm_castsi128_ps(cond2));
					int mask3=_mm_movemask_ps(_mm_castsi128_ps(cond3));
					__m128i idx0=_mm_load_si128((__m128i*)ans_permute_enc+mask0);
					__m128i idx1=_mm_load_si128((__m128i*)ans_permute_enc+mask1);
					__m128i idx2=_mm_load_si128((__m128i*)ans_permute_enc+mask2);
					__m128i idx3=_mm_load_si128((__m128i*)ans_permute_enc+mask3);
					__m128i emit0=_mm_shuffle_epi8(mstate[0], idx0);
					__m128i emit1=_mm_shuffle_epi8(mstate[1], idx1);
					__m128i emit2=_mm_shuffle_epi8(mstate[2], idx2);
					__m128i emit3=_mm_shuffle_epi8(mstate[3], idx3);
					mask0=popcnt_table[mask0];
					mask1=popcnt_table[mask1];
					mask2=popcnt_table[mask2];
					mask3=popcnt_table[mask3];
#ifdef _DEBUG
					if(streamptr-(2*((ptrdiff_t)mask0+mask1+mask2+mask3)+sizeof(__m128i))<=image)
						LOG_ERROR("OOB ptr %016zX <= %016zX", streamptr, image);
#endif
					_mm_storeu_si128((__m128i*)streamptr-1, emit3); streamptr-=mask3;
					_mm_storeu_si128((__m128i*)streamptr-1, emit2); streamptr-=mask2;
					_mm_storeu_si128((__m128i*)streamptr-1, emit1); streamptr-=mask1;
					_mm_storeu_si128((__m128i*)streamptr-1, emit0); streamptr-=mask0;
					{
						__m128i state0=_mm_srli_epi32(mstate[0], 16);
						__m128i state1=_mm_srli_epi32(mstate[1], 16);
						__m128i state2=_mm_srli_epi32(mstate[2], 16);
						__m128i state3=_mm_srli_epi32(mstate[3], 16);
						mstate[0]=_mm_blendv_epi8(mstate[0], state0, cond0);
						mstate[1]=_mm_blendv_epi8(mstate[1], state1, cond1);
						mstate[2]=_mm_blendv_epi8(mstate[2], state2, cond2);
						mstate[3]=_mm_blendv_epi8(mstate[3], state3, cond3);
					}
				}
#ifdef ANS_VAL
				ansval_push(mstate, sizeof(int), NCODERS);
#endif
				//enc update		state += (state*invf>>32>>sh)*negf+cdf		state = state/freq<<12|(cdf+state%freq)
#if 1
				{
					//state += ((state*invf>>32)*(1<<(11-sh))>>11)*negf+cdf
					__m128i lo0=_mm_mul_epu32(mstate[0], minvf[0]);//q = mulhi32(state, invf)
					__m128i lo1=_mm_mul_epu32(mstate[1], minvf[1]);
					__m128i lo2=_mm_mul_epu32(mstate[2], minvf[2]);
					__m128i lo3=_mm_mul_epu32(mstate[3], minvf[3]);
					__m128i hi0=_mm_mul_epu32(_mm_srli_epi64(mstate[0], 32), _mm_srli_epi64(minvf[0], 32));
					__m128i hi1=_mm_mul_epu32(_mm_srli_epi64(mstate[1], 32), _mm_srli_epi64(minvf[1], 32));
					__m128i hi2=_mm_mul_epu32(_mm_srli_epi64(mstate[2], 32), _mm_srli_epi64(minvf[2], 32));
					__m128i hi3=_mm_mul_epu32(_mm_srli_epi64(mstate[3], 32), _mm_srli_epi64(minvf[3], 32));
					__m128i sh0=_mm_srli_epi32(mnegf_sh[0], 16);
					__m128i sh1=_mm_srli_epi32(mnegf_sh[1], 16);
					__m128i sh2=_mm_srli_epi32(mnegf_sh[2], 16);
					__m128i sh3=_mm_srli_epi32(mnegf_sh[3], 16);
					lo0=_mm_srli_epi64(lo0, 32);
					lo1=_mm_srli_epi64(lo1, 32);
					lo2=_mm_srli_epi64(lo2, 32);
					lo3=_mm_srli_epi64(lo3, 32);
					hi0=_mm_srli_epi64(hi0, 32);
					hi1=_mm_srli_epi64(hi1, 32);
					hi2=_mm_srli_epi64(hi2, 32);
					hi3=_mm_srli_epi64(hi3, 32);
					lo0=_mm_mul_epu32(lo0, sh0);
					lo1=_mm_mul_epu32(lo1, sh1);
					lo2=_mm_mul_epu32(lo2, sh2);
					lo3=_mm_mul_epu32(lo3, sh3);
					sh0=_mm_srli_epi64(sh0, 32);
					sh1=_mm_srli_epi64(sh1, 32);
					sh2=_mm_srli_epi64(sh2, 32);
					sh3=_mm_srli_epi64(sh3, 32);
					hi0=_mm_mul_epu32(hi0, sh0);
					hi1=_mm_mul_epu32(hi1, sh1);
					hi2=_mm_mul_epu32(hi2, sh2);
					hi3=_mm_mul_epu32(hi3, sh3);
					lo0=_mm_srli_epi64(lo0, PROBBITS-1);
					lo1=_mm_srli_epi64(lo1, PROBBITS-1);
					lo2=_mm_srli_epi64(lo2, PROBBITS-1);
					lo3=_mm_srli_epi64(lo3, PROBBITS-1);
					hi0=_mm_slli_epi64(hi0, 32-(PROBBITS-1));
					hi1=_mm_slli_epi64(hi1, 32-(PROBBITS-1));
					hi2=_mm_slli_epi64(hi2, 32-(PROBBITS-1));
					hi3=_mm_slli_epi64(hi3, 32-(PROBBITS-1));
					minvf[0]=_mm_blend_epi16(lo0, hi0, 0xCC);
					minvf[1]=_mm_blend_epi16(lo1, hi1, 0xCC);
					minvf[2]=_mm_blend_epi16(lo2, hi2, 0xCC);
					minvf[3]=_mm_blend_epi16(lo3, hi3, 0xCC);
				}
				mstate[0]=_mm_add_epi32(mstate[0], mcdf[0]);
				mstate[1]=_mm_add_epi32(mstate[1], mcdf[1]);
				mstate[2]=_mm_add_epi32(mstate[2], mcdf[2]);
				mstate[3]=_mm_add_epi32(mstate[3], mcdf[3]);
				{
					__m128i lomask=_mm_set1_epi32(0xFFFF);
					__m128i negf0=_mm_and_si128(mnegf_sh[0], lomask);
					__m128i negf1=_mm_and_si128(mnegf_sh[1], lomask);
					__m128i negf2=_mm_and_si128(mnegf_sh[2], lomask);
					__m128i negf3=_mm_and_si128(mnegf_sh[3], lomask);
					minvf[0]=_mm_mullo_epi32(minvf[0], negf0);
					minvf[1]=_mm_mullo_epi32(minvf[1], negf1);
					minvf[2]=_mm_mullo_epi32(minvf[2], negf2);
					minvf[3]=_mm_mullo_epi32(minvf[3], negf3);
				}
#ifdef ANS_VAL
				{
					__m128i one=_mm_set1_epi32(1);
					mmax[0]=_mm_add_epi32(mmax[0], one);
					mmax[1]=_mm_add_epi32(mmax[1], one);
					mmax[2]=_mm_add_epi32(mmax[2], one);
					mmax[3]=_mm_add_epi32(mmax[3], one);
					mmax[0]=_mm_srli_epi32(mmax[0], RANS_STATE_BITS-PROBBITS);
					mmax[1]=_mm_srli_epi32(mmax[1], RANS_STATE_BITS-PROBBITS);
					mmax[2]=_mm_srli_epi32(mmax[2], RANS_STATE_BITS-PROBBITS);
					mmax[3]=_mm_srli_epi32(mmax[3], RANS_STATE_BITS-PROBBITS);
					mmax[0]=_mm_packus_epi32(mmax[0], mmax[1]);
					mmax[2]=_mm_packus_epi32(mmax[2], mmax[3]);
					mmax[1]=mmax[2];
					ansval_push(mmax, sizeof(short), NCODERS);
				}
#endif
#else
				{
					__m128i lo0=_mm_mul_epu32(mstate[0], minvf[0]);//q = mulhi32(state, invf)
					__m128i lo1=_mm_mul_epu32(mstate[1], minvf[1]);
					__m128i lo2=_mm_mul_epu32(mstate[2], minvf[2]);
					__m128i lo3=_mm_mul_epu32(mstate[3], minvf[3]);
					__m128i hi0=_mm_mul_epu32(_mm_srli_epi64(mstate[0], 32), _mm_srli_epi64(minvf[0], 32));
					__m128i hi1=_mm_mul_epu32(_mm_srli_epi64(mstate[1], 32), _mm_srli_epi64(minvf[1], 32));
					__m128i hi2=_mm_mul_epu32(_mm_srli_epi64(mstate[2], 32), _mm_srli_epi64(minvf[2], 32));
					__m128i hi3=_mm_mul_epu32(_mm_srli_epi64(mstate[3], 32), _mm_srli_epi64(minvf[3], 32));
					minvf[0]=_mm_blend_epi16(_mm_srli_epi64(lo0, 32), hi0, 0xCC);
					minvf[1]=_mm_blend_epi16(_mm_srli_epi64(lo1, 32), hi1, 0xCC);
					minvf[2]=_mm_blend_epi16(_mm_srli_epi64(lo2, 32), hi2, 0xCC);
					minvf[3]=_mm_blend_epi16(_mm_srli_epi64(lo3, 32), hi3, 0xCC);
				}
				{
#if 1
					ALIGN(16) volatile short sh[32];
					ALIGN(16) volatile unsigned states[16];

					_mm_store_si128((__m128i*)sh+0, mnegf_sh[0]);
					_mm_store_si128((__m128i*)sh+1, mnegf_sh[1]);
					_mm_store_si128((__m128i*)sh+2, mnegf_sh[2]);
					_mm_store_si128((__m128i*)sh+3, mnegf_sh[3]);
					_mm_store_si128((__m128i*)states+0, minvf[0]);
					_mm_store_si128((__m128i*)states+1, minvf[1]);
					_mm_store_si128((__m128i*)states+2, minvf[2]);
					_mm_store_si128((__m128i*)states+3, minvf[3]);
					states[0x0]>>=sh[0x0*2+1];
					states[0x1]>>=sh[0x1*2+1];
					states[0x2]>>=sh[0x2*2+1];
					states[0x3]>>=sh[0x3*2+1];
					states[0x4]>>=sh[0x4*2+1];
					states[0x5]>>=sh[0x5*2+1];
					states[0x6]>>=sh[0x6*2+1];
					states[0x7]>>=sh[0x7*2+1];
					states[0x8]>>=sh[0x8*2+1];
					states[0x9]>>=sh[0x9*2+1];
					states[0xA]>>=sh[0xA*2+1];
					states[0xB]>>=sh[0xB*2+1];
					states[0xC]>>=sh[0xC*2+1];
					states[0xD]>>=sh[0xD*2+1];
					states[0xE]>>=sh[0xE*2+1];
					states[0xF]>>=sh[0xF*2+1];
					minvf[0]=_mm_load_si128((__m128i*)states+0);
					minvf[1]=_mm_load_si128((__m128i*)states+1);
					minvf[2]=_mm_load_si128((__m128i*)states+2);
					minvf[3]=_mm_load_si128((__m128i*)states+3);
#else
					__m128i sh0=_mm_srli_epi32(mnegf_sh[0], 16);
					__m128i sh1=_mm_srli_epi32(mnegf_sh[1], 16);
					__m128i sh2=_mm_srli_epi32(mnegf_sh[2], 16);
					__m128i sh3=_mm_srli_epi32(mnegf_sh[3], 16);
					minvf[0]=_mm_srlv_epi32(minvf[0], sh0);
					minvf[1]=_mm_srlv_epi32(minvf[1], sh1);
					minvf[2]=_mm_srlv_epi32(minvf[2], sh2);
					minvf[3]=_mm_srlv_epi32(minvf[3], sh3);
#endif
				}
				mstate[0]=_mm_add_epi32(mstate[0], mcdf[0]);
				mstate[1]=_mm_add_epi32(mstate[1], mcdf[1]);
				mstate[2]=_mm_add_epi32(mstate[2], mcdf[2]);
				mstate[3]=_mm_add_epi32(mstate[3], mcdf[3]);
				{
					__m128i lomask=_mm_set1_epi32(0xFFFF);
					__m128i negf0=_mm_and_si128(mnegf_sh[0], lomask);
					__m128i negf1=_mm_and_si128(mnegf_sh[1], lomask);
					__m128i negf2=_mm_and_si128(mnegf_sh[2], lomask);
					__m128i negf3=_mm_and_si128(mnegf_sh[3], lomask);
					minvf[0]=_mm_mullo_epi32(minvf[0], negf0);
					minvf[1]=_mm_mullo_epi32(minvf[1], negf1);
					minvf[2]=_mm_mullo_epi32(minvf[2], negf2);
					minvf[3]=_mm_mullo_epi32(minvf[3], negf3);
#ifdef ANS_VAL
					__m128i mM=_mm_set1_epi32(1<<PROBBITS);
					negf0=_mm_sub_epi16(mM, negf0);
					negf1=_mm_sub_epi16(mM, negf1);
					negf2=_mm_sub_epi16(mM, negf2);
					negf3=_mm_sub_epi16(mM, negf3);
					negf0=_mm_packus_epi32(negf0, negf1);
					negf2=_mm_packus_epi32(negf2, negf3);
					{
						ALIGN(32) unsigned short freqs[NCODERS];
						_mm_store_si128((__m128i*)freqs+0, negf0);
						_mm_store_si128((__m128i*)freqs+1, negf2);
						ansval_push(freqs, sizeof(short), NCODERS);
					}
#endif
				}
#endif
				mstate[0]=_mm_add_epi32(mstate[0], minvf[0]);
				mstate[1]=_mm_add_epi32(mstate[1], minvf[1]);
				mstate[2]=_mm_add_epi32(mstate[2], minvf[2]);
				mstate[3]=_mm_add_epi32(mstate[3], minvf[3]);
#ifdef ANS_VAL
				ansval_push(mstate, sizeof(int), NCODERS);
#endif
			}
		}
		//flush
		streamptr-=sizeof(mstate);
#ifdef _DEBUG
		if(streamptr<=image)
			LOG_ERROR("OOB ptr %016zX <= %016zX", streamptr, image);
#endif
		memcpy(streamptr, mstate, sizeof(mstate));
		prof_checkpoint(isize, "encode main");
		profile_size(streamptr, "/ %9d bytes main", (int)isize);

		//pack hists
		{
			BitPackerLIFO ec;
			bitpacker_enc_init(&ec, image, streamptr);
			if(xremw||yremh)
			{
				enc_packhist(&ec, rhist+256*2, ctxmask, 3*NCTX+2);
				enc_packhist(&ec, rhist+256*1, ctxmask, 3*NCTX+1);
				enc_packhist(&ec, rhist+256*0, ctxmask, 3*NCTX+0);
			}
			{
				int kc;
				for(kc=3*NCTX-1;kc>=0;--kc)
					enc_packhist(&ec, hists+(ptrdiff_t)256*kc, ctxmask, kc);
			}
			bitpacker_enc_flush(&ec);
			streamptr=ec.dstbwdptr;

			streamptr-=8;
			*(unsigned long long*)streamptr=ctxmask;
		}
		prof_checkpoint(((ptrdiff_t)3*NCTX+((xremw||yremh)?3:0))*256, "pack histograms");
		profile_size(streamptr, "/ %9d bytes overhead", (3*NCTX+3)*12<<8>>3);

		//save compressed file
		{
			ptrdiff_t csize2=0;
			FILE *fdst=fopen(dstfn, "wb");
			if(!fdst)
			{
				LOG_ERROR("Cannot open \"%s\" for writing", fdst);
				return 1;
			}
			csize2+=fwrite("32", 1, 2, fdst);
			csize2+=fwrite(&iw, 1, 4, fdst);
			csize2+=fwrite(&ih, 1, 4, fdst);
			{
				int flags=bestrct<<2|(effort&3);
				csize2+=fwrite(&flags, 1, 1, fdst);
			}
			csize2+=fwrite(&dist, 1, 1, fdst);
#ifdef _DEBUG
			if(streamptr>streamstart)
				LOG_ERROR("OOB ptr %016zX > %016zX", streamptr, streamstart);
			if(streamptr<image)
				LOG_ERROR("OOB ptr %016zX < %016zX", streamptr, image);
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
			{
				ptrdiff_t usize2=get_filesize(srcfn);
				printf("%s  WH %d*%d\n", srcfn, iw, ih);
				printf("%8d/%8d bytes\n", (int)csize2, (int)usize2);
			}
#endif
			(void)csize2;
			prof_checkpoint(csize2, "fwrite");
		}
		free(hists);
		free(rhist);
	}
	else
	{
		//deinterleave
		interleave_blocks_inv(interleaved, iw, ih, image);
		prof_checkpoint(usize, "deinterleave");

		if(xremw||yremh)
		{
#ifdef _DEBUG
			if(streamptr>streamend)
				LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
			{
				int kx, ky;
				unsigned state=*(unsigned*)streamptr;
				streamptr+=4;
				for(ky=0;ky<yremh;++ky)
					decode1d(image+rowstride*(blockh*YCODERS+ky), iw, 3, bestrct, &state, (const unsigned char**)&streamptr, streamend, rCDF2syms);
				for(kx=0;kx<xremw;++kx)
					decode1d(image+qxbytes+3*kx, blockh*YCODERS, rowstride, bestrct, &state, (const unsigned char**)&streamptr, streamend, rCDF2syms);
			}
			prof_checkpoint(usize-isize, "remainder");
		}

		//save PPM file
		save_ppm(dstfn, image, iw, ih);
		prof_checkpoint(usize, "fwrite");
	}
	_mm_free(pixels);
	_mm_free(CDF2syms);
	_mm_free(rCDF2syms);
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
#ifdef PRINT_SHIFTBOUNDS
	if(fwd)
		printf("sh %d~%d\n", minsh, maxsh);
#endif
	(void)och_names;
	(void)rct_names;
	return 0;
}