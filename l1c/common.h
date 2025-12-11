//include after config macros
#pragma once
#ifndef INC_COMMON_H
#define INC_COMMON_H
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdarg.h>
#include<time.h>
#if defined _WIN32 || defined WIN32
#define WIN32_LEAN_AND_MEAN
#include<Windows.h>
#endif


#ifdef _MSC_VER
	#define ENABLE_GUIDE		//DEBUG		checks interleaved pixels

//	#define ANS_VAL			//DEBUG
#endif


//macros
#ifdef _MSC_VER
#	define ALIGN(N) __declspec(align(N))
#	define AWM_INLINE __forceinline static
#if _MSC_VER<1900
#define snprintf sprintf_s
#endif
#else
#	define ALIGN(N) __attribute__((aligned(N)))
#	define AWM_INLINE __attribute__((always_inline)) inline static
#	ifndef _countof
#		define _countof(A) (sizeof(A)/sizeof(*(A)))
#	endif
#endif

#define CLAMP2(X, LO, HI)\
	do\
	{\
		if((X)<(LO))X=LO;\
		if((X)>(HI))X=HI;\
	}while(0)

#if defined _M_X64 || defined __x86_64__
#define FLOOR_LOG2(X)\
	(sizeof(X)==8?63-(int32_t)_lzcnt_u64(X):31-_lzcnt_u32((uint32_t)(X)))
#else
AWM_INLINE int floor_log2_64(uint64_t n)
{
	int	logn=-!n;
	int	sh=(n>=1ULL<<32)<<5;	logn+=sh, n>>=sh;
		sh=(n>=1<<16)<<4;	logn+=sh, n>>=sh;
		sh=(n>=1<< 8)<<3;	logn+=sh, n>>=sh;
		sh=(n>=1<< 4)<<2;	logn+=sh, n>>=sh;
		sh=(n>=1<< 2)<<1;	logn+=sh, n>>=sh;
		sh= n>=1<< 1;		logn+=sh;
	return logn;
}
AWM_INLINE int floor_log2_32(uint32_t n)
{
	int	logn=-!n;
	int	sh=(n>=1<<16)<<4;	logn+=sh, n>>=sh;
		sh=(n>=1<< 8)<<3;	logn+=sh, n>>=sh;
		sh=(n>=1<< 4)<<2;	logn+=sh, n>>=sh;
		sh=(n>=1<< 2)<<1;	logn+=sh, n>>=sh;
		sh= n>=1<< 1;		logn+=sh;
	return logn;
}
#define FLOOR_LOG2(X)\
	(sizeof(X)==8?floor_log2_64(X):floor_log2_32((uint32_t)(X)))
#endif


//runtime
static void crash(const char *file, int line, const char *format, ...)
{
	printf("%s(%d):\n", file, line);
	if(format)
	{
		va_list args;
		va_start(args, format);
		vprintf(format, args);
		va_end(args);
	}
	printf("\n");
	exit(1);
}
#define CRASH(FORMAT, ...) crash(__FILE__, __LINE__, FORMAT,##__VA_ARGS__)
static void memfill(void *dst, const void *src, size_t dstbytes, size_t srcbytes)
{
	size_t copied;
	char *d=(char*)dst;
	const char *s=(const char*)src;
	//printf("MEMFILL  %016zX %016zX %016zX %016zX\n", (size_t)dst, (size_t)src, dstbytes, srcbytes);//
#ifdef _DEBUG
	if(!dstbytes||!srcbytes)
	{
		//CRASH("memfill:  dstbytes %zd  srcbytes %zd", dstbytes, srcbytes);
		return;
	}
#endif
	if(dstbytes<srcbytes)
	{
		memcpy(dst, src, dstbytes);
		return;
	}
	copied=srcbytes;
	memcpy(d, s, srcbytes);
	while((copied<<1)<=dstbytes)
	{
		memcpy(d+copied, d, copied);
		copied<<=1;
	}
	if(copied<dstbytes)
		memcpy(d+copied, d, dstbytes-copied);
}
#define FILLMEM(PTR, DATA, ASIZE, ESIZE)\
	do\
	{\
		*(PTR)=(DATA);\
		memfill((PTR)+1, PTR, (ASIZE)-(ESIZE), ESIZE);\
	}while(0)
static double time_sec(void)
{
#if defined _WIN32 || defined WIN32
	static long long t0=0;
	LARGE_INTEGER li;
	double t;
	QueryPerformanceCounter(&li);
	if(!t0)
		t0=li.QuadPart;
	t=(double)(li.QuadPart-t0);
	QueryPerformanceFrequency(&li);
	t/=(double)li.QuadPart;
	return t;
#else
	struct timespec t;
	clock_gettime(CLOCK_REALTIME, &t);//<time.h>
	return t.tv_sec+t.tv_nsec*1e-9;
#endif
}
static int print_timestamp(const char *format)//"%Y-%m-%d_%H%M%S"
{
	char buf[1024]={0};
	time_t tstamp=time(0);
	struct tm *tformat=localtime(&tstamp);
	return (int)strftime(buf, sizeof(buf)-1, format, tformat);
}
static void colorgen(int *colors, int count, int minbrightness, int maxbrightness, int maxtrials)
{
	int brightnessrange=0, dmax0=0, k;
	static const char increments[]=
	{
		//r, g, b
		+1, -1, +0,
		-1, +1, +0,
		+1, +0, -1,
		-1, +0, +1,
		+0, +1, -1,
		+0, -1, +1,
	};

	unsigned char *data=(unsigned char*)colors;
	CLAMP2(minbrightness, 0, 765-1);
	if(maxbrightness<minbrightness+1)
		maxbrightness=minbrightness+1;
	brightnessrange=maxbrightness-minbrightness;
	dmax0=brightnessrange*brightnessrange/3>>1;
	for(k=0;k<count;++k)
	{
		int reject=0, ntrials=0;
		int bestdist=0, bestcolor=0;
		do
		{
			int rem, dmin;
			int r0=rand();//at least 15-bit
			int r1=rand();
			int r2=rand();
			int r=r0>>7;
			int g=r1>>7;
			int b=r2>>7;
			++ntrials;
			if((unsigned)(r+g+b-minbrightness)>=(unsigned)brightnessrange)
			{
				int t0;

				if(ntrials<maxtrials)
				{
					reject=1;
					continue;
				}
				r=g=b=(minbrightness+maxbrightness+3)/6;
				t0=rand();
				while(t0)
				{
					const char *inc;
					int rem=t0;
					t0/=6;
					rem-=t0*6;
					inc=increments+rem*3;
					if(!(((r+inc[0])|(g+inc[1])|(b+inc[2]))>>8))
					{
						r+=inc[0];
						g+=inc[1];
						b+=inc[2];
					}
				}
			}
			rem=(r2&127)<<14|(r1&127)<<7|(r0&127);//21 bit
			dmin=0xFFFFFF;// > 3*255*255
			{
				int k2;
				const unsigned char *p2=(const unsigned char*)data;
				for(k2=0;k2<k;++k2)
				{
					int dr=p2[0]-r;
					int dg=p2[1]-g;
					int db=p2[2]-b;
					int d=dr*dr+dg*dg+db*db;
					if(!k2||dmin>d)
						dmin=d;
					p2+=4;
				}
			}
			if(bestdist<dmin)
			{
				bestdist=dmin;
				bestcolor=b<<16|g<<8|r;
			}
			reject=((unsigned long long)dmin<<21)<(unsigned long long)rem*dmax0;
			if(ntrials>=maxtrials)
				bestcolor=0x808080;
			//	CRASH("%d trials reached, bestcolor %08X", maxtrials, bestcolor);
			//if(!reject&&(unsigned)(r+g+b-minbrightness)>=(unsigned)brightnessrange)
			//	CRASH("");
		}while(reject&&ntrials<maxtrials);
		colors[k]=bestcolor;
	}
}
#define COLORPRINTF_BK_DEFAULT 0x0C0C0C
#define COLORPRINTF_TXT_DEFAULT 0xF2F2F2
static int colorprintf(int textcolor, int bkcolor, const char *format, ...)//0x00BBGGRR
{
	int printed=0, msg=0;
	va_list args;
	char buf[2048];

	printed+=snprintf(buf+printed, sizeof(buf)-1-printed, "\33[48;2;%d;%d;%d;38;2;%d;%d;%dm",
		bkcolor&255, bkcolor>>8&255, bkcolor>>16&255,
		textcolor&255, textcolor>>8&255, textcolor>>16&255
	);
	va_start(args, format);
	msg=vsnprintf(buf+printed, sizeof(buf)-1-printed, format, args);
	printed+=msg;
	va_end(args);
	printed+=snprintf(buf+printed, sizeof(buf)-1-printed, "\33[0m");

	printf("%s", buf);

	return msg;
}


//application-specific
#ifdef ENABLE_GUIDE
static int g_iw=0, g_ih=0;
extern unsigned char *g_image;//debug.c
extern double g_sqe[3];
static void guide_save(const unsigned char *image, int iw, int ih)
{
	int size=3*iw*ih;
	g_iw=iw;
	g_ih=ih;
	g_image=(unsigned char*)malloc(size);
	if(!g_image)
	{
		CRASH("Alloc error");
		return;
	}
	memcpy(g_image, image, size);
}
static void guide_check(const unsigned char *image, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx);
	if(memcmp(image+idx, g_image+idx, 3))
	{
		CRASH("");
		printf("");
	}
}
#else
#define guide_save(...)
#define guide_check(...)
#endif

#ifdef PROFILE_SIZE
static void profile_size(const unsigned char *dstbwdptr, const char *msg, ...)
{
	static ptrdiff_t size=0;
	static const unsigned char *prev=0;
	if(prev)
	{
		ptrdiff_t diff=prev-dstbwdptr;
		size+=diff;
		printf("%10lld (%+10lld) bytes", (int64_t)size, (int64_t)diff);
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
#define PROF_CAP 64
static double prof_timestamp=0;
static SpeedProfilerInfo prof_data[PROF_CAP]={0};
static int prof_count=0;
static void prof_checkpoint(ptrdiff_t size, const char *msg)
{
	double t2=time_sec();
	SpeedProfilerInfo *info=prof_data+prof_count++;
	if(prof_count>=PROF_CAP)
	{
		CRASH("Profiler OOB");
		return;
	}
	info->t=t2-prof_timestamp;
	info->size=size;
	info->msg=msg;
	prof_timestamp=t2;
}
static void prof_print(ptrdiff_t usize)
{
	const int scale=5;
	static char buf[2048]={0};
	int colors[128]={0};
	double timesum=0, tmax=0;
	int prev=0;
	double csum=0;
	int k;

	for(k=0;k<prof_count;++k)
	{
		double t=prof_data[k].t;
		timesum+=t;
		if(tmax<t)
			tmax=t;
	}
	srand((unsigned)__rdtsc());
	colorgen(colors, prof_count, 64, 300, 100);
	printf("1 char = %d ms\n", scale);
	//{
	//	const char m[]="1 sec.";
	//	const int offset=(1000/scale-(sizeof(m)-1))/2;
	//	memset(buf, '-', 1000/scale);
	//	int printed=snprintf(buf+offset, sizeof(buf)-1-offset, "%s", m);
	//	buf[offset+printed]='-';
	//	printf("|%s|\n", buf);
	//}

	printf("|");
	for(k=0;k<prof_count;++k)
	{
		int curr, space;
		int len=0;
		SpeedProfilerInfo *info=prof_data+k;
		csum+=info->t;
		curr=(int)(csum*1000/scale);//fixed scale
		space=curr-prev;
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
			//printf("%s%s", buf, info->msg);
			colorprintf(colors[k], colors[k], buf);
			colorprintf(COLORPRINTF_TXT_DEFAULT, colors[k], "%s", info->msg);
			memset(buf, '-', (ptrdiff_t)space-labelend);
			buf[space-labelend]=0;
			//printf("%s", buf);
			colorprintf(colors[k], colors[k], buf);
		}
		else
		{
			memset(buf, '-', space);
			buf[space]=0;
			//printf("%s", buf);
			colorprintf(colors[k], colors[k], buf);
		}
		printf("|");
		prev=curr;
	}
	printf("\n");

	for(k=0;k<prof_count;++k)
	{
		SpeedProfilerInfo *info=prof_data+k;
		printf("%16.7lf ms %8.4lf%% ", info->t*1000, 100.*info->t/timesum);
		if(info->size)
			printf(" %16.6lf MB/s %16.6lf ms/MB %10lld bytes "
				, info->size/(info->t*1024*1024)
				, info->t*1024*1024*1000/info->size
				, (uint64_t)info->size
			);
		//if(info->msg)
		//	printf("%s", info->msg);
		//else
		//	printf("%d", k);
		if(info->msg)
			colorprintf(COLORPRINTF_TXT_DEFAULT, colors[k], "%-20s", info->msg);
		else
			colorprintf(COLORPRINTF_TXT_DEFAULT, colors[k], "%-20s", "");
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

//cRCT
#if 1
#ifndef ENABLE_RCT_EXTENSION
#define OCHLIST\
	OCH(YX00) OCH(Y0X0) OCH(Y00X)\
	OCH(CX40) OCH(C0X4) OCH(C40X)
#endif
#ifdef ENABLE_RCT_EXTENSION
#define OCHLIST\
	OCH(YX00) OCH(Y0X0) OCH(Y00X)\
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

	OCH_R=OCH_YX00,
	OCH_G=OCH_Y0X0,
	OCH_B=OCH_Y00X,
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
//	II_COEFF_U_SUB_V2,
//	II_COEFF_V_SUB_U2,

	II_COUNT,
} RCTInfoIdx;
//YUV = RCT * RGB	watch out for permutation in last row
//luma: averaging	chroma: subtraction
//example: _X00_40X_3X1 == [1 0 0; -1 0 1; -3/4 1 -1/4]
#ifndef ENABLE_RCT_EXTENSION
#define RCTLIST\
	RCT(_X00_0X0_00X,	OCH_R,		OCH_G,		OCH_B,		0, 1, 2,	0,  0, 0)\
	RCT(_X00_0X0_04X,	OCH_R,		OCH_G,		OCH_BG,		0, 1, 2,	0,  0, 4)\
	RCT(_X00_0X0_40X,	OCH_R,		OCH_G,		OCH_BR,		0, 1, 2,	0,  4, 0)\
	RCT(_0X0_00X_X40,	OCH_G,		OCH_B,		OCH_RG,		1, 2, 0,	0,  4, 0)\
	RCT(_0X0_00X_X04,	OCH_G,		OCH_B,		OCH_RB,		1, 2, 0,	0,  0, 4)\
	RCT(_00X_X00_4X0,	OCH_B,		OCH_R,		OCH_GR,		2, 0, 1,	0,  0, 4)\
	RCT(_00X_X00_0X4,	OCH_B,		OCH_R,		OCH_GB,		2, 0, 1,	0,  4, 0)\
	RCT(_0X0_04X_X40,	OCH_G,		OCH_BG,		OCH_RG,		1, 2, 0,	4,  4, 0)\
	RCT(_0X0_04X_X04,	OCH_G,		OCH_BG,		OCH_RB,		1, 2, 0,	4,  0, 4)\
	RCT(_0X0_X40_40X,	OCH_G,		OCH_RG,		OCH_BR,		1, 0, 2,	4,  0, 4)\
	RCT(_00X_X04_0X4,	OCH_B,		OCH_RB,		OCH_GB,		2, 0, 1,	4,  4, 0)\
	RCT(_00X_X04_4X0,	OCH_B,		OCH_RB,		OCH_GR,		2, 0, 1,	4,  0, 4)\
	RCT(_00X_0X4_X40,	OCH_B,		OCH_GB,		OCH_RG,		2, 1, 0,	4,  0, 4)\
	RCT(_X00_4X0_40X,	OCH_R,		OCH_GR,		OCH_BR,		0, 1, 2,	4,  4, 0)\
	RCT(_X00_4X0_04X,	OCH_R,		OCH_GR,		OCH_BG,		0, 1, 2,	4,  0, 4)\
	RCT(_X00_40X_0X4,	OCH_R,		OCH_BR,		OCH_GB,		0, 2, 1,	4,  0, 4)
#endif
#ifdef ENABLE_RCT_EXTENSION
#define RCTLIST\
	RCT(_X00_0X0_00X,	OCH_R,		OCH_G,		OCH_B,		0, 1, 2,	0,  0, 0)\
	RCT(_X00_0X0_04X,	OCH_R,		OCH_G,		OCH_BG,		0, 1, 2,	0,  0, 4)\
	RCT(_X00_0X0_40X,	OCH_R,		OCH_G,		OCH_BR,		0, 1, 2,	0,  4, 0)\
	RCT(_0X0_00X_X40,	OCH_G,		OCH_B,		OCH_RG,		1, 2, 0,	0,  4, 0)\
	RCT(_0X0_00X_X04,	OCH_G,		OCH_B,		OCH_RB,		1, 2, 0,	0,  0, 4)\
	RCT(_00X_X00_4X0,	OCH_B,		OCH_R,		OCH_GR,		2, 0, 1,	0,  0, 4)\
	RCT(_00X_X00_0X4,	OCH_B,		OCH_R,		OCH_GB,		2, 0, 1,	0,  4, 0)\
	RCT(_0X0_04X_X40,	OCH_G,		OCH_BG,		OCH_RG,		1, 2, 0,	4,  4, 0)\
	RCT(_0X0_04X_X04,	OCH_G,		OCH_BG,		OCH_RB,		1, 2, 0,	4,  0, 4)\
	RCT(_0X0_X40_40X,	OCH_G,		OCH_RG,		OCH_BR,		1, 0, 2,	4,  0, 4)\
	RCT(_00X_X04_0X4,	OCH_B,		OCH_RB,		OCH_GB,		2, 0, 1,	4,  4, 0)\
	RCT(_00X_X04_4X0,	OCH_B,		OCH_RB,		OCH_GR,		2, 0, 1,	4,  0, 4)\
	RCT(_00X_0X4_X40,	OCH_B,		OCH_GB,		OCH_RG,		2, 1, 0,	4,  0, 4)\
	RCT(_X00_4X0_40X,	OCH_R,		OCH_GR,		OCH_BR,		0, 1, 2,	4,  4, 0)\
	RCT(_X00_4X0_04X,	OCH_R,		OCH_GR,		OCH_BG,		0, 1, 2,	4,  0, 4)\
	RCT(_X00_40X_0X4,	OCH_R,		OCH_BR,		OCH_GB,		0, 2, 1,	4,  0, 4)\
	RCT(_X00_0X0_13X,	OCH_R,		OCH_G,		OCH_B1,		0, 1, 2,	0,  1, 3)\
	RCT(_X00_4X0_13X,	OCH_R,		OCH_GR,		OCH_B1,		0, 1, 2,	4,  1, 3)\
	RCT(_X00_00X_3X1,	OCH_R,		OCH_B,		OCH_G1,		0, 2, 1,	0,  3, 1)\
	RCT(_X00_40X_3X1,	OCH_R,		OCH_BR,		OCH_G1,		0, 2, 1,	4,  3, 1)\
	RCT(_0X0_00X_X13,	OCH_G,		OCH_B,		OCH_R1,		1, 2, 0,	0,  1, 3)\
	RCT(_0X0_04X_X13,	OCH_G,		OCH_BG,		OCH_R1,		1, 2, 0,	4,  1, 3)\
	RCT(_0X0_X40_13X,	OCH_G,		OCH_RG,		OCH_B1,		1, 0, 2,	4,  3, 1)\
	RCT(_00X_X04_3X1,	OCH_B,		OCH_RB,		OCH_G1,		2, 0, 1,	4,  1, 3)\
	RCT(_00X_04X_X13,	OCH_B,		OCH_GB,		OCH_R1,		2, 1, 0,	4,  3, 1)\
	RCT(_X00_0X0_22X,	OCH_R,		OCH_G,		OCH_B2,		0, 1, 2,	0,  2, 2)\
	RCT(_X00_4X0_22X,	OCH_R,		OCH_GR,		OCH_B2,		0, 1, 2,	4,  2, 2)\
	RCT(_X00_00X_2X2,	OCH_R,		OCH_B,		OCH_G2,		0, 2, 1,	0,  2, 2)\
	RCT(_X00_40X_2X2,	OCH_R,		OCH_BR,		OCH_G2,		0, 2, 1,	4,  2, 2)\
	RCT(_0X0_00X_X22,	OCH_G,		OCH_B,		OCH_R2,		1, 2, 0,	0,  2, 2)\
	RCT(_0X0_04X_X22,	OCH_G,		OCH_BG,		OCH_R2,		1, 2, 0,	4,  2, 2)\
	RCT(_0X0_X40_22X,	OCH_G,		OCH_RG,		OCH_B2,		1, 0, 2,	4,  2, 2)\
	RCT(_00X_X04_2X2,	OCH_B,		OCH_RB,		OCH_G2,		2, 0, 1,	4,  2, 2)\
	RCT(_00X_0X4_X22,	OCH_B,		OCH_GB,		OCH_R2,		2, 1, 0,	4,  2, 2)\
	RCT(_X00_0X0_31X,	OCH_R,		OCH_G,		OCH_B3,		0, 1, 2,	0,  3, 1)\
	RCT(_X00_4X0_31X,	OCH_R,		OCH_GR,		OCH_B3,		0, 1, 2,	4,  3, 1)\
	RCT(_X00_00X_1X3,	OCH_R,		OCH_B,		OCH_G3,		0, 2, 1,	0,  1, 3)\
	RCT(_X00_40X_1X3,	OCH_R,		OCH_BR,		OCH_G3,		0, 2, 1,	4,  1, 3)\
	RCT(_0X0_00X_X31,	OCH_G,		OCH_B,		OCH_R3,		1, 2, 0,	0,  3, 1)\
	RCT(_0X0_04X_X31,	OCH_G,		OCH_BG,		OCH_R3,		1, 2, 0,	4,  3, 1)\
	RCT(_0X0_X40_31X,	OCH_G,		OCH_RG,		OCH_B3,		1, 0, 2,	4,  1, 3)\
	RCT(_00X_X04_1X3,	OCH_B,		OCH_RB,		OCH_G3,		2, 0, 1,	4,  3, 1)\
	RCT(_00X_0X4_X31,	OCH_B,		OCH_GB,		OCH_R3,		2, 1, 0,	4,  1, 3)
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
#endif

//ANS validation
#ifdef ANS_VAL
#define ANS_VAL_HISTSIZE 128
typedef struct _ANSVALHeader
{
	unsigned short esize, count;
	unsigned idx;
	struct _ANSVALHeader *above, *below;
	unsigned char data[];
} ANSVALNode;
extern ANSVALNode *debugstack;
extern int ansvalidx, ansvalmax;
static void ansval_push(const void *data, int esize, int count)
{
	int size=count*esize;
	ANSVALNode *node=(ANSVALNode*)malloc(sizeof(ANSVALNode)+size);
	if(!node)
	{
		printf("Alloc error\n");
		exit(1);
	}
	memset(node, 0, sizeof(ANSVALNode)+size);
	node->esize=esize;
	node->count=count;
	node->idx=ansvalidx++;
	node->above=0;
	node->below=debugstack;
	if(debugstack)
		debugstack->above=node;
	memcpy(node->data, data, size);
	debugstack=node;
	++ansvalmax;
}
static void ansval_printr(const void *data, int esize, int count, const void *xdata)//print elements in reverse because little-endian
{
	const unsigned char *p=(const unsigned char*)data, *p2=(const unsigned char*)xdata;
	int size=count*esize, k;
	for(k=0;k<size;k+=esize)
	{
		int k2=esize-1;
		printf(" ");
		for(;k2>=0;--k2)
		{
			int val=p[k+k2];
			if(p2)
				val^=p2[k+k2];
			if(p2&&!val)
				printf("--");
			else
				printf("%02X", val);
		}
	}
	printf("\n");
}
static void* ansval_ptrguard(const void *start, const void *end, const void *ptr, ptrdiff_t nbytes)
{
	size_t istart=(size_t)start, iend=(size_t)end;
	ptrdiff_t size=iend-istart;
	size_t ip1=(size_t)ptr, ip2=ip1+nbytes;
	int problems[]=
	{
		size<0,
		(size_t)(ip1-istart)>=(size_t)size,
		(size_t)(ip2-istart)>=(size_t)size,
	};
	if(problems[0]||problems[1]||problems[2])
	{
		printf("\nOOB\n");
		printf("  inc     %+16lld bytes\n", (uint64_t)nbytes);
		printf("  start   %016zd  %16d\n", istart, 0);
		if(nbytes<0)
		{
			printf("  after   %016lld  %16lld%s\n", (uint64_t)ip2, (uint64_t)(ip2-istart), problems[2]?"  <-":"");
			printf("  before  %016lld  %16lld%s\n", (uint64_t)ip1, (uint64_t)(ip1-istart), problems[1]?"  <-":"");
		}
		else
		{
			printf("  before  %016lld  %16lld%s\n", (uint64_t)ip1, (uint64_t)(ip1-istart), problems[1]?"  <-":"");
			printf("  after   %016lld  %16lld%s\n", (uint64_t)ip2, (uint64_t)(ip2-istart), problems[2]?"  <-":"");
		}
		printf("  end     %016lld  %16lld%s\n", (uint64_t)iend, (uint64_t)size, problems[0]?"  <-":"");
		CRASH("\n");
		return 0;
	}
	return (void*)(nbytes<0?ip2:ip1);
}
static void ansval_check(const void *data, int esize, int count)
{
	--ansvalidx;
	if(!debugstack)
	{
		printf("Debug stack is empty\n");
		ansval_printr(data, esize, count, 0);
		CRASH("");
	}
	else if(debugstack->esize!=esize||debugstack->count!=count||memcmp(data, debugstack->data, esize*count))
	{
		printf("\n\nValidation Error  [enc ^ | v dec]\n");
		if(debugstack->above)
		{
			ANSVALNode *node=debugstack->above;
			if(node->above)
			{
				ANSVALNode *node2=node->above;
				printf("[%10d] Verified:   ", node2->idx);
				ansval_printr(node2->data, node2->esize, node2->count, 0);
			}
			printf("[%10d] Verified:   ", node->idx);
			ansval_printr(node->data, node->esize, node->count, 0);
			printf("\n");
		}

		printf("[%10d] Original:   ", debugstack->idx);
		ansval_printr(debugstack->data, esize, count, 0);

		printf("[%10d] Corrupt:    ", debugstack->idx);
		ansval_printr(data, esize, count, 0);
		
		if(debugstack->esize==esize&&debugstack->count==count)
		{
			printf("[%10d] XOR:        ", debugstack->idx);
			ansval_printr(debugstack->data, esize, count, data);
		}
		if(debugstack->below)
		{
			ANSVALNode *node=debugstack->below;
			printf("\n");
			printf("[%10d] Below:      ", node->idx);
			ansval_printr(node->data, node->esize, node->count, 0);
			if(node->below)
			{
				node=node->below;
				printf("[%10d] Below:      ", node->idx);
				ansval_printr(node->data, node->esize, node->count, 0);
			}
		}
		printf("\n\n");
		CRASH("");
	}
	if(debugstack->below)
		debugstack=debugstack->below;
}
#else
#define ansval_push(...)
#define ansval_check(...)
#endif

//LIFO Bypass Coder
#if 1
#define BITPACKERMAX 32
typedef struct _BitPackerLIFO//bwd enc / fwd dec
{
	uint64_t state;
	int32_t enc_nwritten, dec_navailable;//bitcounts, only for tracking renorms
	uint8_t *dstbwdptr;
	const uint8_t *srcfwdptr, *streamend;
} BitPackerLIFO;
AWM_INLINE void bitpacker_enc_init(BitPackerLIFO *ec, const uint8_t *bufstart, uint8_t *bufptr0_OOB)
{
	memset(ec, 0, sizeof(*ec));
	ec->state=1ULL<<32;
	ec->enc_nwritten=33;
	ec->streamend=bufstart;
	ec->dstbwdptr=bufptr0_OOB;
}
AWM_INLINE void bitpacker_dec_init(BitPackerLIFO *ec, const uint8_t *bufptr0_start, const uint8_t *bufend)
{
	memset(ec, 0, sizeof(*ec));
	ec->srcfwdptr=bufptr0_start+8;
	ec->streamend=bufend;
	ec->state=*(const uint64_t*)bufptr0_start;
	ec->dec_navailable=FLOOR_LOG2(ec->state)+1;
}
AWM_INLINE void bitpacker_enc_flush(BitPackerLIFO *ec)
{
	ec->dstbwdptr-=8;
#ifdef _DEBUG
	if(ec->dstbwdptr<ec->streamend)
		CRASH("IntPacker Encoder OOB:  dstbwdptr = 0x%016zX < 0x%016zX", ec->dstbwdptr, ec->streamend);
#endif
	*(uint64_t*)ec->dstbwdptr=ec->state;
}
AWM_INLINE void bitpacker_enc(BitPackerLIFO *ec, int inbits, int sym)
{
#ifdef _DEBUG
	if(inbits>BITPACKERMAX)
		CRASH("BitPacker inbits %d", inbits);
#endif
	//renorm then push inbits
	ec->enc_nwritten+=inbits;
	if(ec->enc_nwritten>64)//renorm on overflow
	{
		ec->enc_nwritten-=32;
		ec->dstbwdptr-=4;
#ifdef _DEBUG
		if(ec->dstbwdptr<ec->streamend)
			CRASH("IntPacker OOB:  dstbwdptr = 0x%016zX < 0x%016zX", ec->dstbwdptr, ec->streamend);
#endif
		*(unsigned*)ec->dstbwdptr=(unsigned)ec->state;
		ec->state>>=32;
#ifdef ANS_VAL
		ansval_push(&ec->state, sizeof(ec->state), 1);
#endif
	}
	ec->state=ec->state<<inbits|sym;
#ifdef ANS_VAL
	ansval_push(&inbits, sizeof(inbits), 1);
	ansval_push(&ec->state, sizeof(ec->state), 1);
#endif
}
AWM_INLINE int bitpacker_dec(BitPackerLIFO *ec, int outbits)
{
	int sym;

#ifdef _DEBUG
	if(outbits>BITPACKERMAX)
		CRASH("BitPacker outbits %d", outbits);
#endif
	sym=ec->state&((1ULL<<outbits)-1);

	//pop outbits then renorm
#ifdef ANS_VAL
	ansval_check(&ec->state, sizeof(ec->state), 1);
	ansval_check(&outbits, sizeof(outbits), 1);
#endif
	ec->dec_navailable-=outbits;
	ec->state>>=outbits;
	if(ec->dec_navailable<=32)
	{
#ifdef ANS_VAL
		ansval_check(&ec->state, sizeof(ec->state), 1);
#endif
		ec->dec_navailable+=32;
#ifdef _DEBUG
		if(ec->srcfwdptr>=ec->streamend)
			CRASH("IntPacker OOB:  srcfwdptr = 0x%016zX >= 0x%016zX", ec->srcfwdptr, ec->streamend);
#endif
		ec->state=ec->state<<32|*(const unsigned*)ec->srcfwdptr;
		ec->srcfwdptr+=4;
	}
	return sym;
}
#endif

//SIMD static-o1 rANS	https://github.com/rygorous/ryg_rans	https://github.com/samtools/htscodecs
#ifdef PROBBITS
typedef struct _rANS_SIMD_SymInfo	//16 bytes/level	4KB/ctx = 1<<12 bytes
{
	uint32_t smax, invf, cdf;
	uint16_t negf, sh;
} rANS_SIMD_SymInfo;
static void enc_hist2stats(int *hist, rANS_SIMD_SymInfo *syminfo, uint64_t *bypassmask, int ctxidx, int sse41)
{
#ifdef ESTIMATE_SIZE
	int count0=0, sum0=0;
#endif
	int sum=0, count=0, ks, rare;
	for(ks=0;ks<256;++ks)
	{
		int freq=hist[ks];
		sum+=freq;
		count+=freq!=0;
	}
	rare=sum<12*256/8;
	*bypassmask|=(uint64_t)rare<<ctxidx;
#ifdef ESTIMATE_SIZE
	count0=count; sum0=sum;
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
	{
		double e=sum0;
		if(count==count0)
		{
			double norm=1./0x1000;
			int ks;
			e=0;
			for(ks=0;ks<256;++ks)//estimate
			{
				int freq=(ks<256-1?hist[ks+1]:1<<PROBBITS)-hist[ks];
				if(freq)
				{
					double p=freq*norm;
					e-=p*log(p);
				}
				if(e!=e)
					CRASH("");
			}
			e*=sum/(8.*M_LN2);
		}
		if(ctxidx&&!(ctxidx%NCTX))
			printf("\n");
		printf("%c  ctx %3d  %12.2lf / %9d bytes%10.2lf%%  %3d %s",
			ctxidx<3*NCTX?"YUV"[ctxidx/NCTX]:"yuv"[ctxidx-3*NCTX],
			ctxidx, e, sum0, 100.*e/sum0, count0, count==count0?"levels":"bypass"
		);
		if(count==count0&&count<256)
		{
			int fmax, ks;

			printf(" %3d", count);
			fmax=0;
			for(ks=0;ks<256;++ks)
			{
				int freq=(ks<256-1?hist[ks+1]:1<<PROBBITS)-hist[ks];
				if(fmax<freq)
					fmax=freq;
			}
			for(ks=0;ks<256;++ks)
			{
				int freq, shade;

				freq=(ks<256-1?hist[ks+1]:1<<PROBBITS)-hist[ks];
				if(!(ks&15))
					printf(" ");

				shade=48+freq*(255-48)/fmax;
				colorprintf(shade<<16|shade<<8|shade, freq?0x808080:COLORPRINTF_BK_DEFAULT, "%c", "0123456789ABCDEF"[ks&15]);
				//int shade=freq*255/fmax;
				//colorprintf(freq?0xFFFFFF:0x808080, shade<<16|0<<8|shade, "%c", "0123456789ABCDEF"[ks&15]);

				//printf("%c", freq?"0123456789ABCDEF"[ks&15]:'-');
			}
		}
		printf("\n");
	}
#ifdef DEBUG_HIST
	if(ctxidx==0)
	{
		const int amplitude=512;
		int ks;

		printf("Context %d: (1 star = %d steps)\n", ctxidx, (1<<PROBBITS)/amplitude);
		for(ks=0;ks<256;++ks)
		{
			int freq, nstars, k;

			freq=(ks<256-1?hist[ks+1]:1<<PROBBITS)-hist[ks];
			nstars=freq*amplitude>>PROBBITS;
			printf("%3d %4d ", ks, freq);
			for(k=0;k<nstars;++k)
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
				uint64_t inv;

				info->sh=FLOOR_LOG2(freq);//eg: x/2 = x*0x80000000>>32>>0
				inv=((0x100000000ULL<<info->sh)+freq-1)/freq;
				info->invf=(uint32_t)inv;
				if(inv>0xFFFFFFFF)
				{
					--info->sh;
					info->invf=(uint32_t)(inv>>1);
				}
			}
#ifdef PRINT_SHIFTBOUNDS
			if(minsh>info->sh)
				minsh=info->sh;
			if(maxsh<info->sh)
				maxsh=info->sh;
#endif
			if(sse41)
				info->sh=1<<(PROBBITS-1-info->sh);
		}
	}
}
static void enc_packhist(BitPackerLIFO *ec, const int *hist, uint64_t bypassmask, int ctxidx)//histogram must be normalized to PROBBITS, with spike at 128
{
	uint16_t CDF[257];
	int ks;

	if(bypassmask>>ctxidx&1)
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
static void dec_unpackhist(BitPackerLIFO *ec, uint32_t *CDF2sym, uint64_t bypassmask, int ctxidx)
{
	uint16_t hist[257];
	int ks;

	if(bypassmask>>ctxidx&1)//rare context
	{
		for(ks=0;ks<256;++ks)//bypass
			hist[ks]=(1<<PROBBITS)/256;
	}
	else
	{
		uint16_t CDF[257]={0};
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
					CRASH("CDF unpack error");
#endif
				break;
			}
		}
		if(CDFlevels)
			CRASH("CDF unpack error");
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
#ifdef DEBUG_HIST
	if(ctxidx==DEBUG_HIST)
	{
		const int amplitude=512;
		int ks;

		printf("Context %d: (1 star = %d steps)\n", ctxidx, (1<<PROBBITS)/amplitude);
		for(ks=0;ks<256;++ks)
		{
			int freq, nstars, k;

			freq=(ks<256-1?hist[ks+1]:1<<PROBBITS)-hist[ks], nstars=freq*amplitude>>PROBBITS;
			printf("%3d %4d ", ks, freq);
			for(k=0;k<nstars;++k)
				printf("*");
			printf("\n");
		}
	}
#endif
}
static void save_ppm(const char *fn, const uint8_t *image, int iw, int ih)
{
	FILE *fdst=fopen(fn, "wb");
	if(!fdst)
	{
		CRASH("Cannot open \"%s\" for writing", fn);
		return;
	}
	fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
	fwrite(image, 1, (ptrdiff_t)3*iw*ih, fdst);
	fclose(fdst);
}
static void decorr1d(uint8_t *data, int count, int bytestride, int bestrct, int *rhist)
{
	const uint8_t *combination=rct_combinations[bestrct];
	int yidx=combination[II_PERM_Y];
	int uidx=combination[II_PERM_U];
	int vidx=combination[II_PERM_V];
	int ufromy=-(combination[II_COEFF_U_SUB_Y]!=0);
	int vc0=combination[II_COEFF_V_SUB_Y];
	int vc1=combination[II_COEFF_V_SUB_U];

	uint8_t *ptr=data;
	int prevy=0, prevu=0, prevv=0, offset=0;
	int k, vpred;
	for(k=0;k<count;++k)
	{
		int y=ptr[yidx]-128;
		int u=ptr[uidx]-128;
		int v=ptr[vidx]-128;
		int sym;
		ptr[0]=sym=(uint8_t)(y-prevy+128);
		++rhist[256*0+sym];
		prevy=y;

		offset=y&ufromy;
		prevu+=offset;
		CLAMP2(prevu, -128, 127);
		ptr[1]=sym=(uint8_t)(u-prevu+128);
		++rhist[256*1+sym];
		prevu=u-offset;

		offset=vc0*y+vc1*u;
		vpred=(prevv+offset)>>2;
		CLAMP2(vpred, -128, 127);
		ptr[2]=sym=(uint8_t)(v-vpred+128);
		++rhist[256*2+sym];
		prevv=4*v-offset;
		ptr+=bytestride;
	}
}
static void encode1d_sse41(uint8_t *data, int count, int bytestride, unsigned *pstate, uint8_t **pstreamptr, const uint8_t *streamend, const rANS_SIMD_SymInfo *rsyminfo)
{
	uint8_t *streamptr=*pstreamptr;
	unsigned state=*pstate;
	uint8_t *ptr=data+(count-(ptrdiff_t)1)*bytestride;
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
				CRASH("OOB ptr %016zX <= %016zX", streamptr, streamend);
#endif
			*(uint16_t*)streamptr=(uint16_t)state;
			state>>=RANS_RENORM_BITS;
		}
		//state += ((state*invf>>32)*(1<<(11-sh))>>11)*negf+cdf
		state+=(((uint64_t)state*info->invf>>32)*info->sh>>(PROBBITS-1))*info->negf+info->cdf;
		//state+=((uint64_t)state*info->invf>>32>>info->sh)*info->negf+info->cdf;
#ifdef ANS_VAL
		ansval_push(&state, sizeof(state), 1);
#endif

		info=rsyminfo+ptr[1]+256*1;
		if(state>info->smax)
		{
			streamptr-=2;
#ifdef _DEBUG
			if(streamptr<=streamend)
				CRASH("OOB ptr %016zX <= %016zX", streamptr, streamend);
#endif
			*(uint16_t*)streamptr=(uint16_t)state;
			state>>=RANS_RENORM_BITS;
		}
		state+=(((uint64_t)state*info->invf>>32)*info->sh>>(PROBBITS-1))*info->negf+info->cdf;
#ifdef ANS_VAL
		ansval_push(&state, sizeof(state), 1);
#endif

		info=rsyminfo+ptr[0]+256*0;
		if(state>info->smax)
		{
			streamptr-=2;
#ifdef _DEBUG
			if(streamptr<=streamend)
				CRASH("OOB ptr %016zX <= %016zX", streamptr, streamend);
#endif
			*(uint16_t*)streamptr=(uint16_t)state;
			state>>=RANS_RENORM_BITS;
		}
		state+=(((uint64_t)state*info->invf>>32)*info->sh>>(PROBBITS-1))*info->negf+info->cdf;
#ifdef ANS_VAL
		ansval_push(&state, sizeof(state), 1);
#endif
		ptr-=bytestride;
	}
	*pstreamptr=streamptr;
	*pstate=state;
}
static void encode1d(uint8_t *data, int count, int bytestride, unsigned *pstate, uint8_t **pstreamptr, const uint8_t *streamend, const rANS_SIMD_SymInfo *rsyminfo)
{
	uint8_t *streamptr=*pstreamptr;
	unsigned state=*pstate;
	uint8_t *ptr=data+(count-(ptrdiff_t)1)*bytestride;
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
				CRASH("OOB ptr %016zX <= %016zX", streamptr, streamend);
#endif
			*(uint16_t*)streamptr=(uint16_t)state;
			state>>=RANS_RENORM_BITS;
		}
		//state += ((state*invf>>32)*(1<<(11-sh))>>11)*negf+cdf
		state+=(uint32_t)((uint64_t)state*info->invf>>32>>info->sh)*info->negf+info->cdf;
#ifdef ANS_VAL
		ansval_push(&state, sizeof(state), 1);
#endif

		info=rsyminfo+ptr[1]+256*1;
		if(state>info->smax)
		{
			streamptr-=2;
#ifdef _DEBUG
			if(streamptr<=streamend)
				CRASH("OOB ptr %016zX <= %016zX", streamptr, streamend);
#endif
			*(uint16_t*)streamptr=(uint16_t)state;
			state>>=RANS_RENORM_BITS;
		}
		state+=(uint32_t)((uint64_t)state*info->invf>>32>>info->sh)*info->negf+info->cdf;
#ifdef ANS_VAL
		ansval_push(&state, sizeof(state), 1);
#endif

		info=rsyminfo+ptr[0]+256*0;
		if(state>info->smax)
		{
			streamptr-=2;
#ifdef _DEBUG
			if(streamptr<=streamend)
				CRASH("OOB ptr %016zX <= %016zX", streamptr, streamend);
#endif
			*(uint16_t*)streamptr=(uint16_t)state;
			state>>=RANS_RENORM_BITS;
		}
		state+=(uint32_t)((uint64_t)state*info->invf>>32>>info->sh)*info->negf+info->cdf;
#ifdef ANS_VAL
		ansval_push(&state, sizeof(state), 1);
#endif
		ptr-=bytestride;
	}
	*pstreamptr=streamptr;
	*pstate=state;
}
static void decode1d(uint8_t *data, int count, int bytestride, int bestrct, unsigned *pstate, const uint8_t **pstreamptr, const uint8_t *streamend, unsigned *rCDF2syms)
{
	const uint8_t *combination=rct_combinations[bestrct];
	int yidx=combination[II_PERM_Y];
	int uidx=combination[II_PERM_U];
	int vidx=combination[II_PERM_V];
	int ufromy=-(combination[II_COEFF_U_SUB_Y]!=0);
	int vc0=combination[II_COEFF_V_SUB_Y];
	int vc1=combination[II_COEFF_V_SUB_U];

	const uint8_t *streamptr=*pstreamptr;
	unsigned state=*pstate;
	uint8_t *ptr=data;
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
				CRASH("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
			state=state<<16|*(uint16_t*)streamptr;
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
				CRASH("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
			state=state<<16|*(uint16_t*)streamptr;
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
				CRASH("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
			state=state<<16|*(uint16_t*)streamptr;
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
#endif

#endif//INC_COMMON_H
