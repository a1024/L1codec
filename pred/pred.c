#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#define WIN32_LEAN_AND_MEAN
#include<Windows.h>


	#define LOUD


#define L1SH 20
#if 0
#define NPREDS 10		//up to 11, otherwise slow
#define PREDLIST\
	PRED(N)\
	PRED(W)\
	PRED(NNN)\
	PRED(WWW)\
	PRED(3*(N-NN)+NNN)\
	PRED(3*(W-WW)+WWW)\
	PRED(W+NE-N)\
	PRED(N+W-NW)\
	PRED(N+NE-NNE)\
	PRED(NEEE)
#endif
#if 1
#define NPREDS 8
#define PREDLIST\
	PRED(N)\
	PRED(W)\
	PRED(3*(N-NN)+NNN)\
	PRED(3*(W-WW)+WWW)\
	PRED(W+NE-N)\
	PRED((WWWWW+WW-W+NNN+N+NEEEEE)>>2)\
	PRED(N+W-NW)\
	PRED(N+NE-NNE)
#endif


static void crash(const char *file, int line, const char *format, ...)
{
	printf("%s(%d):\n", file, line);
	if(format)
	{
		va_list args;
		printf("\n");
		va_start(args, format);
		vprintf(format, args);
		va_end(args);
		printf("\n");
	}
	printf("CRASH\n");
	exit(1);
}
#define CRASH(FORMAT, ...) crash(__FILE__, __LINE__, FORMAT, ##__VA_ARGS__)
static double time_sec(void)
{
#ifdef _WIN32
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
#define CLAMP2(X, LO, HI)\
	do\
	{\
		if((X)<(LO))X=LO;\
		if((X)>(HI))X=HI;\
	}while(0)
static void memfill(void *dst, const void *src, size_t dstbytes, size_t srcbytes)
{
	size_t copied;
	char *d=(char*)dst;
	const char *s=(const char*)src;
#ifdef _DEBUG
	if(!dstbytes||!srcbytes)
		return;
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
typedef enum _OCHIndex
{
#define OCH(X) OCH_##X,
	OCHLIST
#undef  OCH
	OCH_COUNT,
	OCH_C4X0=OCH_CX40,
	OCH_C04X=OCH_C0X4,
	OCH_CX04=OCH_C40X,
	OCH_R=OCH_Y400,
	OCH_G=OCH_Y040,
	OCH_B=OCH_Y004,
	OCH_BG=OCH_C04X,
	OCH_BR=OCH_C40X,
	OCH_RG=OCH_CX40,
	OCH_RB=OCH_CX04,
	OCH_GB=OCH_C0X4,
	OCH_GR=OCH_C4X0,
	OCH_R1=OCH_CX13,
	OCH_G1=OCH_C3X1,
	OCH_B1=OCH_C13X,
	OCH_R2=OCH_CX22,
	OCH_G2=OCH_C2X2,
	OCH_B2=OCH_C22X,
	OCH_R3=OCH_CX31,
	OCH_G3=OCH_C1X3,
	OCH_B3=OCH_C31X,
} OCHIndex;
//static const char *och_names[]=
//{
//#define OCH(X) #X,
//	OCHLIST
//#undef  OCH
//};
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
//	II_COEFF_U_SUB_V_2,
//	II_COEFF_V_SUB_U_2,

	II_COUNT,
} RCTInfoIdx;
//YUV = RCT * RGB	watch out for permutation in last row
//luma: averaging	chroma: subtraction
//example: _400_40X_3X1 == [1 0 0; -1 0 1; -3/4 1 -1/4]
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
#ifdef LOUD
static const char *rct_names[RCT_COUNT]=
{
#define RCT(LABEL, ...) #LABEL,
	RCTLIST
#undef  RCT
};
#endif
int analysis(const unsigned char *image, int iw, int ih)
{
	int bestrct=0;
	int prev[OCH_COUNT]={0};
	long long counters[OCH_COUNT]={0};
	ptrdiff_t size=(ptrdiff_t)3*iw*ih;
	const unsigned char *ptr=image, *end=image+size;
	while(ptr<end)
	{
		int
			r=ptr[0]<<2,
			g=ptr[1]<<2,
			b=ptr[2]<<2,
			rg=r-g,
			gb=g-b,
			br=b-r;
		ptr+=3;
#define UPDATE(IDXA, IDXB, IDXC, A0, B0, C0)\
	do\
	{\
		int a0=A0, b0=B0, c0=C0;\
		counters[IDXA]+=abs(a0-prev[IDXA]);\
		counters[IDXB]+=abs(b0-prev[IDXB]);\
		counters[IDXC]+=abs(c0-prev[IDXC]);\
		prev[IDXA]=a0;\
		prev[IDXB]=b0;\
		prev[IDXC]=c0;\
	}while(0)
		UPDATE(OCH_Y400, OCH_Y040, OCH_Y004, r, g, b);
		UPDATE(OCH_CX40, OCH_C0X4, OCH_C40X, rg, gb, br);

		//r-(3*g+b)/4 = r-g-(b-g)/4
		//g-(3*r+b)/4 = g-r-(b-r)/4
		//b-(3*r+g)/4 = b-r-(g-r)/4
		UPDATE(OCH_CX31, OCH_C3X1, OCH_C31X, rg+(gb>>2), rg+(br>>2), br+(rg>>2));

		//r-(g+3*b)/4 = r-b-(g-b)/4
		//g-(r+3*b)/4 = g-b-(r-b)/4
		//b-(r+3*g)/4 = b-g-(r-g)/4
		UPDATE(OCH_CX13, OCH_C1X3, OCH_C13X, br+(gb>>2), gb+(br>>2), gb+(rg>>2));

		//r-(g+b)/2 = (r-g + r-b)/2
		//g-(r+b)/2 = (g-r + g-b)/2
		//b-(r+g)/2 = (b-r + b-g)/2
		UPDATE(OCH_CX22, OCH_C2X2, OCH_C22X, (rg-br)>>1, (gb-rg)>>1, (br-gb)>>1);
#undef  UPDATE
	}
	{
		long long minerr=0;
		for(int kt=0;kt<RCT_COUNT;++kt)
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
#ifdef LOUD
	printf("\n%-14s\n", rct_names[bestrct]);
#endif
	return bestrct;
}
void predict(unsigned char *image, int iw, int ih, int *prct, int fwd)
{
	double t=time_sec();
	unsigned char *imptr=image;
	int weights[3][NPREDS]={0};
	int paddedwidth=iw+16;
	int psize=(int)sizeof(short[4*3])*paddedwidth;//4 padded rows * 3 channels
	short *pixels=(short*)malloc(psize);
	int rct=fwd?analysis(image, iw, ih):*prct;
	double t2=time_sec();
	double tanalysis=t2-t;
	t=t2;
	if(rct<0||rct>=RCT_COUNT)
	{
		CRASH("Invalid RCT");
		return;
	}
	const unsigned char *combination=rct_combinations[rct];
	int
		yidx=combination[II_PERM_Y],
		uidx=combination[II_PERM_U],
		vidx=combination[II_PERM_V];
	int vfromy=-(combination[II_COEFF_U_SUB_Y]!=0);
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	if(fwd)
		*prct=rct;
	memset(pixels, 0, psize);
	FILLMEM((int*)weights, (1<<L1SH)/NPREDS, sizeof(weights), sizeof(int));
	for(int ky=0;ky<ih;++ky)
	{
		short
			*NNNptr=	pixels+(paddedwidth*((ky-3LL+4)%4)+8)*3,
			*NNptr=		pixels+(paddedwidth*((ky-2LL+4)%4)+8)*3,
			*Nptr=		pixels+(paddedwidth*((ky-1LL+4)%4)+8)*3,
			*currptr=	pixels+(paddedwidth*((ky-0LL+4)%4)+8)*3;
		int yuv[3]={0};
		for(int kx=0;kx<iw;++kx)
		{
			int offset=0;
			if(fwd)
			{
				yuv[0]=imptr[yidx]-128;
				yuv[1]=imptr[uidx]-128;
				yuv[2]=imptr[vidx]-128;
			}
			for(int kc=0;kc<3;++kc)
			{
				int
					NNNWWW		=NNNptr	[-3*3],
					NNNW		=NNNptr	[-1*3],
					NNN		=NNNptr	[+0*3],
					NNNE		=NNNptr	[+1*3],
					NNNEE		=NNNptr	[+2*3],
					NNNEEE		=NNNptr	[+3*3],
					NNNEEEE		=NNNptr	[+4*3],
					NNWWWW		=NNptr	[-4*3],
					NNWWW		=NNptr	[-3*3],
					NNWW		=NNptr	[-2*3],
					NNW		=NNptr	[-1*3],
					NN		=NNptr	[+0*3],
					NNE		=NNptr	[+1*3],
					NNEE		=NNptr	[+2*3],
					NNEEE		=NNptr	[+3*3],
					NNEEEE		=NNptr	[+4*3],
					NWWWW		=Nptr	[-4*3],
					NWWW		=Nptr	[-3*3],
					NWW		=Nptr	[-2*3],
					NW		=Nptr	[-1*3],
					N		=Nptr	[+0*3],
					NE		=Nptr	[+1*3],
					NEE		=Nptr	[+2*3],
					NEEE		=Nptr	[+3*3],
					NEEEE		=Nptr	[+4*3],
					NEEEEE		=Nptr	[+5*3],
					NEEEEEE		=Nptr	[+6*3],
					NEEEEEEE	=Nptr	[+7*3],
					NEEEEEEEE	=Nptr	[+8*3],
					WWWWWWWWW	=currptr[-9*3],
					WWWWWWWW	=currptr[-8*3],
					WWWWWWW		=currptr[-7*3],
					WWWWWW		=currptr[-6*3],
					WWWWW		=currptr[-5*3],
					WWWW		=currptr[-4*3],
					WWW		=currptr[-3*3],
					WW		=currptr[-2*3],
					W		=currptr[-1*3];
#if 1
				(void)NNNWWW	;
				(void)NNNW	;
				(void)NNN	;
				(void)NNNE	;
				(void)NNNEE	;
				(void)NNNEEE	;
				(void)NNNEEEE	;
				(void)NNWWWW	;
				(void)NNWWW	;
				(void)NNWW	;
				(void)NNW	;
				(void)NN	;
				(void)NNE	;
				(void)NNEE	;
				(void)NNEEE	;
				(void)NNEEEE	;
				(void)NWWWW	;
				(void)NWWW	;
				(void)NWW	;
				(void)NW	;
				(void)N		;
				(void)NE	;
				(void)NEE	;
				(void)NEEE	;
				(void)NEEEE	;
				(void)NEEEEE	;
				(void)NEEEEEE	;
				(void)NEEEEEEE	;
				(void)NEEEEEEEE	;
				(void)WWWWWWWWW	;
				(void)WWWWWWWW	;
				(void)WWWWWWW	;
				(void)WWWWWW	;
				(void)WWWWW	;
				(void)WWWW	;
				(void)WWW	;
				(void)WW	;
				(void)W		;
#endif
				int preds[]=
				{
#define PRED(EXPR) EXPR,
					PREDLIST
#undef  PRED
				};
				int *currw=weights[kc];
				int p0=1<<L1SH>>1;
				for(int k=0;k<NPREDS;++k)
					p0+=currw[k]*preds[k];
				p0>>=L1SH;
				int predc=p0;
				int vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;
				CLAMP2(predc, vmin, vmax);
				if(kc)
				{
					predc+=offset;
					CLAMP2(predc, -128, 127);
				}

				int curr;
				if(fwd)
				{
					curr=yuv[kc];
					imptr[kc]=(unsigned char)(curr-predc+128);
				}
				else
					yuv[kc]=curr=(char)(imptr[kc]-128+predc);

				curr-=offset;
				currptr[0]=curr;

				//update
				int e=(curr>p0)-(curr<p0);//L1
				for(int k=0;k<NPREDS;++k)
					currw[k]+=e*preds[k];

				offset=kc?(combination[II_COEFF_V_SUB_Y]*yuv[0]+combination[II_COEFF_V_SUB_U]*yuv[1])>>2:yuv[0]&vfromy;

				++NNNptr;
				++NNptr;
				++Nptr;
				++currptr;
			}
			if(!fwd)
			{
				imptr[yidx]=yuv[0]+128;
				imptr[uidx]=yuv[1]+128;
				imptr[vidx]=yuv[2]+128;
			}
			imptr+=3;
		}
	}
	free(pixels);
	double tpred=time_sec()-t;
	ptrdiff_t usize=(ptrdiff_t)3*iw*ih;
	if(fwd)
		printf("Analysis        %12.6lf sec %12.6lf MB/s\n", tanalysis, usize/(tanalysis*1024*1024));
	printf("Decorrelation   %12.6lf sec %12.6lf MB/s\n", tpred, usize/(tpred*1024*1024));
}
static unsigned char* load_ppm(const char *srcfn, int *ret_iw, int *ret_ih, int *ret_rct)
{
	int rct=-1, iw=0, ih=0, c, nread;
	FILE *fsrc=fopen(srcfn, "rb");
	if(!fsrc)
	{
		CRASH("Cannot open \"%s\"", srcfn);
		return 0;
	}
	fread(&c, 1, 2, fsrc);
	if(c!=('P'|'6'<<8))
	{
		CRASH("Unsupported file \"%s\"", srcfn);
		return 0;
	}
	c=fgetc(fsrc);
	if(c!='\n')
	{
		CRASH("Invalid PPM file");
		return 0;
	}
	c=fgetc(fsrc);
	while(c=='#')
	{
		char comment[1024]={0};
		int idx=0;
		comment[idx]=c=fgetc(fsrc);//skip '#'
		while(c!='\n')
		{
			++idx;
			comment[idx]=c=fgetc(fsrc);
		}
		++idx;
		if(rct==-1)
			rct=atoi(comment);
		c=fgetc(fsrc);//skip newline
	}
	iw=0;
	while((unsigned)(c-'0')<10)iw=10*iw+c-'0', c=fgetc(fsrc);
	while((unsigned)c<=' ')c=fgetc(fsrc);//skip whitespace
	ih=0;
	while((unsigned)(c-'0')<10)ih=10*ih+c-'0', c=fgetc(fsrc);
	while(c=='#')//strip comments
	{
		c=fgetc(fsrc);
		while(c!='\n')c=fgetc(fsrc);
		c=fgetc(fsrc);
	}
	nread=(unsigned)fread(&c, 1, 4, fsrc);
	if(nread!=4||c!=('2'|'5'<<8|'5'<<16|'\n'<<24))
	{
		CRASH("Unsupported file");
		return 0;
	}
	if(iw<1||ih<1)
	{
		CRASH("Unsupported file");
		return 0;
	}
	ptrdiff_t size=(ptrdiff_t)3*iw*ih;
	unsigned char *image=(unsigned char*)malloc(size+32);
	if(!image)
	{
		CRASH("Alloc error");
		return 0;
	}
	fread(image, 1, size, fsrc);//read image
	fclose(fsrc);
	if(ret_iw)*ret_iw=iw;
	if(ret_ih)*ret_ih=ih;
	if(ret_rct)*ret_rct=rct;
	return image;
}
static void print_usage(const char *argv0)
{
	printf(
		"Usage:  \"%s\"  p|r  input.ppm  output.ppm\n"
		"  p:   Predict\n"
		"  r:   Reconstruct\n"
		, argv0
	);
}
int main(int argc, char **argv)
{
	const char *command, *srcfn, *dstfn;
	unsigned char *image;
	int fwd, iw, ih, rct;

	if(argc!=4)
	{
		print_usage(argv[0]);
		return 1;
	}
	command=argv[1];
	srcfn=argv[2];
	dstfn=argv[3];

	fwd=(command[0]&0xDF)=='P'&&!command[1];
	if(!fwd&&!((command[0]&0xDF)=='R'&&!command[1]))
	{
		print_usage(argv[0]);
		CRASH("Invalid command");
		return 1;
	}

	image=load_ppm(srcfn, &iw, &ih, &rct);
	if(!image)
	{
		CRASH("Could not load \"%s\"", srcfn);
		return 1;
	}
	predict(image, iw, ih, &rct, fwd);

	{
		FILE *fdst=fopen(dstfn, "wb");
		if(!fdst)
		{
			CRASH("Cannot open \"%s\" for writing\n", dstfn);
			return 1;
		}
		if(fwd)
			fprintf(fdst, "P6\n# %d\n%d %d\n255\n", rct, iw, ih);
		else
			fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
		fwrite(image, 1, (ptrdiff_t)3*iw*ih, fdst);
		fclose(fdst);
	}
	return 0;
}