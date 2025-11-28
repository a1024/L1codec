#if defined _MSC_VER && !defined _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include"common.h"

//	#define PROFILER

//#ifdef PROFILER
//#include"util.h"
//#endif
#ifdef _MSC_VER
#include<stdlib.h>
#endif
int codec_l1_ssse3(int argc, char **argv);
int codec_l1_sse41(int argc, char **argv);
int codec_l1_avx2(int argc, char **argv);
int codec_l1_avx512(int argc, char **argv);


//	#define CODEC codec_l1_ssse3
//	#define CODEC codec_l1_sse41
	#define CODEC codec_l1_avx2
//	#define CODEC codec_l1_avx512


#ifdef ENABLE_GUIDE
unsigned char *g_image=0;
double g_sqe[3]={0};
#endif
#ifdef ANS_VAL
ANSVALNode *debugstack=0;
int ansvalidx=0, ansvalmax=0;
#endif
int main(int argc, char **argv)
{
	const char *dstfn=//OVERWRITTEN
		"C:/dataset-a-temp/zzz.ppm"
	//	"D:/ML/zzz_deletethis.ppm"
	//	"C:/Projects/datasets/zzz_deletethis.ppm"
	//	"F:/Projects/zzz.ppm"
	;
	const char *tmpfn=//OVERWRITTEN
		"C:/dataset-a-temp/zzz.l1c"
	//	"D:/ML/zzz_deletethis.l1c"
	//	"C:/Projects/datasets/zzz_deletethis.l1c"
	//	"F:/Projects/zzz.l1c"
	;
	const char *srcfn=
	//	"F:/Projects/dataset-GDCC2020-ppm/photo-01.ppm"

	//	"C:/dataset-a70-ppm/20240816_113656_966.ppm"
	//	"C:/dataset-CLIC303-ppm/2048x1320_abigail-keenan-27293.ppm"
	//	"C:/dataset-CLIC303-ppm/2048x1320_alberto-restifo-4549.ppm"
	//	"C:/dataset-CLIC303-ppm/2048x1320_cosmic-timetraveler-29758.ppm"
	//	"C:/dataset-CLIC303-ppm/2048x1320_rosan-harmens-18703.ppm"
	//	"C:/dataset-CLIC303-ppm/2048x1320_zugr-108.ppm"
	//	"C:/dataset-CLIC303-ppm/2048x1320_zugr-108.ppm"
	//	"C:/dataset-DIV2K-ppm"
	//	"C:/dataset-DIV2K-ppm/0801.ppm"
	//	"C:/dataset-DIV2K-ppm/0801.ppm"
	//	"C:/dataset-DIV2K-ppm/0801.ppm"
	//	"C:/dataset-DIV2K-ppm/0805.ppm"
	//	"C:/dataset-DIV2K-ppm/0807.ppm"
	//	"C:/dataset-DIV2K-ppm/0823.ppm"
	//	"C:/dataset-DIV2K-ppm/0843.ppm"
	//	"C:/dataset-DIV2K-ppm/0859.ppm"
	//	"C:/dataset-DIV2K-ppm/0864.ppm"
	//	"C:/dataset-DIV2K-ppm/0880.ppm"
	//	"C:/dataset-DSLR2x4-ppm/DSC_0133.ppm"
	//	"C:/dataset-GDCC2020-ppm/astro-01.ppm"
	//	"C:/dataset-GDCC2020-ppm/astro-01.ppm"
	//	"C:/dataset-GDCC2020-ppm/astro-01.ppm"
	//	"C:/dataset-GDCC2020-ppm/astro-02.ppm"
	//	"C:/dataset-GDCC2020-ppm/astro-06.ppm"
	//	"C:/dataset-GDCC2020-ppm/astro-14.ppm"
	//	"C:/dataset-GDCC2020-ppm/astro-20.ppm"
	//	"C:/dataset-GDCC2020-ppm/astro-30.ppm"
	//	"C:/dataset-GDCC2020-ppm/astro-43.ppm"
		"C:/dataset-GDCC2020-ppm/photo-01.ppm"
	//	"C:/dataset-GDCC2020-ppm/photo-03.ppm"
	//	"C:/dataset-GDCC2020-ppm/photo-05.ppm"
	//	"C:/dataset-GDCC2020-ppm/photo-49.ppm"
	//	"C:/dataset-GDCC2020-ppm/photo-52.ppm"
	//	"C:/dataset-GDCC2020-ppm/photo-67.ppm"
	//	"C:/dataset-HUGE2-ppm/andromeda.ppm"
	//	"C:/dataset-HUGE-ppm/blackmarble.ppm"
	//	"C:/dataset-HUGE-ppm/chaos1.ppm"
	//	"C:/dataset-HUGE-ppm/diagram.ppm"
	//	"C:/dataset-HUGE-ppm/gaia.ppm"
	//	"C:/dataset-HUGE-ppm/jwst.ppm"
	//	"C:/dataset-HUGE-ppm/jwst.ppm"
	//	"C:/dataset-HUGE-ppm/jwst.ppm"
	//	"C:/dataset-HUGE-ppm/kodak.PPM"
	//	"C:/dataset-HUGE-ppm/space_huge.ppm"
	//	"C:/dataset-LPCB-ppm/canon_eos_1100d_01.ppm"
	//	"C:/dataset-LPCB-ppm/canon_eos_1100d_02.ppm"
	//	"C:/dataset-LPCB-ppm/PIA13757.ppm"
	//	"C:/dataset-LPCB-ppm/PIA13803.ppm"
	//	"C:/dataset-LPCB-ppm/PIA13833.ppm"
	//	"C:/dataset-LPCB-ppm/PIA13912.ppm"
	//	"C:/dataset-LPCB-ppm/PIA13915.ppm"	//false color terrain
	//	"C:/dataset-LPCB-ppm/STA13843.ppm"	//space clouds
	//	"C:/dataset-LPCB-ppm/STA13844.ppm"	//space clouds
	//	"C:/dataset-LPCB-ppm/STA13845.ppm"	//space clouds
	//	"C:/dataset-meme-ppm/emoji_u1f628.ppm"
	//	"C:/dataset-memes-ppm/usa.ppm"
	//	"C:/dataset-RAW-ppm/a0014-WP_CRW_6320.ppm"
	//	"C:/dataset-sony-ppm/DSC00315.ppm"
	//	"C:/dataset-synth2-ppm/20240405 1 CPU-load.ppm"
	//	"C:/dataset-synth2-ppm/20240405 1 CPU-load.ppm"
	//	"C:/dataset-synth2-ppm/20240405 1 CPU-load.ppm"
	//	"C:/dataset-synth2-ppm/20240407 blank.ppm"
	//	"C:/dataset-synth2-ppm/20240409 1 LPCB.ppm"
	//	"C:/dataset-synth2-ppm/20240409 1 LPCB.ppm"
	//	"C:/dataset-synth2-ppm/20240412 2 gralic enc.ppm"
	//	"C:/dataset-synth2-ppm/20240419 1 speed for efficiency.ppm"
	//	"C:/dataset-synth2-ppm/20240419 1 speed for efficiency.ppm"
	//	"C:/dataset-synth2-ppm/20240419 3.ppm"
	//	"C:/dataset-synth2-ppm/20240422 1.PPM"
	//	"C:/dataset-synth2-ppm/20240524 numbers.ppm"
	//	"C:/dataset-synth2-ppm/20241006 linux is cursed.ppm"
	//	"C:/dataset-synth2-ppm/art.ppm"
	//	"C:/dataset-synthetic-ppm/20240409 1 LPCB.ppm"
	//	"C:/dataset-synth-ppm/20240421 1 the front.ppm"
	//	"C:/dataset-synth-ppm/20240516 4 DSC_0054.ppm"
	//	"C:/Projects/datasets/0801-cg.ppm"
	//	"C:/Projects/datasets/0868-ecrop.ppm"
	//	"C:/Projects/datasets/20240414-noise.LSIM"
	//	"C:/Projects/datasets/20240414-noise.PPM"
	//	"C:/Projects/datasets/20240513 screenshot.PPM"
	//	"C:/Projects/datasets/20240806 6 why me.PPM"
	//	"C:/Projects/datasets/big_building.PPM"
	//	"C:/Projects/datasets/big_building.PPM"
	//	"C:/Projects/datasets/dataset-CLIC303-ppm/2048x1320_abigail-keenan-27293.ppm"
	//	"C:/Projects/datasets/dataset-CLIC30-ppm/03.ppm"
	//	"C:/Projects/datasets/dataset-DIV2K-ppm/0801.ppm"
	//	"C:/Projects/datasets/dataset-DIV2K-ppm/0801.ppm"
	//	"C:/Projects/datasets/dataset-GDCC2020-ppm/astro-01.ppm"
	//	"C:/Projects/datasets/dataset-GDCC2020-ppm/astro-01.ppm"
	//	"C:/Projects/datasets/dataset-GDCC2020-ppm/photo-01.ppm"
	//	"C:/Projects/datasets/dataset-GDCC2020-ppm/photo-06.ppm"
	//	"C:/Projects/datasets/dataset-kodak-ppm/kodim23.ppm"
	//	"C:/Projects/datasets/dataset-LPCB-ppm/canon_eos_1100d_01.ppm"
	//	"C:/Projects/datasets/dataset-LPCB-ppm/canon_eos_1100d_01.ppm"
	//	"C:/Projects/datasets/dataset-LPCB-ppm/PIA12811.ppm"
	//	"C:/Projects/datasets/dataset-LPCB-ppm/PIA13882.ppm"
	//	"C:/Projects/datasets/dataset-LPCB-ppm/STA13843.ppm"
	//	"C:/Projects/datasets/dataset-LPCB-ppm/STA13843.ppm"	//large
	//	"C:/Projects/datasets/kodim13.ppm"
	//	"C:/Projects/datasets/kodim13.ppm"
	//	"C:/Projects/datasets/kodim13-small16.ppm"
	//	"C:/Projects/datasets/kodim13-small16.PPM"
	//	"C:/Projects/datasets/kodim13-small4.PPM"
	//	"C:/Projects/datasets/kodim24.ppm"
	//	"C:/Projects/datasets/space_huge.ppm"
	//	"C:/Projects/datasets/space_huge.ppm"
	//	"C:/Projects/datasets/temp.c18"
	//	"D:/ML/big_building.PPM"
	//	"D:/ML/big_building.PPM"
	//	"D:/ML/checkboard.PPM"
	//	"D:/ML/dataset-CID22-ppm/pexels-photo-1933873.PPM"
	//	"D:/ML/dataset-CLIC303-ppm/2048x1320_lucas-lof-388.ppm"
	//	"D:/ML/dataset-kodak-ppm/kodim13.c01"
	//	"D:/ML/dataset-kodak-ppm/kodim13.ppm"
	//	"D:/ML/kodim13.ppm"
	//	"D:/ML/kodim24.ppm"
	//	"D:/ML/nice_clock_face.ppm"
	//	"D:/ML/zzz_halfbright.PPM"
	//	"D:/Programs/c29/song.ppm"
	;
	const char *encargs[]=
	{
		argv[0],
		srcfn,
		tmpfn,
		"2",
	//	"17",//near
	};
	const char *decargs[]=
	{
		argv[0],
		tmpfn,
		dstfn,
	};
	if(CODEC(_countof(encargs), (char**)encargs))
		return 1;
	if(CODEC(_countof(decargs), (char**)decargs))
		return 1;
#if defined _MSC_VER && _MSC_VER<1900
	{
		int k;
		printf("Enter 0 to continue: ");
		while(!scanf(" %d", &k));
	}
#endif
	return 0;
}
