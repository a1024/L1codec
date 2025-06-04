#if defined _MSC_VER && !defined _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

//	#define PROFILER

#ifdef PROFILER
#include"util.h"
#endif
#include<stdio.h>
#include<stdlib.h>
int l1_codec(int argc, char **argv);//C32: like C29 but 16 coders

int main(int argc, char **argv)
{
	int retcode=0;
#ifdef PROFILER
	void *prof_ctx=prof_start();
#endif
#ifdef __GNUC__
	retcode=l1_codec(argc, argv);
#else
	const char *dstfn=//OVERWRITTEN
		"C:/Projects/datasets/zzz_deletethis.ppm"

	//	"D:/ML/zzz_deletethis.ppm"
	;
	const char *tmpfn=//OVERWRITE
		"C:/Projects/datasets/zzz_deletethis.lsim"

	//	"D:/ML/zzz_deletethis.lsim"
	;
	const char *srcfn=
	//	"C:/Projects/datasets/0868-ecrop.ppm"
	//	"C:/Projects/datasets/20240806 6 why me.PPM"
	//	"C:/Projects/datasets/big_building.PPM"
	//	"C:/Projects/datasets/dataset-CLIC303-ppm/2048x1320_abigail-keenan-27293.ppm"
	//	"C:/Projects/datasets/dataset-CLIC30-ppm/03.ppm"
		"C:/Projects/datasets/dataset-DIV2K-ppm/0801.ppm"
	//	"C:/Projects/datasets/0801-cg.ppm"
	//	"C:/Projects/datasets/dataset-GDCC2020-ppm/astro-01.ppm"
	//	"C:/Projects/datasets/dataset-kodak-ppm/kodim23.ppm"
	//	"C:/Projects/datasets/dataset-LPCB-ppm/canon_eos_1100d_01.ppm"
	//	"C:/Projects/datasets/dataset-LPCB-ppm/PIA13882.ppm"
	//	"C:/Projects/datasets/dataset-LPCB-ppm/STA13843.ppm"	//large
	//	"C:/Projects/datasets/kodim13.ppm"
	//	"C:/Projects/datasets/kodim13-small16.PPM"
	//	"C:/Projects/datasets/kodim13-small4.PPM"
	//	"C:/Projects/datasets/temp.c18"
		

	//	"C:/dataset-CLIC303-ppm/2048x1320_abigail-keenan-27293.ppm"
	//	"C:/dataset-CLIC303-ppm/2048x1320_cosmic-timetraveler-29758.ppm"
	//	"C:/dataset-CLIC303-ppm/2048x1320_rosan-harmens-18703.ppm"
	//	"C:/dataset-CLIC303-ppm/2048x1320_zugr-108.ppm"
	//	"C:/dataset-DIV2K-ppm/0801.ppm"
	//	"C:/dataset-DSLR2x4-ppm/DSC_0133.ppm"
	//	"C:/dataset-GDCC2020-ppm/astro-01.ppm"
	//	"C:/dataset-HUGE-ppm/kodak.PPM"
	//	"C:/dataset-HUGE-ppm/space_huge.ppm"
	//	"C:/dataset-LPCB-ppm/canon_eos_1100d_01.ppm"
	//	"C:/dataset-LPCB-ppm/PIA13803.ppm"
	//	"C:/dataset-LPCB-ppm/PIA13833.ppm"
	//	"C:/dataset-LPCB-ppm/PIA13912.ppm"
	//	"C:/dataset-LPCB-ppm/PIA13915.ppm"	//false color terrain
	//	"C:/dataset-LPCB-ppm/STA13843.ppm"	//space clouds
	//	"C:/dataset-LPCB-ppm/STA13844.ppm"	//space clouds
	//	"C:/dataset-LPCB-ppm/STA13845.ppm"	//space clouds
	//	"C:/dataset-synthetic-ppm/20240409 1 LPCB.ppm"
	//	"D:/ML/big_building.PPM"
	//	"D:/ML/dataset-CLIC303-ppm/2048x1320_lucas-lof-388.ppm"
	//	"D:/ML/dataset-kodak-ppm/kodim13.c01"
	//	"D:/ML/dataset-kodak-ppm/kodim13.ppm"


	//	"C:/Projects/datasets/20240414-noise.LSIM"
	//	"C:/Projects/datasets/20240414-noise.PPM"
	//	"C:/Projects/datasets/20240513 screenshot.PPM"
	//	"C:/Projects/datasets/big_building.PPM"
	//	"C:/Projects/datasets/dataset-DIV2K-ppm/0801.ppm"
	//	"C:/Projects/datasets/dataset-GDCC2020-ppm/astro-01.ppm"
	//	"C:/Projects/datasets/dataset-LPCB-ppm/PIA12811.ppm"
	//	"C:/Projects/datasets/dataset-LPCB-ppm/STA13843.ppm"
	//	"C:/Projects/datasets/dataset-LPCB-ppm/canon_eos_1100d_01.ppm"
	//	"C:/Projects/datasets/kodim13-small16.ppm"
	//	"C:/Projects/datasets/kodim13.ppm"
	//	"C:/Projects/datasets/kodim24.ppm"
	//	"C:/Projects/datasets/space_huge.ppm"
	//	"C:/Projects/datasets/space_huge.ppm"
	//	"C:/dataset-CLIC303-ppm/2048x1320_alberto-restifo-4549.ppm"
	//	"C:/dataset-CLIC303-ppm/2048x1320_zugr-108.ppm"
	//	"C:/dataset-DIV2K-ppm"
	//	"C:/dataset-DIV2K-ppm/0801.ppm"
	//	"C:/dataset-DIV2K-ppm/0805.ppm"
	//	"C:/dataset-DIV2K-ppm/0807.ppm"
	//	"C:/dataset-DIV2K-ppm/0823.ppm"
	//	"C:/dataset-DIV2K-ppm/0843.ppm"
	//	"C:/dataset-DIV2K-ppm/0859.ppm"
	//	"C:/dataset-DIV2K-ppm/0864.ppm"
	//	"C:/dataset-DIV2K-ppm/0880.ppm"
	//	"C:/dataset-GDCC2020-ppm/astro-01.ppm"
	//	"C:/dataset-GDCC2020-ppm/astro-01.ppm"
	//	"C:/dataset-GDCC2020-ppm/astro-02.ppm"
	//	"C:/dataset-GDCC2020-ppm/astro-06.ppm"
	//	"C:/dataset-GDCC2020-ppm/astro-14.ppm"
	//	"C:/dataset-GDCC2020-ppm/astro-20.ppm"
	//	"C:/dataset-GDCC2020-ppm/astro-30.ppm"
	//	"C:/dataset-GDCC2020-ppm/astro-43.ppm"
	//	"C:/dataset-GDCC2020-ppm/photo-03.ppm"
	//	"C:/dataset-GDCC2020-ppm/photo-05.ppm"
	//	"C:/dataset-GDCC2020-ppm/photo-49.ppm"
	//	"C:/dataset-GDCC2020-ppm/photo-52.ppm"
	//	"C:/dataset-GDCC2020-ppm/photo-67.ppm"
	//	"C:/dataset-HUGE-ppm/blackmarble.ppm"
	//	"C:/dataset-HUGE-ppm/chaos1.ppm"
	//	"C:/dataset-HUGE-ppm/diagram.ppm"
	//	"C:/dataset-HUGE-ppm/gaia.ppm"
	//	"C:/dataset-HUGE-ppm/jwst.ppm"
	//	"C:/dataset-HUGE-ppm/jwst.ppm"
	//	"C:/dataset-HUGE-ppm/jwst.ppm"
	//	"C:/dataset-HUGE2-ppm/andromeda.ppm"
	//	"C:/dataset-LPCB-ppm/PIA13757.ppm"
	//	"C:/dataset-LPCB-ppm/canon_eos_1100d_02.ppm"
	//	"C:/dataset-RAW-ppm/a0014-WP_CRW_6320.ppm"
	//	"C:/dataset-a70-ppm/20240816_113656_966.ppm"
	//	"C:/dataset-meme-ppm/emoji_u1f628.ppm"
	//	"C:/dataset-memes-ppm/usa.ppm"
	//	"C:/dataset-synth-ppm/20240421 1 the front.ppm"
	//	"C:/dataset-synth-ppm/20240516 4 DSC_0054.ppm"
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
	//	"D:/ML/big_building.PPM"
	//	"D:/ML/checkboard.PPM"
	//	"D:/ML/dataset-CID22-ppm/pexels-photo-1933873.PPM"
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
		"1",
	//	"11",//near
	};
	const char *decargs[]=
	{
		argv[0],
		tmpfn,
		dstfn,
	};
	if(l1_codec(_countof(encargs), (char**)encargs))
		return 1;
	if(l1_codec(_countof(decargs), (char**)decargs))
		return 1;
#endif
#ifdef PROFILER
	prof_end(prof_ctx);
#endif
	return retcode;
}
