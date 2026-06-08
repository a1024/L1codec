

#use clang	mingw32-make x64
#use gcc	mingw32-make x64 USEGCC=1


ANDROID_NDK = E:/Programs/AndroidSDK/ndk/30.0.14904198/toolchains/llvm/prebuilt/windows-x86_64/bin

ifeq ($(USEGCC), 1)
CC = gcc
else
CC = clang
endif

CFLAGS = -Wall -O3 -s

ifeq ($(OS), Windows_NT)
LIBS += -static
else
LIBS += -lm
endif


x64:	l1c-avx2  l1c-avx512 l1c-sse41-x64  l1c-port-x64
#l1c-sse41-linux  l1c-avx2-linux
x86:	l1c-port-x86 l1c-ssse3-x86  l1c-sse41-x86


l1c-avx2:	l1c-avx2.c
	$(CC) $(CFLAGS) -mavx2 -mlzcnt -Dcodec_l1_avx2=main l1c-avx2.c -o l1c-avx2 $(LIBS)

#l1c-avx2-prof:	l1c-avx2.c
#	clang -Wall -O2 -g -fprofile-instr-generate -fcoverage-mapping -mavx2 -mlzcnt -Dcodec_l1_avx2=main l1c-avx2.c -o l1c-avx2-prof $(LIBS)

l1c-port-x64:
	$(CC) $(CFLAGS) -Dcodec_l1_port=main l1c-port.c -o l1c-port-x64

l1c-port-x86:
	$(CC) $(CFLAGS) -m32 -march=pentiumpro -fno-vectorize -fno-slp-vectorize -Dcodec_l1_port=main l1c-port.c -o l1c-port-x86

l1c-sse41-x64:	l1c-sse41.c
	$(CC) $(CFLAGS) -msse4.1 -Dcodec_l1_sse41=main l1c-sse41.c -o l1c-sse41-x64 $(LIBS)

l1c-sse41-x86:	l1c-sse41.c
	gcc $(CFLAGS) -m32 -msse4.1 -Dcodec_l1_sse41=main l1c-sse41.c -o l1c-sse41-x86 $(LIBS)


l1c-sse41-linux:	l1c-sse41.c
	$(CC) $(CFLAGS) --target=x86_64-pc-linux-gnu -msse4.1 -Dcodec_l1_sse41=main l1c-sse41.c -o l1c-sse41 $(LIBS)

l1c-avx2-linux:		l1c-avx2.c
	$(CC) $(CFLAGS) --target=x86_64-pc-linux-gnu -mavx2 -mlzcnt -Dcodec_l1_avx2=main l1c-avx2.c -o l1c-avx2-linux $(LIBS)


#incompatible with the rest	TODO
l1c-ssse3-x86:	l1c-ssse3.c
	gcc $(CFLAGS) -m32 -mssse3 -Dcodec_l1_ssse3=main l1c-ssse3.c -o l1c-ssse3-x86 $(LIBS)


#incompatible with the rest	TODO
l1c-avx512:	l1c-avx512.c
	$(CC) $(CFLAGS) -mavx512f -mavx512bw -mavx512cd -mavx512dq -mavx512vl -mavx512vbmi2 -mbmi -mlzcnt -Dcodec_l1_avx512=main l1c-avx512.c -o l1c-avx512 $(LIBS)




l1c-port-aarch64:
	$(CC) $(CFLAGS) -Dcodec_l1_port=main l1c-port.c -o l1c-port-aarch64


l1c-arm-debug:
	$(ANDROID_NDK)/clang --target=aarch64-linux-android33 -Wall -g -Dcodec_l1_arm=main l1c-arm.c -o l1c-arm
	adb push l1c-arm /data/local/tmp/
	adb shell chmod +x /data/local/tmp/l1c-arm

l1c-arm:
	$(ANDROID_NDK)/clang --target=aarch64-linux-android33 -Wall -O3 -Dcodec_l1_arm=main l1c-arm.c -o l1c-arm

l1c-android-dump:
	$(ANDROID_NDK)/llvm-objdump -M intel -d l1c-arm > l1c-android-asm.txt

#l1c-android-test:
#	adb push E:/Programs/AndroidSDK/ndk/30.0.14904198/toolchains/llvm/prebuilt/windows-x86_64/lib/clang/21/lib/linux/aarch64/lldb-server /data/local/tmp/
#	adb push ../canon_eos_1100d_01.ppm /data/local/tmp/
#	adb push l1c-arm /data/local/tmp/
#	adb shell chmod +x /data/local/tmp/l1c-arm
#
#	adb shell
#	cd /data/local/tmp/
#	./l1c-arm canon_eos_1100d_01.ppm zzz.l1c 6
#	exit
#
#	adb pull /data/local/tmp/zzz.l1c
#
#
#	adb shell /data/local/tmp/lldb-server platform --server --listen "*:1234"
#
#	adb forward tcp:1234 tcp:1234
#	lldb
#	platform select remote-android
#	platform connect connect://localhost:1234
#	target create l1c-arm

#l1c-android-profile:
#	adb shell simpleperf stat /data/local/tmp/l1c-arm /data/local/tmp/canon_eos_1100d_01.ppm /data/local/tmp/zzz.l1c 2
#	adb shell simpleperf record /data/local/tmp/l1c-arm /data/local/tmp/canon_eos_1100d_01.ppm /data/local/tmp/zzz.l1c 2



