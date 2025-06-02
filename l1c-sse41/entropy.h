#pragma once
#ifndef INC_AC_H
#define INC_AC_H
#include"util.h"
#include<stdio.h>
#include<string.h>
#include<immintrin.h>
#ifdef __cplusplus
extern "C"
{
#endif


#ifdef ANS_VAL
#define ANS_VAL_HISTSIZE 128
typedef struct _ANSVALHeader
{
	unsigned short esize, count;
	unsigned idx;
} ANSVALHeader;
static ArrayHandle debugstack=0;
static void ansval_push(const void *data, int esize, int count)
{
	static int idx=0;
	ANSVALHeader header={esize, count, idx};
	++idx;
	if(!debugstack)
		ARRAY_ALLOC(char, debugstack, 0, 0, 1024, 0);
	ARRAY_APPEND(debugstack, &header, sizeof(header), 1, 0);//lo header
	ARRAY_APPEND(debugstack, data, (ptrdiff_t)count*esize, 1, 0);
	ARRAY_APPEND(debugstack, &header, sizeof(header), 1, 0);//hi header
}
static void ansval_printr(const void *data, int esize, int count, const void *xdata)//print elements in reverse because little-endian
{
	const unsigned char *p=(const unsigned char*)data;
	const unsigned char *p2=(const unsigned char*)xdata;
	int size=count*esize;
	for(int k=0;k<size;k+=esize)
	{
		printf(" ");
		for(int k2=esize-1;k2>=0;--k2)
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
	//for(int k=size-1;k>=0;--k)
	//{
	//	if((k&3)==3)
	//		printf(" ");
	//	printf("%02X", ((unsigned char*)data)[k]);
	//}
	//printf("\n");
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
		printf("  inc     %+16td bytes\n", nbytes);
		printf("  start   %016zd  %16d\n", istart, 0);
		if(nbytes<0)
		{
			printf("  after   %016zd  %16td%s\n", ip2, ip2-istart, problems[2]?"  <-":"");
			printf("  before  %016zd  %16td%s\n", ip1, ip1-istart, problems[1]?"  <-":"");
		}
		else
		{
			printf("  before  %016zd  %16td%s\n", ip1, ip1-istart, problems[1]?"  <-":"");
			printf("  after   %016zd  %16td%s\n", ip2, ip2-istart, problems[2]?"  <-":"");
		}
		printf("  end     %016zd  %16td%s\n", iend, size, problems[0]?"  <-":"");
		LOG_ERROR("\n");
		return 0;
	}
	return (void*)(nbytes<0?ip2:ip1);
}
static void ansval_check(const void *data, int esize, int count)
{
	static const unsigned char *endptr=0;
	static int totalcount=0, popcount=0;
	int firstpop=!endptr;
	if(firstpop)
		endptr=debugstack->data+debugstack->count;
	const ANSVALHeader *hiheader=(const ANSVALHeader*)(debugstack->data+debugstack->count)-1;
	debugstack->count-=sizeof(ANSVALHeader);
	const unsigned char *data0=debugstack->data+debugstack->count-hiheader->count*hiheader->esize;
	debugstack->count-=hiheader->count*hiheader->esize;
	const ANSVALHeader *loheader=(const ANSVALHeader*)(debugstack->data+debugstack->count)-1;
	debugstack->count-=sizeof(ANSVALHeader);
	if(firstpop)
		totalcount=hiheader->idx+1;
	if(memcmp(hiheader, loheader, sizeof(*hiheader)))
	{
		printf("\n");
		printf("Validation Header Mismatch  idx,esize,count: loheader %d,%d,%d != hiheader %d,%d,%d\n",
			loheader->idx, loheader->esize, loheader->count,
			hiheader->idx, hiheader->esize, hiheader->count
		);
		LOG_ERROR("\n");
	}
	if(esize!=loheader->esize||count!=loheader->count||memcmp(data, data0, count*esize))
	{
		printf("\n");
		printf("Validation Error    pop #%d / total %d,  remaining %d,  using %8.2lf/%8.2lf MB\n",
			popcount,
			totalcount,
			totalcount-popcount-1,
			(double)debugstack->count/(1024*1024),
			(double)debugstack->cap/(1024*1024)
		);
		printf("\n");

		const unsigned char *verptr=debugstack->data+debugstack->count+loheader->count*loheader->esize+sizeof(ANSVALHeader[2]);
		const unsigned char *ptrstack[ANS_VAL_HISTSIZE]={0};
		const ANSVALHeader *loheader2=0, *hiheader2=0;
		const unsigned char *verdata=0, *unverdata=0;
		int nptrs=0;
		printf("Verified pops:\n");
		for(int k=0;k<ANS_VAL_HISTSIZE;++k)
		{
			if(verptr>=endptr)
				break;

			loheader2=(const ANSVALHeader*)ansval_ptrguard(debugstack->data, endptr, verptr, +sizeof(ANSVALHeader));
			verptr+=sizeof(ANSVALHeader);
			verdata=(const unsigned char*)ansval_ptrguard(debugstack->data, endptr, verptr, +loheader2->count*loheader2->esize);
			verptr+=loheader2->count*loheader2->esize;
			hiheader2=(const ANSVALHeader*)verptr;
			verptr+=sizeof(ANSVALHeader);

			ptrstack[nptrs++]=(const unsigned char*)loheader2;
			(void)verdata;
			(void)hiheader2;
		}
		for(int k=nptrs-1;k>=0;--k)
		{
			const unsigned char *ptr=ptrstack[k];
			loheader2=(const ANSVALHeader*)ptr;
			ptr+=sizeof(ANSVALHeader);
			verdata=ptr;
			ptr+=loheader2->count*loheader2->esize;
			hiheader2=(const ANSVALHeader*)ptr;
			ptr+=sizeof(ANSVALHeader);

			printf("  [%7d] %7d B    ", loheader2->idx, loheader2->count*loheader2->esize);
			ansval_printr(verdata, loheader2->esize, loheader2->count, 0);
			(void)hiheader2;
		}
		if(!nptrs)
			printf("  No data\n");
		printf("\n");

		printf("The error:\n");
		printf("  [%7d] Original %7d B    ", loheader->idx, loheader->count*loheader->esize);
		ansval_printr(data0, loheader->esize, loheader->count, 0);
		printf("  [%7d] Corrupt  %7d B    ", loheader->idx, count*esize);
		ansval_printr(data, esize, count, 0);
		printf("  [%7d] XOR      %7s      ", loheader->idx, "");
		ansval_printr(data0, esize, count, data);
		printf("\n");
		
		const unsigned char *unverptr=debugstack->data+debugstack->count;
		printf("Remaining pops:\n");
		for(int k=0;k<ANS_VAL_HISTSIZE;++k)
		{
			if(unverptr<=debugstack->data)
			{
				printf("  No data\n");
				break;
			}

			unverptr=(const unsigned char*)ansval_ptrguard(debugstack->data, endptr, unverptr, -(ptrdiff_t)sizeof(ANSVALHeader));
			hiheader2=(const ANSVALHeader*)unverptr;
			unverptr=(const unsigned char*)ansval_ptrguard(debugstack->data, endptr, unverptr, -(ptrdiff_t)hiheader2->count*hiheader2->esize);
			unverdata=unverptr;
			unverptr=(const unsigned char*)ansval_ptrguard(debugstack->data, endptr, unverptr, -(ptrdiff_t)sizeof(ANSVALHeader));
			loheader2=(const ANSVALHeader*)unverptr;
			
			if(memcmp(hiheader2, loheader2, sizeof(*hiheader2)))
			{
				printf("  Header Mismatch:  (idx,esize,count)\n");
				printf("    hiheader  0x%08X, 0x%08X, 0x%08X    %d, %d, %d\n",
					hiheader2->idx, hiheader2->esize, hiheader2->count,
					hiheader2->idx, hiheader2->esize, hiheader2->count
				);
				printf("    loheader  0x%08X, 0x%08X, 0x%08X    %d, %d, %d\n",
					loheader2->idx, loheader2->esize, loheader2->count,
					loheader2->idx, loheader2->esize, loheader2->count
				);
				break;
			}
			printf("  [%7d] %7d B    ", hiheader2->idx, hiheader2->count*hiheader2->esize);
			ansval_printr(unverdata, hiheader2->esize, hiheader2->count, 0);
		}
		printf("\n");
		LOG_ERROR("");
	}
	++popcount;
}
#endif


//LIFO Bypass Coder
typedef struct _BitPackerLIFO//bwd enc / fwd dec
{
	unsigned long long state;
	int enc_nwritten, dec_navailable;//bitcounts, only for tracking renorms
	unsigned char *dstbwdptr;
	const unsigned char *srcfwdptr, *streamend;
} BitPackerLIFO;
AWM_INLINE void bitpacker_enc_init(BitPackerLIFO *ec, const unsigned char *bufstart, unsigned char *bufptr0_OOB)
{
	memset(ec, 0, sizeof(*ec));
	ec->state=1ULL<<32;
	ec->enc_nwritten=33;
	ec->streamend=bufstart;
	ec->dstbwdptr=bufptr0_OOB;
}
AWM_INLINE void bitpacker_dec_init(BitPackerLIFO *ec, const unsigned char *bufptr0_start, const unsigned char *bufend)
{
	memset(ec, 0, sizeof(*ec));
	ec->srcfwdptr=bufptr0_start+8;
	ec->streamend=bufend;
	ec->state=*(const unsigned long long*)bufptr0_start;
	ec->dec_navailable=FLOOR_LOG2_P1(ec->state);
}
AWM_INLINE void bitpacker_enc_flush(BitPackerLIFO *ec)
{
	ec->dstbwdptr-=8;
#ifdef _DEBUG
	if(ec->dstbwdptr<ec->streamend)
		LOG_ERROR("IntPacker Encoder OOB:  dstbwdptr = 0x%016zX < 0x%016zX", ec->dstbwdptr, ec->streamend);
#endif
	*(unsigned long long*)ec->dstbwdptr=ec->state;
}
AWM_INLINE void bitpacker_enc(BitPackerLIFO *ec, int inbits, int sym)
{
#ifdef _DEBUG
	if(!inbits)
		LOG_ERROR("BitPacker inbits=0");
#endif
	//renorm then push inbits
	ec->enc_nwritten+=inbits;
	if(ec->enc_nwritten>64)//renorm on overflow
	{
		ec->enc_nwritten-=32;
		ec->dstbwdptr-=4;
#ifdef _DEBUG
		if(ec->dstbwdptr<ec->streamend)
			LOG_ERROR("IntPacker OOB:  dstbwdptr = 0x%016zX < 0x%016zX", ec->dstbwdptr, ec->streamend);
#endif
		*(unsigned*)ec->dstbwdptr=(unsigned)ec->state;
		ec->state>>=32;
#ifdef ANS_VAL
		ansval_push(&ec->state, sizeof(ec->state), 1);
#endif
	}
	ec->state=ec->state<<inbits|sym;
#ifdef ANS_VAL
	ansval_push(&ec->state, sizeof(ec->state), 1);
#endif
}
AWM_INLINE int bitpacker_dec(BitPackerLIFO *ec, int outbits)
{
	int sym;
#ifdef _DEBUG
	if(!outbits)
		LOG_ERROR("BitPacker outbits=0");
#endif
	sym=ec->state&((1ULL<<outbits)-1);

	//pop outbits then renorm
#ifdef ANS_VAL
	ansval_check(&ec->state, sizeof(ec->state), 1);
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
			LOG_ERROR("IntPacker OOB:  srcfwdptr = 0x%016zX >= 0x%016zX", ec->srcfwdptr, ec->streamend);
#endif
		ec->state=ec->state<<32|*(const unsigned*)ec->srcfwdptr;
		ec->srcfwdptr+=4;
	}
	return sym;
}


#ifdef __cplusplus
}
#endif
#endif
