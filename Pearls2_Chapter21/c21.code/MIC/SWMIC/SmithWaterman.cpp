#include <stdio.h>
#include <iostream>
#include <malloc.h>
#include <immintrin.h>
#include <time.h>
#include <stdlib.h> 
#include "DataLoader.h"
#include "Blosum62.h"
#include "SmithWaterman.h"
#include <pthread.h>
#include <string.h>
#include <sys/time.h>
#include <scif.h>

using namespace std;

#define MAX_RES 10
#define NUM_THD 240

__declspec(target(mic)) __declspec(align(64)) volatile int g_nThdIndex = 0;

extern "C" struct SWContext
{ 
     char *            m_query_sequence;
     char *            m_query_profile_word;
     int               m_query_length;
     char *            m_db_sequence;
     BufferInfo *      m_pBufferInfo;
     int               m_nAminoCount;
     unsigned short    m_gap_open;
     unsigned short    m_gap_extend;
     void *            m_pWorkspace;
};

struct SWRetVal
{    
	SWRetVal(unsigned short score, unsigned int seq_pos) : m_score(score), m_seq_pos(seq_pos){}
	SWRetVal() :  m_score(0), m_seq_pos(0){}
    unsigned short m_score;
    unsigned int   m_seq_pos;
};

__declspec(target(mic)) void FillProfile(char* query_profile_word, char * query_sequence, unsigned short query_sequence_len, unsigned char vecLen)
{
    unsigned short  segLen = (query_sequence_len + vecLen - 1)/vecLen;
    unsigned short nSize = segLen * vecLen;
    int LenAmino = 25;
    int data = 0;
        
    for(int i = 0; i < LenAmino; i++)
    {
        for(int j = 0; j < segLen; ++j)
        {
            for(int k = j; k < nSize; k += segLen)
            {
                if(k > query_sequence_len)
                    *query_profile_word++  = 0;
                else
                    *query_profile_word++ = OrderedBlosum62[i][query_sequence[k]];
            }
        }
    }
    return;
}


 __declspec(target(mic)) unsigned int SmithWatermanMIC(char *     query_sequence,
                         char *    query_profile_word,
                         const int                 query_length,
                         char *     db_sequence,
                         const int                 db_length,
                         unsigned short      gap_open,
                         unsigned short      gap_extend,
                         void *   pWorkspace)
{
    int     i, j, k;
    unsigned short score = 0;

    int     cmp;
    int     iter = (query_length + 15) / 16;
    int     max_val = 0xFFFF;
    
#ifdef __MIC__
    __m128i *pScore;
    
    __m256i *p;
    __m256i *workspace = (__m256i *) pWorkspace;
    __m256i *pHLoad, *pHStore;
    __m256i *pE;
    
    __m512i E, F, H;

    __m512i v_maxscore;
    __m512i v_gapopen;
    __m512i v_gapextend;

    __m512i v_zero;
    __m512i v_max;
    __m512i v_max_comp;
    __m512i v_temp;
    
    __m512i vec_score;
    
    
    
    // Load gap open penalty to all elements of a constant 
    v_gapopen  = _mm512_extload_epi32((void *)&gap_open  , _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST_1X16, _MM_HINT_NONE);
     
    // Load gap extension penalty to all elements of a constant 
    v_gapextend =_mm512_extload_epi32((void *)&gap_extend  , _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST_1X16, _MM_HINT_NONE);

    v_max =_mm512_extload_epi32((void *)&max_val, _MM_UPCONV_EPI32_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
    v_max_comp = v_max;

    // Load v_maxscore with the zeros. 
    v_maxscore = _mm512_xor_epi32(v_maxscore, v_maxscore);
    
    v_zero = v_maxscore;
    
    // Zero out the storage vector 
    k = 2 * iter;
    int counter = 0;

    p = workspace;
    for (i = 0; i < k; i++)
    {
        _mm512_extstore_epi32(p++, v_maxscore, _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);
    }
    
    pE = workspace;
    pHStore = pE + iter;
    pHLoad = pHStore + iter;

    for (i = 0; i < db_length; ++i)
    {
        pScore = (__m128i *) query_profile_word + db_sequence[i] * iter;
    
        F = _mm512_xor_epi32(F, F);
    
        //load the next h value 
        H = _mm512_extload_epi32 (pHStore + iter - 1, _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST_16X16, _MM_HINT_NONE );
        H = _mm512_alignr_epi32(H, v_zero, 15);
    
        p = pHLoad;
        pHLoad = pHStore;
        pHStore = p;
        
        for (j = 0; j < iter; j++)
        {
            // Load E values 
            E = _mm512_extload_epi32 (pE + j, _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST_16X16, _MM_HINT_NONE );
            
            // Add score to H     
            vec_score = _mm512_extload_epi32 (pScore++, _MM_UPCONV_EPI32_SINT8, _MM_BROADCAST_16X16, _MM_HINT_NONE );
            H = _mm512_add_epi32(H, vec_score);
            v_max = _mm512_max_epi32 (v_max, H);
            H = _mm512_max_epi32(H, v_zero);
    
            // Update highest score encountered this far 
            v_maxscore = _mm512_max_epi32 (v_maxscore, H);
            
              //get max from H, E and F 
            H = _mm512_max_epi32 (H, E);
            H = _mm512_max_epi32 (H, F);
            
            // save H values 
            _mm512_extstore_epi32(pHStore + j, H, _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE );
        
            // subtract the gap open penalty from H 
            H = _mm512_sub_epi32 (H, v_gapopen);
            H = _mm512_max_epi32(H, v_zero);

            // update E value 
            E = _mm512_sub_epi32 (E, v_gapextend);
            E = _mm512_max_epi32 (E, H);
        
            // update F value 
            F = _mm512_sub_epi32 (F, v_gapextend);
            F = _mm512_max_epi32 (F, H);
            
            // save E values 
            _mm512_extstore_epi32(pE + j, E, _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE );
 
            // load the next h value 
            H = _mm512_extload_epi32 (pHLoad + j, _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST_16X16, _MM_HINT_NONE );
        }
    
        // reset pointers to the start of the saved data 
        j = 0;
        H = _mm512_extload_epi32 (pHStore + j, _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST_16X16, _MM_HINT_NONE );

        //  the computed F value is for the given column.  since 
        //  we are at the end, we need to shift the F value over 
        //  to the next column. 
    
        F = _mm512_alignr_epi32(F, v_zero, 15);
      
        v_temp = _mm512_sub_epi32 (H, v_gapopen);
        v_temp = _mm512_max_epi32 (v_temp, v_zero);
        
        cmp =  _mm512_cmp_epi32_mask (F, v_temp, _MM_CMPINT_GT);
 
        while (cmp != 0x0000) 
        {
            E = _mm512_extload_epi32 (pE + j, _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST_16X16, _MM_HINT_NONE );

            H = _mm512_max_epi32 (H, F);

            //save H values 
            _mm512_extstore_epi32(pHStore + j, H, _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE );
          
            // update E in case the new H value would change it
            H = _mm512_sub_epi32 (H, v_gapopen);
            H = _mm512_max_epi32(H, v_zero);
            
            E = _mm512_max_epi32 (E, H);
            _mm512_extstore_epi32 (pE + j, E, _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE );

            // update F value 
            F = _mm512_sub_epi32 (F, v_gapextend);
            F = _mm512_max_epi32 (F, v_zero);
        
            j++;
            if (j >= iter)
            {
                j = 0;
                F = _mm512_alignr_epi32(F, v_zero, 15);
            }
            H = _mm512_extload_epi32 (pHStore + j, _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST_16X16, _MM_HINT_NONE );

            v_temp = _mm512_sub_epi32 (H, v_gapopen);
            v_temp = _mm512_max_epi32 (v_temp, v_zero);
            cmp =  _mm512_cmp_epi32_mask (F, v_temp, _MM_CMPINT_GT);
        }
    }
    
    //Treat overflow case
    cmp = _mm512_cmp_epi32_mask(v_max, v_max_comp, _MM_CMPINT_GT);
    if(cmp)
    {
        return 0xFFFF;
    }
    
    // extract the largest score 
    score = _mm512_reduce_max_epi32(v_maxscore);
#endif
    return score;

}

__declspec(target(mic)) void * fnThreadFunct(void * context)
{
        SWContext * pContext = (SWContext *)context;
        int nThdIndex;
        unsigned short res;
        
        while((nThdIndex = _InterlockedIncrement((void *)&g_nThdIndex)) - 1 < pContext->m_nAminoCount)
        {
            res = SmithWatermanMIC( pContext->m_query_sequence,
                                            pContext->m_query_profile_word,
                                            pContext->m_query_length,
                                            pContext->m_db_sequence + pContext->m_pBufferInfo[nThdIndex].m_dwBufferOffset,
                                            pContext->m_pBufferInfo[nThdIndex].m_dwBufferLength,
                                            pContext->m_gap_open,
                                            pContext->m_gap_extend,
                                            pContext->m_pWorkspace);

            pContext->m_pBufferInfo[nThdIndex].m_dwSWScore = res;
                    
        }
        delete pContext;
        return NULL;
}

__declspec(target(mic)) void DatabaseSearch( char *     query_sequence,
                                             char *    query_profile_word,
                                             const int                 query_length,
                                             char *                 db_sequence,
                                             const int                 db_length,
                                             unsigned short      gap_open,
                                             unsigned short      gap_extend,
                                             BufferInfo *   pBufferInfo,
                                             unsigned int   dwDatabeseSeqNum)
{
    pthread_t thread_vec[NUM_THD];
    void * retval = NULL;

    //Compute query profile
    FillProfile(query_profile_word, query_sequence, query_length, 16);
    
    //Allocate workspace
    unsigned int * pWorkspace = (unsigned int *) _mm_malloc (NUM_THD * 3 * sizeof(__m512i) * ((query_length + 15) / 16) , 64);

    for(int j = 0; j < NUM_THD; j++)
    {
        SWContext * pContext = new SWContext;
        pContext->m_query_sequence = query_sequence;
        pContext->m_query_profile_word = query_profile_word;
        pContext->m_query_length = query_length;
        pContext->m_db_sequence = db_sequence;
        pContext->m_pBufferInfo = pBufferInfo;
        pContext->m_gap_open = gap_open;
        pContext->m_gap_extend = gap_extend;
        pContext->m_nAminoCount = dwDatabeseSeqNum;
        pContext->m_pWorkspace = pWorkspace + (j *  3 * sizeof(__m512i) / sizeof(int) * ((query_length + 15) / 16));
        pthread_create(&thread_vec[j], NULL, fnThreadFunct, (void *)pContext);
    }
    for(int j = 0; j < NUM_THD; j++)
    {
        pthread_join(thread_vec[j], &retval);
    }
    
    _mm_free(pWorkspace);
}

void MergeResult(SWRetVal wRes, SWRetVal * vecRes)
{
    for(int i = MAX_RES - 1; i >= 0; i--)
        if(wRes.m_score > vecRes[i].m_score)
        {
            vecRes[i + 1] = vecRes[i];
            vecRes[i] = wRes;
        }
        else
            return;
}

void PrintResults(BufferInfo * pBufferInfo, unsigned int dwNumSeq, char * szFastaFile)
{
	SWRetVal vecRes[MAX_RES];// = {0};
	char     szDesc[MAX_DESC_LEN];
	
	for(int i = 0; i < dwNumSeq; i++)
		MergeResult(SWRetVal(pBufferInfo[i].m_dwSWScore, i), vecRes);
		
	for(int i = 0; i < MAX_RES; i++)
	{
		if(vecRes[i].m_score)
		{
			GetSequenceDesc(szFastaFile, szDesc, pBufferInfo[vecRes[i].m_seq_pos].m_dwDescOffset);
			printf("Score %d : %s\n", vecRes[i].m_score, szDesc);
		}
	}	
}
int main(int argc, char * argv[])
{
    char             *pDBBuffer = 0L;
    char             *pQueryBuffer = 0L;
    unsigned int    dwDatabaseLen = 0;
    unsigned int    dwDatabeseSeqNum = 0;
    unsigned int    dwQueryLen = 0;
    BufferInfo         *pBufferInfo = 0L; 
    char            szQueryDesc[MAX_DESC_LEN] = {0};
    
    if(argc < 3)
        printf("Usage:\n%s <Database File> <Query File> [Gap Open Cost] [Gap Extend Cost]\n", argv[0]);

    if(!LoadFastaDatabase(argv[1], pDBBuffer, dwDatabaseLen, pBufferInfo, dwDatabeseSeqNum))
    {
        printf("Wrong database file\n");
        return 0;
    }
    printf("Database loaded successfully.\n");
    
    if(!LoadFastaSequence(argv[2], pQueryBuffer, dwQueryLen, szQueryDesc))
    {
        printf("Wrong query file\n");
        return 0;
    }
    printf("Query loaded successfully.\n");
    
    Normalize(pDBBuffer, dwDatabaseLen);
    Normalize(pQueryBuffer, dwQueryLen);
    
    
    unsigned short segLen = (dwQueryLen + 16 - 1)/16;
    unsigned short nSize = segLen * 16;
    unsigned char  GStart = 8;
    unsigned char  GCont = 2;
    char *            query_profile_word = (char *)_mm_malloc(25 * nSize * sizeof(char), 64);
    unsigned int   res = 0;
    
    if(argc > 3)
        GStart = atoi(argv[3]) > 0 ? atoi(argv[3]) : 8;
    
    if(argc > 4)
        GCont= atoi(argv[4]) > 0 ? atoi(argv[4]) : 2;
    
    printf("Database contains %d entries.\n", dwDatabeseSeqNum);
    printf("Starting database search for sequence\n%s\nof length %d with open score %d and extend score %d.\n", szQueryDesc, dwQueryLen, GStart, GCont);
    
    timespec tpStart;
    timespec tpEnd;
    
    clock_gettime(CLOCK_MONOTONIC, &tpStart);
    #pragma offload target(mic:0) \
    in(pQueryBuffer:length(dwQueryLen) alloc_if(1) free_if(0)) \
    in(pDBBuffer:length(dwDatabaseLen) alloc_if(1) free_if(0)) \
    in(query_profile_word:length(25 * nSize) alloc_if(1) free_if(0)) 
    {}

    clock_gettime(CLOCK_MONOTONIC, &tpEnd);
    double time_mem = ((double)tpEnd.tv_sec + 1.0e-9 * (double)tpEnd.tv_nsec) - ((double)tpStart.tv_sec + 1.0e-9 * (double)tpStart.tv_nsec);
    
    clock_gettime(CLOCK_MONOTONIC, &tpStart);
    
    #pragma offload target(mic:0) \
    nocopy(pQueryBuffer:length(dwQueryLen) alloc_if(0) free_if(0)) \
    nocopy(pDBBuffer:length(dwDatabaseLen) alloc_if(0) free_if(0)) \
    nocopy(query_profile_word:length(25 * nSize) alloc_if(0) free_if(0)) \
    inout(pBufferInfo:length(dwDatabeseSeqNum) alloc_if(1) free_if(1))
    DatabaseSearch( pQueryBuffer, query_profile_word, dwQueryLen, pDBBuffer, dwDatabaseLen, GStart, GCont, pBufferInfo, dwDatabeseSeqNum);
    
    clock_gettime(CLOCK_MONOTONIC, &tpEnd);
    
    double time_mic = ((double)tpEnd.tv_sec + 1.0e-9 * (double)tpEnd.tv_nsec) - ((double)tpStart.tv_sec + 1.0e-9 * (double)tpStart.tv_nsec);
    clock_gettime(CLOCK_MONOTONIC, &tpStart);
    
    #pragma offload target(mic:0) \
    nocopy(pQueryBuffer:length(dwQueryLen) alloc_if(0) free_if(1)) \
    nocopy(pDBBuffer:length(dwDatabaseLen) alloc_if(0) free_if(1)) \
    nocopy(query_profile_word:length(25 * nSize) alloc_if(0) free_if(1)) 
    {}

    clock_gettime(CLOCK_MONOTONIC, &tpEnd);
    double time_free = ((double)tpEnd.tv_sec + 1.0e-9 * (double)tpEnd.tv_nsec) - ((double)tpStart.tv_sec + 1.0e-9 * (double)tpStart.tv_nsec);
    printf("Time_mem = %f Time_mic = %f Time_free=%f\n", time_mem, time_mic, time_free);
    
	PrintResults(pBufferInfo, dwDatabeseSeqNum, argv[1]);

    _mm_free(pQueryBuffer);
    _mm_free(pDBBuffer);
    _mm_free(query_profile_word);
    _mm_free(pBufferInfo);
}

