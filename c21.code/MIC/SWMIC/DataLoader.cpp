#include <iostream>
#include <fstream>
#include <memory>
#include <string.h>
#include <emmintrin.h>
#include <map>
#include <vector>
#include "SmithWaterman.h"
#include "DataLoader.h"
using namespace std;

bool LoadFastaSequence(char * FastaFile, char *& bBuffer, unsigned int & dwBufferLen, char * szQueryDesc)
{
        char szBuffer[MAX_DESC_LEN];
        unsigned long dwFileLength = 0;
        ifstream inFastaFile(FastaFile, ios::in);
        
        dwBufferLen = 0;
        
        if(inFastaFile.fail())
        {    
            return false;
        }
        inFastaFile.seekg (0, ios::end);
        dwFileLength = inFastaFile.tellg();
        inFastaFile.seekg(0);
        bBuffer = (char *)_mm_malloc(dwFileLength, 64);
        inFastaFile.getline(szBuffer, sizeof(szBuffer) - 1, 0x0A);
        
        strncpy(szQueryDesc, szBuffer, strlen(szBuffer));
        
        while(!inFastaFile.eof())
        {
            inFastaFile.getline(szBuffer, sizeof(szBuffer) - 1, 0x0A);
            memcpy(bBuffer + dwBufferLen, szBuffer, strlen(szBuffer) ? strlen(szBuffer) - 1 : 0 );
            dwBufferLen += strlen(szBuffer) ? strlen(szBuffer) - 1 : 0;
        }
        return true;
}

bool LoadFastaDatabase(char * FastaFile, char * &bBuffer, unsigned int & dwBufferLen, BufferInfo *& pBufferInfo, unsigned int & dwInfoLen)
{
    char szBuffer[1024];
    char szSequence[256 * 1024] = {0};
    vector<BufferInfo>   vecBufferInfo;
    ifstream inFastaFile(FastaFile, ios::in);
    unsigned long qwFileSize     = 0L;

    if(inFastaFile.fail())
        return false;
    
    inFastaFile.seekg (0, ios::end); 
    qwFileSize = inFastaFile.tellg();
    inFastaFile.seekg(0, ios::beg);
    
    bBuffer = (char *)_mm_malloc(qwFileSize, 64);
    if(!bBuffer)
        return false;
    
	BufferInfo buffInfo = {0};
    while(!inFastaFile.eof())
    {
        inFastaFile.getline(szBuffer, sizeof(szBuffer) - 1, 0x0A);
        if(szBuffer[0] == '>')
        {
            vecBufferInfo.push_back(buffInfo);
            buffInfo.m_dwFileOffset = inFastaFile.tellg();
            buffInfo.m_dwDescOffset = buffInfo.m_dwFileOffset - strlen(szBuffer);
            buffInfo.m_dwBufferOffset = ALIGN_LENGTH(buffInfo.m_dwBufferOffset + buffInfo.m_dwBufferLength);
            buffInfo.m_dwBufferLength = 0;
        }
        else
        {
             memcpy(bBuffer + buffInfo.m_dwBufferOffset + buffInfo.m_dwBufferLength, szBuffer, strlen(szBuffer));
             buffInfo.m_dwBufferLength += strlen(szBuffer);
        }
    }
	pBufferInfo = (BufferInfo *)_mm_malloc(vecBufferInfo.size() * sizeof(BufferInfo), 64);
    dwInfoLen = vecBufferInfo.size();
    dwBufferLen = qwFileSize;
	int i = 0;
    for( vector<BufferInfo>::iterator it = vecBufferInfo.begin(); it != vecBufferInfo.end(); it++)
        memcpy(&pBufferInfo[i++], &(*it), sizeof(BufferInfo));
        
    return true;
}
void Normalize(char * c, unsigned long dwLen)
{
    __m128i * pVector = (__m128i *)c;
    __m128i   offset;
    offset = _mm_insert_epi16(offset, 0x4141, 0); //Offset buffer by A
    offset = _mm_shufflelo_epi16(offset, 0);
    offset = _mm_shuffle_epi32(offset, 0);

    while((void *)pVector < (void *)(c + dwLen))
    {
        *pVector = _mm_sub_epi8(*pVector, offset);
        pVector++;
    }
}

bool GetSequenceDesc(char * FastaFile, char * szSeqDesc, unsigned int nDescOffset)
{
    ifstream inFastaFile(FastaFile, ios::in);
    if(inFastaFile.fail())
        return false;
    inFastaFile.seekg(nDescOffset);
    inFastaFile.getline(szSeqDesc, MAX_DESC_LEN - 1);
    return true;
}
