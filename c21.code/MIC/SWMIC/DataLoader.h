
#include "SmithWaterman.h"

//bool LoadDataFromFasta(char * FastaFile, char * IndexFile, unsigned int dwDesiredLength, char *& bBuffer, unsigned int & nBufferLen);
void Normalize(char * c, unsigned long dwLen);
bool LoadFastaDatabase(char * FastaFile, char * &bBuffer, unsigned int & dwBufferLen, BufferInfo *& pBufferInfo, unsigned int & dwInfoLen);
bool LoadFastaSequence(char * FastaFile, char *& bBuffer, unsigned int & dwBufferLen, char * szQueryDesc);
bool GetSequenceDesc(char * FastaFile, char * szSeqDesc, unsigned int nDescOffset);

#define ALIGN_LENGTH(x)((x) + VECTOR_SIZE - (x) % VECTOR_SIZE)

struct FastaIdx
{
    FastaIdx(unsigned int dwRows, unsigned int dwLength, unsigned int dwOffset, unsigned int dwDescOffset) : m_dwRows(dwRows), m_dwLength(dwLength), m_dwOffset(dwOffset), m_dwDescOffset(dwDescOffset){}
    unsigned int m_dwRows;
    unsigned int m_dwLength;
    unsigned int m_dwOffset;
    unsigned int m_dwDescOffset;
};