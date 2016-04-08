#ifndef SMITH_WATERMAN_H
#define SMITH_WATERMAN_H

extern "C" struct BufferInfo
{
    unsigned int m_dwBufferOffset;
    unsigned int m_dwBufferLength;
    unsigned int m_dwFileOffset;
    unsigned int m_dwDescOffset;
    unsigned int m_dwSWScore;
};

#define VECTOR_SIZE  64
#define MAX_DESC_LEN 256
#endif