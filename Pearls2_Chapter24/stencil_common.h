/*
    The MIT License (MIT)
    
    Copyright (c) 2015 OpenVec
    
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
    
    
    Authors:
    Paulo Souza
    Leonardo Borges
    Cedric Andreolli
    Philippe Thierry

*/
#define UNROLL 4

#define NPOP 8


#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))


#define W0 -9.164532312924f
#define W1 +1.777777777777f
#define W2 -3.111111111111f
#define W3 +7.542087542087e-1f
#define W4 -1.767676767676e-1f 
#define W5 +3.480963480963e-2f
#define W6 -5.180005180005e-3f
#define W7 +5.074290788576e-4f
#define W8 -2.428127428127e-5f


double wallt()
{
  struct timeval tv;
  /*
  struct timezone tz;
  gettimeofday(&tv, &tz);
  */
  gettimeofday(&tv, (struct timezone*)0);
  return tv.tv_sec + tv.tv_usec*1e-6;
}


/* To get perfect alignement n1 must be multiple of OV_FLOAT_WIDTH */

/* Array sizes for big memory machine (requires 3.9 GB) */
/*
#define BSIZEZ 140
#define BSIZEY 24


#define ARRAY_SIZES \
  size_t const n1=704;\
  size_t const n2=670;\
  size_t const n3=714;
*/

/* Array sizes for small memory machine (requires 400 MB) */
#define BSIZEZ 13
#define BSIZEY 9

#define ARRAY_SIZES \
  size_t const n1=320;\
  size_t const n2=320;\
  size_t const n3=320;
