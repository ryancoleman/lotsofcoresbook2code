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
#ifndef _OV_TEST_UTILS_H_
#define _OV_TEST_UTILS_H_


#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))


#ifdef _OV_NEON
/* That's for ARM low prec fast SQRT */
#define MAX_ERR 0.003f

#else
#define MAX_ERR 0.0005f
#endif


#define FAIL(a,b,funcname) do{\
fprintf(stderr, "Error on function " funcname "\n");\
fprintf(stderr, "a = %f != b %f\n", a, b);\
exit(1);\
}while(0)


#define TEST(a, b, funcname) \
do { \
if ((fabs(a)<MAX_ERR) && (fabs(a-b) > MAX_ERR)) FAIL(a,b,funcname); \
else if (fabs(1.0f - (b/a)) > MAX_ERR) FAIL(a,b,funcname);\
}while(0)


#endif
