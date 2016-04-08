# /******************************************************************************
# * Copyright 2014 Intel Corporation All Rights Reserved.                      *
# *                                                                            *
# * The source code, information and material ("Material") contained herein    *
# * is owned by Intel Corporation or its suppliers or licensors, and title to  *
# * such Material remains with Intel Corporation or its suppliers or           *
# * licensors.                                                                 *
# * The Material contains proprietary information of Intel or its suppliers    *
# * and licensors. The Material is protected by worldwide copyright laws and   *
# * treaty provisions. No part of the Material may be used, copied,            *
# * reproduced, modified, published, uploaded, posted, transmitted,            *
# * distributed or disclosed in any way without Intel's prior express written  *
# * permission. No license under any patent, copyright or other intellectual   *
# * property rights in the material is granted to or conferred upon you,       *
# * either expressly, by implication, inducement, estoppel or otherwise. Any   *
# * license under such intellectual property rights must be express and        *
# * approved by Intel in writing.                                              *
# ******************************************************************************/

import numpy as np
import pymic as mic
from gpaw.mic import stream

# load the .so on the MIC
library = None
if stream is not None:
    library = stream.get_device().load_library("libmicblas.so")

_data_type_map = {
    # Python types
    int : 0,
    float : 1,
    complex : 2,

    # Numpy dtypes
    np.dtype(np.int64) : 0,
    np.dtype(np.float64) : 1,
    np.dtype(np.complex128) : 2
}
    
def map_dtype(dtype):
    return _data_type_map[dtype]
    
def gemm(alpha, a, b, beta, c, transa='n'):
    # we want to make sure that we only use OffloadArrays here
    assert isinstance(a, mic.OffloadArray)
    assert isinstance(b, mic.OffloadArray)
    assert isinstance(c, mic.OffloadArray)
    
    # determine the datatype and map it to int
    dt = map_dtype(a.dtype)
    
    # determine sizes of the matrices
    am = a.shape[0]
    ak = np.prod(a.shape[1:])
    bk = b.shape[0]
    bn = np.prod(b.shape[1:])
    cm = c.shape[0]
    cn = np.prod(c.shape[1:])
    
    # just some safety checks
    if transa == 'n':
        assert am == bn
        assert ak == cn
        assert bk == cm
        trans = 0
        m, n, k = ak, bk, am
        lda = a.array.strides[0] / a.array.strides[-1]
        ldb = b.array.strides[0] / b.array.strides[1]
        ldc = c.array.strides[0] / c.array.strides[-1]
    else:
        assert am == cn
        assert ak == bn
        assert bk == cm
        trans = 1
        m, n, k = am, bk, ak
        lda = k
        ldb = b.array.strides[0] / b.array.strides[-1]
        ldc = c.array.strides[0] / c.array.strides[1]

    if a.dtype in [np.complex]:
        alpha = complex(alpha)
        beta = complex(beta)
            
    # perform the offload
    stream.invoke(library.mic_gemm, dt, 
                  a, b, c, 
                  m, n, k, 
                  lda, ldb, ldc, 
                  alpha, beta, trans)
    stream.sync()
                         
def rk(alpha, a, beta, c, trans='c'):
    """Rank-k update of a matrix."""

    assert isinstance(a, mic.OffloadArray)
    assert isinstance(c, mic.OffloadArray)

    dt = map_dtype(a.dtype)
    
    # determine sizes of the matrices
    am = a.shape[0]
    ak = np.prod(a.shape[1:])
    ck = c.shape[0]
    cn = np.prod(c.shape[1:])

    n, k = am, ak
    ldc = c.array.strides[0] / c.array.strides[1]

    
    if a.dtype in [np.complex]:
        alpha = complex(alpha)
        beta = complex(beta)
        
    # perform the offload
    stream.invoke(library.mic_syrk, dt, a, c, n, k, 
                  ldc, alpha, beta)
    stream.sync()
                             
    
def r2k(alpha, a, b, beta, c):
    """Rank-2k update of a matrix."""

    assert isinstance(a, mic.OffloadArray)
    assert isinstance(b, mic.OffloadArray)
    assert isinstance(c, mic.OffloadArray)

    assert(map_dtype(a.dtype) != 2)
    
    # determine sizes of the matrices
    am = a.shape[0]
    ak = np.prod(a.shape[1:])
    bm = b.shape[0]
    bk = np.prod(b.shape[1:])
    ck = c.shape[0]
    cn = np.prod(c.shape[1:])

    n, k = am, ak
    ldc = c.array.strides[0] / c.array.strides[1]

    stream.invoke(library.mic_dsyr2k, a, b, c, n, k, 
                  ldc, alpha, beta)
    stream.sync()
    
# if __name__ == '__main__':    
    # m = 10
    # n = 10
    # k = 10
    # alpha = 1.0
    # beta = 0.0


    # # construct some matrices
    # a = np.random.random(m*k).reshape((m, k))
    # b = np.random.random(k*n).reshape((k, n))
    # c = np.random.random(k*n).reshape((k, n))

    # # print the input of numpy's MxM if it is small enough
    # Am = np.matrix(a)
    # Bm = np.matrix(b)
    # Cm = Am*Bm
    # print "numpy gives us: "
    # print "--------------------------------------"
    # print Cm
    # print 
    # print   

    # oa_a = mic.OffloadArray(a)
    # oa_b = mic.OffloadArray(b)
    # oa_c = mic.OffloadArray(c)

    # oa_a.update_device()
    # oa_b.update_device()
    # oa_c.update_device()
    # gemm(alpha, oa_a, oa_b, beta, oa_c)
    # oa_c.update_host()
    # print "offloading gives us: "
    # print "--------------------------------------"
    # print oa_c
