# /******************************************************************************
#  * Copyright 2014 Intel Corporation All Rights Reserved.                      *
#  *                                                                            *
#  * The source code, information and material ("Material") contained herein    *
#  * is owned by Intel Corporation or its suppliers or licensors, and title to  *
#  * such Material remains with Intel Corporation or its suppliers or           *
#  * licensors.                                                                 *
#  * The Material contains proprietary information of Intel or its suppliers    *
#  * and licensors. The Material is protected by worldwide copyright laws and   *
#  * treaty provisions. No part of the Material may be used, copied,            *
#  * reproduced, modified, published, uploaded, posted, transmitted,            *
#  * distributed or disclosed in any way without Intel's prior express written  *
#  * permission. No license under any patent, copyright or other intellectual   *
#  * property rights in the material is granted to or conferred upon you,       *
#  * either expressly, by implication, inducement, estoppel or otherwise. Any   *
#  * license under such intellectual property rights must be express and        *
#  * approved by Intel in writing.                                              *
#  ******************************************************************************/

all: libmicblas.so

libmicblas.so: micblas.c
	icc -I $(PYMIC_BASE)/include -fPIC -shared -mmic -g -O2 -L$(MKLROOT)/lib/mic -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -o $@ $< -lm 
	chmod a+rX libmicblas.so

clean:
	rm -f libmicblas.so
	
realclean: clean
	rm -f *~

push_all: push push_mkl push_comp push_others
	
push_mkl:	
	sudo scp $(MKLROOT)/lib/mic/*.so* mic0:/lib64
	
push_comp:	
	sudo scp $(TBBROOT)/../compiler/lib/mic/*.so* mic0:/lib64
	
