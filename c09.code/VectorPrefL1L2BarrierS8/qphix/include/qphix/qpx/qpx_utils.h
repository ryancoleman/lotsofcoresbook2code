#ifndef __QPHIX_QPX_UTILS_H__
#define __QPHIX_QPX_UTILS_H__

#if !defined(__xlc__) && !defined(__xlC__)
#include <qpxintrin.h>
#endif

inline vector4double _v4d_int2mask(unsigned int msk) {

    vector4double ret = vec_gpci(00123);
    msk = msk & 0x0F;
    
    switch (msk) {
      /*        
                case 0: ret = vec_gpci(00123); break;
		case 1: ret = vec_gpci(01237); break;
		case 2: ret = vec_gpci(00163); break;
		case 3: ret = vec_gpci(00167); break;
		case 4: ret = vec_gpci(00523); break;
		case 5: ret = vec_gpci(00527); break;
		case 6: ret = vec_gpci(00563); break;
		case 7: ret = vec_gpci(00567); break;
		case 8: ret = vec_gpci(04123); break;
		case 9: ret = vec_gpci(04127); break;
		case 10: ret = vec_gpci(04163); break;
		case 11: ret = vec_gpci(04167); break;
		case 12: ret = vec_gpci(04523); break;
		case 13: ret = vec_gpci(04527); break;
		case 14: ret = vec_gpci(04563); break;
		case 15: ret = vec_gpci(04567); break;
      */
                case 0: ret = vec_gpci(00123); break;
	        case 1: ret = vec_gpci(04123); break;
		case 7: ret = vec_gpci(04563); break;
		case 8: ret = vec_gpci(00127); break;  
		case 14: ret = vec_gpci(00567); break;
		default:
			masterPrintf("Mask is out of range\n");
			ret = vec_gpci(00123);
 			break;
	}
	return ret;
}

#endif
