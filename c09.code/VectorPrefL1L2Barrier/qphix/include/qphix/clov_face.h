#ifndef QPHIX_CLOV_FACE_H
#define QPHIX_CLOV_FACE_H

namespace QPhiX {

  // File scope...
  // (1 + gamma_T) dagger psi
  template<typename FT, int veclen, int soalen, bool compress>
    void ClovDslash<FT,veclen,soalen,compress>::packFaceDir(int tid,
						     const FourSpinorBlock *psi,
						     FT *res,
						     int cb, int dir, int fb, int isPlus)
  {
    int Nxh = s->Nxh();
    int Ny = s->Ny();
    int Nz = s->Nz();
    int Nt = s->Nt();
    int ngy = s->nGY();
    int n_cores = s->getNumCores();
    int Pxyz = s->getPxyz();
    int Pxy = s->getPxy();
    
    int ngy_pack = ngy;
    // If packing x face and ngy is 1, we need to pick alternate rows in Y, so use ngy = 2
    if(dir == 0 && ngy_pack == 1) ngy_pack = 2;
    int n_soa_x = Nxh / soalen;
    int n_blocks_y = Ny/ngy_pack;
    int lens1[4] = {n_soa_x, n_blocks_y, Nz, Nt};
    int lens[4] = {n_soa_x, n_blocks_y, Nz, Nt};
    
    lens[dir] = 1;
    
    int npkts = lens[0] * lens[1] * lens[2] * lens[3];
    int pktsize = veclen;
    unsigned int mask = -1;
    unsigned int mask_xodd[2] = {0, 0};
    if(dir == 0) {
      pktsize = ngy_pack / 2;
      if(ngy == 1) 
	mask = (fb == 0 ? 1 : (1 << (soalen-1)));
      else {
	for(int i = 0; i < veclen; i+=2*soalen)	mask_xodd[0] |= ((fb == 0 ? 1 : (1 << (soalen - 1))) << i);
	mask_xodd[1] = mask_xodd[0] << soalen;
	if(fb == 1) {
	  unsigned int tmp = mask_xodd[0];
	  mask_xodd[0] = mask_xodd[1];
	  mask_xodd[1] = tmp;
	}
      }
    }
    else if(dir == 1) {
      pktsize = soalen;
      mask = (1 << soalen) - 1;
      if(fb == 1) mask <<= ((ngy-1)*soalen);
	}
    
#if defined (__GNUG__) && !defined (__INTEL_COMPILER)
    int* tmpspc __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)))  =&(tmpspc_all[veclen*16*tid]);
#else
    __declspec(align(QPHIX_LLC_CACHE_ALIGN)) int* tmpspc=&(tmpspc_all[veclen*16*tid]);
#endif    
    int *atmp = (int*)((((unsigned long long)tmpspc)+0x3F) & ~0x3F);
    int *offs=&atmp[0];
    int pkts_per_core = (npkts + n_cores - 1) / n_cores;
    
    // Assume thread indexing runs SMT-ID fastest
    // smtid + threads_per_core*cid
    int cid = tid/n_threads_per_core;
    int smtid = tid - n_threads_per_core * cid;
    
    int low_pkt = cid*pkts_per_core;
    int high_pkt = (cid+1)*pkts_per_core;
    if ( high_pkt > npkts ) high_pkt = npkts;
    
    // OK Each core can now work on its pkts:
    for(int pkt = low_pkt + smtid; pkt < high_pkt; pkt+=n_threads_per_core) {
      int coords[4];
      int tmp = pkt;
      for(int j = 0; j < 4; j++) {
	if(j != dir) {
	  int tmp1 = tmp / lens[j];
	  coords[j] = tmp - tmp1*lens[j];
	  tmp = tmp1;
	}
      }
      
      coords[dir] = (fb == 0 ? 0 : lens1[dir] - 1);
      int xblock = coords[0];
      int yblock = coords[1];
      int z = coords[2];
      int t = coords[3];
      int yi = yblock * ngy_pack;
      // yi is always going to be even for dir==0
      int xodd = (z + t + cb) & 1;
      
      if(dir == 0) {
	if(ngy == 1)
	  yi += (fb == 0 ? 1-xodd : xodd);
	else
	  mask = mask_xodd[1-xodd];
      }
      
      int pkt_next = pkt + n_threads_per_core < high_pkt ? pkt + n_threads_per_core : low_pkt + smtid;
      
      tmp = pkt_next;
      for(int j = 0; j < 4; j++) {
	if(j != dir) {
	  int tmp1 = tmp / lens[j];
	  coords[j] = tmp - tmp1*lens[j];
	  tmp = tmp1;
	}
      }
      
      coords[dir] = (fb == 0 ? 0 : lens1[dir] - 1);
      int xblock_next = coords[0];
      int yblock_next = coords[1];
      int z_next = coords[2];
      int t_next = coords[3];
      int yi_next = yblock_next * ngy_pack;
      // yi is always going to be even for dir==0
      int xodd_next = (z_next + t_next + cb) & 1;
      
      if(dir == 0 && ngy == 1) yi_next +=(fb == 0 ? 1-xodd_next : xodd_next);
      // printf("Pack %d %d %d %d\n", t, z, yi, xblock*soalen);
      
      // Now we have x,y,z coordinates, we need the base address of the spinor
      const FourSpinorBlock *xyBase = &psi[t*Pxyz + z*Pxy + yi*n_soa_x + xblock];
      // Offset to next spinor. T-values are the same.
      int off_next = (t_next-t)*Pxyz+(z_next-z)*Pxy+(yi_next-yi)*n_soa_x+(xblock_next-xblock);
      
      // Convert to prefetch distance in floats
      int si_offset = off_next*sizeof(FourSpinorBlock)/sizeof(FT);
      int hsprefdist = (pkt_next - pkt)*sizeof(FourSpinorBlock)/sizeof(FT)/2;
      
      // We are streaming out in sequence
      FT *outbuf = &res[12*pktsize*pkt];
      
      //printf("rank = %d, pkt = %d, outbuf=%p (%lld)\n", myRank, pkt, outbuf, outbuf-res);
      // OK: now we have xyBase, offs, and oubuf -- we should call the kernel.
      if(isPlus)
	face_proj_dir_plus<FT,veclen,soalen,compress>(xyBase, offs, si_offset, outbuf, hsprefdist, mask, dir*2+fb);
      else
	face_proj_dir_minus<FT,veclen,soalen,compress>(xyBase, offs, si_offset, outbuf, hsprefdist, mask, dir*2+fb);
    }
  }
  
  
  //  RECEIVE AND COMPLETE T-FACE FROM FORWARD
  template<typename FT, int veclen,int soalen, bool compress>
    void ClovDslash<FT,veclen,soalen,compress>::completeFaceDir(int tid,
								const FT* psi,
								FourSpinorBlock* res,
								const SU3MatrixBlock* u,
								const CloverBlock* invclov,
								const double beta, 
								int cb, int dir, int fb, int isPlus)
    {
      // This is the total number of veclen in the face. 
      // Guaranteed to be good, since s->Nxh()*s->Ny() is a multiple
      // of VECLEN
      int Nxh = s->Nxh();
      int Ny = s->Ny();
      int Nz = s->Nz();
      int Nt = s->Nt();
      int ngy = s->nGY();
      int n_cores = s->getNumCores();
      int Pxyz = s->getPxyz();
      int Pxy = s->getPxy();
      
      int ngy_pack = ngy;
      // If packing x face and ngy is 1, we need to pick alternate rows in Y, so use ngy = 2
      if(dir == 0 && ngy_pack == 1) ngy_pack = 2;
      int n_soa_x = Nxh / soalen;
      int n_blocks_y = Ny/ngy_pack;
      int lens1[4] = {n_soa_x, n_blocks_y, Nz, Nt};
      int lens[4] = {n_soa_x, n_blocks_y, Nz, Nt};
      
      lens[dir] = 1;
      
      int npkts = lens[0] * lens[1] * lens[2] * lens[3];
      int pktsize = veclen;
      unsigned int mask = -1;
      unsigned int mask_xodd[2] = {0, 0};
      
      if(dir == 0) {
	pktsize = ngy_pack / 2;
	if(ngy == 1) 
	  mask = (fb == 0 ? 1 : (1 << (soalen-1)));
	else {
	  for(int i = 0; i < veclen; i+=2*soalen)	mask_xodd[0] |= ((fb == 0 ? 1 : (1 << (soalen - 1))) << i);
	  mask_xodd[1] = mask_xodd[0] << soalen;
	  if(fb == 1) {
	    unsigned int tmp = mask_xodd[0];
	    mask_xodd[0] = mask_xodd[1];
	    mask_xodd[1] = tmp;
	  }
	}
      }
      else if(dir == 1) {
	pktsize = soalen;
	mask = (1 << soalen) - 1;
	if(fb == 1) mask <<= ((ngy-1)*soalen);
      }
      
#if defined (__GNUG__) && !defined (__INTEL_COMPILER)
      int* tmpspc __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)))  =&(tmpspc_all[veclen*16*tid]);
#else
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) int* tmpspc=&(tmpspc_all[veclen*16*tid]);
#endif    
      int *atmp = (int*)((((unsigned long long)tmpspc)+0x3F) & ~0x3F);
      int *offs = &atmp[0];
      int *gOffs = &atmp[veclen*13];
      
      int pkts_per_core = (npkts + n_cores - 1) / n_cores;
      
      
      // Assume thread indexing runs SMT-ID fastest
      // smtid + threads_per_core*cid
      int cid = tid/n_threads_per_core;
      int smtid = tid - n_threads_per_core * cid;
      
      int low_pkt = cid*pkts_per_core;
      int high_pkt = (cid+1)*pkts_per_core;
      if ( high_pkt > npkts ) high_pkt = npkts;
      
      // OK Each core can now work on its pkts:
      for(int pkt = low_pkt + smtid; pkt < high_pkt; pkt+=n_threads_per_core) {
	int coords[4];
	int tmp = pkt;
	for(int j = 0; j < 4; j++) {
	  if(j != dir) {
	    int tmp1 = tmp / lens[j];
	    coords[j] = tmp - tmp1*lens[j];
	    tmp = tmp1;
	  }
	}
	
	coords[dir] = (fb == 0 ? 0 : lens1[dir] - 1);
	int xblock = coords[0];
	int yblock = coords[1];
	int z = coords[2];
	int t = coords[3];
	int yi = yblock * ngy_pack;
	// yi is always going to be even for dir==0
	int xodd = (z + t + cb) & 1;
	
	if(dir == 0) {
	  if(ngy == 1) 
	    yi += (fb == 0 ? xodd : 1-xodd);
	  else
	    mask = mask_xodd[xodd];
	}
	
	int pkt_next = pkt + n_threads_per_core < high_pkt ? pkt + n_threads_per_core : low_pkt + smtid;
	
	tmp = pkt_next;
	for(int j = 0; j < 4; j++) {
	  if(j != dir) {
	    int tmp1 = tmp / lens[j];
	    coords[j] = tmp - tmp1*lens[j];
	    tmp = tmp1;
	  }
	}
	
	coords[dir] = (fb == 0 ? 0 : lens1[dir] - 1);
	int xblock_next = coords[0];
	int yblock_next = coords[1];
	int z_next = coords[2];
	int t_next = coords[3];
	int yi_next = yblock_next * ngy_pack;
	// yi is always going to be even for dir==0
	int xodd_next = (z_next + t_next + cb) & 1;
	if(dir == 0 && ngy == 1) yi_next += (fb == 0 ? xodd_next : 1-xodd_next);
	// printf("Pack %d %d %d %d\n", t, z, yi, xblock*soalen);
	
	// Now we have x,y,z coordinates, we need the base address of the spinor
	// Now we have x,y,z coordinates, we need the base address of the spinor
	FourSpinorBlock *oBase = &res[t*Pxyz + z*Pxy + yi*n_soa_x + xblock];
	// Offset to next spinor. T-values are the same.
	int off_next = (t_next-t)*Pxyz+(z_next-z)*Pxy+(yi_next-yi)*n_soa_x+(xblock_next-xblock);
	
	int soprefdist = off_next*sizeof(FourSpinorBlock)/sizeof(FT);
	int hsprefdist = (pkt_next - pkt)*sizeof(FourSpinorBlock)/sizeof(FT)/2;
	int goff_next  = ((t_next-t)*Pxyz+(z_next-z)*Pxy+(yi_next-yi)*n_soa_x)/ngy + (xblock_next-xblock);
	int gprefdist = goff_next*sizeof(SU3MatrixBlock)/sizeof(FT);
	
	const SU3MatrixBlock *gBase = &u[(t*Pxyz+z*Pxy+yi*n_soa_x)/ngy+xblock];
	
	const CloverBlock *clBase = &invclov[(t*Pxyz+z*Pxy+yi*n_soa_x)/ngy+xblock];
	const int clov_line_in_floats = sizeof(CloverBlock)/sizeof(FT); // One gauge scanline, in floats
	int clprefdist =goff_next * clov_line_in_floats;
	
	// We are streaming out in sequence
	const FT *inbuf = &psi[12*pktsize*pkt];
	// OK: now we have xyBase, offs, and oubuf -- we should call the kernel.
	FT beta_T = rep<FT,double>(beta);
	
	if(isPlus)
	  face_clov_finish_dir_plus<FT,veclen,soalen,compress>(inbuf, gBase, oBase, clBase, gOffs, offs, hsprefdist, gprefdist, soprefdist, clprefdist, beta_T, mask, dir*2+fb );
	else
	  face_clov_finish_dir_minus<FT,veclen,soalen,compress>(inbuf, gBase, oBase, clBase, gOffs, offs, hsprefdist, gprefdist, soprefdist, clprefdist, beta_T, mask, dir*2+fb );
      }
    } // Function
  
  // Accumulate received back T face (made by packTFaceForwPlus) 
  // Recons_add ( 1 + gamma_T )
  template<typename FT, int veclen,int soalen, bool compress>
    void ClovDslash<FT,veclen,soalen,compress>::completeFaceDirAChiMBDPsi(int tid,
									  const FT* psi,
									  FourSpinorBlock* res,
									  const SU3MatrixBlock* u,
									  const double beta, 
									  int cb, int dir, int fb, int isPlus)
    {
    // This is the total number of veclen in the face. 
    // Guaranteed to be good, since s->Nxh()*s->Ny() is a multiple
    // of VECLEN
    int Nxh = s->Nxh();
    int Ny = s->Ny();
    int Nz = s->Nz();
    int Nt = s->Nt();
    int ngy = s->nGY();
    int n_cores = s->getNumCores();
    int Pxyz = s->getPxyz();
    int Pxy = s->getPxy();

	int ngy_pack = ngy;
	// If packing x face and ngy is 1, we need to pick alternate rows in Y, so use ngy = 2
	if(dir == 0 && ngy_pack == 1) ngy_pack = 2;
    int n_soa_x = Nxh / soalen;
	int n_blocks_y = Ny/ngy_pack;
	int lens1[4] = {n_soa_x, n_blocks_y, Nz, Nt};
	int lens[4] = {n_soa_x, n_blocks_y, Nz, Nt};

	lens[dir] = 1;

	int npkts = lens[0] * lens[1] * lens[2] * lens[3];
	int pktsize = veclen;
	unsigned int mask = -1;
	unsigned int mask_xodd[2] = {0, 0};

	if(dir == 0) {
		pktsize = ngy_pack / 2;
		if(ngy == 1) 
			mask = (fb == 0 ? 1 : (1 << (soalen-1)));
		else {
			for(int i = 0; i < veclen; i+=2*soalen)	mask_xodd[0] |= ((fb == 0 ? 1 : (1 << (soalen - 1))) << i);
			mask_xodd[1] = mask_xodd[0] << soalen;
			if(fb == 1) {
				unsigned int tmp = mask_xodd[0];
				mask_xodd[0] = mask_xodd[1];
				mask_xodd[1] = tmp;
			}
		}
	}
	else if(dir == 1) {
		pktsize = soalen;
		mask = (1 << soalen) - 1;
		if(fb == 1) mask <<= ((ngy-1)*soalen);
	}

      
#if defined (__GNUG__) && !defined (__INTEL_COMPILER)

      int* tmpspc __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)))  =&(tmpspc_all[veclen*16*tid]);
#else
      __declspec(align(QPHIX_LLC_CACHE_ALIGN)) int* tmpspc=&(tmpspc_all[veclen*16*tid]);
#endif    
      int *atmp = (int*)((((unsigned long long)tmpspc)+0x3F) & ~0x3F);
      int *offs = &atmp[0];
      int *gOffs = &atmp[veclen*13];

      int pkts_per_core = (npkts + n_cores - 1) / n_cores;
      
      
      // Assume thread indexing runs SMT-ID fastest
      // smtid + threads_per_core*cid
      int cid = tid/n_threads_per_core;
      int smtid = tid - n_threads_per_core * cid;
      
      int low_pkt = cid*pkts_per_core;
      int high_pkt = (cid+1)*pkts_per_core;
      if ( high_pkt > npkts ) high_pkt = npkts;
      
      // OK Each core can now work on its pkts:
      for(int pkt = low_pkt + smtid; pkt < high_pkt; pkt+=n_threads_per_core) {
	int coords[4];
	int tmp = pkt;
	for(int j = 0; j < 4; j++) {
	  if(j != dir) {
	    int tmp1 = tmp / lens[j];
	    coords[j] = tmp - tmp1*lens[j];
	    tmp = tmp1;
	  }
	}
	
	coords[dir] = (fb == 0 ? 0 : lens1[dir] - 1);
	int xblock = coords[0];
	int yblock = coords[1];
	int z = coords[2];
	int t = coords[3];
	int yi = yblock * ngy_pack;
	// yi is always going to be even for dir==0
	int xodd = (z + t + cb) & 1;
	
	if(dir == 0) {
	  if(ngy == 1) 
	    yi += (fb == 0 ? xodd : 1-xodd);
	  else
	    mask = mask_xodd[xodd];
	}
	
	int pkt_next = pkt + n_threads_per_core < high_pkt ? pkt + n_threads_per_core : low_pkt + smtid;
	
	tmp = pkt_next;
	for(int j = 0; j < 4; j++) {
	  if(j != dir) {
	    int tmp1 = tmp / lens[j];
	    coords[j] = tmp - tmp1*lens[j];
	    tmp = tmp1;
	  }
	}
	
	coords[dir] = (fb == 0 ? 0 : lens1[dir] - 1);
	int xblock_next = coords[0];
	int yblock_next = coords[1];
	int z_next = coords[2];
	int t_next = coords[3];
	int yi_next = yblock_next * ngy_pack;
	// yi is always going to be even for dir==0
	int xodd_next = (z_next + t_next + cb) & 1;
	if(dir == 0 && ngy == 1) yi_next += (fb == 0 ? xodd_next : 1-xodd_next);
	// printf("Pack %d %d %d %d\n", t, z, yi, xblock*soalen);
	
	// Now we have x,y,z coordinates, we need the base address of the spinor
	// Now we have x,y,z coordinates, we need the base address of the spinor
	FourSpinorBlock *oBase = &res[t*Pxyz + z*Pxy + yi*n_soa_x + xblock];
	// Offset to next spinor. T-values are the same.
	int off_next = (t_next-t)*Pxyz+(z_next-z)*Pxy+(yi_next-yi)*n_soa_x+(xblock_next-xblock);
	
	int soprefdist = off_next*sizeof(FourSpinorBlock)/sizeof(FT);
	int hsprefdist = (pkt_next - pkt)*sizeof(FourSpinorBlock)/sizeof(FT)/2;
	int goff_next  = ((t_next-t)*Pxyz+(z_next-z)*Pxy+(yi_next-yi)*n_soa_x)/ngy + (xblock_next-xblock);
	int gprefdist = goff_next*sizeof(SU3MatrixBlock)/sizeof(FT);
	
	const SU3MatrixBlock *gBase = &u[(t*Pxyz+z*Pxy+yi*n_soa_x)/ngy+xblock];
	// We are streaming out in sequence
	const FT *inbuf = &psi[12*pktsize*pkt];
	// OK: now we have xyBase, offs, and oubuf -- we should call the kernel.
	FT beta_T = rep<FT,double>(beta);
	
	if(isPlus)
	  face_finish_dir_plus<FT,veclen,soalen,compress>(inbuf, gBase, oBase, gOffs, offs, hsprefdist, gprefdist, soprefdist, beta_T, mask, dir*2+fb );
	else
	  face_finish_dir_minus<FT,veclen,soalen,compress>(inbuf, gBase, oBase, gOffs, offs, hsprefdist, gprefdist, soprefdist, beta_T, mask, dir*2+fb );
      }
    } // Function 
  
} // Namespace

#endif
