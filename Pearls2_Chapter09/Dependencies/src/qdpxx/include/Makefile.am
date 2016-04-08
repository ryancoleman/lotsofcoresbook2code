PETE_HDRS =  PETE/Combiners.h \
             PETE/Functors.h \
             PETE/PETE.h \
             PETE/TreeNodes.h \
             PETE/CreateLeaf.h \
             PETE/ForEach.h \
             PETE/ForEachInOrder.h \
             PETE/OperatorTags.h \
             PETE/Scalar.h \
             PETE/TypeComputations.h

# Generic optimization headers
genericdir = scalarsite_generic

GENERIC_HDRS = qdp_scalarsite_generic.h \
	$(genericdir)/generic_addvec.h \
	$(genericdir)/generic_adj_mat_vec.h \
	$(genericdir)/generic_mat_vec.h \
	$(genericdir)/generic_mult_aa.h \
	$(genericdir)/generic_mult_an.h \
	$(genericdir)/generic_mult_na.h \
	$(genericdir)/generic_mult_nn.h \
	$(genericdir)/generic_blas_local_sumsq.h \
	$(genericdir)/generic_blas_vadd.h  \
	$(genericdir)/generic_blas_vaxmby3.h \
	$(genericdir)/generic_blas_vaxmby3_norm.h \
	$(genericdir)/generic_blas_vaxmy3.h \
	$(genericdir)/generic_blas_vaxpby3.h \
	$(genericdir)/generic_blas_vaxpy3.h \
	$(genericdir)/generic_blas_vaxpy3_norm.h \
	$(genericdir)/generic_blas_vscal.h \
	$(genericdir)/generic_blas_vsub.h \
	$(genericdir)/generic_blas_vcscal.h \
	$(genericdir)/generic_blas_vcaxpy3.h \
	$(genericdir)/generic_blas_vcaxmy3.h \
	$(genericdir)/generic_blas_vcaxpby3.h \
	$(genericdir)/generic_blas_vcaxmby3.h \
	$(genericdir)/generic_blas_local_vcdot.h \
	$(genericdir)/generic_blas_local_vcdot_real.h \
	$(genericdir)/qdp_scalarsite_generic_linalg.h \
	$(genericdir)/qdp_scalarsite_generic_blas.h \
	$(genericdir)/qdp_scalarsite_generic_cblas.h \
	$(genericdir)/qdp_scalarsite_generic_blas_g5.h \
	$(genericdir)/generic_blas_vadd3_g5.h \
	$(genericdir)/generic_blas_vscal_g5.h \
	$(genericdir)/generic_blas_vaxpy3_g5.h \
	$(genericdir)/generic_blas_vaypx3_g5.h \
	$(genericdir)/generic_blas_vaxpby3_g5.h \
	$(genericdir)/generic_blas_g5.h \
	$(genericdir)/generic_mv_switchbox.h \
	$(genericdir)/generic_spin_proj_inlines.h \
	$(genericdir)/generic_spin_recon_inlines.h \
	$(genericdir)/generic_spin_proj.h \
	$(genericdir)/generic_spin_recon.h \
	$(genericdir)/generic_fused_spin_proj.h \
	$(genericdir)/generic_fused_spin_recon.h \
	$(genericdir)/qdp_generic_spin_project_evaluates.h \
	$(genericdir)/qdp_generic_fused_spin_proj_evaluates.h \
	$(genericdir)/qdp_generic_fused_spin_recon_evaluates.h \
	$(genericdir)/generic_spin_aggregate.h \
	$(genericdir)/qdp_scalarsite_generic_blas_wrapper.h \
	$(genericdir)/qdp_scalarsite_generic_cblas_wrapper.h \
	$(genericdir)/qdp_scalarsite_generic_blas_g5_wrapper.h \
	$(genericdir)/qdp_generic_spin_project_evaluates_wrapper.h \
	$(genericdir)/qdp_generic_fused_spin_proj_evaluates_wrapper.h \
	$(genericdir)/qdp_generic_fused_spin_recon_evaluates_wrapper.h \
	$(genericdir)/qdp_scalarsite_generic_linalg_wrapper.h

# If using SSE, make sure those files are installed
ssesitedir = scalarsite_sse

SSESITE_HDRS = qdp_scalarsite_sse.h \
	qdp_sse_intrin.h \
	$(ssesitedir)/qdp_scalarsite_sse_linalg.h \
	$(ssesitedir)/qdp_scalarsite_sse_blas.h \
	$(ssesitedir)/qdp_scalarsite_sse_blas_g5.h \
	$(ssesitedir)/qdp_scalarsite_sse_blas_g5_includes.h \
	$(ssesitedir)/sse_blas_vadd3_g5.h \
	$(ssesitedir)/sse_blas_vscal3_g5.h \
	$(ssesitedir)/sse_blas_vaxpy3_g5.h \
	$(ssesitedir)/sse_blas_vaypx3_g5.h \
	$(ssesitedir)/sse_blas_vaxpby3_g5.h \
	$(ssesitedir)/sse_spin_aggregate.h \
	$(ssesitedir)/sse_spin_proj_inlines.h \
	$(ssesitedir)/sse_spin_recon_inlines.h \
	$(ssesitedir)/sse_spin_proj.h \
	$(ssesitedir)/sse_spin_recon.h \
	$(ssesitedir)/sse_fused_spin_proj.h \
	$(ssesitedir)/sse_fused_spin_recon.h \
	$(ssesitedir)/qdp_sse_spin_evaluates.h \
	$(ssesitedir)/qdp_sse_fused_spin_proj_evaluates.h \
	$(ssesitedir)/qdp_sse_fused_spin_recon_evaluates.h \
	$(ssesitedir)/qdp_scalarsite_sse_blas_double.h \
	$(ssesitedir)/sse_blas_local_sumsq_double.h \
	$(ssesitedir)/sse_blas_local_vcdot_double.h \
	$(ssesitedir)/sse_blas_local_vcdot_real_double.h \
	$(ssesitedir)/sse_blas_vaxmbyz4_double.h \
	$(ssesitedir)/sse_blas_vaxpbyz4_double.h \
	$(ssesitedir)/sse_blas_vaxmyz4_double.h \
	$(ssesitedir)/sse_blas_vaxpy4_double.h \
	$(ssesitedir)/sse_blas_vaypx4_double.h \
	$(ssesitedir)/sse_blas_vscal4_double.h \
	$(ssesitedir)/sse_linalg_mm_su3_double.h \
	$(ssesitedir)/qdp_scalarsite_sse_linalg_double.h \
	$(ssesitedir)/sse_dcomplex_mult_macros.h \
	$(ssesitedir)/sse_prefetch.h \
	$(ssesitedir)/qdp_scalarsite_sse_blas_wrapper.h \
	$(ssesitedir)/qdp_scalarsite_sse_blas_double_wrapper.h \
	$(ssesitedir)/qdp_sse_spin_evaluates_wrapper.h \
	$(ssesitedir)/qdp_sse_fused_spin_proj_evaluates_wrapper.h \
	$(ssesitedir)/qdp_sse_fused_spin_recon_evaluates_wrapper.h

# Scalarvectsite SSE extensions
ssevecdir = scalarvecsite_sse

SSEVEC_HDRS = qdp_scalarvecsite_sse.h \
	$(ssevecdir)/ssevec_mult_nn.h \
	$(ssevecdir)/qdp_scalarvecsite_sse_linalg.h \
	$(ssevecdir)/qdp_scalarvecsite_sse_blas.h


# QCDOC NONSPECIFIC BAGEL stuff
bageldir = scalarsite_bagel_qdp
BAGEL_HDRS = qdp_scalarsite_bagel_qdp.h \
	$(bageldir)/qdp_scalarsite_bagel_qdp_linalg.h \
	$(bageldir)/qdp_scalarsite_bagel_qdp_blas.h \
	$(bageldir)/qdp_scalarsite_bagel_qdp_blas_g5.h

# HEADERS for the memory allocator
MEMORY_HDRS = qdp_allocator.h \
	      qdp_singleton.h \
	      qdp_default_allocator.h

# All the include files - avoid flattening of dirs by using nobase
nobase_include_HEADERS = \
		qdp_config.h \
		qdp_forward.h \
		qdp_globalfuncs.h \
		qdp_globalfuncs_subtype.h \
		qdp_inner.h \
		qdp_init.h \
		qdp_io.h \
		qdp_stdio.h \
		qdp_layout.h \
		qdp_map.h \
		qdp_multi.h \
		qdp_arrays.h \
		qdp_newops.h \
		qdp_optops.h \
		qdp_outer.h \
		qdp_outersubtype.h \
		qdp_params.h \
		qdp_precision.h \
		qdp_primcolormat.h \
		qdp_primcolorvec.h \
		qdp_primgamma.h \
		qdp_primitive.h \
		qdp_primmatrix.h \
		qdp_primscalar.h \
		qdp_primseed.h \
		qdp_primspinmat.h \
		qdp_primspinvec.h \
		qdp_primvector.h \
		qdp_profile.h \
                qdp_stopwatch.h \
		qdp_flopcount.h \
		qdp_iogauge.h \
		qdp_crc32.h \
		qdp_byteorder.h \
		qdp_util.h \
		qdp_xmlio.h \
		qdp_qdpexpr.h \
		qdp.h \
		qdp_qdpio.h \
		QDPOperators.h \
		qdp_qdpsubtype.h \
		qdp_qdptype.h \
		qdp_random.h \
		qdp_rannyu.h \
		qdp_reality.h \
		qdp_simpleword.h \
		qdp_specializations.h \
		qdp_subset.h \
		qdp_traits.h \
		qdp_word.h \
		qdp_dispatch.h \
	        qdp_scalar_specific.h \
	        qdp_parscalar_specific.h \
	        qdp_scalarvec_specific.h \
	        qdp_parscalarvec_specific.h \
	 	qdp_defs.h \
		qdp_scalarsite_defs.h \
		qdp_scalarvecsite_defs.h \
		qdp_strnlen.h \
		qdp_db.h \
		qdp_db_imp.h \
		qdp_db_stub.h \
		qdp_map_obj.h \
		qdp_map_obj_null.h \
		qdp_map_obj_memory.h \
		qdp_map_obj_disk.h \
		qdp_map_obj_disk_multiple.h \
		qdp_disk_map_slice.h \
		qdp_hdf5.h \
		qdp_threadbind.h \
                $(PETE_HDRS) \
                $(GENERIC_HDRS) \
		$(MEMORY_HDRS) \
		$(SSESITE_HDRS) \
		$(SSEVEC_HDRS) \
	        $(BAGEL_HDRS)

nodist_include_HEADERS = qdp_config_internal.h


