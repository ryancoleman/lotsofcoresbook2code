// -*- C++ -*-

#ifndef __clover_fermact_params_w_h__
#define __clover_fermact_params_w_h__

// These are parameters for setting up a clover term
// I nicked them from chroma, and cut out all the XML 
// IO stuff

//! Parameters for anisotropy
struct AnisoParam_t
{
  AnisoParam_t()
  {
    anisoP = false; // Isotropic by default
    t_dir = 3;      // direction 3 is time (or anisotropic direction anyhow)
    xi_0 = Real(1); // No anisotropy
    nu = Real(1);   // No anisotropy
  }
  ~AnisoParam_t() {}
  
  bool       anisoP;
  int        t_dir;
  Real       xi_0;
  Real       nu;
};


//! Params for clover ferm acts
/*! \ingroup fermacts */
struct CloverFermActParams
{
  CloverFermActParams()
  {
    Mass=Real(0); // Default zero mass
    clovCoeffR=Real(0);
    clovCoeffT=Real(0);
    u0=Real(1);
  }
  
  Real Mass;
  Real clovCoeffR;
  Real clovCoeffT;
  Real u0;
  
  // Optional Anisotropy
  AnisoParam_t anisoParam;
};

#endif
