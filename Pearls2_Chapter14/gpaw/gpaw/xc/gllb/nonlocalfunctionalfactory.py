class NonLocalFunctionalFactory:
    """Factory class.

    NonLocalFunctionalFactory is used by xc_functional.py, when the
    functional starts with words GLLB.
    
    It contains a method called get_functional_by_name, which takes
    the xc-name for non-local functional and returns the corresponding
    NonLocalFunctional object. 

    * GLLB
    * GLLBC (GLLB with screening part from PBE + PBE Correlation)
    * GLLBSC (GLLB with screening part from PBE_SOL + PBE Correlation)

    
    * GLLBNORESP (Just GLLB Screening)
    * GLLBLDA (A test functional, which is just LDA but via
               NonLocalFunctional framework)
    * GLLBGGA (A test functional, which is just GGA but via
               NonLocalFunctional framework)
    """

    def get_functional_by_name(self, name):
        print("Functional name", name)

        from gpaw.xc.gllb.nonlocalfunctional import NonLocalFunctional
        functional = NonLocalFunctional(name)
        
        if name == 'GLLB':
            # Functional GLLB
            # Contains screening part from GGA functional
            # And response part based on simple square root expection
            # of orbital energy differences.
            from gpaw.xc.gllb.c_gllbscr import C_GLLBScr
            from gpaw.xc.gllb.c_response import C_Response
            C_Response(functional, 1.0,
                       C_GLLBScr(functional, 1.0).get_coefficient_calculator())
            return functional
        elif name == 'GLLBM':
            from gpaw.xc.gllb.c_gllbscr import C_GLLBScr
            from gpaw.xc.gllb.c_response import C_Response
            from gpaw.xc.gllb.c_xc import C_XC
            C_Response(functional, 1.0, C_GLLBScr(
                functional, 1.0, metallic=True).get_coefficient_calculator())
            return functional
        elif name == 'GLLBSC':
            from gpaw.xc.gllb.c_gllbscr import C_GLLBScr
            from gpaw.xc.gllb.c_response import C_Response
            from gpaw.xc.gllb.c_xc import C_XC
            C_Response(functional, 1.0, C_GLLBScr(
                functional,
                1.0,
                'GGA_X_PBE_SOL').get_coefficient_calculator())
            C_XC(functional, 1.0, 'GGA_C_PBE_SOL')
            return functional
        elif name == 'GLLBSCM':
            from gpaw.xc.gllb.c_gllbscr import C_GLLBScr
            from gpaw.xc.gllb.c_response import C_Response
            from gpaw.xc.gllb.c_xc import C_XC
            C_Response(functional, 1.0, C_GLLBScr(
                functional,
                1.0,
                'GGA_X_PBE_SOL', metallic=True).get_coefficient_calculator())
            C_XC(functional, 1.0, 'GGA_C_PBE_SOL')
            return functional
        
        elif name == 'GLLBC':
            from gpaw.xc.gllb.c_gllbscr import C_GLLBScr
            from gpaw.xc.gllb.c_response import C_Response
            from gpaw.xc.gllb.c_xc import C_XC
            C_Response(functional, 1.0,
                       C_GLLBScr(functional, 1.0,'X_PBE-None').get_coefficient_calculator())
            C_XC(functional, 1.0, 'None-C_PBE')
            return functional
        elif name == 'GLLBCP86':
            from gpaw.xc.gllb.c_gllbscr import C_GLLBScr
            from gpaw.xc.gllb.c_response import C_Response
            from gpaw.xc.gllb.c_xc import C_XC
            C_Response(functional, 1.0,
                       C_GLLBScr(functional, 1.0).get_coefficient_calculator())
            C_XC(functional, 1.0, 'None-C_P86')
            return functional
        elif name == 'GLLBLDA':
            from gpaw.xc.gllb.c_xc import C_XC
            C_XC(functional, 1.0,'LDA')
            return functional
        elif name == 'GLLBGGA':
            from gpaw.xc.gllb.c_xc import C_XC
            C_XC(functional, 1.0,'PBE')
            return functional
        elif name == 'GLLBSLATER':
            from gpaw.xc.gllb.c_slater import C_Slater
            C_Slater(functional, 1.0)
            return functional
        elif name == 'GLLBNORESP':
            from gpaw.xc.gllb.c_gllbscr import C_GLLBScr
            C_GLLBScr(functional, 1.0)
            return functional
        elif name == 'KLI':
            raise RuntimeError('KLI functional not implemented')
            from gpaw.xc.gllb.c_slater import C_Slater
            from gpaw.xc.gllb.c_response import C_Response
            C_Response(functional, 1.0,
                       C_Slater(functional, 1.0).get_coefficient_calculator())
            return functional
        else:
            raise RuntimeError('Unkown NonLocal density functional: ' + name)

