!
! Copyright (c) Stanford University, The Regents of the University of
!               California, and others.
!
! All Rights Reserved.
!
! See Copyright-SimVascular.txt for additional details.
!
! Permission is hereby granted, free of charge, to any person obtaining
! a copy of this software and associated documentation files (the
! "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish,
! distribute, sublicense, and/or sell copies of the Software, and to
! permit persons to whom the Software is furnished to do so, subject
! to the following conditions:
!
! The above copyright notice and this permission notice shall be included
! in all copies or substantial portions of the Software.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!--------------------------------------------------------------------
!
!     Parameters used for Aliev-Panfilov Ventricular Myocyte Model.
!
!--------------------------------------------------------------------

!     Maximum possible nsd
      INTEGER(KIND=IKIND), PARAMETER :: maxnsd = 3
!--------------------------------------------------------------------
!     Maximum number of properties
      INTEGER(KIND=IKIND), PARAMETER :: maxNProp = 20
!--------------------------------------------------------------------
!     Gauss point coordinates, upto 5 points
      REAL(KIND=RKIND), PARAMETER :: gW(5,5) =
     2   RESHAPE( (/2._RKIND, 0._RKIND, 0._RKIND, 0._RKIND, 0._RKIND,
     3   1._RKIND, 1._RKIND, 0._RKIND, 0._RKIND, 0._RKIND,
     4   0.5555555555555556_RKIND, 0.8888888888888889_RKIND,
     4   0.5555555555555556_RKIND, 0._RKIND, 0._RKIND,
     5   0.3478548451374538_RKIND, 0.6521451548625462_RKIND,
     5   0.6521451548625462_RKIND, 0.3478548451374538_RKIND, 0._RKIND,
     6   0.2369268850561890_RKIND, 0.4786286704993665_RKIND,
     6   0.5688888888888889_RKIND, 0.4786286704993665_RKIND,
     6   0.2369268850561890_RKIND/) , (/5, 5/) )
!--------------------------------------------------------------------
!     Gauss point weights, upto 5 points
      REAL(KIND=RKIND), PARAMETER :: gXi(5,5) =
     2   RESHAPE( (/0._RKIND, 0._RKIND, 0._RKIND, 0._RKIND, 0._RKIND,
     3   -0.57735026918962584_RKIND, 0.57735026918962584_RKIND,
     3   0._RKIND, 0._RKIND, 0._RKIND,
     4   -0.7745966692414834_RKIND, 0._RKIND, 0.7745966692414834_RKIND,
     4   0._RKIND, 0._RKIND,
     5   -0.86113631159405257_RKIND, -0.33998104358485631_RKIND,
     5   0.33998104358485631_RKIND, 0.86113631159405257_RKIND, 0._RKIND,
     6   -0.90617984593866396_RKIND, -0.53846931010568311_RKIND,
     6   0._RKIND, 0.53846931010568311_RKIND,
     6   0.90617984593866396_RKIND/) , (/5, 5/) )
!--------------------------------------------------------------------
      CHARACTER, PARAMETER :: delimiter = "/"
!--------------------------------------------------------------------
!     Possible physical properties. Current maxNPror is 20.
!     When adding more properties, remember to increase maxNProp
!     Density of fluid, density of solid, viscosity of solid, elasticity
!     modulus, Poisson's ratio, conductivity, internal force (X, Y, Z),
!     permeability, stabilization coeff. for backflow divergence,
!     external source, damping, shell thickness, stabilization coeffs.
!     for USTRUCT (mom., cont.)
      INTEGER(KIND=IKIND), PARAMETER :: prop_NA = 0, fluid_density = 1,
     2   solid_density = 2, solid_viscosity = 3, elasticity_modulus = 4,
     3   poisson_ratio = 5, conductivity = 6, f_x = 7, f_y = 8, f_z = 9,
     4   permeability = 10, backflow_stab = 11, source_term = 12,
     5   damping = 13, shell_thickness = 14, ctau_M = 15, ctau_C = 16
!--------------------------------------------------------------------
!     Types of accepted elements
!     Linear (Q1-1D), triangle (P1-2D), tetrahedral (P1-3D), bilinear
!     (Q1-2D), quad (Q2-1D), biquad (Q2-2D), brick (Q1-3D),
!     general NURBS (1-3D), wedge, point (0D), quadratic triangle
!     (P2-2D), quadratic tetrahedron (P2-3D)
      INTEGER(KIND=IKIND), PARAMETER :: eType_NA = 100, eType_LIN = 101,
     2   eType_TRI = 102, eType_TET = 103, eType_BIL = 104,
     3   eType_QUD = 105, eType_BIQ = 106, eType_BRK = 107,
     4   eType_NRB = 108, eType_WDG = 109, eType_PNT = 110,
     5   eType_QTR = 111, eType_QTE = 112
!--------------------------------------------------------------------
!     Types of equations that are included in this solver
!     Fluid equation (Navier-Stokes), nonlinear structure (pure d), heat
!     equation, linear elasticity, heat in fluid (advection-diffusion),
!     fluid-structure-interaction, mesh motion, Shell mechanics,
!     Coupled-Momentum-Method, Cardiac Electro-Physiology,
!     Nonlinear structure (v-p), Stokes equations
      INTEGER(KIND=IKIND), PARAMETER :: phys_NA = 200, phys_fluid = 201,
     2   phys_struct = 202, phys_heatS = 203, phys_lElas = 204,
     3   phys_heatF = 205, phys_FSI = 206, phys_mesh = 207,
     4   phys_shell = 208, phys_CMM = 209, phys_CEP = 210,
     5   phys_ustruct = 211, phys_stokes = 212
!--------------------------------------------------------------------
!     Differenty type of coupling for cplBC
!     Not-available, implicit, semi-implicit, and explicit
      INTEGER(KIND=IKIND), PARAMETER :: cplBC_NA = 400, cplBC_I = 401,
     2   cplBC_SI = 402, cplBC_E = 403
!--------------------------------------------------------------------
!     cplBC type of coupling to between 3D and OD-LPN models:
!     Dirichlet type coupling, Neumann type coupling
      INTEGER(KIND=IKIND), PARAMETER :: cplBC_Dir = 66112,
     2   cplBC_Neu = 66113
!--------------------------------------------------------------------
!     boundary conditions types. Items of this list can be combined
!     BCs from imposing perspective can be Neu/Dir/per
!     BCs time dependence can be std/ustd/cpl/gen/res
!     BCs spatial distribution can be flat/para/ud
!     Beside these nodes at the boundary perimeter can be set to
!     zero and flux through surface can be assigned instead of nodal
!     values.
!     Dirichlet, Neumann, Traction, CMM, Robin, steady, unsteady,
!     coupled, general (combination of ud/ustd), resistance, imposed
!     flux, zero out perimeter, impose BC on the integral of state
!     variable or D (instead of Y), flat profile, parabolic profile,
!     user defined profile, backflow stabilization, BCs for shells
!     (fixed, hinged, free, symmetric), undeforming Neu, RCR-Neu
      INTEGER(KIND=IKIND), PARAMETER :: bType_Dir = 0, bType_Neu = 1,
     2   bType_trac = 2, bType_CMM = 3, bType_Robin = 4, bType_std = 5,
     3   bType_ustd = 6, bType_cpl = 7, bType_gen = 8, bType_res = 9,
     4   bType_flx = 10, bType_zp = 11, bType_impD = 12, bType_flat =13,
     5   bType_para = 14, bType_ud = 15, bType_bfs = 16, bType_fix = 17,
     6   bType_hing = 18, bType_free = 19, bType_symm = 20,
     7   bType_undefNeu = 21, bType_RCR = 22
!--------------------------------------------------------------------
!     Body force types: volumetric (default), traction, Neumann
!     (pressure based), time dependence (steady, unsteady, spatially
!     varying, general)
      INTEGER(KIND=IKIND), PARAMETER :: bfType_vol = 0, bfType_trac = 1,
     2   bfType_Neu = 2, bfType_std = 3, bfType_ustd = 4,
     3   bfType_spl = 5, bfType_gen = 6
!--------------------------------------------------------------------
!     Possible senarios for the output, followed by the possible outputs
!     Undefined output, extract it from A, extract it from Y, extract it
!     from D, calculate WSS, calculate vorticity, energy flux, heat
!     flux, absolute velocity (for FSI)
      INTEGER(KIND=IKIND), PARAMETER :: outGrp_NA = 500, outGrp_A = 501,
     2   outGrp_Y = 502, outGrp_D = 503, outGrp_I = 504,
     3   outGrp_WSS = 505, outGrp_vort = 506, outGrp_eFlx = 507,
     4   outGrp_hFlx = 508, outGrp_absV = 509, outGrp_stInv = 510,
     5   outGrp_vortex = 511, outGrp_trac = 512, outGrp_stress = 513,
     6   outGrp_fN = 514, outGrp_fA = 515, outGrp_J = 516, outGrp_F=517,
     7   outGrp_divV = 518, outGrp_Visc = 519
!--------------------------------------------------------------------
      INTEGER(KIND=IKIND), PARAMETER :: out_velocity = 599,
     2   out_pressure = 598, out_acceleration = 597,
     3   out_temperature = 596, out_WSS = 595, out_vorticity = 594,
     4   out_displacement = 593, out_integ = 592, out_energyFlux = 591,
     5   out_heatFlux = 590, out_strainInv = 589, out_absVelocity = 588,
     6   out_vortex = 587, out_traction = 586, out_stress = 585,
     7   out_fibDir = 584, out_fibAlign = 583, out_actionPotential =582,
     8   out_jacobian = 581, out_defGrad=580, out_divergence = 579,
     9   out_viscosity = 578
!--------------------------------------------------------------------
!     Mesher choice for remeshing for moving wall problems
      INTEGER(KIND=IKIND), PARAMETER :: RMSH_TETGEN = 1,
     2   RMSH_MESHSIM = 2
!--------------------------------------------------------------------
!     Type of constitutive model (isochoric) for structure equation:
!     St.Venant-Kirchhoff, modified St.Venant-Kirchhoff, NeoHookean,
!     Mooney-Rivlin, modified Holzapfel-Gasser-Ogden with dispersion,
!     Linear model (S = mu*I), Guccione (1995), Holzapfel & Ogden model
!     for myocardium (2009)
      INTEGER(KIND=IKIND), PARAMETER :: stIso_NA = 600,
     2   stIso_StVK = 601, stIso_mStVK = 602, stIso_nHook = 603,
     3   stIso_MR = 604, stIso_HGO = 605, stIso_lin = 606,
     4   stIso_Gucci = 607, stIso_HO = 608
!--------------------------------------------------------------------
!     Type of constitutive model (volumetric) for structure eqn:
!     Quadratic, Simo-Taylor91, Miehe94
      INTEGER(KIND=IKIND), PARAMETER :: stVol_NA = 650,
     2   stVol_Quad = 651, stVol_ST91 = 652, stVol_M94 = 653
!--------------------------------------------------------------------
!     Type of fluid viscosity: constant, Carreau-Yasuda shear-thinning
!     model, Cassons non-Newtonian model
      INTEGER(KIND=IKIND), PARAMETER :: viscType_NA = 699,
     2   viscType_Const = 698, viscType_CY = 697, viscType_Cass = 696
!--------------------------------------------------------------------
!     Preconditioner definitions
      INTEGER(KIND=IKIND), PARAMETER :: PREC_NONE = 700,
     2   PREC_FSILS = 701, PREC_TRILINOS_DIAGONAL = 702,
     3   PREC_TRILINOS_BLOCK_JACOBI = 703, PREC_TRILINOS_ILU = 704,
     4   PREC_TRILINOS_ILUT = 705, PREC_TRILINOS_IC = 706,
     5   PREC_TRILINOS_ICT = 707, PREC_TRILINOS_ML = 708,
     6   PREC_RCS = 709
!--------------------------------------------------------------------
!     Solver definitions
      INTEGER(KIND=IKIND), PARAMETER :: lSolver_NA = 799,
     2   lSolver_CG = 798, lSolver_GMRES=797, lSolver_NS=796,
     3   lSolver_BICGS = 795
!--------------------------------------------------------------------
!     Contact model
      INTEGER(KIND=IKIND), PARAMETER :: cntctM_NA = 800,
     2   cntctM_penalty = 801
!--------------------------------------------------------------------
!     Immersed Boundary (IB) treatment
      INTEGER(KIND=IKIND), PARAMETER :: ibMthd_NA = 850,
     2   ibMthd_SSM = 851, ibMthd_IFEM = 852
!--------------------------------------------------------------------
!#######################################################################
