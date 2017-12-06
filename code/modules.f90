!!! Time-stamp: <modules.f90 17:41, 31 Jan 2004 by P Sunthar>

!________________________________________________________________
!   Global modules and interface headers
!________________________________________________________________

!!! $Log: modules.f90,v $
!!! Revision 1.1  2004/01/29 22:13:11  pxs565
!!! Initial revision
!!!


Module Bspglocons ! Bead Spring simulation, Global constant (parameters)
  Save

  Integer, Parameter :: Ndim = 3 ! Dimension of simulation
  Integer, Parameter ::  MaxXdist = 101
  !!! choose precision
  Integer, Parameter :: k4b = Selected_int_kind(9)
  Integer, Parameter :: SNGL = Selected_real_kind(4)
  Integer, Parameter :: DOBL = Selected_real_kind(8)
 ! Integer, Parameter :: DBprec = SNGL
  Integer, Parameter :: DBprec = DOBL


  Integer, Parameter :: NProps = 17
  Integer, Parameter :: MAXCHEB = 500
  Real (DBprec), Parameter :: PI = 3.14159265358979323846
  Real (DBprec), Parameter :: TINI = 1e-25
  Real (DBprec), Parameter :: MYEPS = 1e-6
  Integer, Parameter :: HOOK = 1, FENE = 2, ILC = 3, WLC = 4, Fraenkel = 5
End Module Bspglocons


Module Bspglovars ! Bead Spring simulation, Global variables
  Use Bspglocons


!!$  Real Cubsoln_lu(0:1000), Gama_inc, Gama_max
!!$  Real, Allocatable ::Cubsoln_lu_2d(:,:), Gama_inc_2d(:), Gama_max_2d(:)

  Character(10) :: Prop_names(NProps) 

  Character(10) , Parameter :: Correl_names(1) = (/"Sxy"/)

  Real (DBprec) Imploop_tol, Tstep_conv_tol
End Module Bspglovars

Module Flowvars
  Use bspglocons
  Integer, Parameter :: EQ = 0, SH = 1, UA = 2, PL = 3,  UR = 4, PU = 5, PP = 6
  Integer FlowType
       ! EQ equilibrium, no flow
       ! SH Planar shear
       ! UA Uniaxial Elongational
       ! PL Planar Elongation
       ! UR Uniaxial extension followed by relaxation
       ! PU Periodic uniaxial extension 
  Real (DBprec) gdots
end Module Flowvars



Module Bspintfacs
  Use Bspglocons

  !---------------------------------------------------------------------c	
  !     Driver subroutines                                              c
  !---------------------------------------------------------------------c

  Interface
     Subroutine Initial_position(Stype,N,L0s,R,seed)
       Use bspglocons
       Integer, Intent (in) :: Stype
       Integer(k4b), Intent (in) :: N
       Real (DBprec), intent (in) :: L0s
       Integer(k4b), Intent (inout) :: seed
       Real (DBprec), Intent (out) :: R(:,:)
       !Real, intent (out) :: R(Ndim,N)
     End Subroutine Initial_position
  End Interface

  Interface
     Subroutine GaussRand(N,GX,seed)
       Use bspglocons
       Integer(k4b), Intent (in) :: N
       Integer(k4b), Intent (inout) :: seed
       Real (DBprec), Intent (out), Dimension(N) :: GX
     End Subroutine GaussRand
  End Interface


  Interface
     Subroutine Time_Integrate_Chain(NBeads, R_Bead, spring_type, &
          tcur, tmax, Delts, &
          Hstar, Zstar, Dstar, L0s, Q0s, &
          seed1, Nsamples, times, samples, phi)

       Use bspglocons
       Integer, Intent(in) :: NBeads
       Real (DBprec), Intent (inout), Dimension(:,:) :: R_Bead ! Pos. vector of Beads
       
       Integer, Intent (in) :: spring_type ! Spring force law type
       ! Hookean = 1, FENE = 2, ILC = 3, WLC = 4
       
       
       Real (DBprec), Intent (in) :: tcur,tmax,Delts  ! Integration time interval specs
       
       Real (DBprec),  Intent (in) :: Hstar         ! HI parameter (non-dim)
       Real (DBprec),  Intent (in) :: Zstar,Dstar   ! EV parameters (non-dim)
       Real (DBprec),  Intent (in) :: L0s           ! finite ext., param sqrt(b)
       Real (DBprec),  Intent (in) :: Q0s
       
       Integer (k4b), Intent (inout) :: seed1 ! seed for rnd number
       
       Integer, Intent(in) :: Nsamples     ! Number of sampling points
       Real (DBprec),  Intent (in), Dimension(:) :: times    ! Sampling instances
       Real (DBprec), Intent (inout), Dimension(:,:) :: samples
       Real (DBprec), Intent (in), Dimension(:,:) :: phi
     End Subroutine Time_Integrate_Chain
  End Interface
  

  Interface
     Subroutine polish_poly_root(c,xin,atol)
       Use bspglocons
       Implicit None
       ! use newton raphson to polish the root of a polynomial
       Integer, Parameter :: n=4
       Real (DBprec), Intent (in) :: c(n)
       Real (DBprec), Intent (inout) :: xin
       Real (DBprec), Intent (in) :: atol
     End Subroutine polish_poly_root
  End Interface

  !---------------------------------------------------------------------c	
  !     Utility subroutines                                             c
  !---------------------------------------------------------------------c

  Interface
     Subroutine ran_1(n, R, idum)
       Use bspglocons
       Integer(k4b), Intent(inout) ::idum
       Integer(k4b), Intent(in) :: n
       Real (DBprec), Intent(inout) :: R(n)
     End Subroutine ran_1
  End Interface

  Interface
     Subroutine maxminev_fi(n, A, maxev, minev)
       Use bspglocons
       Integer n
       Real (DBprec) A(:,:,:,:), maxev, minev
     End Subroutine maxminev_fi
  End Interface

  Interface
     Subroutine chbyshv (L, d_a, d_b, a)
       Use bspglocons
       Integer L
       Real (DBprec) d_a, d_b, a(0:MAXCHEB)
     End Subroutine chbyshv
  End Interface

  Interface
     Subroutine cublu
       Use bspglocons
     End Subroutine cublu
  End Interface
  !---------------------------------------------------------------------c
  !     Properties estimators                                           c
  !---------------------------------------------------------------------c	

  Interface
     Subroutine chain_props(N, Rbead, F, props)
       Use bspglocons
       Implicit None
       Integer, intent (in) :: N
       Real (DBprec), Intent (in), Dimension(:,:)  :: Rbead, F
       Real (DBprec), Intent (out),Dimension(:) ::props
     End Subroutine chain_props
  End Interface

  Interface
     Subroutine time_correl(N, R, F, t0, t, correl)
       Use bspglocons
       Use bspglovars
       Implicit None
       Integer, intent (in) :: N
       Real (DBprec), Intent (in), Dimension(:,:)  :: R, F
       Real (DBprec), Intent (out)  :: correl
       Real (DBprec), Intent (in) :: t0, t
     End Subroutine time_correl
  End Interface

  Interface
     Subroutine numint(f,dx,nmin,nmax,nord,sumf)
       Use bspglocons
       Implicit None

       Real (DBprec), Intent(in), Dimension(:) :: f
       Real (DBprec), Intent(in) :: dx
       Real (DBprec), Intent(out) :: sumf
       Integer, Intent (in) :: nmin, nmax,nord
     End Subroutine numint
  End Interface

  Interface
     Subroutine meanerr(vec,mean,err)
       Use bspglocons
       Real (DBprec), Intent(in) :: vec(:)
       Real (DBprec) , Intent(out) :: mean
       Real (DBprec), Intent(out) :: err
     End Subroutine meanerr
  End Interface

  Interface
     Function normeqpr(Stype, r, b)
       Use bspglocons
       Implicit None
       Integer, Intent (in) :: Stype
       Real (DBprec) , Intent (in) :: b, r
       Real (DBprec) normeqpr
     End Function normeqpr
  End Interface

  Interface
     Function  lam1_th(hs, N)
       use bspglocons
       Real (DBprec) lam1_th
       Real (DBprec) ,  Intent (in) :: hs
       Integer ,  Intent (in) :: N
     End Function lam1_th
  End Interface

  Interface
     subroutine get_kappa(t,K)
       use bspglocons
       Implicit none
       Real (DBprec), intent (in) :: t
       Real (DBprec), intent (in), dimension(:,:) :: K
     end subroutine get_kappa
  end Interface
  
  
    
  End Module Bspintfacs
