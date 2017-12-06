!!! Time-stamp: <gsipc.f90 12:45, 06 Apr 2004 by P Sunthar>

!________________________________________________
!   Bead Spring Configuration Space Utilities
!________________________________________________

!!! $Log: gsipc.f90,v $
!!! Revision 1.3  2004/02/09 05:38:36  pxs565
!!! pcerr made relative
!!!
!!! Revision 1.2  2004/02/03 02:02:35  pxs565
!!! cofm initialised correctly
!!!
!!! Revision 1.1  2004/01/29 22:11:08  pxs565
!!! Initial revision
!!!


Module csputls 
!!! This module contains subroutines used to manipulate the 
!!! configuration variable, R and its derivatives
  Use bspglocons

Contains

  Subroutine b2bvector_sym_up(N,R,b2b)
    Implicit None
    ! Calling arguments
    Integer, Intent(in)  :: N
    Real (DBprec), Intent(in)  :: R(:,:)
    Real (DBprec), Intent(out) :: b2b(:,:,:)

!!!_________________________________________________________|
!!! This routine returns the vector one bead to another bead
!!! as a symmetric matrix but only the upper elements filled
!!!_________________________________________________________|

    Integer mu,nu

    ! according to /usr/bin/prof, this takes max time. so commenting out
    ! since only upper diagonal is anyway reqd
    ! b2b = 0

    ! note only upper half (nu > mu) of the matrix is filled
    ! and the convention is R_12 = R_2 - R_1
    ! where R_12 is the vector from 1 to 2
    Do nu = 2,N
       Do mu = 1,nu-1
          b2b(:,mu,nu) = R(:,nu) - R(:,mu)
       End Do
    End Do
  End Subroutine b2bvector_sym_up

  Subroutine modr_sym_up(N,b2bvec,deltaR) 
    Implicit None
    ! calling arguments
    Integer, Intent (in)  :: N
    Real (DBprec), Intent (in)  :: b2bvec(:,:,:)
    Real (DBprec), Intent (out) :: deltaR(:,:) 
!!!_________________________________________________________|
!!! This subroutine returns the magnitude of each of the bead 
!!! to bead vector
!!!_________________________________________________________|

    Integer mu,nu,i
    Real (DBprec) r12(Ndim),modr


    deltaR = 0.d0

    ! note that we cant use a forall since dot_product is a
    ! transformational function
    Do nu = 2,N
       Do mu = 1,nu-1
          r12 = b2bvec(:,mu,nu)
          modr = 0.d0
          Do i=1,Ndim
             modr = modr + r12(i) * r12(i)
          End Do
          If (modr < 1.0e-12) modr = 1.0e-12   !so that initial dr is != 0
          deltaR(mu,nu) = Sqrt(modr)
       End Do
    End Do
  End Subroutine modr_sym_up


  Subroutine tensor_prod_of_vec(alpha,vec,tensor)
    Implicit None
    Real (DBprec), Intent (in) :: alpha
    Real (DBprec), Intent (in) :: vec(:)
    Real (DBprec), Intent (out) ::  tensor(:,:)
    Integer i,j
    Real (DBprec) temp

    Do j = 1,Ndim
       Do i = 1,j
          temp = alpha * vec(i) * vec(j)
          tensor(i,j) = temp
          tensor(j,i) = temp
       End Do
    End Do
  End Subroutine tensor_prod_of_vec
End Module csputls


Module BSpModel
!!! This module contains the constants and functions which are specific
!!! to the model of spring and EV.
  Use bspglocons
  Use csputls
  Implicit None

  ! These are paramters that will be used  by the following procedures
  Real (DBprec) sqrtb

  ! paramters for EV
  Real (DBprec) zsbyds5,pt5bydssq
  Real (DBprec) LJa,LJb,LJpa,LJpb,Rmin,Rcutg,Rcutp, Rlist
  Integer dlist, maxnab
 ! Integer list(:), point(:,:)


Contains 

  Subroutine force_sans_hookean(sptype,sr,ff,Q0s)
    Integer, Intent (in) :: sptype 
    Real (DBprec), Intent (in) :: sr, Q0s
    Real (DOBL), Intent (out) :: ff

    Real (DOBL) r

    r = sr 
    ! quick way to do double prec calculation.  otherwise, the complete 
    ! code needs to be Dprec to incorporate accuracy at every level

    ! compute the force factor, ff
    Select Case (sptype)
    Case (HOOK) ! Hookean
       ff = 1.d0

    Case (FENE) ! Warner Spring
       ff = 1.d0/(1.d0 - r*r)
     
    Case (ILC)  ! Inverse Langevin, Pade approx 
       ff = (3-r*r)/(1 - r*r)/3.d0

    Case (WLC) ! Worm-like Chain, Marko-Siggia interpolation
       ff = 1./(6*r) * ( 4*r+ 1./(1-r)/(1-r) - 1)
    Case (Fraenkel)! Fraenkel spring
       ff = (1.d0 -(Q0s/r))
        
   
    End Select

  End Subroutine force_sans_hookean

  Subroutine Spring_force(sptype,N,R_Bead,b2bvec,dr,Fs,Q0s)
    Integer, Intent (in) :: sptype
    Integer, Intent (in)  :: N
   ! Real (DBprec), Intent (in)  :: b2bvec(:,:,:)! Ndim x Nbeads x Nbeads
    Real (DBprec), Intent (in)  :: R_Bead(:,:)      ! Nbeads x Nbeads
    Real (DBprec), Intent (inout) :: b2bvec(:,:,:)
    Real (DBprec), Intent (inout) :: dr(:,:)
    Real (DBprec), Intent (out) :: Fs(:,:)      ! Ndim x Nbeads
    Real (DBprec), Intent (in) :: Q0s 
!!! Purpose:
    ! computes the net spring force on a bead due to 
    ! connector springs of neighbouring beads

    Integer nu, mu, i
    Real (DBprec) r , Q0sp! = |Q| / sqrtb
    Real (DOBL) Fcnu(Ndim), Fcnu1(Ndim), nonhook
    Real (DBprec) r12(Ndim), modr
    Fs  = 0.d0
    b2bvec = 0.d0  !! check as its not present earlier
    Fcnu1 = 0.d0 ! Fc of nu-1

    Do nu = 1,N-1 ! over each of the connector
       ! Calculating bead to bead distance
       b2bvec(:,nu,nu+1) = R_Bead(:,nu+1) - R_Bead(:,nu)
       r12 = b2bvec(:,nu,nu+1)
       ! Calculating modr
       modr = 0.d0
       Do i=1,Ndim
          modr = modr + r12(i) * r12(i)
       End Do
       If (modr < 1.0e-12) modr = 1.0e-12
       dr(nu,nu+1) =Sqrt(modr)
      ! r = dr(nu,nu+1)/sqrtb ! only upper diagonal
       r = dr(nu,nu+1)/sqrtb
       !write(*,*) "dr", dr
       Q0sp = Q0s/sqrtb
       If (Abs(r - 1.d0) .Lt. MYEPS) r = 1 - MYEPS

       Call force_sans_hookean(sptype,r,nonhook,Q0sp)

       ! connector force from 1 to 2 (nu to nu+1)
       Fcnu = b2bvec(:,nu,nu+1) * nonhook

       ! total spring force on a bead
       Fs(:,nu) = Fcnu - Fcnu1
       Fcnu1 = Fcnu ! save for the next bead
    End Do
      
    ! and the last bead
    Fs(:,N) = - Fcnu1

  End Subroutine Spring_force


  Subroutine solve_implicit_r(sptype,dtby4,gama,natscl,r,ff)
    Integer, Intent (in) :: sptype
    Real (DBprec), Intent (in) :: dtby4, gama, natscl
    Real (DBprec), Intent (out) :: r
    Real (DOBL), Intent (out) :: ff

    !! Purpose
    ! solves the equation r ( 1 + dt/4 ff ) == gama
    ! where ff depends on spring force law F = H Q ff,
    ! for r and returns, r and ff(r)

    Real (DBprec) coeff(4),denom

    ! set up the polynomial equation (cubic), and the guess for
    ! gama >> 1, obtained by the asymptotic behaviour

    coeff(4) = 1.d0

    Select Case (sptype)
    Case (HOOK) ! Hookean
       r = gama/(1.d0 +dtby4)
       ff = 1.d0
       Return

    Case (FENE) ! FENE
       coeff(1) = gama
       coeff(2) = -(1.d0 + dtby4)
       coeff(3) = -gama
       r = 1.d0 - dtby4/2.d0/gama 

    Case (ILC) ! ILC
       denom = 3. + dtby4
       coeff(1) = 3. * gama / denom
       coeff(2) = -(3. + 3.* dtby4) / denom
       coeff(3) = -3. * gama / denom
       r = 1 - dtby4/3/gama

    Case (WLC) ! WLC
       denom = 2*(1.5 + dtby4)
       coeff(1) = -3 *gama/denom
       coeff(2) =  3 * (1. + dtby4 + 2.*gama) / denom
       coeff(3) = -1.5 * (4. + 3.*dtby4 + 2.*gama) / denom
       r = 1 - Sqrt(dtby4/6/gama)

    Case (Fraenkel) ! Fraenkel spring
       r = (gama + (dtby4*natscl))/(1.+dtby4)
       ff = 1.d0-(natscl/r)
       Return

    End Select

    ! the Hookean guess is same for all the force laws
    If (gama < 1.0) Then 
       r = gama/(1.d0 +dtby4)  
    End If

    ! all the common forces laws yeild a cubic in the implicit form
    ! polish the guessed root by a Newt-Raph for polynomials

    ! for WLC, the newt raph, does not converge to the correct root,
    ! for gama > 100 , therefore simply ignore polishing
    If (sptype .Eq. WLC .And. gama > 100) Then
       r = r
    Else 
       Call polish_poly_root(coeff,r,1e-6)
    End If


    Call force_sans_hookean(sptype,r,ff,natscl)

  End Subroutine solve_implicit_r


  Subroutine Excluded_volume_force(N,R_Bead,b2bvec,dr,R_list_save, &
                         R0,point,list,Fev,phi)
    Use bspglocons

!!! This routine returns excluded volume force
    Integer, intent (in) :: N
   ! Real (DBprec), Intent (in)  :: b2bvec(:,:,:)! Ndim x Nbeads x Nbeads
    Real (DBprec), Intent (in)  :: R_Bead(:,:)      ! Nbeads x Nbeads
    Real (DBprec), Intent (inout) :: b2bvec(:,:,:) 
    Real (DBprec), Intent (inout) :: dr(:,:)
    Real (DBprec), Intent (in) :: R_list_save(:,:) ! absolute position
    Real (DBprec), Intent (inout) :: R0(:,:) ! last neighbour list update
    Integer, Intent (inout) :: point(:)
    Integer, Intent (inout) :: list(:)
    Real (DBprec), Intent (out) :: Fev(:,:)     ! Ndim x Nbeads 
    Real (DBprec), Intent (in) :: phi(:,:)
    Integer nu,mu, UPDATE
    Real (DBprec) :: Fpair12(Ndim), r12(Ndim), modr
    Real (DBprec) :: alpha, beta 
!!! Parameters for neighbour list
    Integer nlist
   ! Integer list(maxnab), point(N*N)
    Integer Jbeg, Jend, Jnab, i
    
    dr = 0.d0
    Fev  = 0.d0
    alpha = 3.1730728678d0
    beta = -0.856228645d0
    If (LJa .Gt. 0.d0) Then
       ! caution donot use forall here, i means F12 will b evaluated for 
       ! all the values of mu,nu and then assigned, which is not what we 
       ! want.  The rule to use forall statements is that lhs must be 
       ! unique for-all the values of loop index, and RHS must not depend 
       ! on any other index of lhs, other than the lhs itself
       ! ie, we can have  foo(mu) = foo(mu) + bar(nu)
       ! but not  foo(mu) = foo(mu+3) + bar(nu)
          ! NEIGHBOUR LIST
       Call CheckNeighbour(N,R_list_save,R0,UPDATE)
       !UPDATE = 1 
       Select case (UPDATE)
       Case(1)
             Call SaveConfiguration(N,R_list_save,R0)
             write(*,*) 'dlist_on update', dlist
            nlist = 0
            list = 0        
            Do nu = 2,N
               point(nu) = nlist + 1 
               write(*,*) "point", point(nu), 'nu', nu
               Do mu  = 1,nu-1 
               !!! call bead to bead amd 
               b2bvec(:,mu,nu) = R_Bead(:,nu) - R_Bead(:,mu)
               r12 = b2bvec(:,mu,nu)
               ! Calculating modr
               modr = 0.d0
               Do i=1,Ndim
                   modr = modr + r12(i) * r12(i)
               End Do
               If (modr < 1.0e-12) modr = 1.0e-12
               dr(mu,nu) =Sqrt(modr)
                   
               If (dr(mu,nu) .le. Rlist) Then
                   nlist = nlist + 1
                   list(nlist) = mu 
           
                  ! Time to calculate force
                  ! force between the pair mu,nu, since the EV force is
                  ! repulsive, we take it to be in the opposite direction
                  ! of the connector vector mu -> nu
                  ! convention followed: force is positive for attraction
               If (dr(mu,nu) .le. Rcutg .and. dr(mu,nu) .ge. Rmin) Then
                   Fpair12 = - b2bvec(:,mu,nu) &
                *(LJa/(dr(mu,nu)**(LJpa+2.d0))-LJb/(dr(mu,nu)**(LJpb+2.d0)))
          
                    Fev(:,mu) = Fev(:,mu) + Fpair12
                    Fev(:,nu) = Fev(:,nu) - Fpair12

               Else If (dr(mu,nu) .lt. Rmin) Then
                    Fpair12 = - b2bvec(:,mu,nu) &
                 *(LJa/(Rmin**(LJpa+2.d0))-LJb/(Rmin**(LJpb+2.d0)))

                    Fev(:,mu) = Fev(:,mu) + Fpair12
                    Fev(:,nu) = Fev(:,nu) - Fpair12

               Else If (dr(mu,nu) .le. Rcutp .and. dr(mu,nu) .gt. Rcutg) Then 
                    Fpair12= -b2bvec(:,mu,nu) &
                 *alpha*phi(mu,nu)*sin(alpha*dr(mu,nu)*dr(mu,nu) + beta)
                     Fev(:,mu) = Fev(:,mu) + Fpair12
                     Fev(:,nu) = Fev(:,nu) - Fpair12
 
              End If 
               End If  ! Rlist
               End Do !mu
             End Do !nu
             point(N+1) = nlist + 1 
 !            write(*,*) "b2bvec", b2bvec
 !            write(*,*) "dr", dr
 !            write(*,*) "R_bead", R_bead
 !            write(*,*) "Fev", Fev
       Case(2)
            write(*,*) "dlist without update", dlist
            Do nu = 2,N
              Jbeg = point(nu)
              Jend = point (nu + 1) - 1 
              If (Jbeg .le. Jend) Then
                 write(*,*) "jbeg", jbeg, 'jend', jend, 'fornu', nu
                 Do Jnab = Jbeg, Jend
                   mu = list(Jnab)
                   write(*,*) "mu", mu
                   b2bvec(:,mu,nu) = R_Bead(:,nu) - R_Bead(:,mu)
                   r12 = b2bvec(:,mu,nu)
                   modr = 0.d0
                   Do i=1,Ndim
                       modr = modr + r12(i) * r12(i)
                   End Do
                   If (modr < 1.0e-12) modr = 1.0e-12 
                  !so that initial dr is != 0
                   dr(mu,nu) = Sqrt(modr)

                   ! Force calculation
                   If (dr(mu,nu) .le. Rcutg .and. dr(mu,nu) .ge. Rmin) Then
                      Fpair12 = - b2bvec(:,mu,nu) &
                     *(LJa/(dr(mu,nu)**(LJpa+2.d0))-LJb/(dr(mu,nu)**(LJpb+2.d0)))
                      
                        Fev(:,mu) = Fev(:,mu) + Fpair12
                        Fev(:,nu) = Fev(:,nu) - Fpair12
      

                   Else If (dr(mu,nu) .lt. Rmin) Then
                      Fpair12 = - b2bvec(:,mu,nu) &
                    *(LJa/(Rmin**(LJpa+2.d0))-LJb/(Rmin**(LJpb+2.d0)))

                        Fev(:,mu) = Fev(:,mu) + Fpair12
                        Fev(:,nu) = Fev(:,nu) - Fpair12

                   Else If (dr(mu,nu) .le. Rcutp .and. dr(mu,nu) .gt. Rcutg) Then
                      Fpair12= -b2bvec(:,mu,nu) &
                    *alpha*phi(mu,nu)*sin(alpha*dr(mu,nu)*dr(mu,nu) + beta)                   

                      Fev(:,mu) = Fev(:,mu) + Fpair12
                      Fev(:,nu) = Fev(:,nu) - Fpair12
  
                 End If
!                      Fev(:,mu) = Fev(:,mu) + Fpair12
!                      Fev(:,nu) = Fev(:,nu) - Fpair12
              
               End Do  ! Jnab
               End If ! Jbeg
      !         write(*,*) "outside jnab"
           End Do  ! nu
      !      write(*,*) "b2bvec", b2bvec
      !       write(*,*) "dr", dr
      !       write(*,*) "R_bead", R_bead
      !       write(*,*) "Fev", Fev

       End Select
   End If ! LJa
  End Subroutine Excluded_volume_force

  Subroutine CheckNeighbour(N,R_Bead,R0,UPDATE)
    Integer, Intent(in) :: N
    Real (DBprec), Intent (in) :: R_Bead(:,:)
    Real (DBprec), Intent (in) :: R0(:,:)
    Integer, Intent(inout) :: UPDATE
    Real (DBprec) :: dismax, b2b(Ndim,N), r11(Ndim)
    Real(DBprec) :: modr
    Integer :: i, j
    dismax = 0.d0
    Do i = 1,N
      modr=0.d0
      b2b(:,i) = R_Bead(:,i) -R0(:,i)
      r11= b2b(:,i)
      Do j = 1, Ndim
      
      ! dismax = max(abs(R_Bead(j,i)-R0(j,i)), dismax)
      modr =  modr + r11(j)*r11(j)
      dismax = max(sqrt(modr),dismax)
      End Do
    End Do
   ! dismax = 2.d0 * Sqrt(3.d0*(dismax*dismax))

    If (dismax .gt. (Rlist-Rcutp)) Then
        UPDATE = 1
    Else 
        UPDATE = 2
    End If

    If (dlist .eq. 1) UPDATE =1
       
  End Subroutine CheckNeighbour

  Subroutine SaveConfiguration(N,R,R0)
    Integer , Intent(in) :: N
    Real(DBprec), Intent (in) :: R(:,:)
    Real(DBprec), Intent (inout) :: R0(:,:) 
     R0 = R
  End Subroutine SaveConfiguration
End Module BSpModel


Subroutine Time_Integrate_Chain(NBeads, R_Bead, spring_type, &
     tcur, tmax, Delts, &
     Hstar, Zstar, Dstar, L0s, Q0s, &
     seed1, Nsamples, times, samples, phi)
  Use bspglocons
  Use bspglovars
  Use Flowvars
  Use bspintfacs, TIC_from_bspintfacs => Time_Integrate_Chain
  Use BSpModel
  Use csputls

  !Use blas_interfaces unable to obtain proper interface 
  !which allows matrices of arbitrary rank to b passed
  Implicit None

  Integer, Intent(in) :: NBeads
  Real (DBprec), Intent (inout), Dimension(:,:) :: R_Bead ! Initial positions of Beads

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


!!!_________________________________________________________|
  !     This driver uses:                                              
  !     Predictor-corrector        scheme for excluded volume          
  !     Fully implicit scheme for the spring force                     
  !           based on Ottinger's and Somasi et al.'s suggestions      
  !     Warner spring                                                  
  !     Rotne-Prager-Yamakawa HI tensor with Chebyshev polynomial     
  !         approx. using Fixman's approximation to get conditioner    
  !     Narrow Gaussian EV potential                                   
!!!_________________________________________________________|

  ! external BLAS functions
   !  Real snrm2,sdot
  Real (SNGL) snrm2, sdot
  Real (DOBL) dnrm2, ddot

  Real (DBprec) time, fd_err, fd_err_max, om1, om2, rs,&
       rsbyhs, hsbyrs, hsbyrs2, L_max, L_min, d_a, d_b,&
       gama_mag, lt_err, dtsby4, pcerr, delx2
  Real (DBprec)&
       ! various scratch storages for position vector
       cofm(Ndim), & ! center of mass
       R_pred(Ndim,NBeads), &
       DR_pred(Ndim,NBeads), &
       Ups_pred(Ndim,NBeads), &
       R_pred_old(Ndim,NBeads), &
       R_corr(Ndim,NBeads), &
       Rtemp(Ndim,NBeads), &
       delta_R(NBeads,NBeads),&
       b2b_sup(Ndim,Nbeads,Nbeads), & ! bead to bead vector (super diagonal)
       deltaR_sup(Nbeads,Nbeads), & ! bead to bead dist (super diagonal)
       ! forces
       F_ev(Ndim,NBeads), &     ! excluded volume forces
       F_spring(Ndim,NBeads), & ! spring forces
       F_tot(Ndim,Nbeads), &    ! total forces
       ! symmetric diffusion tensor, only upper diagonal elements stored
       Diffusion_sup(Ndim,NBeads,Ndim,NBeads), &
       ! another scratch storage for diffusion tensor
       Dp_sym_up(Ndim,NBeads,Ndim,NBeads), &
       DelS(Ndim,NBeads),&
       kappa(Ndim,Ndim),&  !flow field divergence
       ! for Chebyshev polynomials
       X_0(Ndim,NBeads), X_l_1(Ndim,NBeads), &
       X_l(Ndim,NBeads), X_lp1(Ndim,NBeads), &
       
       cheba(0:MAXCHEB), gama_mu(Ndim), &
       diffD(Ndim,NBeads-1,Ndim,NBeads), &
       diffUps(Ndim,NBeads-1)
  
  ! some connector forces and related requiring higher precision
  Real (DOBL) ff,  &
       F_con_mu(Ndim), &    ! connector force between beads
       F_con_mu1(Ndim)    ! connector force, save for mu-1 

  Logical fd_check, ncheb_flag

  !     Other definitions...
  Integer i,  l, mu, nu, ncheb, ncheb_old, isample, lt_count
  Real (DBprec) temp1, sqrt2inv, TPI, RPI, TRPI, &
       C1, C2, C3, C4, C5, C6, C7, RPI3by4 
  Real (DBprec) r 


  !     Definitions for the BLAS-2 routine "ssymv"
  Integer Ndof, lda, ldadiff, incx, incy
  Real (DBprec) alpha, beta
  Real (DBprec) :: R0(Ndim,NBeads), R_list_save(Ndim,Nbeads)
  Integer list(Nbeads*Nbeads), point(Nbeads*Nbeads)

  sqrtb = L0s ! sqrtb will be needed by other modules
    
  ! for EV
  zsbyds5 = Zstar/(Dstar**5.0)
  pt5bydssq = 0.5/(Dstar*Dstar)
    
  LJpa = 12.D0  !LJ purely repulsive interaction potential power.
  LJpb = 6.D0  !LJ purely attractive interaction potential power.
  LJa  = LJpa*Zstar*(Dstar**LJpa)
  LJb  = LJpb*Zstar*(Dstar**LJpb)
  Rmin = 0.7D0*Dstar
  Rcutg = Dstar*(2.D0**(1.D0/6.D0))
  Rcutp = 1.5D0*Dstar
  Rlist = 2.d0 * Dstar
  dlist = 0
 ! maxnab = NBeads*NBeads
  sqrt2inv = 1.d0/(2.d0**0.5d0)
  TPI = 2.d0*PI
  RPI = Sqrt(PI)
  TRPI = 2.d0*RPI
  C1 = TPI/3.d0
  C2 = 1.d0
  C3 = 0.75d0/RPI*0.375d0 !C3 = 9./32/RPI
  C4 = 8.d0*PI
  C5 = 14.14855378d0
  C6 = 1.21569221d0
  C7 = 0.09375d0/RPI
  RPI3by4 = RPI*0.75d0


  incx = 1
  incy =1
  Ndof = Ndim*NBeads  ! degrees of freedom
  lda = Ndof

  fd_err_max = 0.0025d0 
  pcerr = 0.0d0


  delta_R = 0.0d0


  Diffusion_sup = 0.0d0
  ! diagonal elements
  Forall (mu = 1:Nbeads, i = 1:Ndim )
     Diffusion_sup(i,mu,i,mu) = 1.d0
  End Forall


  !     Get initial estimates of L_max and L_min using
  !     Kroeger et al's approx. of Zimm theory's preds.;
  !     the number of Chebyshev terms, and the Chebyshev coefficients
  L_max = 2.d0*(1.d0 + PI*Sqrt(Dble(NBeads))*Hstar)
  L_min = (1.d0 - 1.71d0*Hstar)/2     ! Note Hstar < 0.58
  ncheb = Int(Sqrt(L_max/L_min)+0.5d0)+1
  ncheb_flag = .False.
  d_a = 2.d0/(L_max-L_min)
  d_b = -(L_max+L_min)/(L_max-L_min)
  Call chbyshv(ncheb, d_a, d_b, cheba)
  DelS = 0.d0      

  isample = 1
  delx2 = 0

  time = tcur

  !___________________________________________________________________|
  !                The Time integration loop begins here              |
  !___________________________________________________________________|



  Overtime: Do While (time.Le.Tmax+Delts/2.d0)
     dlist = dlist+1
     ! find the center of mass
     cofm = 0.d0
     ! Saving absolute position for Neighbour list
     R_list_save = R_bead
     Do mu = 1, NBeads
        cofm = cofm + R_bead(:,mu)
     End Do

     cofm = cofm/NBeads

     ! shift the origin to the current center of mass
     Forall (mu = 1: NBeads)
        R_Bead(:,mu) = R_Bead(:,mu) - cofm
     End Forall

     ! calculate the bead to bead vector and its magnitude
     ! for use in the forces
     ! Call b2bvector_sym_up(NBeads,R_Bead,b2b_sup)
     ! Call modr_sym_up(NBeads,b2b_sup,deltaR_sup)


     ! change the nature of the spring, in the specific subroutine
     ! in Spring_force(), in module Bead_spring_model
     Call Spring_force(spring_type,NBeads,R_Bead,b2b_sup,deltaR_sup, &
                             F_spring,Q0s)
     Call Excluded_volume_force(NBeads,R_Bead,b2b_sup,deltaR_sup, &
                  R_list_save, R0,point,list,F_ev,phi)
     dlist=dlist+1

     ! Therefore, total force on each bead...    
     F_tot = F_spring + F_ev


     !        Calculate the RPY diffusion tensor
     If (Hstar.Gt.0) Then
        Do nu = 2,NBeads
           Do mu = 1,nu-1
              rs = deltaR_sup(mu,nu)
              hsbyrs = Hstar/rs
              rsbyhs = rs/Hstar
              hsbyrs2 = hsbyrs*hsbyrs
              If (rsbyhs.Ge.TRPI) Then
                 om1 = RPI3by4*hsbyrs*(1.0+C1*hsbyrs2)
                 om2 = RPI3by4*(hsbyrs/(rs*rs))*(1-TPI*hsbyrs2)
              Else
                 om1 = 1-C3*rsbyhs
                 om2 = C7*rsbyhs/(rs*rs)
              End If

              ! store only the upper triangular in the mu,nu index.
              Call tensor_prod_of_vec(om2,b2b_sup(:,mu,nu), &
                   Diffusion_sup(:,mu,:,nu) )

              ! additional factor for diagonal elements, in each mu,nu
              Forall (i = 1:Ndim )
                 Diffusion_sup(i,mu,i,nu) = Diffusion_sup(i,mu,i,nu) + om1
              End Forall
           End Do
        End Do
     End If ! Hstar

     ! mean square displacement of cofm
     ! since the bead positions are centered w.r.t. cofm at every
     ! instant, cofm is the incremental displacement 
     if (time .gt. tcur) then 
        delx2 = delx2 + Sum(cofm * cofm)
     end if
     

     !___________________________________________________________________|
     !       Take samples when required                                  |
     !___________________________________________________________________|

     If (Nsamples.Gt.0) Then
        ! ideally it shud b delts/2, but owing to precision errors
        ! we keep it slighly greater than 0.5
        If(Abs(time-times(isample)).Le.Delts*0.51d0) Then
           Call chain_props (NBeads, R_Bead, F_tot, samples(1:14,isample))

           ! Save information on the error between the predictor
           ! and the final corrector
           samples(15,isample) = pcerr

           !! obtain time correlations
           If (gdots .Eq. 0) Then
              Call time_correl(NBeads, R_Bead, F_tot, tcur, time, &
                   samples(16,isample)) 
           End If

           ! Diffusivity
           If (time > 0.0d0) samples(17,isample) = delx2/2/Ndim/time


           isample = isample + 1
        End If
     End If




     !____________________________________________________________________|
     !    Chebyshev polynomial approximation of DelS begins here          |
     !____________________________________________________________________|

     !Generate the random vector X_0
     X_0 = 0.d0
     Call ran_1(Ndof, X_0, seed1)
     X_0 = X_0 - 0.5d0
     X_0 = (X_0*X_0*C5 + C6)*X_0! Element-wise multiplications

     ! Has check for deviation from fluctuation-dissipation theorem
     ! been performed?
     fd_check = .False.

     If (Hstar.Gt.0) Then          


        FDloop: Do
           ! Update DelS vector
           DelS = cheba(0) * X_0

           ! Shift the D matrix
           Dp_sym_up = d_a * Diffusion_sup

           ! diagonal elements
           Forall (mu = 1:Nbeads, i = 1:Ndim )
              Dp_sym_up(i,mu,i,mu) = Dp_sym_up(i,mu,i,mu) + d_b
           End Forall

           ! Calculate the second Chebyshev vector            
           X_l_1 = X_0
           alpha = 1.D0
           beta = 0.D0

           ! BLAS2 symmetric matrix-vector multiplication
           Call dsymv('U', Ndof, alpha,Dp_sym_up,lda,  X_l_1,incx, &
                beta,X_l,incy)  

           ! Update DelS vector            
           DelS = DelS + cheba(1)*X_l

           Do l = 2,ncheb
              alpha = 2.D0
              beta = 0.D0
              Call dsymv('U', Ndof, alpha,Dp_sym_up,lda, X_l,incx,  &
                   beta,X_lp1,incy)
              X_lp1 = X_lp1-X_l_1
              X_l_1 = X_l

              ! The l-th Chebyshev vector
              X_l = X_lp1

              ! Update DelS vector               
              DelS = DelS + cheba(l)*X_l

           End Do

           ! Calculate the deviation from 
           ! the fluctuation-dissipation theorem

           If (.Not.fd_check) Then
              fd_err = ddot(Ndof, DelS, 1, DelS, 1) ! BLAS-1 function
              alpha = 1.D0
              beta = 0.D0

              !Use D:X_0X_0 = X_0.D.X_0 
              !Get D.X_0 first, and then get the dot product of X_0 with the
              ! resulting vector. X_l is reused 
              Call dsymv('U', Ndof,alpha,Diffusion_sup,lda, &
                   X_0,incx, beta,X_l,incy)       
              temp1 = ddot(Ndof, X_0, 1, X_l, 1)
              fd_err = Abs((fd_err - temp1)/temp1)
           End If

           If ((fd_err.Le.fd_err_max).Or.fd_check) Then
              Exit FDloop
           Else 
              !  If the fd_check has been performed and deviation is large
              !  recaclulate the number of Chebyshev terms required.
              !  First, get the maximum and minimum eigen values of D
              !  using Fixman's suggestion.
              fd_check = .True.
              Call maxminev_fi(Ndof, Diffusion_sup, L_max, L_min)
              ncheb_old = ncheb
              ncheb = Int(Sqrt(L_max/L_min)+0.5)+1
              If ((ncheb/ncheb_old).Gt.2) ncheb_flag = .True.
              If (ncheb.Gt.500) ncheb = 500
              d_a = 2/(L_max-L_min)
              d_b = -(L_max+L_min)/(L_max-L_min)
              Call chbyshv(ncheb, d_a, d_b, cheba)
           End If
        End Do FDloop
     Else ! Hstar == 0, the free-draining case              
        DelS = X_0                                                
     End If

     DelS = DelS * sqrt(delts)


     !____________________________________________________________________|
     !        The predictor step                                          |
     !____________________________________________________________________|


     If (Hstar.Gt.0) Then
        alpha = (2.5D0 - 1)*Delts
        beta = 0.D0
        ! Assigns DR_pred <- 0.25*Delts* D.F     
       ! Select Case (DBprec)
      !  Case(SNGL)
        !   Call ssymv('U', Ndof,alpha,Diffusion_sup,lda, F_tot,incx, &
         !       beta,DR_pred,incy)
       ! Case(DOBL)
           Call dsymv('U', Ndof,alpha,Diffusion_sup,lda, F_tot,incx, &
                beta,DR_pred,incy)
        !End Select
 
       ! Call dsymv('U', Ndof,alpha,Diffusion_sup,lda, F_tot,incx, &
        !     beta,DR_pred,incy)   
     Else
        DR_pred = 0.25d0 * Delts * F_tot
     End If


     ! Add the K.R vector
     call get_kappa(time,kappa)
     DR_pred = DR_pred + Delts * Matmul(kappa,R_Bead) 

     ! Eq (14)
     R_pred = R_Bead + DR_pred + sqrt2inv*DelS

     R_pred_old = R_pred


     !____________________________________________________________________|
     !        The first semi-implicit corrector step                      |
     !____________________________________________________________________|

     ! initialise with part of Eq (18)
     Ups_pred = R_Bead + 0.5d0*DR_pred + sqrt2inv*DelS

     If (Zstar > 0) Then
        ! Calculate distances between beads using R_pred
        ! reuse variables
       ! Call b2bvector_sym_up(NBeads,R_pred,b2b_sup)
       ! Call modr_sym_up(NBeads,b2b_sup,deltaR_sup)
       Call Excluded_volume_force(NBeads,R_pred,b2b_sup,deltaR_sup, &
                   R_list_save,R0,point,list,F_ev,phi)

        ! Calculate D.FEV using R_pred and update
        If (Hstar.Gt.0) Then
           alpha = 0.125D0*Delts      ! The prefactor is 1/8 and not 1/4
           beta = 1.D0               ! Add to existing
         ! Select Case (DBprec)
        !  Case(SNGL)
          !     Call ssymv('U', Ndof, alpha,Diffusion_sup,lda, F_ev,incx,  &
         !       beta,Ups_pred,incy)
         ! Case(DOBL)
          Call dsymv('U', Ndof, alpha,Diffusion_sup,lda, F_ev,incx,  &
                beta,Ups_pred,incy) 
          !End Select

        Else
           Ups_pred = Ups_pred + 0.125d0 * Delts * F_ev 
        End If
     End If

     !        Calculate the 0.5*Delts*K.R_pred vector 

     ! Add the K.R vector
        call get_kappa(time+delts,kappa)
        Ups_pred = Ups_pred + 0.5 * Delts * Matmul(kappa,R_pred) 

     ! Eq (18) is completely assembled now

     ! Generating the matrix obtained by using the D_nu
     ! operator on the D super-matrix

     ! (in the following comments indices i,j are omitted for clarity)
     !    diffD(mu,nu) = D(mu+1,nu) - D(mu,nu)
     ! this is true only for nu > mu, since D's elements are computed only
     ! for nu >= mu, and the RHS depends on mu+1.  for nu <= mu+1, the RHS
     ! is rewritten in terms of the symmetric matrix D
     !    diffD(mu,nu) = D(nu,mu+1) - D(nu,mu)
     ! so that all the elements of the RHS are the computed ones

     Forall (mu=1:NBeads-1)
        diffUps(:,mu) = Ups_pred(:,mu+1) - Ups_pred(:,mu) 
     End Forall

     Forall (mu=1:NBeads-1, nu=1:Nbeads, nu > mu)
        diffD(:,mu,:,nu) = Diffusion_sup(:,mu+1,:,nu) - Diffusion_sup(:,mu,:,nu)
     End Forall
     Forall (mu=1:NBeads-1, nu=1:Nbeads, nu <= mu)
        diffD(:,mu,:,nu) = Diffusion_sup(:,nu,:,mu+1) - Diffusion_sup(:,nu,:,mu)
     End Forall


     ! Start the loop for solving for connector vectors
     R_corr = R_pred

     ldadiff = Ndim*(Nbeads-1)

     lt_count = 0

     dtsby4 = Delts/4.0D0 


     Keepdoing: Do

        F_con_mu1 = 0.d0

        oversprings: Do mu = 1,NBeads-1


           ! the connector force for this spring is obtained from
           ! Fs(mu) = Fc(mu) - Fc(mu-1)

           F_con_mu = F_con_mu1 + F_spring(:,mu) 
           ! note:  F_spring contains forces evaluated with the 
           !          predictor Q for all beads < mu
           !          corrector Q for all beads > mu

           If ( Hstar.Gt.0 ) Then
              gama_mu = 0.125d0*Delts* &
                   Matmul( &
                   Reshape(diffD(:,mu,:,:),(/ Ndim, Ndim*Nbeads /) ),  &
                   Reshape(F_spring,             (/ Ndim*Nbeads /)) &
                   )
           Else
              gama_mu  = 0.125d0*Delts*(F_spring(:,mu+1) - F_spring(:,mu))
           End If

           ! the remaining terms on the RHS of Eq.(20)
           gama_mu = gama_mu + diffUps(:,mu) + 0.25d0 * F_con_mu * Delts 

           gama_mag = dnrm2(Ndim,gama_mu,1) ! BLAS single normal two

           ! r = Q_nu_mag/sqrt(b) varies from (0,1)
           Call solve_implicit_r(spring_type,dtsby4,gama_mag/sqrtb,Q0s/sqrtb,r,ff)
           ! implicit solution of Eq.(21)
           !   r ( 1 + deltat/4 ff) -  |gama_mu|/sqrtb == 0
           ! and returns r and ff, the factor in the spring force other 
           ! than hookean F = H Q ff

           If (Abs(r-1.) .Lt. 1e-6) r = 1 - 1e-6 ! TODO, y this

           ! unit vector for Q_mu same as for Gama_mu
           gama_mu = gama_mu/gama_mag 
           gama_mu = gama_mu*r*sqrtb; ! connector vector update

           ! corrector positions vector update
           R_corr(:,mu+1) = R_corr(:,mu) + gama_mu

           ! Connector force 
           F_con_mu1 = gama_mu * ff ! to b used for next spring

           ! update the spring forces, and connector force for this spring
           ! which will b used at the start of the loop
           F_spring(:,mu)   = F_spring(:,mu)   + (F_con_mu1 - F_con_mu)
           F_spring(:,mu+1) = F_spring(:,mu+1) - (F_con_mu1 - F_con_mu)


        End Do oversprings

        Rtemp = R_corr - R_pred
        !lt_err = dnrm2(Ndof,Rtemp,1)/NBeads
        lt_err = dnrm2(Ndof,Rtemp,1)/dnrm2(Ndof,R_Bead,1)

        lt_count = lt_count + 1

        !If ((lt_err.Lt.imploop_tol).Or.(lt_count.Gt.2*Nbeads*Nbeads)) Exit
        !If ((lt_err.Lt.imploop_tol)) Exit
        If ((lt_err.Lt.imploop_tol).Or.(lt_count > 10*Nbeads)) Exit

        R_pred = R_corr
     End Do Keepdoing

     if (lt_count > 10*NBeads) write (*,1024) lt_count, Delts , Gdots*time
1024 format ('Loop exceeded ', I4, ' for dt = ', F11.4, ' at strain ', G11.4)


     ! Final Major updates, new position vector and time
     R_Bead = R_corr
     time = time + Delts
  ! write(*,*) 'after crossing loop', R_bead
     Rtemp = R_Bead - R_pred_old
     pcerr = dnrm2(Ndof,Rtemp,1)/NBeads
     pcerr = dnrm2(Ndof,Rtemp,1)/dnrm2(Ndof,R_Bead,1)
  ! write(*,*) "AFter 1 time step rtemp", Rtemp 
     If (ncheb_flag) Then
        ncheb = ncheb_old
        ncheb_flag = .False.
     End If
  End Do    Overtime


End Subroutine Time_Integrate_Chain

