! Time-stamp: <sensemble.f90 18:43, 14 Oct 2007 by P Sunthar>


!!! $Log: sensemble.f90,v $
!!! Revision 1.3  2004/03/17 02:05:24  pxs565
!!! Working version of a sample demonstration
!!!
!!! Revision 1.2  2004/03/16 21:48:12  pxs565
!!! minor variations, tried using an ntdone file,
!!!
!!! Revision 1.1  2004/03/16 07:01:14  pxs565
!!! Initial revision
!!!


Program chainsim_p
  Use bspglocons
  Use bspglovars
  Use bspintfacs
  Use Flowvars
  Use Flock_utils
  Implicit None


  !_____________________________________________________________
  !        Other Declarations                                   
  !_____________________________________________________________


  Character (len=12) clk(3)
  Integer (k4b) nseed

  ! File i/o
  Character(10), parameter :: FormatVersion = "GAVG-1.0"
  !Character (10) :: fver
  Character (len=20) infile, outfile, xdfile, gavgfile, ntfile, contfile,fver
  
  Integer, Parameter ::   inunit=10, outunit=11, gavunit=12, &
       resunit=15, cunit=16, inphi=17
  Logical :: Filexists, Flocked

  ! for gfortran
  Character (len=80) :: Format3, Format31, Format4, Format42, &
       Format43, Format5, Format51, Format6

  ! Trajectories related
  Integer, Parameter :: MaxNDT = 10
  Integer i, j, clok(8), k
  Integer nthis, nblock, iblock, ntrajdone, ntraj, Nsamples, &
       idelts, ndelts,  ntrajvals(MaxNDT), nthisvals(MaxNDT)

  Real (DBprec) t1zimm, emax, tol(MaxNDT)
  Real (DBprec) deltseq, deltsne, deltseqvals(10), tlongest, &
            deltsvals(MaxNDT), tmax, teqbm, t1rouse 
            

  Real (DBprec), Allocatable :: PosVecR(:,:), times(:), samples(:,:), &
       avgs(:,:), errs(:,:), &
       global_avgs(:,:), global_errs(:,:) 
  !!!!! phi parameter declared here
  Real (DBprec), Allocatable :: phi(:,:)


  Integer SType,mmult, nsact, cidx, ord, eqprops, neqprops

  Real (DBprec) hstar, zstar, dstar, sqrtb, Q0s
  Real (DBprec) binwidth, intCss, csserr, l1rouse
  Character (len=10), parameter :: ErString = "Error"
  Real (DBprec) rems, reerr, rgms, rgerr, xms, xerr, dfvty, derr, ntot, Conf_t


  Integer :: NBeads
  Real (DBprec) :: Delts

  !_____________________________________________________________
  !        Get input data                                       c
  !_____________________________________________________________

  infile="inputc.dat"


  Open (unit=inunit,file=infile,status="old")
  Read (inunit,*)
  Read (inunit,*) SType, NBeads, hstar, zstar, dstar, sqrtb, Q0s, gdots
  Read (inunit,*)
  Read (inunit,*) emax, Nsamples
  Read (inunit,*)
  Read (inunit,*) ndelts
  Read (inunit,*)

  If (ndelts > MaxNDT) Then
     Write(*,*) 'Number of delt exceeded ', MaxNDT
     Stop
  End If

  Do i = 1, ndelts
     Read(inunit,*) deltseqvals(i), deltsvals(i), nthisvals(i), &
          ntrajvals(i), tol(i)
  End Do
  Close (unit=inunit)
   

  !_____________________________________________________________
  !        Initialization                                       c
  !_____________________________________________________________
  If (Nsamples.Lt.2) Nsamples = 2

  Allocate(PosVecR(Ndim,NBeads))
  Allocate ( &
       times(Nsamples),  &
       samples(NProps,Nsamples), &
       avgs(NProps,Nsamples), &
       errs(NProps,Nsamples))

  Allocate(&
       global_avgs(NProps,Nsamples), &
       global_errs(NProps,Nsamples))
  
 Allocate(phi(NBeads,NBeads))

  global_avgs = 0
  global_errs = 0

    Open (unit=inphi,file="phi.dat",status="old")
    Do k = 1, NBeads
       Read (inphi,*) phi(k,:)
    End Do
   write (*,*) "phi", phi (1,:)
   write (*,*) "Q0 value", Q0s
  !gdots = 0

  times = 0

  ! longest relaxation times Rouse and Zimm
  t1rouse =   0.5/Sin(PI/2/Nbeads)**2
  t1zimm = lam1_th(hstar,NBeads) ! use thurstons formula for Zimm

  tlongest = t1rouse

  teqbm = tlongest ! without EV

  If (zstar .Ne. 0) Then
     teqbm = 30.0*tlongest ! with excluded volume, it takes roughly 3 times
     ! to attain equilibrium even with init dist
  End If

  If (gdots == 0) Then
     tmax = emax
  Else
     tmax = emax/gdots
  End If
     tmax =0.1*t1rouse !!!! kiran 
   write (*,*) "trouse time is", t1rouse, "tmax is", tmax
  !_____________________________________________________________
  !        Initialization variable format expressions
  !         <> language extension not available in gfortran
  !_____________________________________________________________

  Write(Format3,"(a,I3,a)") "('#',", 2*NProps+4, "(A11,1X))"
  ! short cut way to get (nearly) left justification
  Write(Format31,"(a,I3,a)"),  "('#',", 2*NProps+4,"(G2.0,10X))"



  If (gdots == 0) Then
     eqprops = 4
  Else
     neqprops = 4
  End if
     
!!$  Write(Format4 ,"(a,I3,a)") "(", 2*NProps+4,"(G11.4,1X))"
  Write(Format42,"(a,I3,a)") "(", 2*eqprops+1, "(G22.17,1X))"
  Write(Format43,"(a,I3,a)") "(", 2*neqprops+2, "(G11.6,1X))"

!!$  Write(Format5, "(a,I3,a)") "('#',A11,1X,", nsact, "(G11.4,2X))"
!!$  Write(Format51, "(a,I3,a)") "('#',", nsact+1, "(I3,10X))"
!!$
!!$  Write(Format6, "(a,I3,a)") "(", nsact+1, "(G11.4,2X))"


  !_____________________________________________________________
  !        Initialization of variables 
  !_____________________________________________________________


  ! the standard error of mean obtained from t-distribution
  ! for sample mean and sample standard deviation
  Conf_t = 2 ! 95% confidence for degrees of freedom > 20



  Prop_names(1) = "R^2"         ! square of end-to-end vector          
  Prop_names(2) = "S11"         ! 1,1 component of shape tensor       
  Prop_names(3) = "S12"         ! 1,2 component of shape tensor       
  Prop_names(4) = "S13"         ! 1,3 component of shape tensor       
  Prop_names(5) = "S22"         ! 2,2 component of shape tensor       
  Prop_names(6) = "S23"         ! 2,3 component of shape tensor       
  Prop_names(7) = "S33"         ! 3,3 component of shape tensor       
  Prop_names(8) = "N1"          ! First normal stress difference      
  Prop_names(9) = "N2"          ! Second normal stress difference     
  Prop_names(10) = "T12"        ! 1,2 component of stress tensor      
  Prop_names(11) = "X1"         ! Stretch in 1 direction              
  Prop_names(12) = "X2"         ! stretch in 2 direction              
  Prop_names(13) = "X3"         ! stretch in 3 direction              
  Prop_names(14) = "Rg^2"       ! sq. Radius of gyration
  Prop_names(15) = "PC error"   ! Predictor-Corrector Error           
  Prop_names(16) = "<SxySxy>"   ! Auto correlation fn for shear stress
  Prop_names(17) = "Diffusvty"  ! Center of mass diffusivity


  If (gdots == 0 ) Then
     Open (unit = resunit, file = "result.dat", status="unknown")
     Write (resunit,200) "NBeads      ", &
          "sqrtb       ", &
          "h*          ", &
          "z*          ", &
          "d*          ", &
          "dtsne       ", &
          "<R^2>       ", &
          "<Rg^2>      ", &
          "Diffsvty    ", &
          "l_eta       ", &
          "<x>         " 
     Write (resunit,202) (i, i=1,16)
200  Format ('#',6(G12.0,2X),5(G12.0,2X, 'Error       ',2X))
202  Format ('#',16(I2,12X)) !approx LJustification
  end If



  !_____________________________________________________________
  !     Generate the seeds based on the current time
  !_____________________________________________________________

  Call Date_and_time(clk(1), clk(2), clk(3), clok) 
  !       ms              sec         min        hr
  nseed = clok(8)*100000+clok(7)*1000+clok(6)*10+clok(5)

  !-- an unique seed incase clok is the same
  nseed = (nseed + 201271)
  nseed = 201271


  !_____________________________________________________________
  !    Begin the loop for time step sizes                       c
  !_____________________________________________________________

  timesteps: Do  idelts = 1,ndelts

     deltseq = deltseqvals(idelts)
     deltsne = deltsvals(idelts)
     Imploop_tol = tol(idelts)

     mmult = Int(tmax/(Nsamples-1)/deltsne) + 1

     ! actual no of samples to b taken
     nsact = tmax/deltsne/mmult + 1
    ! nsact = Nsamples
     Do i = 1,nsact
        times(i) = (i-1)*mmult*deltsne
     End Do

     ! Obtain the number of traj completed from the disk, if present

     Write (gavgfile, '("gavgs.",I2.2)') idelts 

     Inquire (file=gavgfile, exist=Filexists)

     If (Filexists) Then
        open (unit=gavunit,file=gavgfile,status='old')

        Read (gavunit,*) fver  ! Format version
        if (fver /= FormatVersion) Then
           Write (*,*) 'Incompatible Format version in ', gavgfile
           close(gavunit)
           Go to 99999
        end if

        Read (gavunit,*) ntrajdone
        close(gavunit)


     Else
        ntrajdone = 0
     End If
     
     !--when several procs start at the same time it
     !  can lead to the total trajs done being > ntrajvals
     ntraj = Max(0,ntrajvals(idelts)-ntrajdone)
     
     !--number of trajectories for this run
     nthis = Min(ntraj, nthisvals(idelts))

     nblock = nthis
    ! nsact = Nsamples


     avgs = 0.d0
     errs = 0.d0

     !_____________________________________________________________
     !    Begin the loop for the blocks                            c
     !_____________________________________________________________

     trajectories: Do iblock = 1, nblock
        samples = 0.d0
        
        PosVecR = 0.d0
        Call Initial_position(SType,Nbeads,sqrtb,PosVecR,nseed)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!


        FlowType = EQ
        Call Time_Integrate_Chain(NBeads, PosVecR, SType,  &
             0._DBprec,tmax, deltseq, &
             hstar , zstar, dstar, sqrtb, Q0s,  &
             nseed, nsact, times, samples, phi )

        ! If doing equilibrium studies, use a large deltseq=1 first
        ! and later use deltsne for gathering data
       ! If (zstar .Ne. 0) Then
        !   Call Time_Integrate_Chain(NBeads, PosVecR, SType,  &
         !       0., 4*tlongest, deltsne, &
          !      0. , zstar, dstar, sqrtb, Q0s, &
           !     nseed, 1, times, samples )
       ! End If


       ! FowType = SH
       ! Call Time_Integrate_Chain(NBeads, PosVecR, SType,  &
       !      0., tmax, deltsne, &
       !      0., zstar, dstar, sqrtb, Q0s, &
       !      nseed, nsact, times, samples , phi)


        avgs = avgs + samples
        errs = errs + samples*samples
     End Do trajectories


     !_____________________________________________________________
     !    Consolidate and save final results                       c
     !_____________________________________________________________

     ! obtain global values from the disk
     Inquire (file=gavgfile, exist=Filexists)
     If (Filexists) Then
        !-- check if it is locked
        Flocktest: do 
           call islocked(gavgfile,Flocked)
           if(.not. Flocked) Then
              call lockfile(gavunit,gavgfile) 
              !-- to be unlocked after new data is written

              open (unit=gavunit,file=gavgfile,status='old')
              Read (gavunit,*) fver  ! Format version
              if (fver /= FormatVersion) Then
                 Write (*,*) 'Incompatible Format version in ', gavgfile
                 call unlockfile(gavunit,gavgfile)
                 Go to 99999
              end if
           
              Read (gavunit,*) ntrajdone
              Read (gavunit,*) global_avgs(1,:) , global_errs(1,:)
              Read (gavunit,*) global_avgs(8,:) , global_errs(8,:)
              Read (gavunit,*) global_avgs(10,:), global_errs(10,:)
              Read (gavunit,*) global_avgs(11,:), global_errs(11,:)
              Read (gavunit,*) global_avgs(14,:), global_errs(14,:)
           
              Close (gavunit)
              Exit Flocktest
           Else  ! when File is locked
              ! Sleep for a while and try again until it is unlocked,  
              ! requires -Vaxlib in ifc
              call sleep(1)
           end if
        end do Flocktest
     Else ! if Globalavgs does .not. Filexists

        ! lock the file b4 opening a new one for writing
        call lockfile(gavunit,gavgfile) ! tobe unlocked after new data is written
        global_avgs = 0
        global_errs = 0
        ntrajdone = 0
        
     end If
  

     !-- add the current runs to that in the disk
     global_avgs = global_avgs + avgs
     global_errs = global_errs + errs
     ntrajdone = ntrajdone + nblock

     ! take the sum-total of time-ensemble averages for 
     ! steady equilibrium measurments before averaging
     ! over ensembles

     ntot = nsact * ntrajdone
     rems  = Sum(global_avgs(1,1:nsact))/ntot
     reerr = Conf_t *  &
          ((Sum(global_errs(1,1:nsact)) - &
          ntot * rems*rems)/ntot/(ntot-1))**0.5

     rgms  = Sum(global_avgs(14,1:nsact))/ntot
     rgerr = Conf_t *  &
          ((Sum(global_errs(14,1:nsact)) - & 
          ntot * rgms*rgms )/ntot/(ntot-1))**0.5

     xms  = Sum(global_avgs(11,1:nsact))/ntot
     xerr = Conf_t *  &
          ((Sum(global_errs(11,1:nsact)) - & 
          ntot * xms*xms )/ntot/(ntot-1))**0.5


     ntot = (nsact-1)*ntrajdone

     dfvty  = Sum(global_avgs(17,2:nsact))/ntot
     derr = Conf_t *  &
          ((Sum(global_errs(17,2:nsact)) - & 
          ntot * dfvty*dfvty)/ntot/(ntot-1))**0.5

     ! average and errors over the ensemble
     global_avgs = global_avgs/ntrajdone
     If (ntrajdone.Gt.2) Then
        global_errs = Conf_t*(Abs(global_errs-ntrajdone*global_avgs* &
             global_avgs)/ntrajdone/(ntrajdone-1))**0.5
     End If

     Write (outfile, '("output.",I2.2)') idelts 
     Open (unit = outunit, file = outfile, status="unknown")
     Write (outunit,1) "SpringTyp ", "NBeads    ", "sqrtb     ", &
          "h*        ", "z*        ", "d*        ", &
          "gdot*     ", "dt*eq     ", "dt*ne     ", &
          "imp-tol   ", "TLongest  ", "TEqbm     ", &
          "eMax      ", "Tmax      ", "NTraj     " ,&
          "FlowType  "
     Write (outunit,2) SType, NBeads, sqrtb, hstar, zstar, dstar, &
          gdots, deltseq, deltsne, Imploop_tol, tlongest, teqbm, &
          emax, tmax, ntrajdone, FlowType  

     If (gdots .Eq. 0) Then
        eqprops = 4
        Write (outunit,Format3) "Time       ",  &
             Prop_names(1), Erstring, &
             Prop_names(14), Erstring, &
             Prop_names(17), Erstring, &
             Prop_names(11), Erstring
        Write (outunit,Format31) (i, i=1,1+eqprops*2)
        Write (outunit,Format42) (times(i), &
             global_avgs(1,i), global_errs(1,i), & 
             global_avgs(14,i), global_errs(14,i), & 
             global_avgs(17,i), global_errs(17,i), & 
             global_avgs(11,i), global_errs(11,i), & 
             i = 1,nsact )
      
     Else
        neqprops = 4
        Write (outunit,Format3) &
             "Time      ",  &
             "Strain    ",  &
             Prop_names(1), Erstring, &
             Prop_names(14), Erstring, &
             "psi1     ", Erstring, &
             "etap     ", Erstring 
        Write (outunit,Format31) (i, i=1,2+neqprops*2)
        Write (outunit,Format43) (times(i), times(i)*gdots,  &
             global_avgs(1,i), global_errs(1,i), & 
             global_avgs(14,i), global_errs(14,i), & 
             -global_avgs(8,i)/gdots/gdots, global_errs(8,i)/gdots/gdots, &
             -global_avgs(10,i)/gdots, global_errs(10,i)/gdots, & 
             i = 1,nsact )
     End If


     If (gdots .Eq. 0) Then
        ! The integral of the ensenble averaged Stress-stress correlation 
        ! function gives the intrinsic viscosity at zero shear rate
        cidx = 16

        Do ord=1,3
           Call numint(global_avgs(cidx,:), deltsne*mmult, 1, nsact, &
                ord, intCss)
           Call numint(global_errs(cidx,:), deltsne*mmult, 1, nsact, &
                ord, csserr)
           Write (outunit,101)  ord , intCss, csserr + 1./nsact**ord
101        Format ('# zero sh rate intrinsic visc, ord = ',I2, ':', &
                F10.6,' +/- ',F10.6)
        End Do

        ! lambda_eta from Rouse model
        l1rouse = 0
        Do i=1,Nbeads-1
           l1rouse = l1rouse +  0.5/Sin(i*PI/2/Nbeads)**2
        End Do
        Write (outunit,*) '# from Rouse model = ', l1rouse
        Write (outunit,219) t1rouse,t1zimm
219     Format('#longest Rouse = ',G10.3, ', Zimm = ',G10.3)

        Write (resunit,210) NBeads, sqrtb, hstar, zstar, dstar, deltsne, &
             rems, reerr, rgms, rgerr, dfvty, derr, intCss, csserr, xms, xerr
210     Format (I3,  11X, 16(G12.5,2X))

     end If ! gdots = 0

     Close (unit = outunit)


     !_____________________________________________________________
     ! Save the global averages data on to the disk
     ! But before that revert to the unnormalised values
     ! Note that the gavunit file is still locked
     !_____________________________________________________________
     Open(unit=gavunit,file=gavgfile,status='replace')

     If (ntrajdone.Gt.2) Then
        global_errs = (global_errs/Conf_t)**2 * ntrajdone * (ntrajdone-1)&
             + ntrajdone * global_avgs**2
     end If
     
     global_avgs = global_avgs * ntrajdone

     Write (gavunit,*) FormatVersion
     Write (gavunit,*) ntrajdone
     Write (gavunit,*) global_avgs(1,:) , global_errs(1,:)
     Write (gavunit,*) global_avgs(8,:) , global_errs(8,:)
     Write (gavunit,*) global_avgs(10,:), global_errs(10,:)
     Write (gavunit,*) global_avgs(11,:), global_errs(11,:)
     Write (gavunit,*) global_avgs(14,:), global_errs(14,:)
     Close (gavunit)
     
     !-- global avgs file is released for use with other processors
     call unlockfile(gavunit,gavgfile)


     DeAllocate( global_avgs, global_errs)
     
     !_____________________________________________________________
     !           Write info to disk for continuation
     !_____________________________________________________________
     ntraj = ntrajvals(idelts)-ntrajdone

     Write (contfile, '("continue.",I2.2)') idelts 
     open(unit=cunit,file=contfile,status='replace')

     if (ntraj > 0) Then
        Write (cunit,*) ntraj, '  more trajectories to be completed'
        close(cunit)
     Else 
        close(cunit,status='delete')
     end if
     

  End Do timesteps

  Close(resunit)


  ! Common termination statements
99999   Stop
  !_____________________________________________________________
  !                    Format statements                        c
  !_____________________________________________________________

1  Format ('#',16(A10,2X))
  !2  Format ('#', 2(I10,2X), 12(G10.3,2X), (I10,2X)  )
2  Format ('#', 2(G2.0,10X), 12(G11.4,1X), (G10.3,2X), G2.0,11X  )

!3 Format ('#', <2*NProps+4>(A11,1X))
!3  Format (Format3)
  ! short cut way to get (nearly) left justification
!31 Format ('#', <2*NProps+4>(G2.0,10X))
!31 Format (Format31)

!4 Format (<2*NProps+4>(G11.4,1X))
!42 Format (<2*eqprops+1>(G11.4,1X))
!43 Format (<2*neqprops+2>(G11.4,1X))
!4  Format (Format4)
!42 Format (Format42)
!43 Format (Format43)

!5 Format ('#',A11,1X,<nsact>(G11.4,2X))
!51 Format ('#', <nsact+1>(I3,10X))     
!5  Format (Format5)
!51 Format (Format51)

!6 Format (<nsact+1>(G11.4,2X))     
!6  Format (Format6)

72 Format (5(F10.4,2x))
80 Format (3(F10.4,2x))

End Program chainsim_p


