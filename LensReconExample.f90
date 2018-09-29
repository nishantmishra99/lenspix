!Simple example code to simulate a lensed CMB map, then use it to do temperature quadratic lensing reconstruction
!Not at all realistic (no masks, not even any simulated noise, only isotropic noise in the filter weights)
!AL Apr 2014


!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! Module for the quadratic estimator

    module Recon
    use HealpixObj
    use HealpixVis
    use Random
    use spinalm_tools
    use IniFile
    use AMLUtils
    implicit none

    contains


! initialise the noise power spectrum
    subroutine NoiseInit(Noise, sensitivity, noise_fwhm, lmax)
    Type(HealpixPower):: Noise
    real(dp) sensitivity, noise_fwhm
    real(sp) amp, xlc, sigma2
    integer lmax, L

    call healpixPower_Init(Noise,lmax,.false.)

    ! convert from muK*arcmin to muK*rad
    sensitivity = sensitivity * (3.14159/180.)/60.   ! convert to muK*rad

    ! convert from beam fwhm in arcmin to sigma in rad
    noise_fwhm = noise_fwhm * (3.14159/180.) / 60.   ! convert to rad
    sigma2 = noise_fwhm / sqrt(8.*log(2.))   ! convert from fwhm to sigma
    sigma2 = sigma2**2

    do l=0, lmax
        Noise%Cl(l,1) = sensitivity**2 * exp(l*(l+1)*sigma2)
    end do

    end subroutine NoiseInit



! Compute the lensing noise power spectrum for phi,
! ie the normalization of the phi quadratic estimator.
    subroutine GetA_L(P, Noise, AL, lmax_phi, lmax_est)
    Type(HealpixPower) :: P, Noise
    real(dp) AL(:)
    integer lmax_est, lmax_phi
    integer l1,l2,L
    real(dp), allocatable :: threejs(:)

    allocate(threejs(0:lmax_est + lmax_phi+1))
    AL = 0
    do L=1, lmax_phi
        do l1=2, lmax_est
            call GetThreeJs(threejs(abs(l1-L)),l1,L,0,0)
            do l2 = max(2,abs(l1-L)), min(l1+L,lmax_est)
                if (threejs(l2)==0) cycle
                AL(L) = AL(L)+ (threejs(l2)**2*real((2*l1+1)*(2*l2+1),dp)* &
                ( (L*(L+1) + l2*(l2+1)-l1*(l1+1))*P%Cl(l2,1)+ (L*(L+1) - l2*(l2+1)+l1*(l1+1))*P%Cl(l1,1))**2) / &
                ( 2*(P%Cl(l1,1)+Noise%Cl(l1,1))*(P%Cl(l2,1)+Noise%Cl(l2,1)))
            end do
        end do
        AL(L) = 4*HO_fourpi/ AL(L) *L*(L+1)
        print *, L, AL(L)
    end do

    end subroutine GetA_L



! Non-normalized quadratic estimator for phi
    subroutine QuadraticPhi(H,A, MGradT, LensedCl, Noise, lminest,lmaxest)
    Type(HealpixInfo)  :: H
    Type(HealpixPower) :: LensedCl, Noise
    Type(HealpixAlm) :: A, AOut , PhiSol
    Type(HealpixMap) :: MFilt, MGradT
    integer l,lmaxest, lminest

    print *,'get quadratic', lminest, lmaxest

    ! Create inverse-variance weighted map
    call HealpixAlm_Init(AOut, lmaxest,1)
    AOut%TEB=0
    do l=max(2,lminest),lmaxest
        AOut%TEB(1,l,:) = A%TEB(1,l,:) / (LensedCl%Cl(l,1) + Noise%Cl(l,1))
    end do
    call HealpixAlm2Map(H, AOut,MFilt,H%npix)

    ! Create Wiener-filtered gradient
    do l=max(2,lminest),lmaxest
        AOut%TEB(1,l,:) = Aout%TEB(1,l,:)  * LensedCl%Cl(l,1)
    end do
    call HealpixAlm2GradientMap(H, AOut,MGradT,H%npix,'T')

    ! Do the convolution
    MGradT%SpinField =  MGradT%SpinField * MFilt%TQU(:,1)
    call HealpixMap_Free(Mfilt)
    call HealpixAlm_Free(AOut)

    end subroutine QuadraticPhi

    end module Recon



!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------


    program SimLensReconTT
    use HealpixObj
    use HealpixVis
    use Random
    use spinalm_tools
    use IniFile
    use AMLUtils
    use Recon
    implicit none
    Type(HealpixInfo)  :: H
    Type(HealpixMap)   :: M, GradPhi, MGradT
    Type(HealpixPower) :: UnlensedCl, LensedCl, Noise, P
    Type(HealpixAlm)   :: A,SimAlm,PhiRecon

    integer            :: nside, lmax
    integer(I_NPIX)    :: npix
    character(LEN=1024)  :: w8name = '../Healpix_2.00/data/'
    character(LEN=1024)  :: file_stem, cls_file, out_file_root, cls_lensed_file
    character(LEN=1024) :: healpixloc, aname, in_map
    integer, parameter :: lens_interp =1, lens_exact = 2
    integer :: lens_method = lens_interp
    integer :: mpi_division_method = division_equalrows
    integer ::   rand_seed
    logical :: err, want_pol
    real :: interp_factor
    integer status, i, L
    integer lmax_phi, lmax_est
    real(dp), allocatable :: AL(:)
    real(dp) sensitivity, noise_fwhm
#ifdef MPIPIX

    call mpi_init(i)
#endif

    Ini_Fail_On_Not_Found = .true.
    call Ini_Open(GetParam(1), 3,err)
    if (err) then
#ifdef MPIPIX
        call mpi_finalize(i)
#endif
        stop 'No ini'
    end if



!-----------------------------------------------------------------------------
! Read all parameters from .ini file

    nside  = Ini_Read_Int('nside')
    npix = nside2npix(nside)

    lmax   = Ini_Read_Int('lmax')
    cls_file = Ini_Read_String('cls_file')
    cls_lensed_file = Ini_Read_String('cls_lensed_file')
    out_file_root = Ini_Read_String('out_file_root')

    want_pol = Ini_Read_Logical('want_pol')
    rand_seed = Ini_Read_Int('rand_seed')

    sensitivity = Ini_read_real('sensitivity')  ! in muK*arcmin
    noise_fwhm = Ini_read_real('noise_fwhm') ! in arcmin
    call NoiseInit(Noise, sensitivity, noise_fwhm, lmax)

    !unlensed_power_spectrum_file




    ! What are these lmax exactly?
    !make sure to check this ins the params.ini
    !rewrite
    !
    lmax_phi = Ini_read_int('lmax_phi')
    lmax_est = Ini_read_int('lmax_est',lmax_phi)

    Ini_Fail_On_Not_Found = .false.

    ! unlensed temperature map to be lensed
    in_map = Ini_read_String('input_map')
    if (in_map=='') in_map = 'lensed_sim_cache.fits'


    w8name = Ini_Read_String('w8dir')
    interp_factor=0
    if (lens_method == lens_interp) interp_factor = Ini_Read_Real('interp_factor',3.)
#ifdef MPIPIX
    mpi_division_method = Ini_Read_Int('mpi_division_method',division_balanced);
#endif

    call Ini_Close


!-----------------------------------------------------------------------------
! Create file paths

    file_stem =  trim(out_file_root)//'_lmax'//trim(IntToStr(lmax))//'_nside'//trim(IntTOStr(nside))// &
    '_interp'//trim(RealToStr(interp_factor,3))

    if (want_pol) file_stem=trim(file_stem)//'pol_'
    file_stem = trim(file_stem)//trim(IntToStr(lens_method))

    call SetIdlePriority()


!-----------------------------------------------------------------------------
! Initialize healpix

    if (w8name=='') then
        call get_environment_variable('HEALPIX', healpixloc, status=status)
        if (status==0) then
            w8name = trim(healpixloc)//'/data/'
        end if
    end if

    if (w8name=='') then
        write (*,*) 'Warning: using unit weights as no w8dir found'
        call HealpixInit(H,nside, max(lmax, lmax_phi*2),.true., w8dir='', method= mpi_division_method)
    else
        call HealpixInit(H,nside, max(lmax, lmax_phi*2),.true., w8dir=w8name,method=mpi_division_method)
    end if

!-----------------------------------------------------------------------------

    !All but main thread stay in HealpixInit
    if (H%MpiID ==0) then !if we are main thread

        ! Create empty maps
        ! Rename maps for convenience
        call HealpyAlm_nullify(lensed_Alm)
        call HealpixMap_nullify(Output_lensing_potential_map)
        call HealpixMap_nullify(Output_unlensed_map)

        !call HealpixAlm_nullify(A)
        !call HealpixMap_nullify(GradPhi)
        !call HealpixMap_nullify(M)

        !Read power spectrum of unlensed T, E, B and of phi
        call HealpixPower_ReadFromTextFile(Unlensed_power_spectrum_input, unlensed_power_spectrum_file, lmax, pol=.false., dolens = .false.)
        !call HealpixPower_ReadFromTextFile(UnlensedCl,cls_file,lmax,pol=.true.,dolens = .true.)

        ! Read power spectrum of lensed T, E, B (but not phi)
        call HealpixPower_ReadFromTextFile(Lensed_power_spectrum_input, lensed_power_spectrum_file, lmax, pol=.false., dolens = .false.)
        !call HealpixPower_ReadFromTextFile(LensedCl,cls_lensed_file,lmax,pol=.true.)

        ! Print power spectra to check
        ! rename or remove
        print *,'at L=1500, Lens, Unlens, noise C_l =', Unlensed_power_spectrum_input%Cl(1500,1), Lensed_power_spectrum_input%Cl(1500,1), Noise%Cl(1500,1)
        !print *,'at L=1500, Lens, Unlens, noise C_l =', LensedCl%Cl(1500,1), UnlensedCl%Cl(1500,1), Noise%Cl(1500,1)

        ! Compute normalization for quadratic estimator,
        ! and write to file
        allocate(AL(lmax_phi))
        aname=trim(out_file_root)//'_AL.txt'
        if (.not. FileExists(aname)) then
            call GetA_L(Lensed_power_spectrum_input, Noise, AL, lmax_phi, lmax)
            !call GetA_L(LensedCl, Noise, AL, lmax_phi, lmax)
            call CreateTxtFile(aname,1)
            do i=1, lmax_phi
                write(1,*) i, AL(i)
            end do
            close(1)
        else
            call OpenTxtFile(aname,1)
            do i=1, lmax_phi
                read(1,*) L, AL(i)
            end do
            close(1)
        end if


!-----------------------------------------------------------------------------
! Get the lensed temperature map
! This part doesn't matter because it uses a power spectrum to generate a map when a we are using a map input
        ! If an input lensed map is provided, read it
        if (FileExists(in_map)) then
            print *, 'reading input map', trim(in_map)
            !call HealpixMap_Read(M,in_map)
            call HealpixMap_Read(InputLensedMap, input_lensed_map)
            call HealpixAlm_Nullify(SimAlm)


        ! Otherwise, create a mock lensed map
        else
            ! Generate the GRF unlensed map, and the GRF phi map
            ! write them to SimAlm
            call HealpixAlm_Sim(SimAlm, UnlensedCl, rand_seed, HasPhi=.true., dopol = want_pol)
            ! Measure the unlensed power spectrum, write it to P
            call HealpixAlm2Power(SimAlm,P)
            ! Write the unlensed power spectrum to file
            call HealpixPower_Write(P,trim(file_stem)//'_unlensed_simulated.dat')

            ! Compute the unlensed gradient, save it in "gradPhi"
            call HealpixAlm2GradientMap(H,SimAlm, GradPhi,H%npix,'PHI')
            ! Create lensed temperature map M,
            ! from SimAlm which contains unlensed T and phi,
            ! and from gradPhi which contains the unlensed gradient
            call HealpixInterpLensedMap_GradPhi(H,SimAlm,GradPhi, M, interp_factor, interp_cyl)
            ! Write lensed temperature map M to file
            call HealpixMap_Write(M, in_map)
        end if



!-----------------------------------------------------------------------------
!Test very simple TT reconstruction
!Note no noise added to map, N0=A_L will not be correct noise bias

        ! Get the lensed alm A from the lensed map M
        !call HealpixMap2Alm(H,M, A,lmax)
        call HealpixMap2Alm(H, InputLensedMap, Lensed_Alm, lmax)
        ! Compute non-normalized quadratic estimator for phi,
        ! from lensed alm in A,
        ! write it to MGradT
        call QuadraticPhi(H, Lensed_Alm, Non_normalized_QE, Lensed_power_spectrum_input, Noise, 2, lmax)
        !call QuadraticPhi(H,A, MGradT, LensedCl, Noise, 2,lmax)
        ! convert this map MGradT into alm, written to A
        call HealpixMap2Alm(H, Non-normalized_QE, Lensed_Alm, lmax_est)
        !call HealpixMap2Alm(H, MGradT, A,lmax_est)

        ! compute normalized quadratic estimator for phi,
        ! write it to PhiRecon
        call HealpixAlm_Init(PhiRecon,lmax_phi,npol=0,HasPhi=.true.)
        !call HealpixAlm_Init(PhiRecon,lmax_phi,npol=0,HasPhi=.true.)
        do i=1,lmax_phi
            !PhiRecon
            PhiRecon%Phi(1,i,:) =  A%SpinEB(1,i,:) * AL(i) / sqrt(i*(i+1.))
        end do
        call HealpixAlm2Map(H, PhiRecon ,Output_unlensed_map, nside, DoPhi=.true.)

        ! Compute the power spectrum of the normalized phi quadratic estimator
!        call HealpixAlm2Power(PhiRecon, P)
        ! Write it to file
!        call HealpixPower_write(P, trim(file_stem)//'recon_power.dat')
        ! If the input phi map is known, compute the cross spectrum,
        ! and write it to file
!        if (associated(SimAlm%Phi)) then
!            call HealpixAlm2CrossPhi(PhiRecon, SimAlm, P)
!            call HealpixPower_write(P, trim(file_stem)//'recon_cross_power.dat')
        end if
    end if

#ifdef MPIPIX
    call HealpixFree(H)
    call mpi_finalize(i)
#endif

#ifdef DEBUG
    write (*,*) 'End of program'
    pause
#endif
end program SimLensReconTT
