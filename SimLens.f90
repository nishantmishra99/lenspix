    !Simple program demonstrating how to generate a simulated lensed map
    !AL, Feb 2004; Updated Oct 2007
    program SimLensCMB
    use HealpixObj
    use HealpixVis
    use Random
    use spinalm_tools
    use IniFile
    use AMLUtils
    implicit none
    Type(HealpixInfo)  :: H
    Type(HealpixMap)   :: Unlensed_CMB, Phi, GradPhi, lensed_CMB
    Type(HealpixPower) :: P
    Type(HealpixAlm)   :: Alm_Unlensed_CMB, Alm_Phi

    integer            :: nside, lmax
    integer(I_NPIX)    :: npix
    character(LEN=1024)  :: w8name = '../Healpix_2.00/data/'
    character(LEN=1024)  :: file_stem, cls_file, out_file_root, cls_lensed_file, input_map1, input_lensing, output_map
    character(LEN=1024) :: healpixloc
    integer, parameter :: lens_interp =1, lens_exact = 2
    integer :: lens_method = lens_interp
    integer :: mpi_division_method = division_equalrows
    integer ::  interp_method,  rand_seed, interp_algo
    logical :: err, want_pol
    real :: interp_factor
    integer status
#ifdef MPIPIX
    integer i

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
! Read parameters from input file

    !read the map input maps
    input_map1 = Ini_Read_String('input_map1')
    input_lensing = Ini_Read_String('input_lensing')
    output_map = Ini_Read_String('output_map')


    nside  = Ini_Read_Int('nside')
    npix = nside2npix(nside)

    ! lmax that is computed. Should be larger than desired lmax by +250
    lmax   = Ini_Read_Int('lmax')

    ! power spectrum of unlensed CMB
    cls_file = Ini_Read_String('cls_file')

    out_file_root = Ini_Read_String('out_file_root')

    lens_method = Ini_Read_Int('lens_method')
    want_pol = Ini_Read_Logical('want_pol')
    rand_seed = Ini_Read_Int('rand_seed')

    interp_method = Ini_read_int('interp_method')
    interp_algo = Ini_read_int('interp_algo')

    Ini_Fail_On_Not_Found = .false.

    w8name = Ini_Read_String('w8dir')
    interp_factor=0
    if (lens_method == lens_interp) interp_factor = Ini_Read_Real('interp_factor',3.)
#ifdef MPIPIX
    mpi_division_method = Ini_Read_Int('mpi_division_method',division_balanced);
#endif

    ! close the input file
    call Ini_Close


!-----------------------------------------------------------------------------
! Define the output file name
! for the lensed CMB power spectrum to be computed below

    if (interp_algo == 1) then
        file_stem =  trim(out_file_root)//'_lmax'//trim(IntToStr(lmax))//'_nside'//trim(IntTOStr(nside))// &
        '_interp'//trim(RealToStr(interp_factor,3))//'_method'//trim(IntToStr(interp_method))//'_'
    else
        file_stem =  trim(out_file_root)//'_lmax'//trim(IntToStr(lmax))//'_nside'//trim(IntTOStr(nside))// &
        '_interp'//trim(RealToStr(interp_factor,3))//'_method'//trim(IntToStr(interp_method))//'_algo'//&
        trim(IntToStr(interp_algo))//'_'
    endif

    if (want_pol) file_stem=trim(file_stem)//'pol_'
    file_stem = trim(file_stem)//trim(IntToStr(lens_method))

    cls_lensed_file  = trim(file_stem)//'.dat'

    call SetIdlePriority()


!-----------------------------------------------------------------------------
! Initializing healpix

    if (w8name=='') then
        call get_environment_variable('HEALPIX', healpixloc, status=status)
        if (status==0) then
            w8name = trim(healpixloc)//'/data/'
        end if
    end if

    if (w8name=='') then
        write (*,*) 'Warning: using unit weights as no w8dir found'
        call HealpixInit(H,nside, lmax,.true., w8dir='', method= mpi_division_method)
    else
        call HealpixInit(H,nside, lmax,.true., w8dir=w8name,method=mpi_division_method)
    end if


!-----------------------------------------------------------------------------
! Compute the lensed map and power spectrum


    ! All but the main thread stay in HealpixInit
    ! if we are the main thread
    if (H%MpiID ==0) then
      ! create empty full sky maps
      call HealpixMap_nullify(Unlensed_CMB) !map object
      call HealpixMap_nullify(Phi) !map object
      ! create empty alm for above full sky maps
      call HealpixAlm_nullify(Alm_Unlensed_CMB) !alm object
      call HealpixAlm_nullify(Alm_Phi) !alm object
      ! create empty map for the GradPhi
      call HealpixMap_nullify(GradPhi) !map object
      ! create a final map to put the lensed map on
      call HealpixMap_nullify(lensed_CMB)

      ! Read the .fits map onto the Healpix map object
      call HealpixMap_Read(Unlensed_CMB, input_map1)

      ! Read the .fits map for phi onto the Healpix map object
      call HealpixMap_Read(Phi, input_lensing)

      ! Get unlensed alm from unlensed map
      call HealpixMap2alm(H, Unlensed_CMB, Alm_Unlensed_CMB, lmax, dopol=want_pol)

      ! Get unlensed alm from unlensed map
      call HealpixMap2alm(H, Phi, Alm_Phi, lmax, dopol=want_pol)

      call HealpixAlm2GradientMap(H, Alm_Phi, GradPhi, npix, 'T')

      ! Lens the map: from the unlensed alm (Alm_Unlensed_CMB) and the deflection map (GradPhi),
      ! compute the lensed map M
      if (lens_method == lens_exact) then
          call HealpixExactLensedMap_GradPhi(H,Alm_Unlensed_CMB,GradPhi,lensed_CMB)
      else if (lens_method == lens_interp) then
          call HealpixInterpLensedMap_GradPhi(H,Alm_Unlensed_CMB,GradPhi, lensed_CMB, interp_factor, interp_method, interp_algo)
      else
          stop 'unknown lens_method'
      end if

      !Save lensed map to .fits file
      call HealpixMap_Write(lensed_CMB, output_map, overwrite=.true.)

      ! Convert the lensed map (M) to lensed alm (A),
      ! overwriting the unlensed alm
      ! call HealpixMap2Alm(H,M, A, lmax, dopol = want_pol)

      ! Compute lensed power spectrum
      ! call HealpixAlm2Power(A,P)
      ! Note usually no need to free objects unless memory is short
      ! call HealpixAlm_Free(A)

      ! write lensed power spectrum to file
      ! call HealpixPower_Write(P,cls_lensed_file)

      end if
!-----------------------------------------------------------------------------

#ifdef MPIPIX
    call HealpixFree(H)
    call mpi_finalize(i)
#endif

#ifdef DEBUG
    write (*,*) 'End of program'
    pause
#endif
    end program SimLensCMB
