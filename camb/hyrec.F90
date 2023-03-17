    !---------------------------------------------------------------------------------------------------
    ! Recombination module for CAMB, using HyRec
    ! Author: Antony Lewis
    !---------------------------------------------------------------------------------------------------

    module Recombination
    use constants
    use AMLUtils
    implicit none
    private

    type RecombinationParams
    ! Adding DM params
    real(dl) :: DM_Channel           ! DM decay/annihilation channel
    real(dl) :: Mdm                  ! DM mass in GeV
    real(dl) :: Pann                 ! Dark matter annihilation rate
    real(dl) :: Gamma                ! Dark matter decay width
    real(dl) :: PBH_Model            ! PBH Energy injection mechanism
    real(dl) :: PBH_Distribution     ! PBH_Distribution
    real(dl) :: Mbh                  ! PBH mass
    real(dl) :: fbh                  ! PBH abundance
    real(dl) :: PBH_Lognormal_Sigma  ! PBH Log-normal distribution width
    real(dl) :: PBH_PWL_Mmax         ! PBH_PWL_Mmax
    real(dl) :: PBH_PWL_Gamma        ! PBH_PWL_Gamma
    real(dl) :: PBH_Spin             ! PBH_Spin
    
    end type RecombinationParams

    character(LEN=*), parameter :: Recombination_Name = 'HyRec'

    public RecombinationParams, Recombination_xe, Recombination_tm, Recombination_init,   &
        Recombination_ReadParams, Recombination_SetDefParams, &
        Recombination_Validate, Recombination_Name


    contains

    subroutine Recombination_ReadParams(R, Ini)
    use IniFile
    Type(RecombinationParams) :: R
    Type(TIniFile) :: Ini
    ! Read Params, the 3rd inout of Ini_Read_Double_File is the default
        R%DM_Channel = Ini_Read_Double_File(Ini, 'DM_Channel', 2.0D0) 
        R%Mdm = Ini_Read_Double_File(Ini, 'Mdm', 1.0D0) 
        R%Pann = Ini_Read_Double_File(Ini, 'Pann', 0.0D0)
        R%Gamma = Ini_Read_Double_File(Ini, 'Gamma', 0.0D0)
        R%PBH_Model = Ini_Read_Double_File(Ini, 'PBH_Model', 2.0D0)
        R%PBH_Distribution = Ini_Read_Double_File(Ini, 'PBH_Distribution', 1.0D0)
        R%Mbh = Ini_Read_Double_File(Ini, 'Mbh', 1.0D15)
        R%fbh = Ini_Read_Double_File(Ini, 'fbh', 0.0D0)
        R%PBH_Lognormal_Sigma = Ini_Read_Double_File(Ini, 'PBH_Lognormal_Sigma', 0.5D0)
        R%PBH_PWL_Mmax = Ini_Read_Double_File(Ini, 'PBH_PWL_Mmax', 1.0D17)
        R%PBH_PWL_Gamma = Ini_Read_Double_File(Ini, 'PBH_PWL_Gamma', 0.5D0)
        R%PBH_Spin = Ini_Read_Double_File(Ini, 'PBH_Spin', 0.0D0)

    end subroutine Recombination_ReadParams


    subroutine Recombination_SetDefParams(R)
    type (RecombinationParams) ::R


    end subroutine Recombination_SetDefParams


    subroutine Recombination_Validate(R, OK)
    Type(RecombinationParams), intent(in) :: R
    logical, intent(inout) :: OK


    end subroutine Recombination_Validate



    function Recombination_tm(a)
    real(dl), intent(in) :: a
    real(dl) Recombination_tm, hyrec_tm
    external hyrec_tm

    Recombination_tm =  hyrec_tm(a);

    end function Recombination_tm


    function Recombination_xe(a)
    real(dl), intent(in) :: a
    real(dl) Recombination_xe,hyrec_xe
    external hyrec_xe

    Recombination_xe = hyrec_xe(a);

    end function Recombination_xe


    subroutine Recombination_init(Recomb, OmegaC, OmegaB, OmegaN, Omegav, h0inp, tcmb, yp, num_nu)
    use AMLUtils
    implicit none
    Type (RecombinationParams), intent(in) :: Recomb
    real(dl), intent(in) :: OmegaC, OmegaB, OmegaN, OmegaV, h0inp, tcmb, yp, num_nu
    external rec_build_history_camb

    call rec_build_history_camb(OmegaC, OmegaB, OmegaN, Omegav, h0inp, tcmb, yp, num_nu,&
                               Recomb%DM_Channel, Recomb%Mdm, Recomb%Pann, Recomb%Gamma, Recomb%PBH_Model,&
                               Recomb%PBH_Distribution, Recomb%Mbh, Recomb%fbh, Recomb%PBH_Lognormal_Sigma,&
                               Recomb%PBH_PWL_Mmax, Recomb%PBH_PWL_Gamma, Recomb%PBH_Spin)

    end subroutine Recombination_init




    end module Recombination

