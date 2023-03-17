/* Dark Matter Energy injection module */

// Rho_cr * c^2 /h^2, in eV/cm^3
#define Rhocr_C2_no_h2 1.054041699E4
#define Nearly_Zero 1.0E-50
#define Zmax_Energy_Injection 2000
#define debug_mode 0


void Validate_Inputs(REC_COSMOPARAMS *params)
{
    /* Check input params
    1. Particle Channel
    2. PBH Spin
    3. PBH Model: Currently only support [1 2 3]
    */
    int Particle_Channel, PBH_Model;
    Particle_Channel = Convert_to_Int(params->DM_Channel);
    PBH_Model = Convert_to_Int(params->PBH_Model);
    if ((Particle_Channel < 1) || (Particle_Channel > 12))
    {
        printf("Error from Validate_Inputs@HyRec: Unknown particle channel selected\n");
        exit(1);
    }
    if ((params->PBH_Spin > 0.99999) || (params->PBH_Spin < -1.0E-5))
    {
        printf("Error from Validate_Inputs@HyRec: wrong PBH_Spin.\n");
        exit(1);
    }
    if ((PBH_Model < 1) || (PBH_Model > 3))
    {
        printf("Error from Validate_Inputs@HyRec: Wrong choice of PBH_Model\n");
        exit(1);
    }
}

// ALL dEdVdt are in ev/cm^3/s unit

double dEdVdt_decay_inj(double z, REC_COSMOPARAMS *params)
{
    double Gamma, Omch2, r;
    Gamma = params->Gamma;
    Omch2 = params->odmh2;
    if (Gamma > 0)
    {
        r = Gamma * Omch2 * pow(1 + z, 3.0) * Rhocr_C2_no_h2;
    }
    else
    {
        r = 0.0;
    }
    return r;
}

double dEdVdt_ann_inj(double z, REC_COSMOPARAMS *params)
{
    double Omch2, r;
    Omch2 = params->odmh2;
    // 1E9 converts GeV to eV
    r = params->Pann * square(Rhocr_C2_no_h2 * Omch2 * cube(1. + z)) * 1.0E-9;
    return r;
}


double dEdVdt_Hawking_inj(double z, double Mbh, REC_COSMOPARAMS *params)
{
    // Hawking Radiation injection, Normalised to mbh>10^17 g
    // See Eq (3.10) of arxiv 2108.13256
    double r;
    r = (5.626976744186047e+29 / cube(Mbh) * params->fbh * params->odmh2 * cube(1. + z));
    return r;
}

double dEdVdt_decay_dep(double z, REC_COSMOPARAMS *params, int dep_channel)
{
    double inj, r, EFF, Mdm;
    int DM_Channel;
    DM_Channel = Convert_to_Int(params->DM_Channel);
    inj = dEdVdt_decay_inj(z, params);
    EFF = Interp_EFF_DM_Decay(params->Mdm, z, dep_channel, DM_Channel);
    r = EFF * inj;
    return r;
}


double dEdVdt_Hawking_Mono_dep(double z, double Mbh, REC_COSMOPARAMS *params, int dep_channel)
{
    // Hawking Radiation monochromatic deosition rate
    double inj, r, EFF, Mdm;
    inj = dEdVdt_Hawking_inj(z, Mbh, params);
    EFF = Interp_EFF_Hawking(params->Mbh, z, params->PBH_Spin, params->PBH_Model, dep_channel);
    r = EFF * inj;
    return r;
}

double dEdVdt_Hawking_dep(double z, REC_COSMOPARAMS *params, int dep_channel)
{
    // Hawking Radiation for general mass distributions
    // Currently only allow monochromatic, more on the way
    return dEdVdt_Hawking_Mono_dep(z, params->Mbh, params, dep_channel);
}

double dEdVdt_ann_dep(double z, REC_COSMOPARAMS *params, int dep_channel)
{
    double inj, r, EFF, Mdm;
    int DM_Channel;
    DM_Channel = (int)round(params->DM_Channel);
    inj = dEdVdt_ann_inj(z, params);
    EFF = Interp_EFF_DM_Annihilation(params->Mdm, z, dep_channel, DM_Channel);
    // printf("EFF = %E\n",EFF);
    r = EFF * inj;
    return r;
}

double dEdVdt_deposited(double z, REC_COSMOPARAMS *params, int dep_channel)
{
    /* Energy Deposition Rate in ev/cm^3/s
     -- inputs --
     dep_channel = 1: HIon
                   3: LyA
                   4: Heating
    */

    double r_dec, r_ann, r_Hawking, r;

    // Check params
    if ((dep_channel == 2) || (dep_channel == 5))
    {
        fprintf(stderr, "HeIon and Continnum dep channels not allowed, exitting\n");
        exit(1);
    }
    if (params->Gamma < Nearly_Zero)
    {
        r_dec = 0.0;
    }
    else
    {
        r_dec = dEdVdt_decay_dep(z, params, dep_channel);
    }
    if (params->Pann < Nearly_Zero)
    {
        r_ann = 0.0;
    }
    else
    {
        r_ann = dEdVdt_ann_dep(z, params, dep_channel);
    }
    if (params->fbh < Nearly_Zero)
    {
        r_Hawking = 0.0;
    }
    else
    {
        r_Hawking = dEdVdt_Hawking_dep(z, params, dep_channel);
    }

    r = r_dec + r_ann + r_Hawking;
    if (z > Zmax_Energy_Injection)
    {
        r = 0.0;
    }
    return r;
}

void Update_DarkArray(double z, REC_COSMOPARAMS *params, double *DarkArray)
{
    double dEdVdt_HIon, dEdVdt_LyA, dEdVdt_Heat, nH;
    dEdVdt_HIon = dEdVdt_deposited(z, params, 1);
    dEdVdt_LyA = dEdVdt_deposited(z, params, 3);
    dEdVdt_Heat = dEdVdt_deposited(z, params, 4);
    // Let's fill dEdVdt first and see what we can do with it
    DarkArray[0] = dEdVdt_HIon;
    DarkArray[1] = dEdVdt_LyA;
    DarkArray[2] = dEdVdt_Heat;
    // H number density in cm^3
    nH = params->nH0*cube(1+z)*1.E-6;
    DarkArray[3] = nH;
    // printf("%E  %E  %E\n", DarkArray[0], DarkArray[1], DarkArray[2]);
}

void Check_Error(double xe, double T)
{
    // Check for inifnity and NaN in Xe and T
    if (isfinite(xe) == 0)
    {
        printf("Error from Check_Error@HyRec: xe is NaN or infinite\n");
        exit(1);
    }
    if (isfinite(T) == 0)
    {
        printf("Error from Check_Error@HyRec: T is NaN or infinite\n");
        exit(1);
    }
}
