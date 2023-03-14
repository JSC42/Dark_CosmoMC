cd /Users/cangtao/IHEPBox/Projects/GitHub/HyRec/DarkSide

Input = 'EFF/Energy_1.txt'
Input = 'EFF/PBH_HIon_K1.txt'
Tab_Name = 'x_axis'

clc

OutFile = 'Tables.h';
delete(OutFile)
FileID=fopen(OutFile,'w');

% Define table sizes
fprintf(FileID,'#define Redshift_Size   63\n');
fprintf(FileID,'#define Electron_Mass_eV    510998.946\n');
fprintf(FileID,'#define Kinetic_Energy_Axis_Size	40\n');
fprintf(FileID,'#define Mdm_Axis_Size	62\n\n');

fprintf(FileID,'#define Kerr_PBH_Temperature_Size	200\n\n');

fclose(FileID);

% -- Redshift (1+z)

Input = 'EFF/Redshift_1.txt';
Tab_Name = 'Redshift_Axis';
Printer(Input, Tab_Name)


% -- Kerr PBH Temperature

Input = 'EFF/PBH_Temperature_1.txt';
Tab_Name = 'Kerr_PBH_Temperature_Axis';
Printer(Input, Tab_Name)

% -- Kinetic Energy

Input = 'EFF/Energy_1.txt';
Tab_Name = 'Kinetic_Energy_Axis';
Printer(Input, Tab_Name)

% -- SM Decay Mass Axis

Input = 'EFF/M_DEC_1.txt';
Tab_Name = 'Mdm_Axis_Decay';
Printer(Input, Tab_Name)

% -- SM Annihilation Mass Axis

Input = 'EFF/M_ANN_1.txt';
Tab_Name = 'Mdm_Axis_Annihilation';
Printer(Input, Tab_Name)

% -- Phot Ann ----

% HIon
Input = 'EFF/Phot_HIon_HMG.txt';
Tab_Name = 'EFF_Phot_HIon_Ann_HMG';
Printer(Input, Tab_Name)


% LyA
Input = 'EFF/Phot_LyA_HMG.txt';
Tab_Name = 'EFF_Phot_LyA_Ann_HMG';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/Phot_Heat_HMG.txt';
Tab_Name = 'EFF_Phot_Heat_Ann_HMG';
Printer(Input, Tab_Name)

% -- Phot Decay ----

% HIon
Input = 'EFF/Phot_HIon_DECAY.txt';
Tab_Name = 'EFF_Phot_HIon_Decay';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/Phot_LyA_DECAY.txt';
Tab_Name = 'EFF_Phot_LyA_Decay';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/Phot_Heat_DECAY.txt';
Tab_Name = 'EFF_Phot_Heat_Decay';
Printer(Input, Tab_Name)

% -- Elec Ann ----

% HIon
Input = 'EFF/Elec_HIon_HMG.txt';
Tab_Name = 'EFF_Elec_HIon_Ann_HMG';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/Elec_LyA_HMG.txt';
Tab_Name = 'EFF_Elec_LyA_Ann_HMG';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/Elec_Heat_HMG.txt';
Tab_Name = 'EFF_Elec_Heat_Ann_HMG';
Printer(Input, Tab_Name)

% -- Elec Decay ----

% HIon
Input = 'EFF/Elec_HIon_DECAY.txt';
Tab_Name = 'EFF_Elec_HIon_Decay';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/Elec_LyA_DECAY.txt';
Tab_Name = 'EFF_Elec_LyA_Decay';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/Elec_Heat_DECAY.txt';
Tab_Name = 'EFF_Elec_Heat_Decay';
Printer(Input, Tab_Name)

% -- Higgs Ann ----

% HIon
Input = 'EFF/Higgs_HIon_HMG.txt';
Tab_Name = 'EFF_Higgs_HIon_Ann_HMG';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/Higgs_LyA_HMG.txt';
Tab_Name = 'EFF_Higgs_LyA_Ann_HMG';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/Higgs_Heat_HMG.txt';
Tab_Name = 'EFF_Higgs_Heat_Ann_HMG';
Printer(Input, Tab_Name)

% -- Higgs Decay ----

% HIon
Input = 'EFF/Higgs_HIon_DEC.txt';
Tab_Name = 'EFF_Higgs_HIon_Decay';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/Higgs_LyA_DEC.txt';
Tab_Name = 'EFF_Higgs_LyA_Decay';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/Higgs_Heat_DEC.txt';
Tab_Name = 'EFF_Higgs_Heat_Decay';
Printer(Input, Tab_Name)

% -- Muon Ann ----

% HIon
Input = 'EFF/Muon_HIon_HMG.txt';
Tab_Name = 'EFF_Muon_HIon_Ann_HMG';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/Muon_LyA_HMG.txt';
Tab_Name = 'EFF_Muon_LyA_Ann_HMG';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/Muon_Heat_HMG.txt';
Tab_Name = 'EFF_Muon_Heat_Ann_HMG';
Printer(Input, Tab_Name)

% -- Muon Decay ----

% HIon
Input = 'EFF/Muon_HIon_DEC.txt';
Tab_Name = 'EFF_Muon_HIon_Decay';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/Muon_LyA_DEC.txt';
Tab_Name = 'EFF_Muon_LyA_Decay';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/Muon_Heat_DEC.txt';
Tab_Name = 'EFF_Muon_Heat_Decay';
Printer(Input, Tab_Name)

% -- Tau Ann ----

% HIon
Input = 'EFF/Tau_HIon_HMG.txt';
Tab_Name = 'EFF_Tau_HIon_Ann_HMG';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/Tau_LyA_HMG.txt';
Tab_Name = 'EFF_Tau_LyA_Ann_HMG';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/Tau_Heat_HMG.txt';
Tab_Name = 'EFF_Tau_Heat_Ann_HMG';
Printer(Input, Tab_Name)

% -- Tau Decay ----

% HIon
Input = 'EFF/Tau_HIon_DEC.txt';
Tab_Name = 'EFF_Tau_HIon_Decay';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/Tau_LyA_DEC.txt';
Tab_Name = 'EFF_Tau_LyA_Decay';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/Tau_Heat_DEC.txt';
Tab_Name = 'EFF_Tau_Heat_Decay';
Printer(Input, Tab_Name)

% -- Q Ann ----

% HIon
Input = 'EFF/Q_HIon_HMG.txt';
Tab_Name = 'EFF_Q_HIon_Ann_HMG';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/Q_LyA_HMG.txt';
Tab_Name = 'EFF_Q_LyA_Ann_HMG';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/Q_Heat_HMG.txt';
Tab_Name = 'EFF_Q_Heat_Ann_HMG';
Printer(Input, Tab_Name)

% -- Q Decay ----

% HIon
Input = 'EFF/Q_HIon_DEC.txt';
Tab_Name = 'EFF_Q_HIon_Decay';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/Q_LyA_DEC.txt';
Tab_Name = 'EFF_Q_LyA_Decay';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/Q_Heat_DEC.txt';
Tab_Name = 'EFF_Q_Heat_Decay';
Printer(Input, Tab_Name)

% -- Charm Ann ----

% HIon
Input = 'EFF/Charm_HIon_HMG.txt';
Tab_Name = 'EFF_Charm_HIon_Ann_HMG';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/Charm_LyA_HMG.txt';
Tab_Name = 'EFF_Charm_LyA_Ann_HMG';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/Charm_Heat_HMG.txt';
Tab_Name = 'EFF_Charm_Heat_Ann_HMG';
Printer(Input, Tab_Name)

% -- Charm Decay ----

% HIon
Input = 'EFF/Charm_HIon_DEC.txt';
Tab_Name = 'EFF_Charm_HIon_Decay';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/Charm_LyA_DEC.txt';
Tab_Name = 'EFF_Charm_LyA_Decay';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/Charm_Heat_DEC.txt';
Tab_Name = 'EFF_Charm_Heat_Decay';
Printer(Input, Tab_Name)

% -- Bottom Ann ----

% HIon
Input = 'EFF/Bottom_HIon_HMG.txt';
Tab_Name = 'EFF_Bottom_HIon_Ann_HMG';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/Bottom_LyA_HMG.txt';
Tab_Name = 'EFF_Bottom_LyA_Ann_HMG';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/Bottom_Heat_HMG.txt';
Tab_Name = 'EFF_Bottom_Heat_Ann_HMG';
Printer(Input, Tab_Name)

% -- Bottom Decay ----

% HIon
Input = 'EFF/Bottom_HIon_DEC.txt';
Tab_Name = 'EFF_Bottom_HIon_Decay';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/Bottom_LyA_DEC.txt';
Tab_Name = 'EFF_Bottom_LyA_Decay';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/Bottom_Heat_DEC.txt';
Tab_Name = 'EFF_Bottom_Heat_Decay';
Printer(Input, Tab_Name)

% -- Top Ann ----

% HIon
Input = 'EFF/Top_HIon_HMG.txt';
Tab_Name = 'EFF_Top_HIon_Ann_HMG';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/Top_LyA_HMG.txt';
Tab_Name = 'EFF_Top_LyA_Ann_HMG';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/Top_Heat_HMG.txt';
Tab_Name = 'EFF_Top_Heat_Ann_HMG';
Printer(Input, Tab_Name)

% -- Top Decay ----

% HIon
Input = 'EFF/Top_HIon_DEC.txt';
Tab_Name = 'EFF_Top_HIon_Decay';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/Top_LyA_DEC.txt';
Tab_Name = 'EFF_Top_LyA_Decay';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/Top_Heat_DEC.txt';
Tab_Name = 'EFF_Top_Heat_Decay';
Printer(Input, Tab_Name)

% -- W Ann ----

% HIon
Input = 'EFF/W_HIon_HMG.txt';
Tab_Name = 'EFF_W_HIon_Ann_HMG';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/W_LyA_HMG.txt';
Tab_Name = 'EFF_W_LyA_Ann_HMG';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/W_Heat_HMG.txt';
Tab_Name = 'EFF_W_Heat_Ann_HMG';
Printer(Input, Tab_Name)

% -- W Decay ----

% HIon
Input = 'EFF/W_HIon_DEC.txt';
Tab_Name = 'EFF_W_HIon_Decay';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/W_LyA_DEC.txt';
Tab_Name = 'EFF_W_LyA_Decay';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/W_Heat_DEC.txt';
Tab_Name = 'EFF_W_Heat_Decay';
Printer(Input, Tab_Name)

% -- Z Ann ----

% HIon
Input = 'EFF/Z_HIon_HMG.txt';
Tab_Name = 'EFF_Z_HIon_Ann_HMG';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/Z_LyA_HMG.txt';
Tab_Name = 'EFF_Z_LyA_Ann_HMG';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/Z_Heat_HMG.txt';
Tab_Name = 'EFF_Z_Heat_Ann_HMG';
Printer(Input, Tab_Name)

% -- Z Decay ----

% HIon
Input = 'EFF/Z_HIon_DEC.txt';
Tab_Name = 'EFF_Z_HIon_Decay';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/Z_LyA_DEC.txt';
Tab_Name = 'EFF_Z_LyA_Decay';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/Z_Heat_DEC.txt';
Tab_Name = 'EFF_Z_Heat_Decay';
Printer(Input, Tab_Name)

% -- Gluon Ann ----

% HIon
Input = 'EFF/Gluon_HIon_HMG.txt';
Tab_Name = 'EFF_Gluon_HIon_Ann_HMG';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/Gluon_LyA_HMG.txt';
Tab_Name = 'EFF_Gluon_LyA_Ann_HMG';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/Gluon_Heat_HMG.txt';
Tab_Name = 'EFF_Gluon_Heat_Ann_HMG';
Printer(Input, Tab_Name)

% -- Gluon Decay ----

% HIon
Input = 'EFF/Gluon_HIon_DEC.txt';
Tab_Name = 'EFF_Gluon_HIon_Decay';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/Gluon_LyA_DEC.txt';
Tab_Name = 'EFF_Gluon_LyA_Decay';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/Gluon_Heat_DEC.txt';
Tab_Name = 'EFF_Gluon_Heat_Decay';
Printer(Input, Tab_Name)

% -- Hawking K1 ----

% HIon
Input = 'EFF/PBH_HIon_K1.txt';
Tab_Name = 'EFF_Hawking_HIon_K1';
Printer(Input, Tab_Name)

Input = 'EFF/PBH_HIon_K1B.txt';
Tab_Name = 'EFF_Hawking_HIon_K1B';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/PBH_LyA_K1.txt';
Tab_Name = 'EFF_Hawking_LyA_K1';
Printer(Input, Tab_Name)

Input = 'EFF/PBH_LyA_K1B.txt';
Tab_Name = 'EFF_Hawking_LyA_K1B';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/PBH_Heat_K1.txt';
Tab_Name = 'EFF_Hawking_Heat_K1';
Printer(Input, Tab_Name)

Input = 'EFF/PBH_Heat_K1B.txt';
Tab_Name = 'EFF_Hawking_Heat_K1B';
Printer(Input, Tab_Name)

% -- Hawking K2 ----

% HIon
Input = 'EFF/PBH_HIon_K2.txt';
Tab_Name = 'EFF_Hawking_HIon_K2';
Printer(Input, Tab_Name)

Input = 'EFF/PBH_HIon_K2B.txt';
Tab_Name = 'EFF_Hawking_HIon_K2B';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/PBH_LyA_K2.txt';
Tab_Name = 'EFF_Hawking_LyA_K2';
Printer(Input, Tab_Name)

Input = 'EFF/PBH_LyA_K2B.txt';
Tab_Name = 'EFF_Hawking_LyA_K2B';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/PBH_Heat_K2.txt';
Tab_Name = 'EFF_Hawking_Heat_K2';
Printer(Input, Tab_Name)

Input = 'EFF/PBH_Heat_K2B.txt';
Tab_Name = 'EFF_Hawking_Heat_K2B';
Printer(Input, Tab_Name)

% -- Hawking K3 ----

% HIon
Input = 'EFF/PBH_HIon_K3.txt';
Tab_Name = 'EFF_Hawking_HIon_K3';
Printer(Input, Tab_Name)

Input = 'EFF/PBH_HIon_K3B.txt';
Tab_Name = 'EFF_Hawking_HIon_K3B';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/PBH_LyA_K3.txt';
Tab_Name = 'EFF_Hawking_LyA_K3';
Printer(Input, Tab_Name)

Input = 'EFF/PBH_LyA_K3B.txt';
Tab_Name = 'EFF_Hawking_LyA_K3B';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/PBH_Heat_K3.txt';
Tab_Name = 'EFF_Hawking_Heat_K3';
Printer(Input, Tab_Name)

Input = 'EFF/PBH_Heat_K3B.txt';
Tab_Name = 'EFF_Hawking_Heat_K3B';
Printer(Input, Tab_Name)

% -- Hawking K4 ----

% HIon
Input = 'EFF/PBH_HIon_K4.txt';
Tab_Name = 'EFF_Hawking_HIon_K4';
Printer(Input, Tab_Name)

Input = 'EFF/PBH_HIon_K4B.txt';
Tab_Name = 'EFF_Hawking_HIon_K4B';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/PBH_LyA_K4.txt';
Tab_Name = 'EFF_Hawking_LyA_K4';
Printer(Input, Tab_Name)

Input = 'EFF/PBH_LyA_K4B.txt';
Tab_Name = 'EFF_Hawking_LyA_K4B';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/PBH_Heat_K4.txt';
Tab_Name = 'EFF_Hawking_Heat_K4';
Printer(Input, Tab_Name)

Input = 'EFF/PBH_Heat_K4B.txt';
Tab_Name = 'EFF_Hawking_Heat_K4B';
Printer(Input, Tab_Name)

% -- Hawking K5 ----

% HIon
Input = 'EFF/PBH_HIon_K5.txt';
Tab_Name = 'EFF_Hawking_HIon_K5';
Printer(Input, Tab_Name)

Input = 'EFF/PBH_HIon_K5B.txt';
Tab_Name = 'EFF_Hawking_HIon_K5B';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/PBH_LyA_K5.txt';
Tab_Name = 'EFF_Hawking_LyA_K5';
Printer(Input, Tab_Name)

Input = 'EFF/PBH_LyA_K5B.txt';
Tab_Name = 'EFF_Hawking_LyA_K5B';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/PBH_Heat_K5.txt';
Tab_Name = 'EFF_Hawking_Heat_K5';
Printer(Input, Tab_Name)

Input = 'EFF/PBH_Heat_K5B.txt';
Tab_Name = 'EFF_Hawking_Heat_K5B';
Printer(Input, Tab_Name)

% -- Hawking K6 ----

% HIon
Input = 'EFF/PBH_HIon_K6.txt';
Tab_Name = 'EFF_Hawking_HIon_K6';
Printer(Input, Tab_Name)

Input = 'EFF/PBH_HIon_K6B.txt';
Tab_Name = 'EFF_Hawking_HIon_K6B';
Printer(Input, Tab_Name)

% LyA
Input = 'EFF/PBH_LyA_K6.txt';
Tab_Name = 'EFF_Hawking_LyA_K6';
Printer(Input, Tab_Name)

Input = 'EFF/PBH_LyA_K6B.txt';
Tab_Name = 'EFF_Hawking_LyA_K6B';
Printer(Input, Tab_Name)

% Heat
Input = 'EFF/PBH_Heat_K6.txt';
Tab_Name = 'EFF_Hawking_Heat_K6';
Printer(Input, Tab_Name)

Input = 'EFF/PBH_Heat_K6B.txt';
Tab_Name = 'EFF_Hawking_Heat_K6B';
Printer(Input, Tab_Name)

