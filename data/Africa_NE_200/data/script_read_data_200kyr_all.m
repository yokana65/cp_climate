% Script read all data fpr 200 kyr analysis NE Africa
% ML Fischer 2024

warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')


% ODP 967
data_odp_967_22     = readtable('data_odp_967_22_43247_2021_339_MOESM2_ESM_XRF_Ti_Al.txt');
data_odp_967_22     = data_odp_967_22(:,[8 19]);
data_odp_967_22     = table2array(data_odp_967_22);

index           = data_odp_967_22(:, 1) < 200;
data_odp_967_22 = data_odp_967_22(index, :);


% ODP 721/722
data_odp721_722_terr     = readtable('data_odp_721_722_terr.txt');
data_odp721_722_terr     = table2array(data_odp721_722_terr);
data_odp721_722_terr     = data_odp721_722_terr(:,[2 3]);

index                = data_odp721_722_terr(:, 1) < 200;
data_odp721_722_terr = data_odp721_722_terr(index, :);


% ODP 709
data_odp_709  = readtable('ODP709_Ca_K_Ti_Fe_ratio_600kyr.txt');
data_odp_709  = table2array(data_odp_709);

index        = data_odp_709(:, 1) < 200;
data_odp_709 = data_odp_709(index, :);


% ICPD Chew Bahir
data_icdp_chb     = readtable('data_icdp_chb14_2_cb01_cb03_age_rrmarch2021_mht560.txt');
data_icdp_chb     = table2array(data_icdp_chb);

index         = data_icdp_chb(:, 1) < 200;
data_icdp_chb = data_icdp_chb(index, :);


% KL 09
data_kl09     = readtable('data_KL09_41467_2014_BFncomms6076_MOESM1344_ESM.txt');
data_kl09     = table2array(data_kl09);
data_kl09     = data_kl09(:,[1 2]);

index     = data_kl09(:, 1) < 200;
data_kl09 = data_kl09(index, :);


% KL 11
data_kl11        = readtable('data_KL11_geochem.txt');
data_kl11        = table2array(data_kl11);
data_kl11        = data_kl11(:,[2 3]);

index     = data_kl11(:, 1) < 200;
data_kl11 = data_kl11(index, :);


% KL 15
data_kl15       = readtable('data_KL15_XRF_tuned.txt');
datakl15agem    = readtable('data_KL15-2_smooth_spline_ages.txt');

data_kl15       = renamevars(data_kl15, "NEW_depth_cm","depth");
datakl15agem    = renamevars(datakl15agem, "best","age");
T               = outerjoin(datakl15agem,data_kl15);
data_kl15       = table2array(T);
data_kl15       = data_kl15(2:2168,[4 8]);
data_kl15(:,1)  = data_kl15(:,1)./1000;

index     = data_kl15(:, 1) < 200;
data_kl15 = data_kl15(index, :);

data_kl15(isnan(data_kl15(:,2))==1,:) = [];


% Lake Tana
data_lake_tana        = readtable('Lake Tana XRF data file Lamb et al 2018.txt');
data_lake_tana        = table2array(data_lake_tana);
data_lake_tana        = data_lake_tana(:,[2 4]);
data_lake_tana(:,1)  = data_lake_tana(:,1)./1000;

index     = data_lake_tana(:, 1) < 200;
data_lake_tana = data_lake_tana(index, :);






