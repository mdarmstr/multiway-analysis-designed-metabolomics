

clear all %#ok

load wheat
load parafac_met

parafac_met = md_mean;

load parafasca_met

parafasca_met = md_mean;

original_position = 1:1:size(var_l,1);

[~,Iparafac] = sort(parafac_met,"descend");
[~,Iparafasca] = sort(parafasca_met,"descend");

parafac_position = original_position(Iparafac);
parafasca_position = original_position(Iparafasca);

displacement_vec = zeros(size(parafasca_position));

for ii = 1:size(parafasca_position,2)
    I = find(parafac_position(ii) == parafasca_position);
    displacement_vec(ii) = I -ii;
end

displacement_vec = displacement_vec';

disp(displacement_vec)

