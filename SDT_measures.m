function [d_prime,A_prime,nat_lgB,c,GrierB] = SDT_measures (hits, fa, n_sig, n_noise )
%% ----------------------- Script Description -----------------------------
% Function to calculate signal detection theory metrices.
% Takes in input:
% 1) hits: hit rate
% 2) fa: false alarm
% * Both should be expressed between 0 to 1
% 3) n_sig : number of signal trials
% 4) n_noise : number of noise trials
%
% Produces output:
% 1) d_prime : Measures distance between signal and noise in SD
% 2) A_prime : Non-parametric measure of sensitivity
% 3) ln(B) : Natural log of response bias
% 4) c : distance between the criterion and the neutral point
% 5) Grier's B : Non-parametric measure of response bias (different from
% the B in (3))

% All formulas based on Stanislaw & Todorov 1999
% Written on 15/1/2016

%% ----------------------- Starting Script -----------------------------

% Use loglinear if hits or fa is 0 or 1
if hits == 0 || hits == 1 || fa == 0 || fa == 1
    hit_num = hits * n_sig;
    fa_num = fa * n_noise;
    new_hit = hit_num + .5;
    new_fa = fa_num + .5;
    hits = new_hit / (n_sig + 1);
    fa = new_fa / (n_noise + 1);
end


diff = hits - fa;
z_hits = norminv(hits,0,1);
z_fa = norminv(fa,0,1);

% Calculating d'
d_prime = z_hits - z_fa;

% Calculating A'

if diff >= 0
    A_prime =.5 + (((hits - fa)*(1+hits-fa))/(4*hits*(1 - fa)));
elseif diff < 0
    A_prime =.5 - (((fa - hits)*(1+fa-hits))/(4*fa*(1 - hits)));
end

% Calculating ln(B)
nat_lgB = ((z_fa^2) - (z_hits^2))/2;

% Calculating c
c = -(z_hits + z_fa) / 2;

% Calculating Grier's B
if diff >= 0
    GrierB = ((hits*(1 - hits)) - (fa*(1 - fa))) / ...
        ((hits*(1 - hits)) + (fa*(1 - fa)));
elseif diff < 0
    GrierB = ((fa*(1 - fa)) - (hits*(1 - hits))) / ...
        ((fa*(1 - fa)) + (hits*(1 - hits)));
end