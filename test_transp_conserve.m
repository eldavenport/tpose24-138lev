% test_transp_conserve.m
% Quick sanity check: run transport-conserving interpolation for the north
% boundary only, at the first BC time step (iter=72). Prints source TP6
% transport, raw TP24 transport, and corrected TP24 transport.

clear all
addpath /home/mmazloff/ANALYSIS

% ---- TP6 grid ----
tp6_grid = '/data/SO6/TPOSE/tpose6/grid_6';
XC6    = rdmds([tp6_grid '/XC']);  YC6  = rdmds([tp6_grid '/YC']);
hFacS6 = rdmds([tp6_grid '/hFacS']);
DXG6   = rdmds([tp6_grid '/DXG']);
DRF6   = rdmds([tp6_grid '/DRF']);  DRF6 = DRF6(:);
RC6    = rdmds([tp6_grid '/RC']);   RC_tp6 = -RC6(:);
[nx6, ny6, nz6] = size(hFacS6);
xc6 = XC6(:,1);  yc6 = YC6(1,:);

% ---- TP24 grid ----
tp24_grid = '/data/SO3/edavenport/tpose24/oct2012_TP6Vel_3month';
XC24    = rdmds([tp24_grid '/XC']);  YC24  = rdmds([tp24_grid '/YC']);
hFacS24 = rdmds([tp24_grid '/hFacS']);
DXG24   = rdmds([tp24_grid '/DXG']);
DRF24   = rdmds([tp24_grid '/DRF']);  DRF24 = DRF24(:);
[nx24, ny24, nz24] = size(hFacS24);
xc24 = XC24(:,1);  yc24 = YC24(1,:);

% ---- 138-level vertical grid ----
zc_tp24 = zeros(1, nz24);
for k = 1:43
    zc_tp24(2*k-1) = RC_tp6(k) - DRF6(k)/4;
    zc_tp24(2*k)   = RC_tp6(k) + DRF6(k)/4;
end
for k = 1:12; zc_tp24(86+k) = RC_tp6(43+k); end
zc_tp24(99:138) = 1550:100:5450;

% ---- Boundary areas ----
drf24 = DRF24';
areaN = bsxfun(@times, DXG24(:,ny24), squeeze(hFacS24(:,ny24,:))) .* drf24;
wet_n24 = sum(areaN(:));
maskN = areaN > 0;

% ---- TP6 source area (north boundary) ----
[~, j_n6] = min(abs(yc6 - yc24(ny24)));
i_rng6 = find(xc6 >= xc24(1)-1 & xc6 <= xc24(end)+1);
dA_n6 = bsxfun(@times, DXG6(i_rng6, j_n6), ...
               squeeze(hFacS6(i_rng6, j_n6, :))) .* DRF6';

% ---- Read TP6 V at iter=72 ----
tp6_run = '/data/SO3/edavenport/tpose6/sep2012/run_iter14';
iter = 72;
V6 = rdmds([tp6_run '/diag_state'], iter, 'rec', 4);
V6(V6==0) = nan;

% ---- Source transport ----
Tr_n6 = nansum(nansum(squeeze(V6(i_rng6, j_n6, :)) .* dA_n6));
fprintf('TP6 source transport (north):   %+.4e m^3/s\n', Tr_n6);

% ---- Interpolate V onto north boundary ----
nhpts = nx24;
tmp1 = zeros(nhpts, nz6, 'single');
for k = 1:nz6
    col = interp2(yc6, xc6, V6(:,:,k), yc24(ny24)*ones(1,nx24), xc24');
    tmp1(:,k) = col(:);
end
tmp1(isnan(tmp1)) = 0;

% Vertical interp (batched)
tmp_filled = tmp1;
for ip = 1:nhpts
    prof = double(tmp1(ip,:)); prof = prof(:);
    msk = prof==0;
    if all(msk) || ~any(msk); continue; end
    prof(msk) = interp1(RC_tp6(~msk), prof(~msk), RC_tp6(msk), 'linear','extrap');
    tmp_filled(ip,:) = single(prof);
end
wet = max(abs(tmp_filled),[],2) > 0;
Vn = zeros(nhpts, nz24, 'single');
if any(wet)
    out_t = interp1(RC_tp6, tmp_filled(wet,:)', zc_tp24, 'linear','extrap');
    Vn(wet,:) = single(out_t');
end
Vn(isnan(Vn)) = 0;

% ---- Raw TP24 transport ----
Tr_n24_raw = sum(Vn(:) .* areaN(:));
fprintf('TP24 raw transport (north):     %+.4e m^3/s\n', Tr_n24_raw);

% ---- Apply correction ----
delta_n = (Tr_n24_raw - Tr_n6) / wet_n24;
Vn_corr = Vn;
Vn_corr(maskN) = Vn_corr(maskN) - delta_n;
Tr_n24_corr = sum(Vn_corr(:) .* areaN(:));
fprintf('TP24 corrected transport (north):%+.4e m^3/s\n', Tr_n24_corr);
fprintf('Barotropic correction applied:  %+.4e m/s\n', delta_n);
fprintf('Residual error:                 %+.4e m^3/s  (%.2f%%)\n', ...
    Tr_n24_corr - Tr_n6, 100*(Tr_n24_corr-Tr_n6)/abs(Tr_n6));
