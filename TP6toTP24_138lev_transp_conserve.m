% TP6toTP24_138lev_transp_conserve.m
%
% Interpolates TP6 model output onto the TP24 138-level grid for initial
% conditions, bathymetry, and open boundary conditions.
%
% Velocity boundary conditions are transport-conserved: after regridding,
% a uniform barotropic correction is applied so that transport through each
% face matches the TP6 source transport. A final residual correction is
% applied to the east boundary to enforce zero net volume flux.
%
% Control flags at the top let you skip ICs/bathy when only BCs are needed.

clear all
addpath /home/mmazloff/ANALYSIS

% =========================================================================
% Control flags
% =========================================================================
do_ics   = 0;   % generate initial conditions
do_bathy = 0;   % generate bathymetry files
do_obcs  = 1;   % generate transport-conserved boundary conditions

out_dir = '/data/SO3/edavenport/tpose24/setup/transp_conserve_BC/';
if ~exist(out_dir,'dir'); mkdir(out_dir); end

% =========================================================================
% TP6 grid
% =========================================================================
tp6_grid = '/data/SO6/TPOSE/tpose6/grid_6';
XC6    = rdmds([tp6_grid '/XC']);
YC6    = rdmds([tp6_grid '/YC']);
XG6    = rdmds([tp6_grid '/XG']);
YG6    = rdmds([tp6_grid '/YG']);
hFacC6 = rdmds([tp6_grid '/hFacC']);
hFacS6 = rdmds([tp6_grid '/hFacS']);
hFacW6 = rdmds([tp6_grid '/hFacW']);
DXG6   = rdmds([tp6_grid '/DXG']);   % zonal spacing at V-pts  [nx6 x ny6]
DYG6   = rdmds([tp6_grid '/DYG']);   % meridional spacing at U-pts [nx6 x ny6]
DRF6   = rdmds([tp6_grid '/DRF']);   % level thicknesses [nz6 x 1]
RC6    = rdmds([tp6_grid '/RC']);     % level centers, negative downward
Depth6 = rdmds([tp6_grid '/Depth']);
[nx6, ny6, nz6] = size(hFacC6);
xc6 = XC6(:,1);   yc6 = YC6(1,:);
RC_tp6 = -RC6(:); % positive downward
DRF6   = DRF6(:);

% =========================================================================
% TP24 grid (from existing run)
% =========================================================================
tp24_grid = '/data/SO3/edavenport/tpose24/oct2012_TP6Vel_3month';
XC24    = rdmds([tp24_grid '/XC']);
YC24    = rdmds([tp24_grid '/YC']);
hFacC24 = rdmds([tp24_grid '/hFacC']);
hFacS24 = rdmds([tp24_grid '/hFacS']);
hFacW24 = rdmds([tp24_grid '/hFacW']);
DXG24   = rdmds([tp24_grid '/DXG']);
DYG24   = rdmds([tp24_grid '/DYG']);
DRF24   = rdmds([tp24_grid '/DRF']);
[nx24, ny24, nz24] = size(hFacC24);
xc24 = XC24(:,1);   yc24 = YC24(1,:);
DRF24 = DRF24(:);

% =========================================================================
% 138-level vertical grid (center depths, positive downward)
% =========================================================================
zc_tp24 = zeros(1, nz24);
for k = 1:43
    zc_tp24(2*k-1) = RC_tp6(k) - DRF6(k)/4;
    zc_tp24(2*k)   = RC_tp6(k) + DRF6(k)/4;
end
for k = 1:12
    zc_tp24(86+k) = RC_tp6(43+k);
end
zc_tp24(99:138) = 1550:100:5450;

% =========================================================================
% Source TP6 run paths
% =========================================================================
tp6_run  = '/data/SO3/edavenport/tpose6/sep2012/run_iter14';
iter_ic  = 72*30;
iters_bc = 72:72:8784;
nt       = length(iters_bc);
nobcs    = 8;

% =========================================================================
% Compute Tref / Sref from IC time step for T/S land fill on boundaries
% (needed whether or not do_ics is set)
% =========================================================================
fprintf('Computing T/S reference profiles for land fill...\n');
Tref = zeros(nz24,1); Sref = zeros(nz24,1);
for iref = 1:2
    Q = rdmds([tp6_run '/diag_state'], iter_ic, 'rec', iref);
    Q(Q==0) = nan;
    tmp1 = zeros(nx24, ny24, nz6, 'single');
    for k = 1:nz6
        tmp1(:,:,k) = interp2(yc6, xc6, Q(:,:,k), yc24, xc24);
    end
    tmp1(isnan(tmp1)) = 0;
    tmp2 = zeros(nx24, ny24, nz24, 'single');
    for ix = 1:nx24
        for iy = 1:ny24
            prof = double(tmp1(ix,iy,:)); prof = prof(:);
            msk = prof==0;
            if all(msk); continue; end
            if any(msk)
                prof(msk) = interp1(RC_tp6(~msk), prof(~msk), RC_tp6(msk), 'linear','extrap');
            end
            tmp2(ix,iy,:) = single(interp1(RC_tp6, prof, zc_tp24, 'linear','extrap'));
        end
    end
    for k = 1:nz24
        sl = tmp2(:,:,k); sl(sl==0) = nan;
        Qref = nanmean(sl(:));
        if isnan(Qref) && k>1; Qref = (iref==1)*Tref(k-1) + (iref==2)*Sref(k-1); end
        if iref==1; Tref(k)=Qref; else; Sref(k)=Qref; end
    end
end

% =========================================================================
% INITIAL CONDITIONS
% =========================================================================
if do_ics
fprintf('--- Initial conditions ---\n');
ic_files = {'Tini_frmTP6Vel_EMv1te.bin', 'Sini_frmTP6Vel_EMv1te.bin', ...
            'Uini_frmTP6Vel_EMv1te.bin', 'Vini_frmTP6Vel_EMv1te.bin', ...
            'ETAini_frmTP6Vel_EMv1.bin'};
for iini = 1:5
    switch iini
        case 1; Q = rdmds([tp6_run '/diag_state'], iter_ic, 'rec', 1);
        case 2; Q = rdmds([tp6_run '/diag_state'], iter_ic, 'rec', 2);
        case 3; Q = rdmds([tp6_run '/diag_state'], iter_ic, 'rec', 3);
        case 4; Q = rdmds([tp6_run '/diag_state'], iter_ic, 'rec', 4);
        case 5; Q = rdmds([tp6_run '/diag_surf'],  iter_ic, 'rec', 1);
    end
    Q(Q==0) = nan;
    if iini < 5
        tmp1 = zeros(nx24, ny24, nz6, 'single');
        for k = 1:nz6
            tmp1(:,:,k) = interp2(yc6, xc6, Q(:,:,k), yc24, xc24);
        end
        tmp1(isnan(tmp1)) = 0;
        tmp2 = batch_vert_interp(tmp1, RC_tp6, zc_tp24, nx24, ny24, nz24);
    else
        tmp2 = single(interp2(yc6, xc6, Q, yc24, xc24));
    end
    tmp2(isnan(tmp2)) = 0;
    if iini < 3
        for k = 1:nz24
            ref = tmp2(:,:,k); ref(ref==0) = nan;
            tmp3 = fillininils(nx24, ny24, tmp2(:,:,k), [27 15], 1.5);
            tmp3(tmp3==0) = (iini==1)*Tref(k) + (iini==2)*Sref(k);
            tmp2(:,:,k) = tmp3;
        end
    end
    for i = 1:nobcs
        tmp2(i,:,:)       = tmp2(nobcs+1,:,:);
        tmp2(end-i+1,:,:) = tmp2(end-nobcs,:,:);
        tmp2(:,i,:)       = tmp2(:,nobcs+1,:);
        tmp2(:,end-i+1,:) = tmp2(:,end-nobcs,:);
    end
    fid = fopen([out_dir ic_files{iini}],'w','b');
    fwrite(fid, tmp2, 'single'); fclose(fid);
    fprintf('  wrote %s\n', ic_files{iini});
end
end % do_ics

% =========================================================================
% BATHYMETRY
% =========================================================================
if do_bathy
fprintf('--- Bathymetry ---\n');
Depth6(Depth6==0) = nan;
tmp = single(interp2(yc6, xc6, -abs(Depth6), yc24, xc24));
tmp(isnan(tmp)) = 0;
for i = 1:nobcs
    tmp(i,:)       = tmp(nobcs+1,:);
    tmp(end-i+1,:) = tmp(end-nobcs,:);
    tmp(:,i)       = tmp(:,nobcs+1);
    tmp(:,end-i+1) = tmp(:,end-nobcs);
end
fid = fopen([out_dir 'Bathy_frmTP6_EMv1.bin'],'w','b'); fwrite(fid,tmp,'single'); fclose(fid);
tmp2 = tmp; tmp2(1,:)=0; tmp2(end,:)=0; tmp2(:,1)=0; tmp2(:,end)=0;
fid = fopen([out_dir 'Bathy_frmTP6_EMv1_CB.bin'],'w','b'); fwrite(fid,tmp2,'single'); fclose(fid);
fprintf('  wrote bathymetry files\n');
end % do_bathy

% =========================================================================
% BOUNDARY CONDITIONS WITH TRANSPORT CONSERVATION
% =========================================================================
if do_obcs
fprintf('--- Boundary conditions ---\n');

% ------------------------------------------------------------------
% Pre-compute TP24 boundary face areas  [nbnd_pts x nz24]
% MITgcm convention:
%   obcn/V: south face of row ny24  → hFacS(:,ny24,:)
%   obcs/V: south face of row 2     → hFacS(:,2,:)
%   obcw/U: west face of col 2      → hFacW(2,:,:)
%   obce/U: west face of col nx24   → hFacW(nx24,:,:)
% ------------------------------------------------------------------
drf = DRF24';   % 1 x nz24 for broadcasting

areaN = bsxfun(@times, DXG24(:,ny24),   squeeze(hFacS24(:,ny24,:)))   .* drf;  % [nx24 x nz24]
areaS = bsxfun(@times, DXG24(:,2),      squeeze(hFacS24(:,2,:)))      .* drf;
areaW = bsxfun(@times, DYG24(2,:)',     squeeze(hFacW24(2,:,:)))       .* drf;  % [ny24 x nz24]
areaE = bsxfun(@times, DYG24(nx24,:)',  squeeze(hFacW24(nx24,:,:)))    .* drf;

wet_n24 = sum(areaN(:));
wet_s24 = sum(areaS(:));
wet_w24 = sum(areaW(:));
wet_e24 = sum(areaE(:));

maskN = areaN > 0;   % logical masks for wet points
maskS = areaS > 0;
maskW = areaW > 0;
maskE = areaE > 0;

% ------------------------------------------------------------------
% Pre-compute TP6 source boundary face areas
% Find TP6 row/col indices at each TP24 boundary location
% ------------------------------------------------------------------
[~, j_n6] = min(abs(yc6 - yc24(ny24)));
[~, j_s6] = min(abs(yc6 - yc24(1)));
[~, i_w6] = min(abs(xc6 - xc24(1)));
[~, i_e6] = min(abs(xc6 - xc24(end)));

% TP6 index ranges covering TP24 domain extent
i_rng6 = find(xc6 >= xc24(1)-1 & xc6 <= xc24(end)+1);
j_rng6 = find(yc6 >= yc24(1)-1 & yc6 <= yc24(end)+1);

drf6 = DRF6';   % 1 x nz6

dA_n6 = bsxfun(@times, DXG6(i_rng6, j_n6), squeeze(hFacS6(i_rng6, j_n6, :))) .* drf6;
dA_s6 = bsxfun(@times, DXG6(i_rng6, j_s6), squeeze(hFacS6(i_rng6, j_s6, :))) .* drf6;
dA_w6 = bsxfun(@times, DYG6(i_w6, j_rng6)', squeeze(hFacW6(i_w6, j_rng6, :))) .* drf6;
dA_e6 = bsxfun(@times, DYG6(i_e6, j_rng6)', squeeze(hFacW6(i_e6, j_rng6, :))) .* drf6;

% ------------------------------------------------------------------
% Open output files (all 4 boundaries x 5 variables)
% ------------------------------------------------------------------
bnd_names = {'obcn','obcs','obcw','obce'};
for ib = 1:4
    b = bnd_names{ib};
    fidT(ib) = fopen([out_dir 'T' b '_frmTP6Vel_EMv1te.bin'],'w','b');
    fidS(ib) = fopen([out_dir 'S' b '_frmTP6Vel_EMv1te.bin'],'w','b');
    fidU(ib) = fopen([out_dir 'U' b '_frmTP6Vel_EMv1te.bin'],'w','b');
    fidV(ib) = fopen([out_dir 'V' b '_frmTP6Vel_EMv1te.bin'],'w','b');
    fidE(ib) = fopen([out_dir 'Eta' b '_frmTP6Vel_EMv1.bin'],'w','b');
end

% Boundary slice dimensions on TP24 grid: [horiz_pts x nz24]
% obcn/obcs: nx24 pts;  obcw/obce: ny24 pts
bnd_i = {1:nx24, 1:nx24, ones(1,ny24),   nx24*ones(1,ny24)};
bnd_j = {ny24*ones(1,nx24), ones(1,nx24), 1:ny24, 1:ny24};
bnd_nhpts = [nx24, nx24, ny24, ny24];

% ------------------------------------------------------------------
% Time loop
% ------------------------------------------------------------------
for t = 1:nt
    iter = iters_bc(t);
    if mod(t,10)==0
        fprintf('  t=%d/%d  iter=%d\n', t, nt, iter);
    end

    % Read full TP6 3-D fields for this time step
    V6 = rdmds([tp6_run '/diag_state'], iter, 'rec', 4);  V6(V6==0)=nan;
    U6 = rdmds([tp6_run '/diag_state'], iter, 'rec', 3);  U6(U6==0)=nan;
    T6 = rdmds([tp6_run '/diag_state'], iter, 'rec', 1);  T6(T6==0)=nan;
    S6 = rdmds([tp6_run '/diag_state'], iter, 'rec', 2);  S6(S6==0)=nan;
    E6 = rdmds([tp6_run '/diag_surf'],  iter, 'rec', 1);  E6(E6==0)=nan;

    % ---- Source transport on TP6 grid ----
    % squeeze removes the singleton j- or i-dimension before element-wise multiply
    Tr_n6 = nansum(nansum(squeeze(V6(i_rng6, j_n6, :)) .* dA_n6));
    Tr_s6 = nansum(nansum(squeeze(V6(i_rng6, j_s6, :)) .* dA_s6));
    Tr_w6 = nansum(nansum(squeeze(U6(i_w6, j_rng6, :)) .* dA_w6));
    Tr_e6 = nansum(nansum(squeeze(U6(i_e6, j_rng6, :)) .* dA_e6));

    % Storage for all 4 boundary slices, all 5 variables
    % Cell arrays indexed by [boundary, variable]
    slc = cell(4,5);  % slc{ib,ivar} = [nhpts x nz24] (or [nhpts x 1] for SSH)

    % ---- Interpolate each variable onto all 4 boundaries ----
    Q6_all = {T6, S6, U6, V6, E6};
    for ivar = 1:5
        Q = Q6_all{ivar};
        for ib = 1:4
            nhpts = bnd_nhpts(ib);
            ii = bnd_i{ib}; jj = bnd_j{ib};

            if ivar < 5
                % Horizontal interp (all TP6 levels) → [nhpts x nz6]
                % Use (:) on both query args to force column vectors and
                % prevent interp2 from creating a 2-D grid.
                tmp1 = zeros(nhpts, nz6, 'single');
                yq = yc24(jj); yq = yq(:);   % force column [nhpts x 1]
                xq = xc24(ii); xq = xq(:);
                for k = 1:nz6
                    col = interp2(yc6, xc6, Q(:,:,k), yq, xq);
                    tmp1(:,k) = col(:);
                end
                tmp1(isnan(tmp1)) = 0;
                % Vertical interp to 138 levels (batched)
                tmp2 = batch_vert_interp_1d(tmp1, RC_tp6, zc_tp24, nhpts, nz24);
            else
                yq = yc24(jj); yq = yq(:);
                xq = xc24(ii); xq = xq(:);
                col = interp2(yc6, xc6, Q, yq, xq);
                tmp2 = single(col(:));
            end
            tmp2(isnan(tmp2)) = 0;

            % T/S land fill along boundary
            if ivar < 3
                ref = (ivar==1)*Tref + (ivar==2)*Sref;
                for k = 1:nz24
                    col = tmp2(:,k);
                    col(col==0) = ref(k);
                    tmp2(:,k) = col;
                end
            end
            if ivar==5
                m = tmp2; m(m==0)=nan;
                tmp2(tmp2==0) = nanmean(m(:));
            end
            tmp2(isnan(tmp2)) = 0;
            slc{ib,ivar} = tmp2;
        end
    end

    % ================================================================
    % Transport conservation for V (boundaries 1=N, 2=S)
    % ================================================================
    % North
    Vn = slc{1,4};   % [nx24 x nz24]
    Tr_n24 = sum(Vn(:) .* areaN(:));
    delta_n = (Tr_n24 - Tr_n6) / wet_n24;
    Vn(maskN) = Vn(maskN) - delta_n;
    slc{1,4} = Vn;

    % South
    Vs = slc{2,4};
    Tr_s24 = sum(Vs(:) .* areaS(:));
    delta_s = (Tr_s24 - Tr_s6) / wet_s24;
    Vs(maskS) = Vs(maskS) - delta_s;
    slc{2,4} = Vs;

    % ================================================================
    % Transport conservation for U (boundaries 3=W, 4=E)
    % ================================================================
    % West
    Uw = slc{3,3};   % [ny24 x nz24]
    Tr_w24 = sum(Uw(:) .* areaW(:));
    delta_w = (Tr_w24 - Tr_w6) / wet_w24;
    Uw(maskW) = Uw(maskW) - delta_w;
    slc{3,3} = Uw;

    % East
    Ue = slc{4,3};
    Tr_e24 = sum(Ue(:) .* areaE(:));
    delta_e = (Tr_e24 - Tr_e6) / wet_e24;
    Ue(maskE) = Ue(maskE) - delta_e;
    slc{4,3} = Ue;

    % ================================================================
    % Net volume flux balance: adjust east boundary
    % Net inflow = Tr_s_corr + Tr_w_corr - Tr_n_corr - Tr_e_corr
    % (Tr values after per-boundary correction should equal TP6 transports,
    %  but TP6 itself may not be exactly balanced → zero it here)
    % ================================================================
    Vn_corr = slc{1,4};  Vs_corr = slc{2,4};
    Uw_corr = slc{3,3};  Ue_corr = slc{4,3};
    Tr_n_c = sum(Vn_corr(:) .* areaN(:));
    Tr_s_c = sum(Vs_corr(:) .* areaS(:));
    Tr_w_c = sum(Uw_corr(:) .* areaW(:));
    Tr_e_c = sum(Ue_corr(:) .* areaE(:));
    net_flux = Tr_s_c + Tr_w_c - Tr_n_c - Tr_e_c;  % positive = net inflow
    % Remove net inflow by increasing eastward outflow on east boundary
    delta_net = net_flux / wet_e24;
    Ue_corr(maskE) = Ue_corr(maskE) + delta_net;
    slc{4,3} = Ue_corr;

    % ================================================================
    % Write to binary files
    % ================================================================
    fids = {fidT, fidS, fidU, fidV, fidE};
    for ib = 1:4
        for ivar = 1:5
            fwrite(fids{ivar}(ib), slc{ib,ivar}, 'single');
        end
    end

end % time loop

for ib = 1:4
    fclose(fidT(ib)); fclose(fidS(ib)); fclose(fidU(ib));
    fclose(fidV(ib)); fclose(fidE(ib));
end
fprintf('Done. Files written to %s\n', out_dir);

end % do_obcs


% =========================================================================
% Local helper: vertical interpolation for a 3D field [nx x ny x nz_in]
% Returns [nx x ny x nz_out]
% =========================================================================
function out = batch_vert_interp(tmp1, z_in, z_out, nx, ny, nz_out)
out = zeros(nx, ny, nz_out, 'single');
for ix = 1:nx
    for iy = 1:ny
        prof = double(tmp1(ix,iy,:)); prof = prof(:);
        msk = prof==0;
        if all(msk); continue; end
        if any(msk)
            prof(msk) = interp1(z_in(~msk), prof(~msk), z_in(msk), 'linear','extrap');
        end
        out(ix,iy,:) = single(interp1(z_in, prof, z_out, 'linear','extrap'));
    end
end
end

% =========================================================================
% Local helper: vertical interpolation for a 2D field [nhpts x nz_in]
% Returns [nhpts x nz_out]
% =========================================================================
function out = batch_vert_interp_1d(tmp1, z_in, z_out, nhpts, nz_out)
out = zeros(nhpts, nz_out, 'single');
% Bulk interpolation for non-land columns (fast path)
% For the boundary slice, most columns are wet, so batch all at once
% then fix any all-zero columns
nz_in = length(z_in);
% Identify wet columns (have at least one non-zero value)
col_max = max(abs(tmp1), [], 2);   % [nhpts x 1]  -- sum over z_in dim
wet = col_max > 0;
if any(wet)
    % Build a padded profile: fill zero (land) values via neighbor interp
    tmp_filled = tmp1;
    for ip = find(wet)'
        prof = double(tmp1(ip,:)); prof = prof(:);
        msk = prof==0;
        if any(msk) && ~all(msk)
            prof(msk) = interp1(z_in(~msk), prof(~msk), z_in(msk), 'linear','extrap');
        end
        tmp_filled(ip,:) = single(prof);
    end
    % interp1 on matrix: interpolates each column of Y
    % tmp_filled is [nhpts x nz_in]; interp1 treats rows as samples → transpose
    out_t = interp1(z_in, tmp_filled(wet,:)', z_out, 'linear','extrap');  % [nz_out x nwet]
    out(wet,:) = single(out_t');
end
end
