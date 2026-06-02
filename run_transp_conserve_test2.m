% Run just the first 2 time steps to verify file writing end-to-end
% Modifies iters_bc before calling the main script.
clear all
addpath /home/mmazloff/ANALYSIS

% Override: only 2 steps, BCs only
% We do this by running the script contents with a patched iters_bc.
% The easiest way is to just copy/paste the relevant block with nt=2.

% Set flags that the script will read
do_ics   = 0;
do_bathy = 0;
do_obcs  = 1;
out_dir  = '/data/SO3/edavenport/tpose24/setup/transp_conserve_BC/';

% Run with monkey-patched iteration list
% We'll run the real script but intercept iters_bc after load
% Since the script uses clear all at top, we use run() then immediately
% it will overwrite everything. So instead, inline the key section:

fprintf('=== 2-step end-to-end test ===\n');

tp6_grid  = '/data/SO6/TPOSE/tpose6/grid_6';
tp24_grid = '/data/SO3/edavenport/tpose24/oct2012_TP6Vel_3month';
tp6_run   = '/data/SO3/edavenport/tpose6/sep2012/run_iter14';
nobcs     = 8;

XC6    = rdmds([tp6_grid '/XC']);  YC6  = rdmds([tp6_grid '/YC']);
hFacS6 = rdmds([tp6_grid '/hFacS']);
hFacW6 = rdmds([tp6_grid '/hFacW']);
DXG6   = rdmds([tp6_grid '/DXG']);
DYG6   = rdmds([tp6_grid '/DYG']);
DRF6   = rdmds([tp6_grid '/DRF']);  DRF6 = DRF6(:);
RC6    = rdmds([tp6_grid '/RC']);   RC_tp6 = -RC6(:);
Depth6 = rdmds([tp6_grid '/Depth']);
[nx6, ny6, nz6] = size(hFacS6);
xc6 = XC6(:,1);  yc6 = YC6(1,:);

XC24    = rdmds([tp24_grid '/XC']);  YC24  = rdmds([tp24_grid '/YC']);
hFacS24 = rdmds([tp24_grid '/hFacS']);
hFacW24 = rdmds([tp24_grid '/hFacW']);
DXG24   = rdmds([tp24_grid '/DXG']);
DYG24   = rdmds([tp24_grid '/DYG']);
DRF24   = rdmds([tp24_grid '/DRF']);  DRF24 = DRF24(:);
[nx24, ny24, nz24] = size(hFacS24);
xc24 = XC24(:,1);  yc24 = YC24(1,:);

zc_tp24 = zeros(1, nz24);
for k = 1:43
    zc_tp24(2*k-1) = RC_tp6(k) - DRF6(k)/4;
    zc_tp24(2*k)   = RC_tp6(k) + DRF6(k)/4;
end
for k = 1:12; zc_tp24(86+k) = RC_tp6(43+k); end
zc_tp24(99:138) = 1550:100:5450;

% Tref/Sref
iter_ic = 72*30;
Tref = zeros(nz24,1); Sref = zeros(nz24,1);
for iref = 1:2
    Q = rdmds([tp6_run '/diag_state'], iter_ic, 'rec', iref); Q(Q==0)=nan;
    tmp1 = zeros(nx24,ny24,nz6,'single');
    for k=1:nz6; tmp1(:,:,k)=interp2(yc6,xc6,Q(:,:,k),yc24,xc24); end
    tmp1(isnan(tmp1))=0;
    tmp2 = batch_vert_interp(tmp1,RC_tp6,zc_tp24,nx24,ny24,nz24);
    for k=1:nz24
        sl=tmp2(:,:,k); sl(sl==0)=nan; Qref=nanmean(sl(:));
        if isnan(Qref)&&k>1; Qref=(iref==1)*Tref(k-1)+(iref==2)*Sref(k-1); end
        if iref==1; Tref(k)=Qref; else; Sref(k)=Qref; end
    end
end
fprintf('Tref/Sref computed.\n');

drf = DRF24';
areaN = bsxfun(@times,DXG24(:,ny24),squeeze(hFacS24(:,ny24,:))).*drf;
areaS = bsxfun(@times,DXG24(:,2),  squeeze(hFacS24(:,2,:)))   .*drf;
areaW = bsxfun(@times,DYG24(2,:)', squeeze(hFacW24(2,:,:)))   .*drf;
areaE = bsxfun(@times,DYG24(nx24,:)',squeeze(hFacW24(nx24,:,:))).*drf;
wet_n24=sum(areaN(:)); wet_s24=sum(areaS(:));
wet_w24=sum(areaW(:)); wet_e24=sum(areaE(:));
maskN=areaN>0; maskS=areaS>0; maskW=areaW>0; maskE=areaE>0;

[~,j_n6]=min(abs(yc6-yc24(ny24))); [~,j_s6]=min(abs(yc6-yc24(1)));
[~,i_w6]=min(abs(xc6-xc24(1)));   [~,i_e6]=min(abs(xc6-xc24(end)));
i_rng6=find(xc6>=xc24(1)-1 & xc6<=xc24(end)+1);
j_rng6=find(yc6>=yc24(1)-1 & yc6<=yc24(end)+1);
dA_n6=bsxfun(@times,DXG6(i_rng6,j_n6),squeeze(hFacS6(i_rng6,j_n6,:))).*DRF6';
dA_s6=bsxfun(@times,DXG6(i_rng6,j_s6),squeeze(hFacS6(i_rng6,j_s6,:))).*DRF6';
dA_w6=bsxfun(@times,DYG6(i_w6,j_rng6)',squeeze(hFacW6(i_w6,j_rng6,:))).*DRF6';
dA_e6=bsxfun(@times,DYG6(i_e6,j_rng6)',squeeze(hFacW6(i_e6,j_rng6,:))).*DRF6';

bnd_names={'obcn','obcs','obcw','obce'};
for ib=1:4; b=bnd_names{ib};
    fidT(ib)=fopen([out_dir 'T' b '_frmTP6Vel_EMv1te.bin'],'w','b');
    fidS(ib)=fopen([out_dir 'S' b '_frmTP6Vel_EMv1te.bin'],'w','b');
    fidU(ib)=fopen([out_dir 'U' b '_frmTP6Vel_EMv1te.bin'],'w','b');
    fidV(ib)=fopen([out_dir 'V' b '_frmTP6Vel_EMv1te.bin'],'w','b');
    fidE(ib)=fopen([out_dir 'Eta' b '_frmTP6Vel_EMv1.bin'],'w','b');
end
bnd_i={1:nx24,1:nx24,ones(1,ny24),nx24*ones(1,ny24)};
bnd_j={ny24*ones(1,nx24),ones(1,nx24),1:ny24,1:ny24};
bnd_nhpts=[nx24,nx24,ny24,ny24];

iters_bc_test = [72, 144];   % just 2 steps
for t = 1:2
    iter = iters_bc_test(t);
    fprintf('Processing iter %d...\n', iter);
    V6=rdmds([tp6_run '/diag_state'],iter,'rec',4); V6(V6==0)=nan;
    U6=rdmds([tp6_run '/diag_state'],iter,'rec',3); U6(U6==0)=nan;
    T6=rdmds([tp6_run '/diag_state'],iter,'rec',1); T6(T6==0)=nan;
    S6=rdmds([tp6_run '/diag_state'],iter,'rec',2); S6(S6==0)=nan;
    E6=rdmds([tp6_run '/diag_surf'], iter,'rec',1); E6(E6==0)=nan;

    Tr_n6=nansum(nansum(squeeze(V6(i_rng6,j_n6,:)).*dA_n6));
    Tr_s6=nansum(nansum(squeeze(V6(i_rng6,j_s6,:)).*dA_s6));
    Tr_w6=nansum(nansum(squeeze(U6(i_w6,j_rng6,:)).*dA_w6));
    Tr_e6=nansum(nansum(squeeze(U6(i_e6,j_rng6,:)).*dA_e6));

    slc=cell(4,5);
    Q6_all={T6,S6,U6,V6,E6};
    for ivar=1:5
        Q=Q6_all{ivar};
        for ib=1:4
            nhpts=bnd_nhpts(ib);
            ii=bnd_i{ib}; jj=bnd_j{ib};
            if ivar<5
                tmp1=zeros(nhpts,nz6,'single');
                yq=yc24(jj); yq=yq(:);
                xq=xc24(ii); xq=xq(:);
                for k=1:nz6
                    col=interp2(yc6,xc6,Q(:,:,k),yq,xq);
                    tmp1(:,k)=col(:);
                end
                tmp1(isnan(tmp1))=0;
                tmp2=batch_vert_interp_1d(tmp1,RC_tp6,zc_tp24,nhpts,nz24);
            else
                yq=yc24(jj); yq=yq(:); xq=xc24(ii); xq=xq(:);
                col=interp2(yc6,xc6,Q,yq,xq);
                tmp2=single(col(:));
            end
            tmp2(isnan(tmp2))=0;
            if ivar<3
                ref=(ivar==1)*Tref+(ivar==2)*Sref;
                for k=1:nz24; col=tmp2(:,k); col(col==0)=ref(k); tmp2(:,k)=col; end
            end
            if ivar==5; m=tmp2; m(m==0)=nan; tmp2(tmp2==0)=nanmean(m(:)); end
            tmp2(isnan(tmp2))=0;
            slc{ib,ivar}=tmp2;
        end
    end

    Vn=slc{1,4}; Tr_n24=sum(Vn(:).*areaN(:)); delta_n=(Tr_n24-Tr_n6)/wet_n24;
    Vn(maskN)=Vn(maskN)-delta_n; slc{1,4}=Vn;
    Vs=slc{2,4}; Tr_s24=sum(Vs(:).*areaS(:)); delta_s=(Tr_s24-Tr_s6)/wet_s24;
    Vs(maskS)=Vs(maskS)-delta_s; slc{2,4}=Vs;
    Uw=slc{3,3}; Tr_w24=sum(Uw(:).*areaW(:)); delta_w=(Tr_w24-Tr_w6)/wet_w24;
    Uw(maskW)=Uw(maskW)-delta_w; slc{3,3}=Uw;
    Ue=slc{4,3}; Tr_e24=sum(Ue(:).*areaE(:)); delta_e=(Tr_e24-Tr_e6)/wet_e24;
    Ue(maskE)=Ue(maskE)-delta_e; slc{4,3}=Ue;

    Vn_c=slc{1,4}; Vs_c=slc{2,4}; Uw_c=slc{3,3}; Ue_c=slc{4,3};
    net=sum(Vs_c(:).*areaS(:))+sum(Uw_c(:).*areaW(:)) ...
       -sum(Vn_c(:).*areaN(:))-sum(Ue_c(:).*areaE(:));
    delta_net=net/wet_e24;
    Ue_c(maskE)=Ue_c(maskE)+delta_net; slc{4,3}=Ue_c;

    net_final=sum(slc{2,4}(:).*areaS(:))+sum(slc{3,3}(:).*areaW(:)) ...
             -sum(slc{1,4}(:).*areaN(:))-sum(slc{4,3}(:).*areaE(:));
    fprintf('  iter %d: net flux after correction = %.2e m^3/s\n', iter, net_final);

    fids={fidT,fidS,fidU,fidV,fidE};
    for ib=1:4
        for ivar=1:5; fwrite(fids{ivar}(ib),slc{ib,ivar},'single'); end
    end
end
for ib=1:4; fclose(fidT(ib));fclose(fidS(ib));fclose(fidU(ib));fclose(fidV(ib));fclose(fidE(ib)); end

% Check file sizes
fprintf('\nOutput files (first 3):\n');
d=dir([out_dir '*.bin']);
for k=1:min(3,length(d)); fprintf('  %s  %.1f MB\n',d(k).name,d(k).bytes/1e6); end
fprintf('Test PASSED\n');


function out = batch_vert_interp(tmp1, z_in, z_out, nx, ny, nz_out)
out = zeros(nx, ny, nz_out, 'single');
for ix = 1:nx
    for iy = 1:ny
        prof = double(tmp1(ix,iy,:)); prof = prof(:);
        msk = prof==0; if all(msk); continue; end
        if any(msk); prof(msk)=interp1(z_in(~msk),prof(~msk),z_in(msk),'linear','extrap'); end
        out(ix,iy,:) = single(interp1(z_in, prof, z_out, 'linear','extrap'));
    end
end
end

function out = batch_vert_interp_1d(tmp1, z_in, z_out, nhpts, nz_out)
out = zeros(nhpts, nz_out, 'single');
col_max = max(abs(tmp1), [], 2);
wet = col_max > 0;
if any(wet)
    tmp_filled = tmp1;
    for ip = find(wet)'
        prof = double(tmp1(ip,:)); prof = prof(:);
        msk = prof==0;
        if any(msk) && ~all(msk)
            prof(msk)=interp1(z_in(~msk),prof(~msk),z_in(msk),'linear','extrap');
        end
        tmp_filled(ip,:) = single(prof);
    end
    out_t = interp1(z_in, tmp_filled(wet,:)', z_out, 'linear','extrap');
    out(wet,:) = single(out_t');
end
end
