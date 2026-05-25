clear all
addpath ~/

dobcs = 1;

load /data/SO6/TPOSE/tpose6/grid_6/grid  XC YC XG YG hFacC DRF RC Depth
[nx ny nz] = size(hFacC);
XC = XC(:,1); YC = YC(1,:);
x1 = 210; x2 = 230; %E
y1 = -5;  y2 = 10;
DRFem(1:2:2*nz) = DRF/2;
DRFem(2:2:2*nz) = DRF/2;
DRFem(end-1:end)=[];

%FOR TELE - 138 LEVEL VERSION
%
% Vertical grid structure:
%
%   Levels 1-86:   TP6 levels 1-43 halved via DRFem (0-500m)    -- UNCHANGED
%   Levels 87-98:  TP6 levels 44-55 direct copy (500-1500m)     -- UNCHANGED
%   Levels 99-138: 40 x 100m boxes (1500-5500m)                 -- NEW
%
% Original levels 99-108 were: 200/300/300/300/400/500x5m (10 levels, 1500-5500m)
% Replaced with 40 x 100m boxes covering the same total depth (1500-5500m).
% New levels 99-138 are vertically interpolated from TP6 levels 55-66
%   (centers 1450-5750m), which safely brackets all new target depths (1550-5450m).

NZ_NEW = 138;

KI1  = 1:86;    % shallow doubled zone (unchanged)
KI2i = 44:55;   % TP6 levels for direct-copy zone (levels 87-98, unchanged)

% TP6 source levels for interpolating new levels 99-138
KI_interp_tp6 = 55:66;   % TP6 levels 55-66, centers 1450-5750m
KI_interp_out = 99:138;  % output levels 99-138, 100m boxes 1500-5500m

% TP6 level center depths (positive downward)
RC_tp6 = -RC;

% Center depths of new levels 99-138 (100m boxes, centers at 1550,1650,...,5450m)
zc_interp_out = 1550:100:5450;   % 40 values

% Sanity check
fprintf('New grid: %d levels\n', NZ_NEW);
fprintf('  Levels 1-86:   TP6 1-43 doubled (0-500m)\n');
fprintf('  Levels 87-98:  TP6 44-55 direct copy (500-1500m)\n');
fprintf('  Levels 99-138: 40 x 100m interpolated (1500-5500m)\n');

I1 = max(find(XC<x1))-1;
I2 = max(find(XC<x2))+2;
J1 = max(find(YC<y1))-1;
J2 = max(find(YC<y2))+2;

%ADD IN OBCS AND ADJUST DOMAIN SIZE
 xgOrigin = XG(I1,J1) - 8/24;
 ygOrigin = YG(I1,J1) - 4/24;

 xgFin = XG(I2+1,J2+1) + 8/24;
 ygFin = YG(I2+1,J2+1) + 4/24;

 nx = length(xgOrigin:1/24:xgFin);
 ny = length(ygOrigin:1/24:ygFin);

xc = [xgOrigin+1/48:1/24:xgFin]';
yc = ygOrigin+1/48:1/24:ygFin;

%ADD OBCS width
nobcs = 8;

% TP6 source level centers for interp1
zc_src = RC_tp6(KI_interp_tp6);

%NOW GET TSUV
for iini = 1:4
  switch iini
    case 1; ivar = 'THETA'; bstr = 'Tini_frmTP6Vel_EMv1te.bin';
	   Q = rdmds('/data/SO3/edavenport/tpose6/sep2012/run_iter14/diag_state',72*30,'rec',1);
    case 2; ivar = 'SALT'; bstr = 'Sini_frmTP6Vel_EMv1te.bin';
           Q = rdmds('/data/SO3/edavenport/tpose6/sep2012/run_iter14/diag_state',72*30,'rec',2);
    case 3; ivar = 'UVEL'; bstr = 'Uini_frmTP6Vel_EMv1te.bin';
           Q = rdmds('/data/SO3/edavenport/tpose6/sep2012/run_iter14/diag_state',72*30,'rec',3);
    case 4; ivar = 'VVEL'; bstr = 'Vini_frmTP6Vel_EMv1te.bin';
           Q = rdmds('/data/SO3/edavenport/tpose6/sep2012/run_iter14/diag_state',72*30,'rec',4);
  end

  % Step 1: horizontal interpolation onto new horizontal grid (all TP6 levels)
  tmp1 = zeros(nx,ny,nz,'single');
  Q(Q==0) = nan;
  for k = 1:nz
    tmp1(:,:,k) = interp2(YC,XC,Q(:,:,k),yc,xc);
  end

  % Step 2: vertical mapping onto 138-level grid
  tmp2 = zeros(nx,ny,NZ_NEW,'single');

  % Levels 1-86: double TP6 levels 1-43 (unchanged)
  for k = 1:43
    tmp2(:,:,k*2-1) = tmp1(:,:,k);
    tmp2(:,:,k*2)   = tmp1(:,:,k);
  end

  % Levels 87-98: direct copy of TP6 levels 44-55 (unchanged)
  for k = 1:length(KI2i)
    tmp2(:,:,86+k) = tmp1(:,:,KI2i(k));
  end

  % Levels 99-138: vertical interpolation from TP6 levels 55-66
  tmp1(isnan(tmp1)) = 0;
  for ix = 1:nx
    for iy = 1:ny
      prof = double(tmp1(ix,iy,KI_interp_tp6));
      mask = prof == 0;
      if all(mask)
        tmp2(ix,iy,KI_interp_out) = 0;
      else
        if any(mask)
          prof(mask) = interp1(zc_src(~mask), prof(~mask), ...
                               zc_src(mask), 'linear', 'extrap');
        end
        tmp2(ix,iy,KI_interp_out) = single(interp1(zc_src, prof, ...
                                            zc_interp_out, 'linear', 'extrap'));
      end
    end
  end

  tmp2(isnan(tmp2)) = 0;

  if iini < 3  % EXTRAP T and S onto land
    for k = 1:NZ_NEW
      ref = tmp2(:,:,k); ref(ref==0) = nan;
      Qref(k) = nanmean(ref(:));
      if isnan(Qref(k)); Qref(k) = Qref(k-1); end
      tmp3 = fillininils(nx,ny,tmp2(:,:,k),[27 15],1.5);
      tmp3(tmp3==0) = Qref(k);
      tmp2(:,:,k) = tmp3;
    end
    if iini==1; Tref=Qref; elseif iini==2; Sref=Qref; end
  end

  for i = 1:nobcs
    tmp2(i,:,:)       = tmp2(nobcs+1,:,:);
    tmp2(end-i+1,:,:) = tmp2(end-nobcs,:,:);
    tmp2(:,i,:)       = tmp2(:,nobcs+1,:);
    tmp2(:,end-i+1,:) = tmp2(:,end-nobcs,:);
  end

  fid = fopen(['/data/SO3/edavenport/tpose24/setup/' bstr],'w','b');
  fwrite(fid,tmp2,'single');
  fclose(fid);
end

%AND SSH
  Q = rdmds('/data/SO3/edavenport/tpose6/sep2012/run_iter14/diag_surf',72*30,'rec',1);
  Q(Q==0) = nan;
  tmp = zeros(nx,ny,'single');
  tmp(:,:) = interp2(YC,XC,Q(:,:),yc,xc);
  tmp(isnan(tmp)) = 0;
  for i = 1:nobcs
    tmp(i,:)       = tmp(nobcs+1,:);
    tmp(end-i+1,:) = tmp(end-nobcs,:);
    tmp(:,i)       = tmp(:,nobcs+1);
    tmp(:,end-i+1) = tmp(:,end-nobcs);
  end
  tmp = fillininils(nx,ny,tmp,[27 15],1.5);
  fid = fopen(['/data/SO3/edavenport/tpose24/setup/ETAini_frmTP6Vel_EMv1.bin'],'w','b');
  fwrite(fid,tmp,'single');
  fclose(fid);

%AND BATHY
  tmp(:,:) = interp2(YC,XC,-1*abs(Depth),yc,xc);
  for i = 1:nobcs
    tmp(i,:)       = tmp(nobcs+1,:);
    tmp(end-i+1,:) = tmp(end-nobcs,:);
    tmp(:,i)       = tmp(:,nobcs+1);
    tmp(:,end-i+1) = tmp(:,end-nobcs);
  end
  fid = fopen(['/data/SO3/edavenport/tpose24/setup/setup/Bathy_frmTP6_EMv1.bin'],'w','b');
  fwrite(fid,tmp,'single');
  fclose(fid);

  tmp(1,:)=0;  tmp(end,:)=0;
  tmp(:,1)=0;  tmp(:,end)=0;
  fid = fopen(['/data/SO3/edavenport/tpose24/setup/Bathy_frmTP6_EMv1_CB.bin'],'w','b');
  fwrite(fid,tmp,'single');
  fclose(fid);

if dobcs  %AND MAKE OBCS
for iobcs = 1:4
  switch iobcs
    case 1; ivar = 'obcn'; i = 1:nx; j = ny-nobcs;
    case 2; ivar = 'obcs'; i = 1:nx; j = 1;
    case 3; ivar = 'obcw'; i = nx-nobcs; j = 1:ny;
    case 4; ivar = 'obce'; i = 1; j = 1:ny;
  end
  fidT = fopen(['/data/SO3/edavenport/tpose24/setup/T' ivar '_frmTP6Vel_EMv1te.bin'],'w','b');
  fidS = fopen(['/data/SO3/edavenport/tpose24/setup/S' ivar '_frmTP6Vel_EMv1te.bin'],'w','b');
  fidU = fopen(['/data/SO3/edavenport/tpose24/setup/U' ivar '_frmTP6Vel_EMv1te.bin'],'w','b');
  fidV = fopen(['/data/SO3/edavenport/tpose24/setup/V' ivar '_frmTP6Vel_EMv1te.bin'],'w','b');
  fidE = fopen(['/data/SO3/edavenport/tpose24/setup/Eta' ivar '_frmTP6Vel_EMv1.bin'],'w','b');

  for iter = 72:72:8784
    iter
    for iini = 1:5
      clear tmp1 tmp2 tmp3 tmp4
      switch iini
        case 1;
             Q = rdmds('/data/SO3/edavenport/tpose6/sep2012/run_iter14/diag_state',iter,'rec',1);
        case 2;
             Q = rdmds('/data/SO3/edavenport/tpose6/sep2012/run_iter14/diag_state',iter,'rec',2);
        case 3;
             Q = rdmds('/data/SO3/edavenport/tpose6/sep2012/run_iter14/diag_state',iter,'rec',3);
        case 4;
             Q = rdmds('/data/SO3/edavenport/tpose6/sep2012/run_iter14/diag_state',iter,'rec',4);
        case 5;
             Q = rdmds('/data/SO3/edavenport/tpose6/sep2012/run_iter14/diag_surf',iter,'rec',1);
      end
      Q(Q==0) = nan;

      tnx = length(i);
      tny = length(j);
      tnz = NZ_NEW; if iini == 5; tnz = 1; end

      tmp2 = zeros(tnx,tny,tnz,'single');

      if iini < 5
        % Horizontal interp onto boundary slice (all TP6 levels)
        tmp1 = zeros(tnx,tny,nz,'single');
        for k = 1:nz
          tmp1(:,:,k) = interp2(YC,XC,Q(:,:,k),yc(j),xc(i));
        end

        % Levels 1-86: double TP6 1-43
        for k = 1:43
          tmp2(:,:,k*2-1) = tmp1(:,:,k);
          tmp2(:,:,k*2)   = tmp1(:,:,k);
        end

        % Levels 87-98: direct copy TP6 44-55
        for k = 1:length(KI2i)
          tmp2(:,:,86+k) = tmp1(:,:,KI2i(k));
        end

        % Levels 99-138: vertical interpolation from TP6 55-66
        tmp1(isnan(tmp1)) = 0;
        for ix = 1:tnx
          for iy = 1:tny
            prof = double(tmp1(ix,iy,KI_interp_tp6));
            mask = prof == 0;
            if all(mask)
              tmp2(ix,iy,KI_interp_out) = 0;
            else
              if any(mask)
                prof(mask) = interp1(zc_src(~mask), prof(~mask), ...
                                     zc_src(mask), 'linear', 'extrap');
              end
              tmp2(ix,iy,KI_interp_out) = single(interp1(zc_src, prof, ...
                                                  zc_interp_out, 'linear', 'extrap'));
            end
          end
        end

      else  % SSH
        tmp2(:,:) = interp2(YC,XC,Q(:,:),yc(j),xc(i));
      end

      tmp2(isnan(tmp2)) = 0;
      tmp2 = squeeze(tmp2);

      if iini < 3  % T/S extrapolation on boundary slice
        tmp3 = tmp2; tmp3(tmp3==0) = nan; tmp3 = nanmean(tmp3,1);
        for k = 1:tnz
          tmp4 = tmp2(:,k);
          tmp4(tmp4==0) = tmp3(k);
          if iini==1
            tmp4(tmp4==0) = Tref(k);
          elseif iini==2
            tmp4(tmp4==0) = Sref(k);
          end
          tmp2(:,k) = tmp4;
        end
      end

      if iini == 5
        tmp3 = tmp2; tmp3(tmp3==0) = nan; tmp3 = nanmean(tmp3);
        tmp2(tmp2==0) = tmp3;
      end

      tmp2(isnan(tmp2)) = 0;

      switch iini
        case 1;  fwrite(fidT,tmp2,'single');
        case 2;  fwrite(fidS,tmp2,'single');
        case 3;  fwrite(fidU,tmp2,'single');
        case 4;  fwrite(fidV,tmp2,'single');
        case 5;  fwrite(fidE,tmp2,'single');
      end
    end %variable
  end %iter

  fclose(fidT);
  fclose(fidS);
  fclose(fidU);
  fclose(fidV);
  fclose(fidE);
end %obcs
end
