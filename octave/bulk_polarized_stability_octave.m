function result = bulk_polarized_stability_octave()
  % GNU Octave/MATLAB cross-check of the Python bulk LSWT calculation.
  % Run from the project root:
  %   addpath('octave'); bulk_polarized_stability_octave();
  % No add-on packages or MEX files are required.

  ws = make_ewald_workspace(0.22, 42.0, 8.0);
  models = {'M0','M1','M2','M3','M4'};
  rows = cell(1, numel(models));
  for i = 1:numel(models)
    rows{i} = locate_mode(models{i}, ws);
    fprintf('%s: Hc=%.12f T, q=(%.9f, %.9f, %.9f)\n', ...
      models{i}, rows{i}.critical_field_T, rows{i}.critical_wavevector_rlu);
  end

  muB = 0.67171381563;
  gz = 4.716;
  mev2k = 11.60451812155;
  analytic_m0 = (0.1225 * mev2k * (3 + 1.5*0.636))/(gz*muB);
  assert(abs(rows{1}.critical_field_T - analytic_m0) < 1e-7);
  result = struct('method','Octave LSWT/Ewald independent cross-check', ...
                  'analytic_M0_field_T',analytic_m0,'rows',{rows});
  if ~exist('results','dir'), mkdir('results'); end
  fid = fopen(fullfile('results','bulk_polarized_stability_octave.json'),'w');
  if fid < 0, error('cannot open Octave output file'); end
  fprintf(fid,'%s',jsonencode(result));
  fclose(fid);
end


function row = locate_mode(model, ws)
  j0 = coupling_fourier([0;0;0], model, ws);
  if any(strcmp(model, {'M0','M1','M2'}))
    starts = [1/3,2/3; 2/3,1/3; 1/3,1/3; 2/3,2/3]';
    dim = 2;
  else
    starts = [1/3,2/3,0.5; 2/3,1/3,0.5; ...
              1/3,2/3,1/3; 2/3,1/3,2/3; ...
              1/3,1/3,1/3; 2/3,2/3,2/3]';
    dim = 3;
  end
  opts = optimset('TolX',1e-10,'TolFun',1e-11,'MaxIter',1500, ...
                  'MaxFunEvals',5000,'Display','off');
  best_h = -Inf; best_q = zeros(3,1);
  for s = 1:size(starts,2)
    f = @(x) -hc_energy(mod(x(:),1), model, ws, j0, dim);
    [x,fval] = fminsearch(f, starts(:,s), opts);
    if -fval > best_h
      best_h = -fval;
      if dim == 2, best_q = [mod(x(:),1);0]; else, best_q = mod(x(:),1); end
    end
  end
  row = struct('model',model,'boundary','tin_foil', ...
               'critical_wavevector_rlu',best_q', ...
               'critical_energy_K',best_h, ...
               'critical_field_T',best_h/(4.716*0.67171381563));
end


function h = hc_energy(x, model, ws, j0, dim)
  if dim == 2, q = [x(1);x(2);0]; else, q = x(1:3); end
  jk = coupling_fourier(q, model, ws);
  h = 0.5*(j0(3,3)-min(eig(jk(1:2,1:2))));
end


function out = coupling_fourier(q, model, ws)
  out = exchange_fourier(q, model);
  if any(strcmp(model, {'M3','M4'}))
    geom = ewald_fourier(q, ws);
    pref = 53.68151120621697; % micro-eV Angstrom^3
    g = diag([4.2,4.2,4.716]);
    out = out + pref*g*geom*g*(11.60451812155/1000);
  end
  out = (out+out')/2;
end


function out = exchange_fourier(q, model)
  mev2k = 11.60451812155;
  jzz = 0.1225; delta = 0.636; jpp = -0.0019;
  kd = 0.00036; gxy = 4.2; gz = 4.716;
  disps = [1,0,0; 0,1,0; 1,-1,0]';
  angles = [0,pi/3,-pi/3];
  out = zeros(3);
  for n = 1:3
    theta = angles(n);
    if strcmp(model,'M0')
      xy = jzz*delta; zz = jzz; pp = 0;
    elseif strcmp(model,'M2')
      xy = jzz*delta-kd*gxy^2/2;
      zz = jzz+kd*gz^2;
      pp = jpp-3*kd*gxy^2/4;
    else
      xy = jzz*delta; zz = jzz; pp = jpp;
    end
    c = cos(2*theta); s = sin(2*theta);
    tensor = mev2k*[xy+2*pp*c,-2*pp*s,0; ...
                    -2*pp*s,xy-2*pp*c,0; 0,0,zz];
    out = out + 2*cos(2*pi*disps(:,n)'*q)*tensor;
  end
  if strcmp(model,'M4')
    inter_d = [0,0,1; 1,0,1; 0,1,1; 1,-1,1; ...
               -1,0,1; 0,-1,1; -1,1,1]';
    inter_j = [2.8,0,0,0,7.9,7.9,7.9];
    for n = 1:numel(inter_j)
      tensor = diag([0.08,0.08,1])*inter_j(n)*(mev2k/1000);
      out = out + 2*cos(2*pi*inter_d(:,n)'*q)*tensor;
    end
  end
end


function ws = make_ewald_workspace(alpha, real_cut, recip_cut)
  a = 5.312; c = 6.983;
  cellm = [a,a/2,0; 0,sqrt(3)*a/2,0; 0,0,c];
  recip = 2*pi*inv(cellm)';
  rv = integer_vectors(cellm, real_cut);
  lens = sqrt(sum(rv.^2,1));
  keep = lens > 1e-13 & lens <= real_cut;
  rv = rv(:,keep); lens = lens(keep);
  rt = zeros(9,numel(lens));
  for n = 1:numel(lens)
    r = rv(:,n); l = lens(n); e = exp(-(alpha*l)^2);
    b = erfc(alpha*l)/l^3 + 2*alpha/sqrt(pi)*e/l^2;
    cc = 3*erfc(alpha*l)/l^5 + 2*alpha/sqrt(pi)*e*(3/l^4+2*alpha^2/l^2);
    t = b*eye(3)-cc*(r*r'); rt(:,n) = t(:);
  end
  maxk = 0.5*sum(sqrt(sum(recip.^2,1)));
  gv = integer_vectors(recip, recip_cut+maxk);
  ws = struct('cell',cellm,'recip',recip,'volume',abs(det(cellm)), ...
              'alpha',alpha,'real_cut',real_cut,'recip_cut',recip_cut, ...
              'rv',rv,'rt',rt,'gv',gv, ...
              'self',4*alpha^3/(3*sqrt(pi))*eye(3));
end


function vecs = integer_vectors(basis, cutoff)
  singular = svd(basis); nmax = ceil(cutoff/min(singular))+1;
  [i,j,k] = ndgrid(-nmax:nmax,-nmax:nmax,-nmax:nmax);
  coeff = [i(:),j(:),k(:)]';
  vecs = basis*coeff;
end


function out = ewald_fourier(frac, ws)
  red = mod(frac(:)+0.5,1)-0.5;
  k = ws.recip*red;
  weights = cos(k'*ws.rv)';
  realpart = reshape(ws.rt*weights,3,3);
  q = ws.gv+k;
  q2 = sum(q.^2,1);
  keep = q2 > 1e-24 & q2 <= ws.recip_cut^2;
  qq = q(:,keep); q2 = q2(keep);
  factors = (4*pi/ws.volume)*exp(-q2/(4*ws.alpha^2))./q2;
  recipart = (qq .* factors)*qq';
  out = realpart+recipart-ws.self;
  out = (out+out')/2;
end
