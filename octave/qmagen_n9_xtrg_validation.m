function payload = qmagen_n9_xtrg_validation(D)
  % Validate the public QMagen XTRG engine against exact N=9 diagonalization.
  % The upstream repository is only read.  Results are written in the project.
  if nargin < 1, D = 8; end

  script_dir = fileparts(mfilename('fullpath'));
  project_root = fileparts(script_dir);
  repo = fullfile(project_root,'sources','raw','QMagen_public','V0.0');
  run_dir = fullfile(repo,'RunScript');
  old_dir = pwd();
  cleanup = onCleanup(@() cd(old_dir));
  cd(run_dir);

  addpath(fullfile(repo,'Class'));
  addpath(fullfile(repo,'UtilityFunc'));
  addpath(genpath(fullfile(repo,'SpinModel')));
  addpath(fullfile(repo,'ManyBodySolver'));
  addpath(fullfile(repo,'ManyBodySolver','TensorFunction'));
  addpath(genpath(fullfile(repo,'ManyBodySolver','XTRG')));
  addpath(genpath(fullfile(repo,'svd_lapack_interface')));
  addpath(fullfile(script_dir,'qmagen_octave_compat'),'-begin');

  global EVOFLAG EVO_check PLOTFLAG SAVEFLAG;
  EVOFLAG = 0; EVO_check = 0; PLOTFLAG = 0; SAVEFLAG = 0;

  Jzz_K = 0.1225 * 11.60451812155;
  Delta = 0.636;
  gz = 4.716;
  muB_over_kB = 0.67171381563;
  field_T = 1.50;

  Para = struct();
  Para.ManyBodySolver = 'XTRG';
  Para.IntrcMap_Name = 'IntrcMap_TLXXZ';
  Para.d = 2; Para.L = 9;
  Para.Geo = struct('Lx',3,'Ly',3,'BCX','OBC','BCY','PBC','L',9);
  Para.Model = struct('J1xy',Delta,'J1z',1.0,'J2xy',0.0,'J2z',0.0, ...
                      'ES',Jzz_K,'g',[4.2,4.2,gz]);
  Para.Field = struct('h',[0,0,gz*muB_over_kB*field_T/Jzz_K], ...
                      'B',[0,0,field_T]);
  Para.ThDQ = 'Cm';
  Para.Parallel = false; Para.RunParallel = false;
  Para.RunParallelNum = 1; Para.RunParallelmaxNumCompThreads = 1;
  Para.tau = 0.00025;
  Para.It = 15;
  Para.D_list = D*ones(1,Para.It);
  Para.MCrit = D;
  Para.Ver = 'Memory';
  Para.VariProd_step_max = 50;
  Para.VariSum_step_max = 50;
  Para.SETTN_init_step_max = 50;
  Para.TStr = sprintf('qmagen_n9_D%d',D);
  Para.TStr_log = 'st0';
  tmp = fullfile(repo,'Tmp',['tmp_',Para.TStr]);
  if ~exist(tmp,'dir'), mkdir(tmp); end

  started = tic();
  [H,Id,Op] = InitHam(Para);
  Hmat = MPO2H(H.A)*exp(H.lgnorm);
  Hmat = (Hmat+Hmat')/2;
  hermiticity = max(max(abs(Hmat-Hmat')));
  evals = eig(Hmat);
  evals = real(evals);

  xtrg = XTRG_Main(Para,H,Id,Op,0,0);
  beta = xtrg.beta(:);
  exact_logz = zeros(size(beta));
  exact_energy = zeros(size(beta));
  exact_entropy = zeros(size(beta));
  e0 = min(evals);
  shifted = evals-e0;
  for i = 1:numel(beta)
    w = exp(-beta(i)*shifted);
    zshift = sum(w);
    exact_logz(i) = (-beta(i)*e0+log(zshift))/Para.L;
    exact_energy(i) = (e0+sum(w.*shifted)/zshift)/Para.L;
    exact_entropy(i) = exact_logz(i)+beta(i)*exact_energy(i);
  end
  xtrg_entropy = xtrg.LnZ(:)+beta.*xtrg.En(:);
  temperature_K = Jzz_K./beta;

  rows = cell(numel(beta),1);
  for i = 1:numel(beta)
    rows{i} = struct('temperature_K',temperature_K(i),'beta_dimensionless',beta(i), ...
      'exact_logZ_per_spin',exact_logz(i),'xtrg_logZ_per_spin',xtrg.LnZ(i), ...
      'logZ_error_per_spin',xtrg.LnZ(i)-exact_logz(i), ...
      'exact_energy_scaled_per_spin',exact_energy(i), ...
      'xtrg_energy_scaled_per_spin',xtrg.En(i), ...
      'energy_error_scaled_per_spin',xtrg.En(i)-exact_energy(i), ...
      'exact_entropy_kB_per_spin',exact_entropy(i), ...
      'xtrg_entropy_kB_per_spin',xtrg_entropy(i), ...
      'entropy_error_kB_per_spin',xtrg_entropy(i)-exact_entropy(i));
  end
  cold = find(temperature_K <= 0.12,1,'first');
  if isempty(cold), cold = numel(beta); end
  payload = struct('engine','public QMagen XTRG','upstream_commit', ...
    'dcf17db27682d4c505155682108c612bdd11274c','upstream_modified',false, ...
    'N',9,'geometry','3x3, OBC x / PBC y as public TLXXZ default', ...
    'model','M0 2026 parameters','field_T',field_T,'bond_dimension',D, ...
    'tau',Para.tau,'iterations',Para.It,'hermiticity_error',hermiticity, ...
    'elapsed_seconds',toc(started),'cold_validation_row',rows{cold}, ...
    'max_abs_entropy_error_kB_per_spin',max(abs(xtrg_entropy-exact_entropy)), ...
    'rows',{rows});

  output = fullfile(project_root,'results',sprintf('qmagen_n9_xtrg_D%d.json',D));
  fid = fopen(output,'w');
  if fid < 0, error('cannot open QMagen validation output'); end
  fprintf(fid,'%s',jsonencode(payload)); fclose(fid);
  fprintf('QMagen N9 D=%d: cold T=%.6f K, dS=%+.6g, elapsed=%.2f s\n', ...
          D,rows{cold}.temperature_K,rows{cold}.entropy_error_kB_per_spin, ...
          payload.elapsed_seconds);
end
