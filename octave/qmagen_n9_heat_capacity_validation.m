function payload = qmagen_n9_heat_capacity_validation(D)
  % Four-stagger XTRG validation of QMagen's heat-capacity postprocessing.
  if nargin < 1, D = 24; end
  script_dir = fileparts(mfilename('fullpath'));
  project_root = fileparts(script_dir);
  repo = fullfile(project_root,'sources','raw','QMagen_public','V0.0');
  old_dir = pwd(); cleanup = onCleanup(@() cd(old_dir));
  cd(fullfile(repo,'RunScript'));
  addpath(fullfile(repo,'UtilityFunc'));
  addpath(genpath(fullfile(repo,'SpinModel')));
  addpath(fullfile(repo,'ManyBodySolver'));
  addpath(fullfile(repo,'ManyBodySolver','TensorFunction'));
  addpath(genpath(fullfile(repo,'ManyBodySolver','XTRG')));
  addpath(genpath(fullfile(repo,'svd_lapack_interface')));
  addpath(fullfile(script_dir,'qmagen_octave_compat'),'-begin');
  global EVOFLAG EVO_check PLOTFLAG SAVEFLAG;
  EVOFLAG=0; EVO_check=0; PLOTFLAG=0; SAVEFLAG=0;

  Jzz_K=0.1225*11.60451812155; gz=4.716; field_T=1.50;
  Para=struct('ManyBodySolver','XTRG','IntrcMap_Name','IntrcMap_TLXXZ', ...
    'd',2,'L',9,'Geo',struct('Lx',3,'Ly',3,'BCX','OBC','BCY','PBC','L',9), ...
    'Model',struct('J1xy',0.636,'J1z',1.0,'J2xy',0.0,'J2z',0.0, ...
                   'ES',Jzz_K,'g',[4.2,4.2,gz]), ...
    'Field',struct('h',[0,0,gz*0.67171381563*field_T/Jzz_K], ...
                   'B',[0,0,field_T]),'ThDQ','Cm','Parallel',false, ...
    'RunParallel',false,'RunParallelNum',1,'RunParallelmaxNumCompThreads',1, ...
    'tau',0.00025,'It',15,'D_list',D*ones(1,15),'MCrit',D, ...
    'Ver','Memory','VariProd_step_max',50,'VariSum_step_max',50, ...
    'SETTN_init_step_max',50,'TStr',sprintf('qmagen_n9_C_D%d',D),'TStr_log','st');
  tmp=fullfile(repo,'Tmp',['tmp_',Para.TStr]);if ~exist(tmp,'dir'),mkdir(tmp);end

  started=tic(); [H,Id,Op]=InitHam(Para);
  Hmat=MPO2H(H.A)*exp(H.lgnorm);Hmat=(Hmat+Hmat')/2;
  evals=real(eig(Hmat));e0=min(evals);shifted=evals-e0;
  beta=[];lnz=[];en=[];entropy=[];
  for st=0:3
    x=XTRG_Main(Para,H,Id,Op,st,0);
    beta=[beta;x.beta(:)];lnz=[lnz;x.LnZ(:)];en=[en;x.En(:)];
  end
  [beta,order]=sort(beta);lnz=lnz(order);en=en(order);
  entropy=lnz+beta.*en;
  rg=struct('beta',beta,'En',en);
  [cx,beta_c]=ThDQ_func(rg,Para);
  exact_c=zeros(size(beta_c));
  for i=1:numel(beta_c)
    w=exp(-beta_c(i)*shifted);z=sum(w);emean=e0+sum(w.*shifted)/z;
    var=sum(w.*((evals-emean).^2))/z;
    exact_c(i)=beta_c(i)^2*var/Para.L;
  end
  rows=cell(numel(beta_c),1);
  for i=1:numel(beta_c)
    rows{i}=struct('temperature_K',Jzz_K/beta_c(i), ...
      'beta_dimensionless',beta_c(i),'exact_C_kB_per_spin',exact_c(i), ...
      'xtrg_C_kB_per_spin',cx(i),'C_error_kB_per_spin',cx(i)-exact_c(i));
  end
  band=find(Jzz_K./beta_c>=0.08 & Jzz_K./beta_c<=0.30);
  payload=struct('engine','public QMagen four-stagger XTRG','N',9, ...
    'geometry','3x3 OBCx/PBCy','model','M0 2026 parameters','field_T',field_T, ...
    'bond_dimension',D,'elapsed_seconds',toc(started), ...
    'max_abs_C_error_80_to_300mK',max(abs(cx(band)-exact_c(band))), ...
    'rmse_C_error_80_to_300mK',sqrt(mean((cx(band)-exact_c(band)).^2)), ...
    'rows',{rows});
  out=fullfile(project_root,'results',sprintf('qmagen_n9_heat_capacity_D%d.json',D));
  fid=fopen(out,'w');fprintf(fid,'%s',jsonencode(payload));fclose(fid);
  fprintf('QMagen heat capacity D=%d: RMSE(80-300mK)=%.6g, max=%.6g, elapsed=%.1fs\n', ...
    D,payload.rmse_C_error_80_to_300mK,payload.max_abs_C_error_80_to_300mK, ...
    payload.elapsed_seconds);
end
