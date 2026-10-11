function [U,S,V,t_error_relative,entropy] = svdT(A,varargin)
  % Compatibility copy of public QMagen svdT.m.
  % The upstream code tests find(...,1) but then calls find(...) without the
  % scalar limit.  MATLAB/Octave consequently turn D_max into a vector when
  % several singular values are below epsilon.  This copy changes only that
  % index selection and guarantees at least one retained singular value.
  try
    [U,S,V] = svd(A,'econ');
  catch
    [U,S,V] = svd_lapack(A,'econ','gesvd');
  end
  if isnan(norm(U)), [U,S,V] = svd_lapack(A,'econ','gesvd'); end
  S = diag(S); ns = norm(S); t_error_relative = 0; entropy = 0;
  if ~isempty(varargin)
    D_max = length(S); epsilon = 0;
    switch varargin{1}
      case 'Nkeep', D_max = min(D_max,varargin{2});
      case 'epsilon', epsilon = varargin{2};
    end
    if length(varargin) == 4
      switch varargin{3}
        case 'Nkeep', D_max = min(D_max,varargin{4});
        case 'epsilon', epsilon = varargin{4};
      end
    end
    below = find(S/ns < epsilon,1);
    if ~isempty(below), D_max = min(below-1,D_max); end
    D_max = max(1,D_max);
    Sp = S/ns;
    S = S(1:D_max);
    t_error_relative = sum(Sp(D_max+1:end).^2);
    entropy = -sum(Sp.^2 .* log(Sp.^2));
    U = U(:,1:D_max); V = V(:,1:D_max);
  end
  V = V';
end
