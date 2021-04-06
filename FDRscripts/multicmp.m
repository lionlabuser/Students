%  Function File: multicmp
%
%          USAGE: [padj,alpha] = multicmp (p,option,alpha)
%
%
%  A set of multiple comparison tools that apply procedures to control
%  either the family-wise error rate (FWER) or the false discovery rate
%  (FDR) when performing post-tests. The FWER is the probability of
%  making at least one false discovery (type 1 error) among the entire
%  set of hypotheses when multiple testing. In contrast, the FDR is the
%  expected proportion of false discoveries as a fraction the number of
%  null hypothesis rejected among a set of hypotheses. For example, for
%  an experimental set of 1000 t-tests, 50 of the tests (5 %) could be
%  false positive (<0.05). Typically, when one conducts a statistical
%  test one expects that a Type I error (one or more false positives)
%  will occur once in every 20 times the test is repeated. FWER control
%  maintains this interpretation for a set of multiple tests, such that
%  if one repeats the 1000-test experiment 20 times, at least one of the
%  repeats could have a false positive. FDR control permits some false
%  positives to occur at a tolerable level (e.g. 0.05). Thus, each time
%  one runs the 1000-test experiment, 5 % of the test results considered
%  significant will be declared false positive.
%
%  The default method for controlling the FWER is the universally
%  applicable step-down Bonferroni procedure of Holm. This function
%  reports multiplicity adjusted p-values [3]. Other method options
%  available are described in detail below. They are summarised in
%  the following list, starting with the most conservative, then
%  becoming more liberal as one descends:
%
%  - Family-wise Error Rate (FWER) control:
%
%    * Holm-Bonferroni (type: 'down')
%    * Hochberg-Bonferroni (type: 'up')
%
%  - False Discovery Rate (FDR) control
%
%    * Benjamini-Hochberg (type: 'fdr')
%
%
%  The input takes the following arguments
%
%   'p'       A column vector of raw p-values to correct for multiple
%             comparisons.
%
%   'type'    Sets the type of stepwise Bonferroni procedure ('down' or
%             'up') to control the FWER, or uses a procedure to control
%             the false discovery rate ('fdr').
%
%             By default, the type is 'down', which implements Holm's
%             universally applicable sequentially rejective step-down
%             Bonferroni procedure [4].
%
%             If type is set to 'up', Hochberg's step-up Bonferroni
%             procedure is implemented [5,6]. This procedure is less
%             conservative than Holm's procedure in controlling
%             the FWER across comparisons of positively correlated
%             measures.
%
%             If type is set to 'fdr', the Benjamini-Hochberg
%             procedure to control the false discovery rate is
%             implemented at a tolerance level set to alpha [7].
%             This procedure is less conservative than the step-wise
%             Bonferroni procedures but will not maintain the FWER.
%             The FDR is maintained where tests are independent
%             or positively correlated. The FDR is not maintained
%             when tests are negatively correlated.
%
%   'alpha'   A scalar value defining the alpha level to obtain
%             100*(1-alpha)% coverage of confidence intervals after
%             multiple testing.
%
%  The output takes the following arguments
%
%   'padj'    Vector of P-values adjusted for multiple comparisons.
%
%   'alpha'   Vector of adjusted alpha levels for setting corrected
%             confidence intervals. Adjusted alpha levels are available
%             for the Holm-Bonferroni [8] and the false discovery rate
%             procedures [9]. For the latter, adjusted confidence intervals
%             can only be made for tests that were significant after
%             controlling the FDR.
%
%
% The algorithms for adjusting p-values are summarised in reference [10]
%
%
%  [1]  Welch (1947) The generalization of "Student's" problem when
%        several different population variances are involved.
%        Biometrika. 34 (1–2): 28–35.
%  [2]  Ruxton (2006) The unequal variance t-test is an underused
%        alternative to Student's t-test and Mann-Whitney U test.
%        Behavioral Ecology. 17(4):688.
%  [3]  Wright (1992) Adjusted P-Values for Simultaneous Inference.
%        Biometrics. 48,1000-1013.
%  [4]  Holm (1979) A simple sequentially rejective multiple test
%        procedure. Scandanavian Journal of Statistics, 6:65-70
%  [5]  Holland and Copenhaver (1987) An improved sequentially
%        rejective Bonferroni test procedure. Biometrics. 43: 417-423
%  [6]  Hochberg (1988) A sharper Bonferroni procedure for multiple
%        tests of significance. Biometrika. 75: 800-802
%  [7]  Benjamini and Hochberg (1995) Controlling the False Discovery
%        Rate: A Practical and Powerful Approach to Multiple Testing.
%        Journal of the Royal Statistical Society. Series B. 57:289-300
%  [8]  Serlin, R. (1993). Confidence intervals and the scientific method:
%        A case for Holm on the range. Journal of Experimental Education.
%        61(4), 350–360.
%  [9]  Benjamini and Yekutieli (1995) False Discovery Rate–Adjusted Multiple
%        Confidence Intervals for Selected Parameters. JASA 100(469):71-81
%  [10] Westfall (1997) Multiple Testing of General Contrasts Using
%        Logical Constraints and Correlations. JASA 92(437):299-306
%
%  multicmp v1.1.0.1 (last updated: 02/01/2020)
%  Author: Andrew Charles Penn
%  https://www.researchgate.net/profile/Andrew_Penn/

%  Example data tested in Westfall (1997), JASA, 92(437): 299-306
%  p=[0.005708;0.023544;0.024193;0.044895;0.048805;0.221227;0.395867;0.693051;0.775755]


function [padj,alpha] = multicmp(p,type,alpha)

  % Allocate input arguments to function variables
  if size(p,2) == 1
    dim = 1;
  elseif size(p,1) == 1
    dim = 2;
  else
    error('Input p-values must be a vector');
  end
  if dim == 2
    p = p.';
  end

  % Unless option is specified, set multicmp to report Holm-Bonferroni adjusted p-values
  if nargin<2
    % Set default type to Holm-Bonferroni
    type = 'down';
  end
  if nargin<3
    % Set default central coverage of two-sided confidence intervals to 95%
    alpha = 0.05;
  end
  if nargout>2
    error('Invalid number of output arguments');
  end

  % Implement stepwise multiple comparison procedure
  switch lower(type)
    case 'down'
      if nargout > 1
        [padj,alpha] = holm(p,alpha);
      else
        padj = holm(p);
      end
    case 'up'
      if nargout > 1
        [padj,alpha] = hochberg(p);
      else
        padj = hochberg(p);
      end
    case 'fdr'
      if nargout > 1
        [padj,alpha] = fdr(p,alpha);
      else
        padj = fdr(p);
      end
  end

  if dim == 2
    padj = padj.';
    if nargout > 1
      alpha = alpha.';
    end
  end

end

%--------------------------------------------------------------------------

function [padj,alpha] = holm(p,alpha)

  % Function File: [padj,alpha] = holm(p,alpha)
  % Holm's step-down procedure
  %
  % References:
  % Serlin, R. (1993). Confidence intervals and the scientific method: A case
  % for Holm on the range. Journal of Experimental Education, 61(4), 350–360

  % Order raw p-values
  [p,idx] = sort(p,'ascend');

  % Initialize
  m = numel(p);
  padj = nan(m,1);
  alpha = alpha*ones(m,1);

  % Holm's step-down procedure
  padj(1) = m*p(1);
  alpha(1) = alpha(1)/m;
  for i = 2:m
    padj(i) = max(padj(i-1),(m-i+1)*p(i));
    alpha(i) = alpha(i)/(m-i+1);
  end
  R = find(padj>0.05,1,'first');
  alpha(R:end) = 0.05/(m-R+1);

  % Truncate adjusted p-values to 1.0
  padj(padj>1) = 1;

  % Return padj and alpha to the same order as p
  [sorted,idx] = sort(idx);
  padj = padj(idx);
  alpha = alpha(idx);

end

%--------------------------------------------------------------------------

function [padj,alpha] = hochberg(p)

  % Function File: [padj] = hochberg (p)
  % Hochberg's step-up procedure.

  % Order raw p-values
  [p,idx] = sort(p,'ascend');

  % Initialize
  m = numel(p);
  padj = nan(m,1);
  alpha = nan(m,1);

  % Hochberg's step-up procedure
  padj(m) = p(m);
  for j = 1:m-1
    i = m-j;
    padj(i) = min(padj(i+1),(m-i+1)*p(i));
  end

  % Truncate adjusted p-values to 1.0
  padj(padj>1) = 1;

  % Return padj and alpha to the same order as p
  [sorted,idx] = sort(idx);
  padj = padj(idx);

end

%--------------------------------------------------------------------------

function [padj,alpha] = fdr(p,alpha)

  % Function File: [padj,alpha] = fdr(p,alpha)
  % Benjamini-Hochberg procedure to control the false discovery rate

  % Order raw p-values
  [p,idx] = sort(p,'ascend');

  % Initialize
  q = alpha;
  m = numel(p);
  padj = nan(m,1);
  alpha = nan(m,1);

  % False discovery rate step-up procedure
  padj(m) = p(m);
  for j = 1:m-1
    i = m-j;
    padj(i) = min(padj(i+1),m/i*p(i));
  end

  % False coverage rate procedure
  i=(1:m).';
  R = find(p<i*q/m,1,'last');
  alpha(1:R) = R*q/m;

  % Truncate adjusted p-values to 1.0
  padj(padj>1) = 1;

  % Return padj and alpha to the same order as the input p values
  [sorted,idx] = sort(idx);
  padj = padj(idx);
  alpha = alpha(idx);

end

