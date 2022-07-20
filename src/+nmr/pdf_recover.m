function [X_optimal] = pdf_recover(pdf,sample_size,weightType,nonincreasing)
% The recovered solution from deautocorrelation may not a true
% pdf, need refine step to move into the possibility space.
% weightType: 
%       'Average': assign same weight on each element of 'pdf'.
%       'Tail': assign more weight on the tails of 'pdf'.
%       'Center': assign more weight on the center part of 'pdf'.
% nonincreasing:
%       1: Set the refined distribution nonincreasing from origin.
%       0: No limitation on the monotonicity.


if(~exist('weightType','var'))
    weightType = 'Average';
end
if(~exist('nonincreasing','var'))
    nonincreasing = 0;
end
if(~exist('sample_size','var'))
    sample_size = 0;
end

%% recover pdf
X = nmr.deautocorrelation(pdf);
%% refine
% X_optimal = refine(pdf,X,weightType,nonincreasing);
X_optimal = nmr.refine_tik(pdf,X,sample_size,weightType,nonincreasing);

end