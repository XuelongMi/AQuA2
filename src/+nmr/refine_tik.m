function [X_optimal] = refine_tik(pdf,X,sample_size,weightType,nonincreasing)
% The recovered solution from deautocorrelation may not a true
% pdf, need refine step to move into the possibility space.
% weightType: 
%       'Average': assign same weight on each element of 'pdf'.
%       'Tail': assign more weight on the tails of 'pdf'.
%       'Center': assign more weight on the center part of 'pdf'.
% nonincreasing:
%       1: Set the refined distribution nonincreasing from origin.
%       0: No limitation on the monotonicity.

% errorbound = 7e-5;
errorbound = 10^(-0.15-0.5*log10(sample_size));
if(~exist('weightType','var'))
    weightType = 'Average';
end
if(~exist('nonincreasing','var'))
    nonincreasing = 0;
end

tolerrance = 1e-8;
if(sum(X(X<0))>=tolerrance && sum(X)<=1+tolerrance && sum(X)>=1-tolerrance)
    X_optimal = X;
else
%% Postprocess - iterative way to optimize the pdf.
N_X = (length(X)-1)/2;
if length(pdf)~=4*N_X+1
   pdf = [0;pdf;0]; 
end

X0 = X(N_X+1:end);
A = [];
b = [];
if(nonincreasing==1)
    A = zeros(length(X0)-1,length(X0)); % no inequality constraints
    for i = 1:length(X0)-1
       A(i,i) = -1;
       A(i,i+1) = 1;
    end
    b = zeros(length(X0)-1,1);
end
Aeq = [1,2*ones(1,N_X)];    % possibility constraint
beq = 1;
lb = zeros(size(X0));   % non-negative constraint
ub = [];
nonlcon = [];

weightMatrix = ones(size(pdf)); % defined by users
weightMatrix2 = 2*ones(1,N_X+1);
weightMatrix2(1) = 1;
switch weightType
    case 'Average'
        weightMatrix = ones(size(pdf));
    case 'Tail'
        weightMatrix = abs([-2*N_X:2*N_X])'+1;
        weightMatrix = weightMatrix/sum(weightMatrix)*(4*N_X+1);
    case 'Center'
        weightMatrix = abs([-2*N_X:2*N_X])';
        weightMatrix = (2*N_X+1-weightMatrix)/sum(weightMatrix)*(4*N_X+1);
end

% find local optimal
options = optimoptions(@fmincon,'Display','off','Algorithm','sqp','MaxIterations',1000);   
X_star = X0;
t1 = 1.01;
t2 = 3;
alpha = 0.5; %(from 0 to 1)
k = 2;
lastError = 1;
while(1)
    fun = @(x) sum(abs(pdf-conv([flipud(x(2:end));x],[flipud(x(2:end));x])).^2.*weightMatrix) + alpha*weightMatrix2*((x-X_star).^2); % objective function
    [X_new,val] = fmincon(fun,X0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    X_optimal = [flipud(X_new(2:end));X_new];
    error = sqrt(sum((pdf - conv(X_optimal,X_optimal)).^2));
    if(alpha>=1 || abs(lastError-error)/error < 1e-2 ||error>=t1*errorbound && error<=t2*errorbound)
        break;
    else
        if(error<t1*errorbound)
            % increase alpha
            alpha = alpha * k;
        else                            % error>t2*errorbound)
            % decrease alpha
            alpha = alpha / k;
        end
    end
    lastError = error;
end

X_optimal(X_optimal<0)=0;   % in case still has negative values
% figure;
% plot([-N_X:N_X],X,'LineWidth',1.5);
% hold on;
% plot([-N_X:N_X],X_optimal,'LineWidth',1.5);
% legend('Before refine','After refine');
% grid on;
end
end
