function [X,pdf] = deautocorrelation(pdf)
% Deautocorrelation by calculating zeros to determine
% the sign of fourier transform

%% 
% regularize pdf. let pdf is a vector N*1
[M,N] = size(pdf);
if(N>1)
   pdf = pdf';
else
   N = M; 
end

% make sure pdf is even, with N = 2N0+1 terms
if(mod(N,2)==0)
    N0 = N/2; 
else
    N0 = (N-1)/2;
end
pdf = [flipud(pdf(N0+2:end));pdf(N0+1:end)];

% Assume X is even, with 2N_X + 1 term, then conv(X,X) is 4N_X + 1 term.
% thus, N0 must be even, let N0 be even, we add 0.
N0 = (length(pdf)-1)/2;
if(mod(N0,2)~=0)
    pdf = [0;pdf;0];
end
N0 = (length(pdf)-1)/2;
pdf = pdf/sum(pdf);


%% Initialization - Determine the sign from fourier
z = roots(pdf);
N = length(pdf);
err_tolerance = 1e-5; % considering quantization error
z = z(abs(z)<(1+err_tolerance) & abs(z)>(1-err_tolerance)); % select
w = length(pdf)/2/pi*(angle(z)-sqrt(-1)*log(abs(z))); % zeros on the real axis
pdf_rotate = [pdf(N0+1:end);pdf(1:N0)]; % shift pdf
Y1 = real(fft(pdf_rotate)); % FFT. Omit quantization error, only contain real value

% Find the non-negative part
w = sort(real(w));
w = w(w>=0);
M = N0; % For the case if no sampling error
for i=1:length(w)
   num_wi = sum(w<w(i)+err_tolerance & w>w(i)-err_tolerance);
   if(mod(num_wi,2)~=0)
       M = floor(w(i));
      break; 
   end
end

Y2 = zeros(size(Y1));   % Construct new vector for Phi_X
Y2(1) = sqrt(Y1(1));
for i = 1:M % Assume on non-negative part
    num_zero = sum(w<i & w>=0); % calculate the number of roots in [0,a]
    if mod(num_zero,4)~=0
        curSign = -1;
    else
        curSign = 1;
    end
    Y2(1+i) = curSign*sqrt(Y1(1+i));
    Y2(1+i) = curSign*sqrt(Y1(1+i));
end
Y2(N0+2:end) = flipud(Y2(2:N0+1));  % Make it symmetric. values in [-M,-1]
N_X = N0/2;
X = real(ifft(Y2)); % may have some numerical problem, cause imaginary part
X = X(1:N_X+1);
X = [flipud(X(2:end));X];   % intercept values in [-N_X,N_X]
% figure;plot([-N_X:N_X]*binSize,X,'Linewidth',1.5);title('Root');
% figure;plot([-N0:N0]*binSize,pdf,'Linewidth',1.5);hold on; plot([-N0:N0]*binSize,conv(X,X),'Linewidth',1.5);legend('Target','Convolution');
end
