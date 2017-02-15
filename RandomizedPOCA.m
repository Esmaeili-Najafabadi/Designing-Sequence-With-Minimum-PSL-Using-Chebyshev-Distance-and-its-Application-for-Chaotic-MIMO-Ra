
%Initializations
clc
clear
close all

N = 100;
epsilon = .2*10^-14;           % do not change this value
% P = 80;
% Q = 36;                  %We just simulate the P = Q value. 
% 
% T = [eye(26), zeros(26,10);...
%      zeros(44,36);...
%      zeros(10,26), eye(10)];
%  ind = [1:26,71:80]  ;        % is determined from the T matrix

P = 40;
Q = P;                  %We just simulate the P = Q value. 

T = eye(P);
 ind = [1:P]  ;        % is determined from the T matrix

%Golomb Siquence
golomb = zeros(N,1);
for n = 1:N
        golomb(n) = exp(1i*pi*(n-1)*n/N);
end
 x = golomb;
% x = rand(N,1); 
NumberOfIterations = 0;

 
 %Constructing Xhat
Xhat = zeros(N+P-1,P);
for i = 1:N+P-1
    for j = 1:P
        if ((i-j+1 <= N) && (i-j+1>0))
            Xhat(i,j) = x(i-j+1);
        end
    end
end

XtildaNew = Xhat*T;
Xtilda = zeros(N+P-1,Q);
x = zeros(N,1);
k = 4;
Omega = randn(N+P-1,k);
 
 while norm(Xtilda(:)-XtildaNew(:),'inf')>epsilon
           
    
      %For debugging
     NumberOfIterations = NumberOfIterations+1
     norm(Xtilda-XtildaNew,'fro')
  % [ISL] = sidelobe(x);
   
     Xtilda = XtildaNew;
     
     
%      [U1, Sigma, U2] = svd(Xtilda','econ');                 %equation 16
  
    A = Xtilda';
    Y = A*Omega;
    [QQ,R] = qr(Y);
 
    B = QQ'*A;
    [Uhat, Sigma, U2] = svd(B,'econ');

    U1 = QQ*Uhat; 
%     [U1, Sigma, U2] = svdfast(Xtilda',4);                 %equation 16
     U = U2*U1';                                            %equation 17
     
     clear mu_k
     %% Step 2 
     for n = 1:N
         S = 0;
         for i = 1:Q
               mu_k(i) = (sqrt(N))*U(n-1+ind(i),i);
               
         end

             x(n) =   (max(mu_k)+min(mu_k))/2;                  %equation 20
     end
   
     % Constructing Xtilda from x
    Xhat = zeros(N+P-1,P);
    for i = 1:N+P-1
        for j = 1:P
            if ((i-j+1 <= N) && (i-j+1>0))
                Xhat(i,j) =x(i-j+1);
            end
        end
    end
    XtildaNew = Xhat*T;
   
 end
 
 %% Computing Autocorrelation function
 AKF = zeros(2*N-1,1);
 for k = 0:N-1
    S = 0;
    for  n = k+1:N
        S = S+x(n)*x(n-k)';
    end
        AKF(N+k) = S;
        AKF(N-k) = S';
 end     
AKF = AKF/AKF(N);

figure
plot([-N+1:N-1],db(abs(AKF)),'LineWidth',1.2);
grid minor;
xlabel('lag k');
ylabel('Autocorrelation (dB)');
print '-depsc2' RPmCAPPSL2TA1000 
vec = AKF(N+ind(1:end-1));
r0 = AKF(N);
MMF = N^2/(2*sum((abs(vec)).^2))
MPCL = max(abs(vec/r0))

 