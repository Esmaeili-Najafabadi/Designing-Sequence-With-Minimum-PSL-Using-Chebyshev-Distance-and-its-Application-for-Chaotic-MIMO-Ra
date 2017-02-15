%New algorithm for designining the Unimodular 
%First Algorithm CAP-PSL  step 2 is different (minimizing PSL is importanat
%here

function AKF = CAPPSL3smallestCircle

%Initializations

warning off

N = 13;
epsilon = 10^-12;
% P = 24;
% Q = 8;                  %We just simulate the P = Q value. 
% 
% T = [eye(4), zeros(4);...
%      zeros(16,8);...
%      zeros(4), eye(4)];
%  ind = [1:4,20:24]  ;        % is determined from the T matrix

P = 12;
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
iii = 1;
 
 while norm(Xtilda(:)-XtildaNew(:),'inf')>epsilon
           
    
      %For debugging
     NumberOfIterations = NumberOfIterations+1;
     norm(Xtilda-XtildaNew,'fro');
  % [ISL] = sidelobe(x);
   
     Xtilda = XtildaNew;
     
     
     [U1, Sigma, U2] = svd(Xtilda','econ');                 %equation 16
     U = U2*U1';                                            %equation 17
     
     clear mu_k
     %% Step 2 
     for n = 1:N
         S = 0;
         for i = 1:Q
               mu_k(i) = (sqrt(N))*U(n-1+ind(i),i);
               
         end

             [x(n) ] =   SmallestCircle(mu_k);               %equation 20
             iii = iii+1;
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

 