% This program will compute the smallest cicle according to quadratic
% programming

function [center] = SmallestCircle(Mu)

delta = 7;

xR0 = real((max(Mu)+min(Mu))./2);
xI0 = imag((max(Mu)+min(Mu)./2));
z = 0;
H = 2*[1 0 0;0 1 0; 0 0 0];
f = [0 0 -1];
A = ones(length(Mu),3);
bb = ones(length(Mu),1);
for i = 1:length(Mu)
   A(i,1) = -2*real(Mu(i));
   A(i,2) = -2*imag(Mu(i));
   bb(i) = -real(Mu(i))^2-imag(Mu(i))^2;
end

xR = xR0; 
xI = xI0;
R = max(abs(xR+1i*xI-Mu));
while abs(xR^2+xI^2-z)>delta
%      abs(xR^2+xI^2-z)
    R = max(abs(xR+1i*xI-Mu));
    
    [x] = quadprog(H,f,A,R+bb);


    xR = x(1);
    xI = x(2);
    z = x(3);
end

center = xR+1i*xI;


end