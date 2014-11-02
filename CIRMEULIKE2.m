function sumll = CIRMEULIKE2(para,Z, tau, nrow, ncol) %calculate log likelihood


% initialize the parameter for CIR model
theta = para(1); kappa = para(2); sigma = para(3); lambda = para(4);
theta2 = para(5); kappa2 = para(6); sigma2 = para(7); lambda2 = para(8);
alfa0 = para(9);
% sigmai = 0.005*ones(1,ncol);%volatility of measurement error
sigmai = para(10:end);
R = eye(ncol);
for i = 1:ncol
    R(i,i) = sigmai(i)^2;
end
dt = 1/12; %monthly data
initx = [theta ,theta2]';
initV = [sigma^2*theta/(2*kappa),0;0,sigma2^2*theta2/(2*kappa2)];

% parameter setting for transition equation
C = [theta*(1-exp(-kappa*dt)),theta2*(1-exp(-kappa2*dt))]' ;
F = [exp(-kappa*dt),0;0,exp(-kappa2*dt)];

% parameter setting for measurement equation
A = zeros(ncol, 1);
H = zeros(ncol,2);

for i = 1:ncol
    
    h1=sqrt(((kappa+lambda)^2)+2*sigma*sigma);
    h2=sqrt(((kappa2+lambda2)^2)+2*sigma2*sigma2);
    
    AffineAlpha = ((2*kappa*theta)/(sigma*sigma))*log((2*h1*exp(0.5*(kappa+lambda+h1)*tau(i)))/(h1-(kappa+lambda)+(kappa+lambda+h1)*exp(h1*tau(i))));
    AffineAlpha2 = ((2*kappa2*theta2)/(sigma2*sigma2))*log((2*h2*exp(0.5*(kappa2+lambda2+h2)*tau(i)))/(h2-(kappa2+lambda2)+(kappa2+lambda2+h2)*exp(h2*tau(i))));
    
    AffineBeta = (2*(exp(h1*tau(i))-1))/(h1-(kappa+lambda)+(kappa+lambda+h1)*exp(h1*tau(i)));
    AffineBeta2 = (2*(exp(h2*tau(i))-1))/(h2-(kappa2+lambda2)+(kappa2+lambda2+h2)*exp(h2*tau(i)));

    A(i) = alfa0-(AffineAlpha/tau(i))-(AffineAlpha2/tau(i));
    H(i,1) = AffineBeta/tau(i);
    H(i,2) = AffineBeta2/tau(i);
end

%now recursive steps
EY = initx;
VARY = initV;
ll = zeros(nrow,1); %log-likelihood


EZ = A+H*EY;
VARZ = H*VARY*H'+R;




Erro = Z(1,:)-EZ';
Erro=Erro';
K = VARY*H'*inv(VARZ);
EY = EY +K*Erro;

if EY(1) < 0
        EY(1) =0; 
        
    else
        EY(1) =EY(1);
end
    
    if EY(2) < 0
        EY(2) =EY(2);%0
        
    else
        EY(2) =EY(2);
    end
VARY = (eye(2)-K*H)*VARY;
ll(1) = -(ncol/2)*log(2*pi)-0.5*log(det(VARZ))-0.5*Erro'*inv(VARZ)*Erro;

estado1(1)=EY(1);
estado2(1)=EY(2);


for i = 2:nrow
    
    FUTUROY = C+F*EY;
    
    if FUTUROY(1) < 0
        FUTUROY(1) =0; 
        
    else
        FUTUROY(1) =FUTUROY(1);
    end
    
    if FUTUROY(2) < 0
        FUTUROY(2) =FUTUROY(2); %0
        
    else
        FUTUROY(2) =FUTUROY(2);
    end
    
    Q1 = ((sigma*sigma)/kappa)*(1-exp(-kappa*dt))*((theta/2)*(1-exp(-kappa*dt))+exp(-kappa*dt)*FUTUROY(1));
    Q2 = ((sigma2*sigma2)/kappa2)*(1-exp(-kappa2*dt))*((theta2/2)*(1-exp(-kappa2*dt))+exp(-kappa2*dt)*FUTUROY(2));
    Q = [Q1,0;0,Q2];
    
%     VARY = VARY-F*VARY*F'+Q;
VARY = F*VARY*F'+Q;
    EZ = A+H*FUTUROY;
    
    VARZ = H*VARY*H'+R;
    Erro = Z(i,:)-EZ';
    Erro=Erro';
    K = VARY*H'*inv(VARZ);
    FUTUROY = FUTUROY +K*Erro;
    
    if FUTUROY(1) < 0
        FUTUROY(1) =0; 
        
    else
        FUTUROY(1) =FUTUROY(1);
    end
    
    if FUTUROY(2) < 0
        FUTUROY(2) =FUTUROY(2); %0
        
    else
        FUTUROY(2) =FUTUROY(2);
    end
    
    VARY = (eye(2)-K*H)*VARY;
    EY = FUTUROY;
    
    estado1(i) = EY(1);
    estado2(i) = EY(2);
    
    
    ll(i) = -(ncol/2)*log(2*pi)-0.5*log(det(VARZ))-0.5*Erro'*inv(VARZ)*Erro;
    
   
end

estado = [estado1' estado2'];
csvwrite('novo.csv',estado)

sumll = -sum(ll)

end
