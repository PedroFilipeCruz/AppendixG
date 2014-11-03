
%kalman filter
function sumll = CIRDefaultent(para,input, datescup,parametros,cupao) %calculate log likelihood

estado1 = input(:,3);
estado2 = input(:,4);

Z = input(:,2)/100;

datesobs = input(:,1);
datesobsmatlab = datesobs + datenum('30-Dec-1899');
datescupmatlab = datescup + datenum('30-Dec-1899');
[nrow,ncol] = size(datesobsmatlab); %numero de observações e de obrigacoes


theta3 = para(1); k3 = para(2); sigma3 = para(3); eta3 = para(4); gama0 = para(5);
beta1 = para(6); beta2 = para(7);
sigmai = para(8);

theta1 = parametros(1,1) ; theta2 = parametros(1,2);
k1 = parametros(2,1) ; k2 = parametros(2,2); 
sigma1 = parametros(3,1) ; sigma2 = parametros(3,2); 
eta1 = parametros(4,1); eta2 = parametros(4,2); 
alfa0 = parametros(1,3);


%n de cupoes
ncupaux = size(datescup);
ncup = ncupaux(1);

h1=sqrt(((k1+eta1)^2)+2*(1+beta1)*sigma1*sigma1);
h2=sqrt(((k2+eta2)^2)+2*(1+beta2)*sigma2*sigma2);
h3=sqrt(((k3+eta3)^2)+2*sigma3*sigma3);


AffineAlpha1 = @(tau) ((2*k1*theta1)/(sigma1*sigma1))*log((2*h1*exp(0.5*(k1+eta1+h1)*tau))/(h1-(k1+eta1)+(k1+eta1+h1)*exp(h1*tau)));
AffineAlpha2 = @(tau) ((2*k2*theta2)/(sigma2*sigma2))*log((2*h2*exp(0.5*(k2+eta2+h2)*tau))/(h2-(k2+eta2)+(k2+eta2+h2)*exp(h2*tau)));
AffineAlpha3 = @(tau) ((2*k3*theta3)/(sigma3*sigma3))*log((2*h3*exp(0.5*(k3+eta3+h3)*tau))/(h3-(k3+eta3)+(k3+eta3+h3)*exp(h3*tau)));


AffineBeta1 = @(tau) (2*(1+beta1)*(exp(h1*tau)-1))/(h1-(k1+eta1)+(k1+eta1+h1)*exp(h1*tau));
AffineBeta2 = @(tau) (2*(1+beta1)*(exp(h2*tau)-1))/(h2-(k2+eta2)+(k2+eta2+h2)*exp(h2*tau));
AffineBeta3 = @(tau) (2*(exp(h3*tau)-1))/(h3-(k3+eta3)+(k3+eta3+h3)*exp(h3*tau));


f = @(t,valor1,valor2) exp(-alfa0*t+AffineAlpha1(t)+AffineAlpha2(t)-AffineBeta1(t)*valor1-AffineBeta2(t)*valor2);
g = @(t,valor3) exp(-gama0*t+AffineAlpha3(t)-AffineBeta3(t)*valor3);


R = eye(ncol);
for i = 1:ncol
    R(i,i) = sigmai(i)^2;
end
dt = 1/12; 
initx = theta3 ;
initV = sigma3^2*theta3/(2*k3);


C = theta3*(1-exp(-k3*dt)) ;
F = exp(-k3*dt);


EY = initx;
VARY = initV;
ll = zeros(nrow,1); %log-likelihood

dates = [];
Aaux = 0;
Haux = 0;
for i=1:ncup
   if datescupmatlab(i) > datesobsmatlab(1)
       dates = [dates,datescupmatlab(i)];
   end

end

numdatesaux = size(dates);
numdates = numdatesaux(2);



for j=1:numdates
Aaux=Aaux+cupao*f((ceil((dates(j)-datesobsmatlab(1))/31))/12,estado1(1),estado2(1))*g((ceil((dates(j)-datesobsmatlab(1))/31))/12,EY)+cupao*EY*AffineBeta3((ceil((dates(j)-datesobsmatlab(1))/31))/12)*f((ceil((dates(j)-datesobsmatlab(1))/31))/12,estado1(1),estado2(1))*g((ceil((dates(j)-datesobsmatlab(1))/31))/12,EY);
Haux = Haux - AffineBeta3((ceil((dates(j)-datesobsmatlab(1))/31))/12)*cupao*f((ceil((dates(j)-datesobsmatlab(1))/31))/12,estado1(1),estado2(1))*g((ceil((dates(j)-datesobsmatlab(1))/31))/12,EY);
end
A(1)=Aaux+(1+AffineBeta3((ceil((dates(end)-datesobsmatlab(1))/31))/12)*EY)*f((ceil((dates(end)-datesobsmatlab(1))/31))/12,estado1(1),estado2(1))*g((ceil((dates(end)-datesobsmatlab(1))/31))/12,EY);
H(1) = Haux - AffineBeta3((ceil((dates(end)-datesobsmatlab(1))/31))/12)*f((ceil((dates(end)-datesobsmatlab(1))/31))/12,estado1(1),estado2(1))*g((ceil((dates(end)-datesobsmatlab(1))/31))/12,EY);
% Aaux = 0;
Haux = 0;
%dates =[];




 EZ = A(1)+H(1)*EY;
% EZ = Aaux + f((ceil((dates(end)-datesobsmatlab(1))/31))/12,estado1(1),estado2(1),estado3(1))*g((ceil((dates(end)-datesobsmatlab(1))/31))/12,EY);
Aaux=0; %retirar
dates =[];%retirar
VARZ = H(1)*VARY*H(1)'+R;
Erro = Z(1,:)-EZ';
Erro=Erro';
K = VARY*H(1)'*inv(VARZ);
EY = EY +K*Erro;

 if EY < 0
        EY =0;
        
    else
        EY = EY;
    end

VARY = (eye(1)-K*H(1))*VARY;
ll(1) = -(ncol/2)*log(2*pi)-0.5*log(det(VARZ))-0.5*Erro'*inv(VARZ)*Erro;

estado3(1)=EY;



for i = 2:nrow
    
    FUTUROY = C+F*EY;
    
     if FUTUROY < 0
        FUTUROY =0;
        
    else
        FUTUROY = FUTUROY;
     end
    
 %%%%%%%%%parte nova    
   for o=1:ncup
   if datescupmatlab(o) > datesobsmatlab(i)
       dates = [dates,datescupmatlab(o)];
   end

end

numdatesaux = size(dates);
numdates = numdatesaux(2);


for j=1:numdates
Aaux=Aaux+cupao*f((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i))*g((ceil((dates(j)-datesobsmatlab(i))/31))/12,FUTUROY)+cupao*FUTUROY*AffineBeta3((ceil((dates(j)-datesobsmatlab(i))/31))/12)*f((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i))*g((ceil((dates(j)-datesobsmatlab(i))/31))/12,FUTUROY);
Haux = Haux - AffineBeta3((ceil((dates(j)-datesobsmatlab(i))/31))/12)*cupao*f((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i))*g((ceil((dates(j)-datesobsmatlab(i))/31))/12,FUTUROY);
end
A(i)=Aaux+(1+AffineBeta3((ceil((dates(end)-datesobsmatlab(i))/31))/12)*FUTUROY)*f((ceil((dates(end)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i))*g((ceil((dates(end)-datesobsmatlab(i))/31))/12,FUTUROY);
H(i) = Haux - AffineBeta3((ceil((dates(end)-datesobsmatlab(i))/31))/12)*f((ceil((dates(end)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i))*g((ceil((dates(end)-datesobsmatlab(i))/31))/12,FUTUROY);
%Aaux = 0;
Haux = 0;
%dates =[];      
     
     
%     Q = (theta*sigma*sigma/(2*kappa))*(1-exp(-kappa*dt))*(1-exp(-kappa*dt)) + ((sigma*sigma)/(kappa))*(exp(-kappa*dt)-exp(-2*kappa*dt))*FUTUROY;
   % Q = ((sigma*sigma)/kappa)*(1-exp(-kappa*dt))*((theta/2)*(1-exp(-kappa*dt))+exp(-kappa*dt)*FUTUROY(1));
    Q = theta3*sigma3*sigma3*(1-exp(-k3*dt))^2/(2*k3)+sigma3*sigma3/k3*(exp(-k3*dt)-exp(-2*k3*dt))*FUTUROY;
    VARY = VARY-F*VARY*F'+Q;
 EZ = A(i)+H(i)*FUTUROY;
% EZ = Aaux + f((ceil((dates(end)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i),estado3(i))*g((ceil((dates(end)-datesobsmatlab(i))/31))/12,FUTUROY);
Aaux =0; %retirar
dates =[]; %retirar
    VARZ = H(i)*VARY*H(i)'+R;
    Erro = Z(i,:)-EZ';
    Erro=Erro';
    K = VARY*H(i)'*inv(VARZ);
    FUTUROY = FUTUROY +K*Erro;
    
     if FUTUROY < 0
        FUTUROY =0;
        
    else
        FUTUROY = FUTUROY;
    end
    
    
    VARY = (eye(1)-K*H(i))*VARY;
    EY = FUTUROY;
    
    
    
    
    estado3(i) = EY;
    
    
    
    ll(i) = -(ncol/2)*log(2*pi)-0.5*log(det(VARZ))-0.5*Erro'*inv(VARZ)*Erro;
    
   
end

estado = estado3';
csvwrite('default.csv',estado)

sumll = -sum(ll)
 
end



