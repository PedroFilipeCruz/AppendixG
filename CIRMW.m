function sumll = CIRMW(para,input, datescup,parametros,cupao,mw) %calculate log likelihood

estado1 = input(:,3);
estado2 = input(:,4);

estado4 = input(:,5); %este é o de default

Z = input(:,2)/100;

datesobs = input(:,1);
datesobsmatlab = datesobs + datenum('30-Dec-1899');
datescupmatlab = datescup + datenum('30-Dec-1899');
[nrow,ncol] = size(datesobsmatlab); %numero de observações e de obrigacoes


theta5 = para(1); k5 = para(2); sigma5 = para(3); eta5 = para(4); epsilon0 = para(5);
epsilon4 = para(6); sigmai = para(7);

theta1 = parametros(1,1) ; theta2 = parametros(1,2);  theta4 = parametros(1,3);
k1 = parametros(2,1) ; k2 = parametros(2,2); k4 = parametros(2,3);
sigma1 = parametros(3,1) ; sigma2 = parametros(3,2);  sigma4 = parametros(3,3);
eta1 = parametros(4,1); eta2 = parametros(4,2);  eta4 = parametros(4,3);
alfa0 = parametros(1,4); gama0 = parametros(2,4); beta1 = parametros(3,4); beta2 = parametros(4,4);


%n de cupoes
ncupaux = size(datescup);
ncup = ncupaux(1);

h1=sqrt(((k1+eta1)^2)+(1+beta1)*sigma1*sigma1);
h2=sqrt(((k2+eta2)^2)+(1+beta2)*sigma2*sigma2);
h5=sqrt(((k5+eta5)^2)+2*sigma5*sigma5);

AffineAlpha1 = @(tau) ((2*k1*theta1)/(sigma1*sigma1))*log((2*h1*exp(0.5*(k1+eta1+h1)*tau))/(h1-(k1+eta1)+(k1+eta1+h1)*exp(h1*tau)));
AffineAlpha2 = @(tau) ((2*k2*theta2)/(sigma2*sigma2))*log((2*h2*exp(0.5*(k2+eta2+h2)*tau))/(h2-(k2+eta2)+(k2+eta2+h2)*exp(h2*tau)));

AffineAlpha5 = @(tau) ((2*k5*theta5)/(sigma5*sigma5))*log((2*h5*exp(0.5*(k5+eta5+h5)*tau))/(h5-(k5+eta5)+(k5+eta5+h5)*exp(h5*tau)));

AffineBeta1 = @(tau) (2*(1+beta1)*(exp(h1*tau)-1))/(h1-(k1+eta1)+(k1+eta1+h1)*exp(h1*tau));
AffineBeta2 = @(tau) (2*(1+beta2)*(exp(h2*tau)-1))/(h2-(k2+eta2)+(k2+eta2+h2)*exp(h2*tau));

AffineBeta5 = @(tau) (2*(exp(h5*tau)-1))/(h5-(k5+eta5)+(k5+eta5+h5)*exp(h5*tau));

f = @(t,valor1,valor2) exp(-(alfa0+gama0)*t+AffineAlpha1(t)+AffineAlpha2(t)-AffineBeta1(t)*valor1-AffineBeta2(t)*valor2);
g = @(t,valor5) exp(-epsilon0*t+AffineAlpha5(t)-AffineBeta5(t)*valor5);

kimmela = (2*(theta4^2)*(k4^2))/(sigma4^4)-(2*theta4*k4)/(sigma4^2)+(4*mw*epsilon4)/(sigma4^2) +(3/8);
kimmelb = 0.5*sqrt(((k4+eta4)^2)+2*(sigma4^2));
kimmeld = (-1)*((theta4*k4*(k4*eta4))/(sigma4^2));
kimmelalfa = ((2*theta4*k4)/(sigma4^2))-0.5;
kimmelgama = 0.5*(1-sqrt(1+8*kimmela));
kimmeldelta = @(t) 1-exp(-2*kimmelb*t);
kimmely = @(valor4) (2*sqrt(valor4))/(sigma4);
kimmelz = @(t,valor4) sqrt(2*kimmelb)*exp(-kimmelb*t)*kimmely(valor4);

auxiliarpi1= @(t,valor4) ((kimmelalfa-kimmelgama)*(kimmelalfa+kimmelgama-1)/(2*(kimmelz(t,valor4))))+((2*kimmelalfa+1)/4)*(1-((k4+eta4)/(2*kimmelb)))+(kimmelz(t,valor4)^2)*(1/8)*((1-((k4+eta4)/(2*kimmelb)))^2);
auxiliarpi2 = @(t,valor4) ((kimmelalfa-kimmelgama)*(kimmelalfa-kimmelgama-2)*(kimmelalfa+kimmelgama-1)*(kimmelalfa+kimmelgama-3))/(4*(kimmelz(t,valor4)^4));
auxiliarpi3 = @(t,valor4) ((2*kimmelalfa-1)*(kimmelalfa-kimmelgama)*(kimmelalfa+kimmelgama-1)*(1-((k4+eta4)/(2*kimmelb))))/(4*kimmelz(t,valor4)^2);
auxiliarpi4 = @(t,valor4) ((1/16)*(2*kimmelalfa+3)*(2*kimmelalfa+1)+(1/8)*(kimmelalfa-kimmelgama)*(kimmelalfa+kimmelgama-1))*((1-((k4+eta4)/(2*kimmelb)))^2);
auxiliarpi5 = @(t,valor4) (1/16)*(2*kimmelalfa+3)*(kimmelz(t,valor4)^2)*((1-((k4+eta4)/(2*kimmelb)))^3)+(1/64)*(kimmelz(t,valor4)^4)*((1-((k4+eta4)/(2*kimmelb)))^4);

kimmelpi = @(t,valor4) 1+kimmeldelta(t)*auxiliarpi1(t,valor4)+0.5*(kimmeldelta(t)^2)*(auxiliarpi2(t,valor4)+auxiliarpi3(t,valor4)+auxiliarpi4(t,valor4)+auxiliarpi5(t,valor4));
kimmelwum = @(t,valor4) ((kimmelz(t,valor4)/(sqrt(2*kimmelb)))^(kimmelalfa-kimmelgama))*exp((1/4)*(1-((k4+eta4)/(2*kimmelb)))*(kimmelz(t,valor4)^2))*kimmelpi(t,valor4);
kimmelh = @(t,valor4) ((kimmelz(t,valor4)/(sqrt(2*kimmelb)))^kimmelgama)*exp(-0.5*kimmelb*(kimmely(valor4)^2)-(0.5*kimmelb+kimmeld)*t)*kimmelwum(t,valor4);

kimmelpifinal = @(t,valor4) ((4*valor4/(sigma4^2))^(0.25-((theta4*k4)/(sigma4^2))))*exp((valor4*(k4+eta4))/(sigma4^2))*kimmelh(t,valor4);


R = eye(ncol);
for i = 1:ncol
    R(i,i) = sigmai(i)^2;
end
dt = 1/12; %monthly data
initx = theta5 ;
initV = sigma5^2*theta5/(2*k5);

% parameter setting for transition equation
C = theta5*(1-exp(-k5*dt)) ;
F = exp(-k5*dt);

% parameter setting for measurement equation

 
%EKF
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


%Agora e ao contrário das treasury yields os A e H's vão ser diferentes
%para cada observação uma vez que estamos a actualizar diferente numero de
%cupões
for j=1:numdates
Aaux=Aaux+cupao*kimmelpifinal((ceil((dates(j)-datesobsmatlab(1))/31))/12,estado4(1))*f((ceil((dates(j)-datesobsmatlab(1))/31))/12,estado1(1),estado2(1))*g((ceil((dates(j)-datesobsmatlab(1))/31))/12,EY)+cupao*EY*kimmelpifinal((ceil((dates(j)-datesobsmatlab(1))/31))/12,estado4(1))*AffineBeta5((ceil((dates(j)-datesobsmatlab(1))/31))/12)*f((ceil((dates(j)-datesobsmatlab(1))/31))/12,estado1(1),estado2(1))*g((ceil((dates(j)-datesobsmatlab(1))/31))/12,EY);
 Haux = Haux - AffineBeta5((ceil((dates(j)-datesobsmatlab(1))/31))/12)*cupao*kimmelpifinal((ceil((dates(j)-datesobsmatlab(1))/31))/12,estado4(1))*f((ceil((dates(j)-datesobsmatlab(1))/31))/12,estado1(1),estado2(1))*g((ceil((dates(j)-datesobsmatlab(1))/31))/12,EY);
end
A(1)=Aaux+(1+AffineBeta5((ceil((dates(end)-datesobsmatlab(1))/31))/12)*EY)*kimmelpifinal((ceil((dates(end)-datesobsmatlab(1))/31))/12,estado4(1))*f((ceil((dates(end)-datesobsmatlab(1))/31))/12,estado1(1),estado2(1))*g((ceil((dates(end)-datesobsmatlab(1))/31))/12,EY);
H(1) = Haux - AffineBeta5((ceil((dates(end)-datesobsmatlab(1))/31))/12)*kimmelpifinal((ceil((dates(end)-datesobsmatlab(1))/31))/12,estado4(1))*f((ceil((dates(end)-datesobsmatlab(1))/31))/12,estado1(1),estado2(1))*g((ceil((dates(end)-datesobsmatlab(1))/31))/12,EY);
% Aaux = 0;
Haux = 0;
%dates =[];




 EZ = A(1)+H(1)*EY;
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

estado5(1)=EY;



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
Aaux=Aaux+cupao*kimmelpifinal((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado4(i))*f((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i))*g((ceil((dates(j)-datesobsmatlab(i))/31))/12,FUTUROY)+cupao*FUTUROY*kimmelpifinal((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado4(i))*AffineBeta5((ceil((dates(j)-datesobsmatlab(i))/31))/12)*f((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i))*g((ceil((dates(j)-datesobsmatlab(i))/31))/12,FUTUROY);
Haux = Haux - AffineBeta5((ceil((dates(j)-datesobsmatlab(i))/31))/12)*cupao*kimmelpifinal((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado4(i))*f((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i))*g((ceil((dates(j)-datesobsmatlab(i))/31))/12,FUTUROY);
end
A(i)=Aaux+(1+AffineBeta5((ceil((dates(end)-datesobsmatlab(i))/31))/12)*FUTUROY)*kimmelpifinal((ceil((dates(end)-datesobsmatlab(i))/31))/12,estado4(i))*f((ceil((dates(end)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i))*g((ceil((dates(end)-datesobsmatlab(i))/31))/12,FUTUROY);
H(i) = Haux - AffineBeta5((ceil((dates(end)-datesobsmatlab(i))/31))/12)*kimmelpifinal((ceil((dates(end)-datesobsmatlab(i))/31))/12,estado4(i))*f((ceil((dates(end)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i))*g((ceil((dates(end)-datesobsmatlab(i))/31))/12,FUTUROY);
Haux = 0;
    
     %kimmelpifinal((ceil((dates(end)-datesobsmatlab(i))/31))/12,estado4(i))*
     

Q = theta5*sigma5*sigma5*(1-exp(-k5*dt))^2/(2*k5)+sigma5*sigma5/k5*(exp(-k5*dt)-exp(-2*k5*dt))*FUTUROY;
VARY = VARY-F*VARY*F'+Q;
EZ = A(i)+H(i)*FUTUROY;
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
    
    
    
    
    estado5(i) = EY;
    
    
    
    ll(i) = -(ncol/2)*log(2*pi)-0.5*log(det(VARZ))-0.5*Erro'*inv(VARZ)*Erro;
    
   
end

estado = estado5';
csvwrite('makewhole.csv',estado)

sumll = -sum(ll)
 
end



