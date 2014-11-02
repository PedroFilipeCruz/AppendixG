function [para, sumll] = CIRMEUtent(a,b,tau)
%a->numero de factores no cir
%b-> ficheiro dados a utilizar

%%EXEMPLO
%%%CIRMEUtent(2,'dadostreasurytese.xlsx',[3/12 6/12 1 2 3 5 7 10])

warning('off')
global u

u = a;

%colocar as yields sem ser em percentagem
 Y = xlsread(b);
 X=Y(:,1:end);
 [nrow, ncol] = size(X);
 Z=[];

Settle=today(); 
CurveDates=Settle+tau*365+1; 
 
%passar da par curve à zero curve
for i=1:nrow
 [zeroyields,ooo]= pyld2zero(X(i,:)',CurveDates',Settle,1,0,2,0);
     for j=1:ncol
      vectorzeros(j) = log(1+zeroyields(j));   
      end
  Z=[Z;vectorzeros];
end 
 
% valor para iniciar minimização
 if a==1
  para0 = [0.0325, 0.4, 0.0586, -0.2537, -0.01, 0.01*rand(1,ncol).*ones(1,ncol)]; 
 elseif a==2
    para0 = [0.3,0.07,0.05,-0.01,0.09,0.04,0.01,-0.008,-0.1, 0.01*rand(1,ncol).*ones(1,ncol)];
     else 
    para0 = [0.1,0.5,0.5,0.18,0.2,0.41,0.25,0.04,0.8,0.005,0.05,-0.001,-0.3, 0.01*rand(1,ncol).*ones(1,ncol)];
 end

    
   

options = struct('MaxFunEvals', 10000);


%minimização
if u==2
[x, fval] = fmincon(@(x) CIRMEULIKE2(x,Z,tau,nrow,ncol), para0,[],[],[],[],[0.00001,0.00001,0.00001, -0.99,0.00001,0.00001,0.00001, -0.99,-1,0.00001*ones(1,ncol)],[1,1,1,0,1,1,1,0,1,1*ones(1,ncol)],@inequa,options);
%[x, fval] = ga(@(x) CIRMEULIKE2(x,Z,tau,nrow,ncol),17,[],[],[],[],[0.00001,0.00001,0.00001, -1,0.00001,0.00001,0.00001, -1,-1,0.00001*ones(1,ncol)],[1,1,1,0,1,1,1,0,1,1*ones(1,ncol)],@inequa,options);

elseif u==1
[x, fval] = fmincon(@(x) CIRMEULIKE1(x,Z,tau,nrow,ncol), para0,[],[],[],[],[0.00001,0.00001,0.00001, -1,-1,0.00001*ones(1,ncol)],[1,1,1, 1,1,1*ones(1,ncol)],@inequa,options);
else
[x, fval] = fmincon(@(x) CIRMEULIKE3(x,Z,tau,nrow,ncol), para0,[],[],[],[],[0.001,0.00001,0.00001, -1,0.0001,0.00001,0.00001, -1,0.00001,0.00001,0.00001, -1,-1,0.00001*ones(1,ncol)],[1,1,1,0,1,1,1,0,1,1,1,0,1,1*ones(1,ncol)],@inequa,options);  
end


para = x;
sumll = fval;

%calcular as yields com o fit

if u==2
%%%%%%%%%%%%%%CASO CIR2%%%%%%%%%%%

M = csvread('novo.csv');
estado1 = M(:,1);
estado2 = M(:,2);

theta = para(1); kappa = para(2); sigma = para(3); lambda = para(4);
theta2 = para(5); kappa2 = para(6); sigma2 = para(7); lambda2 = para(8);
alfa0 = para(9);

for i = 1:ncol
    
    h1=sqrt(((kappa+lambda)^2)+2*sigma*sigma);
    h2=sqrt(((kappa2+lambda2)^2)+2*sigma2*sigma2);
    
    AffineAlpha = ((2*kappa*theta)/(sigma*sigma))*log((2*h1*exp(0.5*(kappa+lambda+h1)*tau(i)))/(h1-(kappa+lambda)+(kappa+lambda+h1)*exp(h1*tau(i))));
    AffineAlpha2 =((2*kappa2*theta2)/(sigma2*sigma2))*log((2*h2*exp(0.5*(kappa2+lambda2+h2)*tau(i)))/(h2-(kappa2+lambda2)+(kappa2+lambda2+h2)*exp(h2*tau(i))));
    
    AffineBeta = (2*(exp(h1*tau(i))-1))/(h1-(kappa+lambda)+(kappa+lambda+h1)*exp(h1*tau(i)));
    AffineBeta2 = (2*(exp(h2*tau(i))-1))/(h2-(kappa2+lambda2)+(kappa2+lambda2+h2)*exp(h2*tau(i)));

    A(i) = alfa0-(AffineAlpha/tau(i))-(AffineAlpha2/tau(i));
    H(i,1) = AffineBeta/tau(i);
    H(i,2) = AffineBeta2/tau(i);
end

dt = 1/12;
deltat = 0:dt:dt*(nrow-1);

fit=[];

for i=1:ncol
   auxiliar = A(i) + H(i,1)*estado1 + H(i,2)*estado2;
   fit=[fit,auxiliar];
    
end

parametros=para(1:9);

printmat(parametros,'parametros','parametros', 'theta1 k1 sigma1 eta1 theta2 k2 sigma2 eta2 alfa0')

elseif u==3
%%%%%%%%%%%%%%CASO CIR3%%%%%%%%%%%

M = csvread('novo.csv');
estado1 = M(:,1);
estado2 = M(:,2);
estado3 = M(:,3);

theta = para(1); kappa = para(2); sigma = para(3); lambda = para(4);
theta2 = para(5); kappa2 = para(6); sigma2 = para(7); lambda2 = para(8);
theta3 = para(9); kappa3 = para(10); sigma3 = para(11); lambda3 = para(12);
alfa0 = para(13);

for i = 1:ncol
    
    h1=sqrt(((kappa+lambda)^2)+2*sigma*sigma);
    h2=sqrt(((kappa2+lambda2)^2)+2*sigma2*sigma2);
    h3=sqrt(((kappa3+lambda3)^2)+2*sigma3*sigma3);
    
    AffineAlpha = ((2*kappa*theta)/(sigma*sigma))*log((2*h1*exp(0.5*(kappa+lambda+h1)*tau(i)))/(h1-(kappa+lambda)+(kappa+lambda+h1)*exp(h1*tau(i))));
    AffineAlpha2 =((2*kappa2*theta2)/(sigma2*sigma2))*log((2*h2*exp(0.5*(kappa2+lambda2+h2)*tau(i)))/(h2-(kappa2+lambda2)+(kappa2+lambda2+h2)*exp(h2*tau(i))));
    AffineAlpha3 =((2*kappa3*theta3)/(sigma3*sigma3))*log((2*h3*exp(0.5*(kappa3+lambda3+h3)*tau(i)))/(h3-(kappa3+lambda3)+(kappa3+lambda3+h3)*exp(h3*tau(i))));
    
    AffineBeta = (2*(exp(h1*tau(i))-1))/(h1-(kappa+lambda)+(kappa+lambda+h1)*exp(h1*tau(i)));
    AffineBeta2 = (2*(exp(h2*tau(i))-1))/(h2-(kappa2+lambda2)+(kappa2+lambda2+h2)*exp(h2*tau(i)));
    AffineBeta3 = (2*(exp(h3*tau(i))-1))/(h3-(kappa3+lambda3)+(kappa3+lambda3+h3)*exp(h3*tau(i)));

    A(i) = alfa0-(AffineAlpha/tau(i))-(AffineAlpha2/tau(i));
    H(i,1) = AffineBeta/tau(i);
    H(i,2) = AffineBeta2/tau(i);
    H(i,3) = AffineBeta3/tau(i);
end

dt = 1/12;
deltat = 0:dt:dt*(nrow-1);

fit=[];

for i=1:ncol
   auxiliar = A(i) + H(i,1)*estado1 + H(i,2)*estado2 + H(i,3)*estado3;
   fit=[fit,auxiliar];
    
end
%Tabela
parametros=para(1:13);

printmat(parametros,'parametros','parametros', 'theta1 k1 sigam1 eta1 theta2 k2 sigam2 eta2 theta3 k3 sigam3 eta3 alfa0')




else
    
%%%%%%%%%%%%%%%%%%%%%%CASO CIR1

M = csvread('novo.csv');
estado1 = M(:);

theta = para(1); kappa = para(2); sigma = para(3); lambda = para(4);
alfa0 = para(5);

for i = 1:ncol
    
    h1=sqrt(((kappa+lambda)^2)+2*sigma*sigma);
    AffineAlpha = ((2*kappa*theta)/(sigma*sigma))*log((2*h1*exp(0.5*(kappa+lambda+h1)*tau(i)))/(h1-(kappa+lambda)+(kappa+lambda+h1)*exp(h1*tau(i))));
    AffineBeta = (2*(exp(h1*tau(i))-1))/(h1-(kappa+lambda)+(kappa+lambda+h1)*exp(h1*tau(i)));
   
    A(i) = alfa0-(AffineAlpha/tau(i));
    H(i) = AffineBeta/tau(i);
    
end

dt = 1/12;
deltat = 0:dt:dt*(nrow-1);


fit=[];

for i=1:ncol
   auxiliar = A(i) + H(i)*estado1;
   fit=[fit,auxiliar];
    
end

parametros=para(1:5);

printmat(parametros,'parametros','parametros', 'theta1 k1 sigam1 eta1 alfa0')

end


%calcular o RMSE
diferencas=[];
for i=1:ncol
    AuxiliarDiferencas = ((Z(:,i))-fit(:,i)).^2;
    diferencas=[diferencas,AuxiliarDiferencas];
end

Error=[];


RMSEtot=10000*sqrt(sum(sum(diferencas))/(nrow*ncol));
     
Error=[Error;RMSEtot];

for j=1:ncol
RMSE = 10000*sqrt(sum(diferencas(:,j))/(nrow));  
Error=[Error;RMSE];
end

Error

Colors=jet(ncol);

%plots dos fits
for j=1:ncol
figure (j)
plot(deltat,Z(:,j),'k')
hold on
plot(deltat,fit(:,j),'Color',Colors(j,:))  
end




end

%impor as condições de feller na minimização
function [c,ceq]=inequa(x)

global u

if u==1 
 c =  x(3)^2 - 2*x(1)*x(2);
 ceq = [];   
elseif u==2
  c(1) = x(3)^2 - 2*x(1)*x(2);
 %c(2) = x(7)^2 - 2*x(5)*x(6);  
   ceq = [];
else
 c(1) = x(3)^2 - 2*x(1)*x(2);
 c(2) = x(7)^2 - 2*x(5)*x(6);  
 c(3) = x(11)^2 - 2*x(9)*x(10);
   ceq = [];
end

end

