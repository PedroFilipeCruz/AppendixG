function yield = ytmpedro(cupao,anos,preco,facevalue)

options = optimset('TolX',1e-12);



equation = @(x) abs(preco-(0.5*(cupao*facevalue)*((1-((1+0.5*x)^(-2*anos)))/(0.5*x)) + (facevalue/((1+0.5*x)^(2*anos)))));

 [a,b]=fminbnd(equation,0,1,options);
 
 if anos > 0 
yield=a;
 else
     yield=0;
 end

end