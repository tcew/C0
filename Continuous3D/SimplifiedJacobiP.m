
function phi = SimplifiedJacobiP(x,alpha,beta,N)

  phi= JacobiP(x,alpha,beta,N);
  return;
  
  if(N==0)
    phi = ones(size(x));
    return;
  end

  if(N==1)
    phi = (1+alpha) + (alpha+beta+2)*(x-1)/2;
    return
  end

  if(N==2)
    phi = 0.5*(alpha+1)*(alpha+2) + (alpha+2)*(alpha+beta+3)*(x-1)/2 + 0.5*(alpha+beta+3)*(alpha+beta+4)*((x-1)/2).^2;
    return
  end

  c1 = (2*N+alpha+beta-1)*((2*N+alpha+beta)*(2*N+alpha+beta-2)*x + alpha^2 - beta^2);
  c2 = -2*(N+alpha-1)*(N+beta-1)*(2*N+alpha+beta);
  c3 = 2*N*(N+alpha+beta)*(2*N+alpha+beta-2);
  phi = (c1.*SimplifiedJacobiP(x, alpha, beta, N-1)+c2*SimplifiedJacobiP(x, alpha, beta, N-2))/c3;
