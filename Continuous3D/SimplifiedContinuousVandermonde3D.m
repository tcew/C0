function V = SimplifiedContinuousVandermonde3D(N, r, s, t)

  [a,b,c] = rsttoabc(r,s,t);

  cnt = 1;

  Np = (N+1)*(N+2)*(N+3)/6;
  
  V = zeros(length(a),Np);
  
  %% vertices
  V(:,cnt) = (1-a).*(1-b).*(1-c)/8; cnt = cnt+1;
  V(:,cnt) = (1+a).*(1-b).*(1-c)/8; cnt = cnt+1;
  V(:,cnt) = (1+b).*(1-c)/4; cnt = cnt+1; %% NOTE FIX
  V(:,cnt) = (1+c)/2; cnt = cnt+1;

  %% bubble functions
  si = 1;
  for i=0:N
    Ba(:,i+1) = 0.25*(1-a).*(1+a).*SimplifiedJacobiP(a,si,si,i);
    Bb(:,i+1) = 0.25*(1-b).*(1+b).*SimplifiedJacobiP(b,si,si,i);
    Bc(:,i+1) = 0.25*(1-c).*(1+c).*SimplifiedJacobiP(c,si,si,i); 
  end

  %% AB
  for i=0:N-2
    V(:,cnt) = Ba(:,i+1).*((0.5*(1-b)).^(i+2)).*((0.5*(1-c)).^(i+2)); cnt=cnt+1;
  end

  %% AC
  for i=0:N-2
    V(:,cnt) = 0.5*(1-a).*Bb(:,i+1).*((0.5*(1-c)).^(i+2)); cnt=cnt+1;
  end

  %% BC
  for i=0:N-2
    V(:,cnt) = 0.5*(1+a).*Bb(:,i+1).*((0.5*(1-c)).^(i+2)); cnt=cnt+1;
  end
    
  %% AD
  for i=0:N-2
    V(:,cnt) = 0.25*(1-a).*(1-b).*Bc(:,i+1); cnt=cnt+1;
  end

  %% BD
  for i=0:N-2
    V(:,cnt) = 0.25*(1+a).*(1-b).*Bc(:,i+1); cnt=cnt+1;
  end


  %% CD
  for i=0:N-2
    V(:,cnt) = 0.25*(1+b).*Bc(:,i+1); cnt=cnt+1;
  end

  %% ABC
  for i=0:N-3
    for j=0:N-3-i
      V(:,cnt) = Ba(:,i+1).*((0.5*(1-b)).^(i+1)).*Bb(:,j+1).*((0.5*(1-c)).^(i+j+3)); cnt=cnt+1;
    end
  end

  %% ABD
  for i=0:N-3
    for j=0:N-3-i
      V(:,cnt) = Ba(:,i+1).*((0.5*(1-b)).^(i+2)).*Bc(:,j+1).*((0.5*(1-c)).^(i+1)); cnt=cnt+1; %% +2 ?
    end
  end
  

  %% ACD
  for i=0:N-3
    for j=0:N-3-i
      V(:,cnt) = 0.5*(1-a).*Bb(:,i+1).*((0.5*(1-c)).^(i+1)).*Bc(:,j+1); cnt=cnt+1; %%  +2 ?
    end
  end


  %% BCD
  for i=0:N-3
    for j=0:N-3-i
      V(:,cnt) = 0.5*(1+a).*Bb(:,i+1).*((0.5*(1-c)).^(i+1)).*Bc(:,j+1); cnt=cnt+1; %%  +2 ?
    end
  end

  %% interior

  for i=0:N-4
    for j=0:N-4-i
      for k=0:N-4-i-j
	V(:,cnt) = Ba(:,i+1).*((0.5*(1-b)).^(i+1)).*Bb(:,j+1).*((0.5*(1-c)).^(i+j+2)).*Bc(:,k+1); cnt = cnt+1;
      end
    end
  end
      
