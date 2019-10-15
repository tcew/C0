
for N=1:10
  
  %% warp blend nodes
  [rwb,swb,twb] = Nodes3D(N);
  [rwb,swb,twb] = xyztorst(rwb,swb,twb);

  %% 1D Gauss-Legendre nodes to build tensor product with
  [glz,glw] = JacobiGQ(0,0,N+1);

  %% 3D tensor product of GL nodes
  Nq = length(glz);
  cnt = 1;
  JW = [];  agl = []; bgl = []; cgl = [];
  for c=1:Nq
    for b=1:Nq
      for a=1:Nq
	agl(cnt,1) = glz(a);
	bgl(cnt,1) = glz(b);
	cgl(cnt,1) = glz(c);
	JW(cnt,1) = glw(a)*glw(b)*glw(c)*0.125*(1-glz(b))*(1-glz(c))^2;
	cnt = cnt+1;
      end
    end
  end

  %% Duffy transform
  rgl = 0.25*(1+agl).*(1-bgl).*(1-cgl)-1;
  sgl = 0.5*(1+bgl).*(1-cgl)-1;
  tgl = cgl;

  %% build Vandermonde matrix using simplified C0 KS basis
  SVgl = SimplifiedContinuousVandermonde3D(N, rgl, sgl, tgl);
  SVwb = SimplifiedContinuousVandermonde3D(N, rwb, swb, twb);

  %% build Vandermonde matrix using original C0 KS basis
  Vgl = ContinuousVandermonde3D(N, rgl, sgl, tgl);
  Vwb = ContinuousVandermonde3D(N, rwb, swb, twb);

  I = Vgl/Vwb;

  cnt = 1;
  maxerr = [];
  for k=0:N
    for j=0:N-k
      for i=0:N-j-k
	u   = (rwb.^i).*(swb.^j).*(twb.^k);
	ugl = (rgl.^i).*(sgl.^j).*(tgl.^k);
	Iu = I*u;

	maxerr(cnt) = max(max(max(abs(ugl-Iu))));
	cnt = cnt+1;
      end
    end
  end

  maxmaxerr = max(maxerr)

%  condVwb = cond(Vwb)
%  condSimplifiedVwb = cond(SVwb)

  MM = transpose(Vgl)*diag(JW)*Vgl;
  SMM = transpose(SVgl)*diag(JW)*SVgl;

  % build Warp Blend Lagrange basis mass matrix two ways
  Iwb = Vgl/Vwb;
  SIwb = SVgl/SVwb;

  MMwb = transpose(Iwb)*diag(JW)*Iwb;
  SMMwb = transpose(SIwb)*diag(JW)*SIwb;

  maxDiffInMMwb = max(max(abs(MMwb-SMMwb)))

end
