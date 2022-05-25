function energy = TCEuler_energy(u,u0,x1,x2,afa,lambda)

%%%%%%%%%%%%%%%%%%%%%  Energy
norm_x = sqrt(x1.^2 + x2.^2);
normx = afa.*norm_x;

energy = 0.5*lambda*norm(u-u0,'fro')^2 + sum(normx(:));
