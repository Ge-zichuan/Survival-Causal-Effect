function xout = gendata_x(n,xdimp,xtype)
% generate data X, n: sample size, xdim: dimension of each observation
% xtype---
%       norm: standard normal ranform variable
%       unif: standard uniform ranform variable
switch xtype
    case 'unif'
        xout = rand(n,xdimp);
    case 'norm'
        xout = normrnd(0,1,n,xdimp);
    case 'pm1a'
        xout(:,1) = binornd(1,0.5,n);
        xout(:,2) = rand(n,1);
    case 'pm2a'
        xout = binornd(1,0.5,n,ximp);
end
end