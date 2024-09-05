function parasout = gendata_beta(xdimp,xdimd,bttype)
% ma71: [0,-1,0,1,0,-1]
% ma72: [1.3,-1.3,1,-0.5,0.5,-0.5]
% ma93: [-0.6,0,-0.3,-0.1,0,0.1,0.3,0,0.6]
% pm1a: 1
% pm1b: 0
% add1: [0.5,0,1,0.5,0,0.5]
% gex1: [-3.75,2,-1,0.75,0,-0.75,1,-2,3.75]/5
% ge52: [2.75,-0.75,-1,2
%        -3.125,-1.125,1,-2]
% posi: [0,1,0,1,0,1]

switch bttype
    case 'self'
    case 'ma71'
        parasout = [0,-1,0,1,0,-1]';
    case 'ma72'
        parasout = [1.3,-1.3,1,-0.5,0.5,-0.5]';
    case 'ma61'
        parasout = [1.3,-1.3,-1,-0.5,0.5]'/2;
    case 'ma62'
        parasout = [-1.3,-1,-0.5,0.5;-1,0.8,-0.8,0.5]'/2;
    case 'ma91'
        parasout = [-0.6,0,-0.3,-0.1,0,0.1,0.3,-0.5]';
    case 'ma21'
        parasout = [-0.6]';
    case 'ma92'
        parasout = [-0.6,-0.3,-0.1,0,0.1,0.3,-0.5;...
            1,-0.5,0.5,-1,0.8,-0.8,0.5]';
    case 'mimic1'
        parasout = [-0.6,-0.3,-0.1,0.2,0.5,-0.2]';
    case 'mimic2'
        parasout = [-0.3,-0.1,0.1,0.3,-0.5,0;...
            -0.5,0.4,-1,0.6,-0.8,-0.2]';
    case 'mimic3'
        parasout = [-0.6]';
    case 'pm1a'
        parasout = 1;
    case 'pm1b'
        parasout = 0;
    case 'add1'
        parasout = [0.5,0,1,0.5,0,0.5]';
    case 'gex1'
        parasout = [-3.75,2,-1,0.75,0,-0.75,1,-2,3.75]';
        parasout = parasout/5;
    case 'ge52'
        parasout = [2.75,-0.75,-1,2;
                    -3.125,-1.125,1,-2]';
    case 'posi'
        parasout = [0,1,0,1,0,1]';
    case 'tang2019'
        parasout = ones(xdimp-1,xdimd);
end
parasout = (parasout);
end
