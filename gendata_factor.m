function dataF = gendata_factor(ftype,samplesize)
switch ftype
    case 'Jiang6'
        matSigma = 0.5.^abs(repmat(0:(2-1),2,1)-repmat(0:(2-1),2,1)');
        dataF12 = (mvnrnd(zeros(2,1),matSigma,samplesize));
        errorxi1 = (normrnd(0,1,samplesize,1));
        dataF3 = abs(sum(dataF12,2))+dataF12(:,1).*errorxi1;
        errorxi2 = (normrnd(0,1,samplesize,1));
        dataF4 = sum(dataF12,2).^2+abs(dataF12(:,2)).*errorxi2;
        dataF5 = (binornd(1,logsig(dataF12(:,2)),samplesize,1));
        dataF6 = (binornd(1,normcdf(dataF12(:,2)),samplesize,1));
        dataF = [dataF12, dataF3, dataF4, dataF5, dataF6];
    case 'Jiang9'
        matSigma = 0.5.^abs(repmat(0:(2-1),2,1)-repmat(0:(2-1),2,1)');
        dataF12 = (mvnrnd(zeros(2,1),matSigma,samplesize));
        errorxi1 = (normrnd(0,1,samplesize,1));
        dataF3 = abs(sum(dataF12,2))+dataF12(:,1).*errorxi1;
        errorxi2 = (normrnd(0,1,samplesize,1));
        dataF4 = sum(dataF12,2).^2+abs(dataF12(:,2)).*errorxi2;
        dataF5 = (binornd(1,logsig(dataF12(:,2)),samplesize,1));
        dataF6 = (binornd(1,normcdf(dataF12(:,2)),samplesize,1));
        dataF7 = (binornd(1,normcdf(dataF3),samplesize,1));
        dataF8 = (normrnd(0,1,samplesize,1));
        dataF9 = (normrnd(0,1,samplesize,1));
        dataF = [dataF12, dataF3, dataF4, dataF5, dataF6, dataF7, dataF8, dataF9];
    case 'tang2019'
        dataF = exprnd(1,samplesize,9);
    case 'Jiang9New'
        matSigma = 0.5.^abs(repmat(0:(2-1),2,1)-repmat(0:(2-1),2,1)');
        dataF12 = (mvnrnd(zeros(2,1),matSigma,samplesize));
        errorxi1 = (normrnd(0,1,samplesize,1));
        dataF3 = abs(sum(dataF12,2))+dataF12(:,1).*errorxi1;
        errorxi2 = (normrnd(0,1,samplesize,1));
        dataF4 = sum(dataF12,2).^2+abs(dataF12(:,2)).*errorxi2;
        dataF5 = (binornd(1,logsig(dataF12(:,2)),samplesize,1));
        dataF6 = (binornd(1,normcdf(dataF12(:,2)),samplesize,1));
        dataF7 = (binornd(1,normcdf(dataF3),samplesize,1));
        dataF8 = (unifrnd(0,1,samplesize,1));
        dataF9 = (unifrnd(0,1,samplesize,1));
        dataF = [dataF12, dataF3, dataF4, dataF5, dataF6, dataF7, dataF8, dataF9];
    case 'Jiang2New'
        matSigma = 0.5.^abs(repmat(0:(2-1),2,1)-repmat(0:(2-1),2,1)');
        dataF12 = (mvnrnd(zeros(2,1),matSigma,samplesize));
        dataF = [dataF12];
    case 'mimic'
        % Disease.stage, 1--4 ordinal, [0.1,0.1,0.5,0.3]
        temp = rand(samplesize,1);
        dataF1 = 1*(temp<=0.1)+2*(temp>0.1 & temp<=0.2)+...
            3*(temp>0.2 & temp<=0.7)+4*(temp>0.7);
        % Double.or.Triple.Hit.Lymphoma, binary
        dataF2 = binornd(1,0.3,samplesize,1);
        % Age, continuous
        dataF3 = gamrnd(3,14,samplesize,1);
        % Gender, binary
        dataF4 = binornd(1,0.45,samplesize,1);
        % Prior.Lines.of.Therapy, 0-3 ordinal, [0.2,0.3,0.3,0.2]
        temp = rand(samplesize,1);
        dataF5 = 0*(temp<=0.2)+1*(temp>0.2 & temp<=0.5)+...
            2*(temp>0.5 & temp<=0.8)+3*(temp>0.8);
        % Prior.Autologous.Stem.Cell.Transplant..ASCT. binary
        dataF6 = binornd(1,0.2,samplesize,1);
        % Started.other.therapy.post.CAR.T.or.post.stem.cell.transplant, binary
        dataF7 = binornd(1,0.6,samplesize,1);
        % Ethnic.Group, binary. 
        dataF8 = binornd(1,0.8,samplesize,1);
        dataF = [dataF1, dataF2, dataF3, dataF4, dataF5, dataF6, dataF7];
    case 'mimic3'
        % Disease.stage, 1--4 ordinal, [0.1,0.1,0.5,0.3]
        temp = rand(samplesize,1);
        dataF1 = 1*(temp<=0.1)+2*(temp>0.1 & temp<=0.2)+...
            3*(temp>0.2 & temp<=0.7)+4*(temp>0.7);
        % Double.or.Triple.Hit.Lymphoma, binary
        dataF2 = gamrnd(3,14,samplesize,1);
        dataF = [dataF1, dataF2];
end
end
