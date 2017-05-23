% funmatdemo    
% by X. Armangue
% (c) Mr3D - University of Girona, September 2002
%
% Demostration of the Fundamental Matrix Estimation Toolbox
clear all; close all;

% List of methods
p7			=	1;
ls			=	2;
eig		    =	3;
Faugeras	=	4;
ilmmdpel	=	5;
nr			=	6;
inmmdpel	=	7;
gradls	    =	8;
gradeig	    =	9;
Mestls	    =	10;
Mesteig	    =	11;
MestTorr	=	12;
LMedSls	    =	13;
LMedSeig	=	14;
RANSAC	    =	15;
mtorr_estf  =   16;                     % copyright Philip Torr and Microsoft Corp 2002
mtorr_estf_bookstein        =   17;     % copyright Philip Torr and Microsoft Corp 2002
mtorr_estf_bookstein_sam    =   18;     % copyright Philip Torr and Microsoft Corp 2002
mtorr_estf_lin_non_lin      =   19;     % copyright Philip Torr and Microsoft Corp 2002
mtorr_F_constrained_fit     =   20;     % copyright Philip Torr and Microsoft Corp 2002
mlesac      =   21;                     % copyright Philip Torr and Microsoft Corp 2002
mapsac      =   22;                     % copyright Philip Torr and Microsoft Corp 2002
fns         =   23;                     % copyright Anton van den Hengel form University of Adelaide
cfns        =   24;                     % copyright Anton van den Hengel form University of Adelaide

image=5;            % choose a number of image 
redraw=0;           % choose draw epipolar geometry
normal=2;           % choose normalization method 

% choose one or more methods
%method=[1,20,2,3,16,4];             % linear methods
%method=[5,17,18,6,19,7,8,9,23,24];  % iterative methods
%method=[10,11,12,13,14,15,21,22];   % robust methods
method=[1,20,2,3,16,4,5,17,18,6,19,7,8,9,23,24,10,11,12,13,14,15,21,22]; % all methods
%method=[13,14,15,21,22]; 

threshold=0.75;

[M,nomia,nomib]=loadtestimages(image); % loading image points 
plotimagepoints(M,nomia,nomib)         % display image points
disp('Display image points. (Press Return)')
pause
disp('*******************************************************************')

% Normalization
Mini=M;
switch normal,
    case 1,
        disp('Normalization [-1 1]');
        [M,T1,T2]=normalonetoone(Mini);
    case 2,
        disp('Hartley Normalization');
        [M,T1,T2]=normalHartley(Mini);
    case 3,
        disp('No Normalization');
        T1=[1 0 0; 0 1 0; 0 0 1];
        T2=[1 0 0; 0 1 0; 0 0 1];
end

bars=[];
bartime=[];
xlabelbars=[];
indexbars=1;
for i=1:size(method,2),
    redrawnow=0;
    switch method(i),
        case 1,
            disp('Method: seven points')
            tic
            F7p=funmat7p(M(:,1:7));
            time7p=toc;
            F7p=T1'*F7p*T2;
            F7p=F7p./norm(F7p);
            [mean7p,stddev7p,minim7p,maxim7p]=funmatError(Mini,F7p);
            disp('Funamental Matrix')
            disp(F7p)
            disp(sprintf('Distance point-epipolar line\nmean: %f  stdev: %f  min: %f  max: %f\n',mean7p,stddev7p,minim7p,maxim7p));
            disp(sprintf('Rank-2: %d',rank(F7p)==2))            
            disp(sprintf('Time: %f',time7p))            
            disp('*******************************************************************')
            if redraw, 
                fummatplot(Mini,F7p);
            end
            bars=[bars;mean7p,stddev7p];
            bartime=[bartime;time7p];
            xlabelbars=[xlabelbars sprintf(' %d.- 7p ',indexbars)];
            if mod(indexbars,8)==0, xlabelbars=[xlabelbars '\n']; end
            indexbars=indexbars+1;
          
        case 2,
            disp('Method: Least-squares')
            tic
            Fls=funmatls(M);
            timels=toc;
            Fls=T1'*Fls*T2;
            Fls=Fls./norm(Fls);
            [meanls,stddevls,minimls,maximls]=funmatError(Mini,Fls);
            disp('Funamental Matrix')
            disp(Fls)
            disp(sprintf('Distance point-epipolar line\nmean: %f  stdev: %f  min: %f  max: %f\n',meanls,stddevls,minimls,maximls));
            disp(sprintf('Rank-2: %d',rank(Fls)==2))            
            disp(sprintf('Time: %f',timels))            
            disp('*******************************************************************')
            if redraw, 
                fummatplot(Mini,Fls)
            end
            bars=[bars;meanls,stddevls];
            bartime=[bartime;timels];
            xlabelbars=[xlabelbars sprintf(' %d.- ls ',indexbars)];
            if mod(indexbars,8)==0, xlabelbars=[xlabelbars '\n']; end
            indexbars=indexbars+1;
            
        case 3,
            disp('Method: Least-squares with eigen analysis')
            tic
            Feig=funmateig(M);
            timeeig=toc;
            Feig=T1'*Feig*T2;
            Feig=Feig./norm(Feig);
            [meaneig,stddeveig,minimeig,maximeig]=funmatError(Mini,Feig);
            disp('Funamental Matrix')
            disp(Feig)
            disp(sprintf('Distance point-epipolar line\nmean: %f  stdev: %f  min: %f  max: %f\n',meaneig,stddeveig,minimeig,maximeig));
            disp(sprintf('Rank-2: %d',rank(Feig)==2))            
            disp(sprintf('Time: %f',timeeig))            
            disp('*******************************************************************')
            if redraw, 
                fummatplot(Mini,Feig)
            end
            bars=[bars;meaneig,stddeveig];
            bartime=[bartime;timeeig];
            xlabelbars=[xlabelbars sprintf(' %d.- lseig ',indexbars)];
            if mod(indexbars,8)==0, xlabelbars=[xlabelbars '\n']; end
            indexbars=indexbars+1;
            
        case 4,
            disp('Method: linear rank-2 constraint (by Faugeras)')
            tic
            FFaugeras=funmatFaugeras(M);
            timeFaugeras=toc;
            FFaugeras=T1'*FFaugeras*T2;
            FFaugeras=FFaugeras./norm(FFaugeras);
            [meanFaugeras,stddevFaugeras,minimFaugeras,maximFaugeras]=funmatError(Mini,FFaugeras);
            disp('Funamental Matrix')
            disp(FFaugeras)
            disp(sprintf('Distance point-epipolar line\nmean: %f  stdev: %f  min: %f  max: %f\n',meanFaugeras,stddevFaugeras,minimFaugeras,maximFaugeras));
            disp(sprintf('Rank-2: %d',rank(FFaugeras)==2))            
            disp(sprintf('Time: %f',timeFaugeras))            
            disp('*******************************************************************')
            if redraw, 
                fummatplot(Mini,FFaugeras)
            end
            bars=[bars;meanFaugeras,stddevFaugeras];
            bartime=[bartime;timeFaugeras];
            xlabelbars=[xlabelbars sprintf(' %d.- Faugeras ',indexbars)];
            if mod(indexbars,8)==0, xlabelbars=[xlabelbars '\n']; end
            indexbars=indexbars+1;
            
        case 5,
            disp('Method: Iterative Linear Method Minimizing Distances of Points to Epipolar Line')
            tic
            Filmmdpel=funmatilmmdpel(M,funmatls(M),100,10^(-5));
            timeilmmdpel=toc;
            Filmmdpel=T1'*Filmmdpel*T2;
            Filmmdpel=Filmmdpel./norm(Filmmdpel);
            [meanilmmdpel,stddevilmmdpel,minimilmmdpel,maximilmmdpel]=funmatError(Mini,Filmmdpel);
            disp('Funamental Matrix')
            disp(Filmmdpel)
            disp(sprintf('Distance point-epipolar line\nmean: %f  stdev: %f  min: %f  max: %f\n',meanilmmdpel,stddevilmmdpel,minimilmmdpel,maximilmmdpel));
            disp(sprintf('Rank-2: %d',rank(Filmmdpel)==2))            
            disp(sprintf('Time: %f',timeilmmdpel))            
            disp('*******************************************************************')
            if redraw, 
                fummatplot(Mini,Filmmdpel)
            end
            bartime=[bartime;timeilmmdpel];
            bars=[bars;meanilmmdpel,stddevilmmdpel];
            xlabelbars=[xlabelbars sprintf(' %d.- ilmmdpel ',indexbars)];
            if mod(indexbars,8)==0, xlabelbars=[xlabelbars '\n']; end
            indexbars=indexbars+1;
            
        case 6,
            disp('Method: iterative Newton-Raphson')
            tic
            Fnr=funmatnr(M,funmatls(M),100,10^(-5));
            timenr=toc;
            Fnr=T1'*Fnr*T2;
            Fnr=Fnr./norm(Fnr);
            [meannr,stddevnr,minimnr,maximnr]=funmatError(Mini,Fnr);
            disp('Funamental Matrix')
            disp(Fnr)
            disp(sprintf('Distance point-epipolar line\nmean: %f  stdev: %f  min: %f  max: %f\n',meannr,stddevnr,minimnr,maximnr));
            disp(sprintf('Rank-2: %d',rank(Fnr)==2))            
            disp(sprintf('Time: %f',timenr))            
            disp('*******************************************************************')
            if redraw, 
                fummatplot(Mini,Fnr)
            end
            bartime=[bartime;timenr];
            bars=[bars;meannr,stddevnr];
            xlabelbars=[xlabelbars sprintf(' %d.- nr ',indexbars)];
            if mod(indexbars,8)==0, xlabelbars=[xlabelbars '\n']; end
            indexbars=indexbars+1;
            
        case 7,
            disp('Method: Iterative Nonlinear Method Minimizing Distances of Points to Epipolar Lines')
            tic
            Finmmdpel=funmatinmmdpel(M,funmatls(M),100);
            timeinmmdpel=toc;
            Finmmdpel=T1'*Finmmdpel*T2;
            Finmmdpel=Finmmdpel./norm(Finmmdpel);
            [meaninmmdpel,stddevinmmdpel,miniminmmdpel,maximinmmdpel]=funmatError(Mini,Finmmdpel);
            disp('Funamental Matrix')
            disp(Finmmdpel)
            disp(sprintf('Distance point-epipolar line\nmean: %f  stdev: %f  min: %f  max: %f\n',meaninmmdpel,stddevinmmdpel,miniminmmdpel,maximinmmdpel));
            disp(sprintf('Rank-2: %d',rank(Finmmdpel)==2))            
            disp(sprintf('Time: %f',timeinmmdpel))            
            disp('*******************************************************************')
            if redraw, 
                fummatplot(Mini,Finmmdpel)
            end
            bars=[bars;meaninmmdpel,stddevinmmdpel];
            bartime=[bartime;timeinmmdpel];
            xlabelbars=[xlabelbars sprintf(' %d.- inmmdpel ',indexbars)];
            if mod(indexbars,8)==0, xlabelbars=[xlabelbars '\n']; end
            indexbars=indexbars+1;
           
        case 8,
            disp('Method: gradient using Least-Squares')
            tic
            Fgradls=funmatgradls(M,10,10^(-5));
            timegradls=toc;
            Fgradls=T1'*Fgradls*T2;
            Fgradls=Fgradls./norm(Fgradls);
            [meangradls,stddevgradls,minimgradls,maximgradls]=funmatError(Mini,Fgradls);
            disp('Funamental Matrix')
            disp(Fgradls)
            disp(sprintf('Distance point-epipolar line\nmean: %f  stdev: %f  min: %f  max: %f\n',meangradls,stddevgradls,minimgradls,maximgradls));
            disp(sprintf('Rank-2: %d',rank(Fgradls)==2))            
            disp(sprintf('Time: %f',timegradls))            
            disp('*******************************************************************')
            if redraw, 
                fummatplot(Mini,Fgradls)
            end
            bars=[bars;meangradls,stddevgradls];
            bartime=[bartime;timegradls];
            xlabelbars=[xlabelbars sprintf(' %d.- gradls ',indexbars)];
            if mod(indexbars,8)==0, xlabelbars=[xlabelbars '\n']; end
            indexbars=indexbars+1;
            
        case 9,
            disp('Method: gradient using Least-Squares with eigen analysis')
            tic
            Ffunmatgradeig=funmatgradeig(M,10,10^(-5));
            timegradeig=toc;
            Ffunmatgradeig=T1'*Ffunmatgradeig*T2;
            Ffunmatgradeig=Ffunmatgradeig./norm(Ffunmatgradeig);
            [meangradeig,stddevgradeig,minimgradeig,maximgradeig]=funmatError(Mini,Ffunmatgradeig);
            disp('Funamental Matrix')
            disp(Ffunmatgradeig)
            disp(sprintf('Distance point-epipolar line\nmean: %f  stdev: %f  min: %f  max: %f\n',meangradeig,stddevgradeig,minimgradeig,maximgradeig));
            disp(sprintf('Rank-2: %d',rank(Ffunmatgradeig)==2))            
            disp(sprintf('Time: %f',timegradeig))            
            disp('*******************************************************************')
            if redraw, 
                fummatplot(Mini,Ffunmatgradeig)
            end
            bars=[bars;meangradeig,stddevgradeig];
            bartime=[bartime;timegradeig];
            xlabelbars=[xlabelbars sprintf(' %d.- gradeig ',indexbars)];
            if mod(indexbars,8)==0, xlabelbars=[xlabelbars '\n']; end
            indexbars=indexbars+1;
          
        case 10,
            disp('Method: M-Estimator using Least-Squares')
            tic
            [FMestls,wMestls]=funmatMestls(M,10,10^(-5));
            timeMestls=toc;
            FMestls=T1'*FMestls*T2;
            FMestls=FMestls./norm(FMestls);
            [meanMestls,stddevMestls,minimMestls,maximMestls]=funmatError(goodpoints(Mini,wMestls,threshold),FMestls);
            disp('Funamental Matrix')
            disp(FMestls)
            disp(sprintf('Distance point-epipolar line\nmean: %f  stdev: %f  min: %f  max: %f\n',meanMestls,stddevMestls,minimMestls,maximMestls));
            disp(sprintf('Rank-2: %d',rank(FMestls)==2))
            w=wMestls;
            redrawnow=1;
            Mfi=goodpoints(Mini,w,threshold);
            disp(sprintf('Outliers: %d/%d',size(M,2)-size(Mfi,2),size(M,2)))
            disp(sprintf('Time: %f',timeMestls))            
            disp('*******************************************************************')
            if redraw, 
                fummatplot(Mfi,FMestls)
            end
            bars=[bars;meanMestls,stddevMestls];
            bartime=[bartime;timeMestls];
            xlabelbars=[xlabelbars sprintf(' %d.- Mestls ',indexbars)];
            if mod(indexbars,8)==0, xlabelbars=[xlabelbars '\n']; end
            indexbars=indexbars+1;
            
        case 11,
            disp('Method: M-Estimator using Least-Squares with eigen analysis')
            tic
            [FMesteig,wMesteig]=funmatMesteig(M,10,10^(-5));
            timeMesteig=toc;
            FMesteig=T1'*FMesteig*T2;
            FMesteig=FMesteig./norm(FMesteig);
            [meanMesteig,stddevMesteig,minimMesteig,maximMesteig]=funmatError(goodpoints(Mini,wMesteig,threshold),FMesteig);
            disp('Funamental Matrix')
            disp(FMesteig)
            disp(sprintf('Distance point-epipolar line\nmean: %f  stdev: %f  min: %f  max: %f\n',meanMesteig,stddevMesteig,minimMesteig,maximMesteig));
            disp(sprintf('Rank-2: %d',rank(FMesteig)==2))            
            w=wMesteig;
            redrawnow=1;
            Mfi=goodpoints(Mini,w,threshold);
            disp(sprintf('Outliers: %d/%d',size(M,2)-size(Mfi,2),size(M,2)))
            disp(sprintf('Time: %f',timeMesteig))            
            disp('*******************************************************************')
            if redraw, 
                fummatplot(Mfi,FMesteig)
            end
            bars=[bars;meanMesteig,stddevMesteig];
            bartime=[bartime;timeMesteig];
            xlabelbars=[xlabelbars sprintf(' %d.- Mesteig ',indexbars)];
            if mod(indexbars,8)==0, xlabelbars=[xlabelbars '\n']; end
            indexbars=indexbars+1;
            
        case 12,
            disp('Method: M-Estimator proposed by Torr (IJCV97)')
            tic
            [FMestTorr,wMestTorr]=funmatMestTorr(M,10,10^(-5));
            timeMestTorr=toc;
            FMestTorr=T1'*FMestTorr*T2;
            FMestTorr=FMestTorr./norm(FMestTorr);
            [meanMestTorr,stddevMestTorr,minimMestTorr,maximMestTorr]=funmatError(goodpoints(Mini,wMestTorr,threshold),FMestTorr);
            disp('Funamental Matrix')
            disp(FMestTorr)
            disp(sprintf('Distance point-epipolar line\nmean: %f  stdev: %f  min: %f  max: %f\n',meanMestTorr,stddevMestTorr,minimMestTorr,maximMestTorr));
            disp(sprintf('Rank-2: %d',rank(FMestTorr)==2))            
            w=wMestTorr;
            redrawnow=1;
            Mfi=goodpoints(Mini,w,threshold);
            disp(sprintf('Outliers: %d/%d',size(M,2)-size(Mfi,2),size(M,2)))
            disp(sprintf('Time: %f',timeMestTorr))            
            disp('*******************************************************************')
            if redraw, 
                fummatplot(Mfi,FMestTorr)
            end
            bars=[bars;meanMestTorr,stddevMestTorr];
            bartime=[bartime;timeMestTorr];
            xlabelbars=[xlabelbars sprintf(' %d.- MestTorr ',indexbars)];
            if mod(indexbars,8)==0, xlabelbars=[xlabelbars '\n']; end
            indexbars=indexbars+1;
           
        case 13,   
            disp('Method: LMedS using least squares')
            tic
            [FLMedSls,wLMedSls]=funmatLMedSls(M,8,0.99,0.25);
            timeLMedSls=toc;
            FLMedSls=T1'*FLMedSls*T2;
            FLMedSls=FLMedSls./norm(FLMedSls);
            [meanLMedSls,stddevLMedSls,minimLMedSls,maximLMedSls]=funmatError(goodpoints(Mini,wLMedSls,threshold),FLMedSls);
            disp('Funamental Matrix')
            disp(FLMedSls)
            disp(sprintf('Distance point-epipolar line\nmean: %f  stdev: %f  min: %f  max: %f\n',meanLMedSls,stddevLMedSls,minimLMedSls,maximLMedSls));
            disp(sprintf('Rank-2: %d',rank(FLMedSls)==2))            
            w=wLMedSls;
            redrawnow=1;
            Mfi=goodpoints(Mini,w,threshold);
            disp(sprintf('Outliers: %d/%d',size(M,2)-size(Mfi,2),size(M,2)))
            disp(sprintf('Time: %f',timeLMedSls))            
            disp('*******************************************************************')
            if redraw, 
                fummatplot(Mfi,FLMedSls)
            end
            bars=[bars;meanLMedSls,stddevLMedSls];
            bartime=[bartime;timeLMedSls];
            xlabelbars=[xlabelbars sprintf(' %d.- LMedSls ',indexbars)];
            if mod(indexbars,8)==0, xlabelbars=[xlabelbars '\n']; end
            indexbars=indexbars+1;
            
        case 14,
            disp('Method: LMedS using least squares with eigen analysis')
            tic
            [FLMedSeig,wLMedSeig]=funmatLMedSeig(M,8,0.99,0.25);
            timeLMedSeig=toc;
            FLMedSeig=T1'*FLMedSeig*T2;
            FLMedSeig=FLMedSeig./norm(FLMedSeig);
            [meanLMedSeig,stddevLMedSeig,minimLMedSeig,maximLMedSeig]=funmatError(goodpoints(Mini,wLMedSeig,threshold),FLMedSeig);
            disp('Funamental Matrix')
            disp(FLMedSeig)
            disp(sprintf('Distance point-epipolar line\nmean: %f  stdev: %f  min: %f  max: %f\n',meanLMedSeig,stddevLMedSeig,minimLMedSeig,maximLMedSeig));
            disp(sprintf('Rank-2: %d',rank(FLMedSeig)==2))            
            w=wLMedSeig;
            redrawnow=1;
            Mfi=goodpoints(Mini,w,threshold);
            disp(sprintf('Outliers: %d/%d',size(M,2)-size(Mfi,2),size(M,2)))
            disp(sprintf('Time: %f',timeLMedSeig))            
            disp('*******************************************************************')
            if redraw, 
                fummatplot(Mfi,FLMedSeig)
            end
            bars=[bars;meanLMedSeig,stddevLMedSeig];
            bartime=[bartime;timeLMedSeig];
            xlabelbars=[xlabelbars sprintf(' %d.- LMedSeig ',indexbars)];
            if mod(indexbars,8)==0, xlabelbars=[xlabelbars '\n']; end
            indexbars=indexbars+1;
            
        case 15,
            disp('Method: RANSAC')
            tic
            [FRANSAC,wRANSAC]=funmatRANSAC(M,8,0.99,0.25);
            timeRANSAC=toc;
            FRANSAC=T1'*FRANSAC*T2;
            FRANSAC=FRANSAC./norm(FRANSAC);
            [meanRANSAC,stddevRANSAC,minimRANSAC,maximRANSAC]=funmatError(goodpoints(Mini,wRANSAC,threshold),FRANSAC);
            disp('Funamental Matrix')
            disp(FRANSAC)
            disp(sprintf('Distance point-epipolar line\nmean: %f  stdev: %f  min: %f  max: %f\n',meanRANSAC,stddevRANSAC,minimRANSAC,maximRANSAC));
            disp(sprintf('Rank-2: %d',rank(FRANSAC)==2))            
            w=wRANSAC;
            redrawnow=1;
            Mfi=goodpoints(Mini,w,threshold);
            disp(sprintf('Outliers: %d/%d',size(M,2)-size(Mfi,2),size(M,2)))
            disp(sprintf('Time: %f',timeRANSAC))            
            disp('*******************************************************************')
            if redraw, 
                fummatplot(Mfi,FRANSAC)
            end
            bars=[bars;meanRANSAC,stddevRANSAC];
            bartime=[bartime;timeRANSAC];
            xlabelbars=[xlabelbars sprintf(' %d.- RANSAC ',indexbars)];
            if mod(indexbars,8)==0, xlabelbars=[xlabelbars '\n']; end
            indexbars=indexbars+1;
            
        case 16,
            disp('Method: linear implemented by Torr (equivalent to least-squares with eigen analysis)')
            tic
            [f, f_sq_errors, n_inliers,inlier_index,Ftorr_estf] = torr_estimateF( M', 1, [], 'linear', 0);
            timetorr_estf=toc;
            Ftorr_estf=T1'*Ftorr_estf'*T2;
            Ftorr_estf=Ftorr_estf./norm(Ftorr_estf);
            [meantorr_estf,stddevtorr_estf,minimtorr_estf,maximtorr_estf]=funmatError(Mini,Ftorr_estf);
            disp('Funamental Matrix')
            disp(Ftorr_estf)
            disp(sprintf('Distance point-epipolar line\nmean: %f  stdev: %f  min: %f  max: %f\n',meantorr_estf,stddevtorr_estf,minimtorr_estf,maximtorr_estf));
            disp(sprintf('Rank-2: %d',rank(Ftorr_estf)==2))            
            disp(sprintf('Time: %f',timetorr_estf))            
            disp('*******************************************************************')
            if redraw, 
                fummatplot(Mini,Ftorr_estf)   
            end
            bars=[bars;meantorr_estf,stddevtorr_estf];
            bartime=[bartime;timetorr_estf];
            xlabelbars=[xlabelbars sprintf(' %d.- torr estf ',indexbars)];
            if mod(indexbars,8)==0, xlabelbars=[xlabelbars '\n']; end
            indexbars=indexbars+1;
            
        case 17,
            disp('Method: bookstein implemented by Torr')
            tic
            [f, f_sq_errors, n_inliers,inlier_index,Ftorr_estf_bookstein] = torr_estimateF( M', 1, [], 'bookstein', 0);
            timetorr_estf_bookstein=toc;
            Ftorr_estf_bookstein=T1'*Ftorr_estf_bookstein'*T2;
            Ftorr_estf_bookstein=Ftorr_estf_bookstein./norm(Ftorr_estf_bookstein);
            [meantorr_estf_bookstein,stddevtorr_estf_bookstein,minimtorr_estf_bookstein,maximtorr_estf_bookstein]=funmatError(Mini,Ftorr_estf_bookstein);
            disp('Funamental Matrix')
            disp(Ftorr_estf_bookstein)
            disp(sprintf('Distance point-epipolar line\nmean: %f  stdev: %f  min: %f  max: %f\n',meantorr_estf_bookstein,stddevtorr_estf_bookstein,minimtorr_estf_bookstein,maximtorr_estf_bookstein));
            disp(sprintf('Rank-2: %d',rank(Ftorr_estf_bookstein)==2))            
            disp(sprintf('Time: %f',timetorr_estf_bookstein))            
            disp('*******************************************************************')
            if redraw, 
                fummatplot(Mini,Ftorr_estf_bookstein)   
            end
            bars=[bars;meantorr_estf_bookstein,stddevtorr_estf_bookstein];
            bartime=[bartime;timetorr_estf_bookstein];
            xlabelbars=[xlabelbars sprintf(' %d.- torr estf bookstein ',indexbars)];
            if mod(indexbars,8)==0, xlabelbars=[xlabelbars '\n']; end
            indexbars=indexbars+1;
            
        case 18,
            disp('Method: bookstein+sampson implemented by Torr')
            tic
            [f, f_sq_errors, n_inliers,inlier_index,Ftorr_estf_bookstein_sam] = torr_estimateF( M', 1, [], 'b+sampson', 0);
            timetorr_estf_bookstein_sam=toc;
            Ftorr_estf_bookstein_sam=T1'*Ftorr_estf_bookstein_sam'*T2;
            Ftorr_estf_bookstein_sam=Ftorr_estf_bookstein_sam./norm(Ftorr_estf_bookstein_sam);
            [meantorr_estf_bookstein_sam,stddevtorr_estf_bookstein_sam,minimtorr_estf_bookstein_sam,maximtorr_estf_bookstein_sam]=funmatError(Mini,Ftorr_estf_bookstein_sam);
            disp('Funamental Matrix')
            disp(Ftorr_estf_bookstein_sam)
            disp(sprintf('Distance point-epipolar line\nmean: %f  stdev: %f  min: %f  max: %f\n',meantorr_estf_bookstein_sam,stddevtorr_estf_bookstein_sam,minimtorr_estf_bookstein_sam,maximtorr_estf_bookstein_sam));
            disp(sprintf('Rank-2: %d',rank(Ftorr_estf_bookstein_sam)==2))            
            disp(sprintf('Time: %f',timetorr_estf_bookstein_sam))            
            disp('*******************************************************************')
            if redraw, 
                fummatplot(Mini,Ftorr_estf_bookstein_sam)   
            end
            bars=[bars;meantorr_estf_bookstein_sam,stddevtorr_estf_bookstein_sam];
            bartime=[bartime;timetorr_estf_bookstein_sam];
            xlabelbars=[xlabelbars sprintf(' %d.- torr estf bookstein sam ',indexbars)];
            if mod(indexbars,8)==0, xlabelbars=[xlabelbars '\n']; end
            indexbars=indexbars+1;
            
        case 19,
            disp('Method: linear+non linear implemented by Torr')
            tic
            [f, f_sq_errors, n_inliers,inlier_index,Ftorr_estf_lin_non_lin] = torr_estimateF( M', 1, [], 'lin+non_lin', 0);
            timetorr_estf_lin_non_lin=toc;
            Ftorr_estf_lin_non_lin=T1'*Ftorr_estf_lin_non_lin'*T2;
            Ftorr_estf_lin_non_lin=Ftorr_estf_lin_non_lin./norm(Ftorr_estf_lin_non_lin);
            [meantorr_estf_lin_non_lin,stddevtorr_estf_lin_non_lin,minimtorr_estf_lin_non_lin,maximtorr_estf_lin_non_lin]=funmatError(Mini,Ftorr_estf_lin_non_lin);
            disp('Funamental Matrix')
            disp(Ftorr_estf_lin_non_lin)
            disp(sprintf('Distance point-epipolar line\nmean: %f  stdev: %f  min: %f  max: %f\n',meantorr_estf_lin_non_lin,stddevtorr_estf_lin_non_lin,minimtorr_estf_lin_non_lin,maximtorr_estf_lin_non_lin));
            disp(sprintf('Rank-2: %d',rank(Ftorr_estf_lin_non_lin)==2))            
            disp(sprintf('Time: %f',timetorr_estf_lin_non_lin))            
            disp('*******************************************************************')
            if redraw, 
                fummatplot(Mini,Ftorr_estf_lin_non_lin)   
            end
            bars=[bars;meantorr_estf_lin_non_lin,stddevtorr_estf_lin_non_lin];
            bartime=[bartime;timetorr_estf_lin_non_lin];
            xlabelbars=[xlabelbars sprintf(' %d.- torr estf lin non lin ',indexbars)];
            if mod(indexbars,8)==0, xlabelbars=[xlabelbars '\n']; end
            indexbars=indexbars+1;
            
        case 20,
            disp('Method: seven points implemented by Torr')
            tic
            [nf,f]=torr_F_constrained_fit(M(1,1:7),M(2,1:7),M(3,1:7),M(4,1:7),1);
            timetorr_F_constrained_fit=toc;
            Ftorr_F_constrained_fit=[f(nf,1) f(nf,2) f(nf,3); f(nf,4) f(nf,5) f(nf,6); f(nf,7) f(nf,8) f(nf,9)];
            Ftorr_F_constrained_fit=T1'*Ftorr_F_constrained_fit'*T2;
            Ftorr_F_constrained_fit=Ftorr_F_constrained_fit./norm(Ftorr_F_constrained_fit);
            [meantorr_F_constrained_fit,stddevttorr_F_constrained_fit,minimtorr_F_constrained_fit,maximtorr_F_constrained_fit]=funmatError(Mini,Ftorr_F_constrained_fit);
            disp('Funamental Matrix')
            disp(Ftorr_F_constrained_fit)
            disp(sprintf('Distance point-epipolar line\nmean: %f  stdev: %f  min: %f  max: %f\n',meantorr_F_constrained_fit,stddevttorr_F_constrained_fit,minimtorr_F_constrained_fit,maximtorr_F_constrained_fit));
            disp(sprintf('Rank-2: %d',rank(Ftorr_F_constrained_fit)==2))            
            disp(sprintf('Time: %f',timetorr_F_constrained_fit))            
            disp('*******************************************************************')
            if redraw, 
                fummatplot(Mini,Ftorr_F_constrained_fit)   
            end
            bars=[bars;meantorr_F_constrained_fit,stddevttorr_F_constrained_fit];
            bartime=[bartime;timetorr_F_constrained_fit];
            xlabelbars=[xlabelbars sprintf(' %d.- torr F constrained fit ',indexbars)];
            if mod(indexbars,8)==0, xlabelbars=[xlabelbars '\n']; end
            indexbars=indexbars+1;
            
        case 21,
            disp('Method: MLESAC implemented by Torr')
            tic
            [f, f_sq_errors, n_inliers,inlier_index,Fmlesac] = torr_estimateF( M', 1, [round(log(1-0.99)/log(1-(1-0.25)^7)) 0.0039], 'mlesac', 0);
            [f,f_sq_errors] = torr_estimateF(M(:,inlier_index)', 1, [], 'non_linear',1,f);
            timeMLESAC=toc;
            Fmlesac=reshape(f,3,3);
            Fmlesac=T1'*Fmlesac*T2;
            Fmlesac=Fmlesac./norm(Fmlesac);
            [meanmlesac,stddevmlesac,minimmlesac,maximmlesac]=funmatError(Mini,Fmlesac);
            disp('Funamental Matrix')
            disp(Fmlesac)
            disp(sprintf('Distance point-epipolar line\nmean: %f  stdev: %f  min: %f  max: %f\n',meanmlesac,stddevmlesac,minimmlesac,maximmlesac));
            disp(sprintf('Rank-2: %d',rank(Fmlesac)==2))            
            disp(sprintf('Outliers: %d/%d',size(M,2)-size(Mini(:,inlier_index),2),size(M,2)))
            disp(sprintf('Time: %f',timeMLESAC))            
            disp('*******************************************************************')
            if redraw, 
                fummatplot(Mini,Fmlesac)   
            end
            redrawnow=1;
            Mfi=Mini(:,inlier_index);
            bars=[bars;meanmlesac,stddevmlesac];
            bartime=[bartime;timeMLESAC];
            xlabelbars=[xlabelbars sprintf(' %d.- MLESAC ',indexbars)];
            if mod(indexbars,8)==0, xlabelbars=[xlabelbars '\n']; end
            indexbars=indexbars+1;
            
        case 22,
            disp('Method: MAPSAC implemented by Torr')
            tic
            [f, f_sq_errors, n_inliers,inlier_index,Fmapsac] = torr_estimateF( M', 1, [round(log(1-0.99)/log(1-(1-0.25)^7)) 0.0039], 'mapsac', 0);
            [f,f_sq_errors] = torr_estimateF(M(:,inlier_index)', 1, [], 'non_linear',1,f);
            timeMAPSAC=toc;
            Fmapsac=reshape(f,3,3);
            Fmapsac=T1'*Fmapsac*T2;
            Fmapsac=Fmapsac./norm(Fmapsac);
            [meanmapsac,stddevmapsac,minimmapsac,maximmapsac]=funmatError(Mini(:,inlier_index),Fmapsac);
            disp('Funamental Matrix')
            disp(Fmapsac)
            disp(sprintf('Distance point-epipolar line\nmean: %f  stdev: %f  min: %f  max: %f\n',meanmapsac,stddevmapsac,minimmapsac,maximmapsac));
            disp(sprintf('Rank-2: %d',rank(Fmapsac)==2))            
            disp(sprintf('Outliers: %d/%d',size(M,2)-size(Mini(:,inlier_index),2),size(M,2)))
            disp(sprintf('Time: %f',timeMAPSAC))            
            disp('*******************************************************************')
            if redraw, 
                fummatplot(Mini(:,inlier_index),Fmapsac);
            end
            redrawnow=1;
            Mfi=Mini(:,inlier_index);
            bars=[bars;meanmapsac,stddevmapsac];
            bartime=[bartime;timeMAPSAC];
            xlabelbars=[xlabelbars sprintf(' %d.- MAPSAC ',indexbars)];
            if mod(indexbars,8)==0, xlabelbars=[xlabelbars '\n']; end
            indexbars=indexbars+1;
            
        case 23,
            disp('Method: FNS implemented by Anton van den Hengel')
            tic
            [f, f_sq_errors, n_inliers,inlier_index,Ffns] = torr_estimateF( M', 1, [500 1], 'fns', 0);         
            timeFNS=toc;
            Ffns=T1'*Ffns'*T2;
            Ffns=Ffns./norm(Ffns);
            [meanfns,stddevfns,minimfns,maximfns]=funmatError(Mini,Ffns);
            disp('Funamental Matrix')
            disp(Ffns)
            disp(sprintf('Distance point-epipolar line\nmean: %f  stdev: %f  min: %f  max: %f\n',meanfns,stddevfns,minimfns,maximfns));
            disp(sprintf('Rank-2: %d',rank(Ffns)==2))            
            disp(sprintf('Time: %f',timeFNS))            
            disp('*******************************************************************')
            if redraw, 
                fummatplot(Mini,Ffns)
            end
            bars=[bars;meanfns,stddevfns];
            bartime=[bartime;timeFNS];
            xlabelbars=[xlabelbars sprintf(' %d.- FNS ',indexbars)];
            if mod(indexbars,8)==0, xlabelbars=[xlabelbars '\n']; end
            indexbars=indexbars+1;
            
        case 24,
            disp('Method: CFNS implemented by Anton van den Hengel')
            tic
            [f, f_sq_errors, n_inliers,inlier_index,Fcfns] = torr_estimateF( M', 1, [500 1], 'cfns', 0);
            timeCFNS=toc;
            Fcfns=T1'*Fcfns'*T2;
            Fcfns=Fcfns./norm(Fcfns);
            [meancfns,stddevcfns,minimcfns,maximcfns]=funmatError(Mini,Fcfns);
            disp('Funamental Matrix')
            disp(Fcfns)
            disp(sprintf('Distance point-epipolar line\nmean: %f  stdev: %f  min: %f  max: %f\n',meancfns,stddevcfns,minimcfns,maximcfns));
            disp(sprintf('Rank-2: %d',rank(Fcfns)==2))
            disp(sprintf('Time: %f',timeCFNS))            
            disp('*******************************************************************')
            if redraw, 
                fummatplot(Mini,Fcfns) 
            end   
            bars=[bars;meancfns,stddevcfns];
            bartime=[bartime;timeCFNS];
            xlabelbars=[xlabelbars sprintf(' %d.- CFNS ',indexbars)];
            if mod(indexbars,8)==0, xlabelbars=[xlabelbars '\n']; end
            indexbars=indexbars+1;
            
    end
        
    if redrawnow & redraw,
        plotimagepoints(Mfi,nomia,nomib)
    end
%    pause
end

figure;
bar(bars);
legend('Mean','Stdev')               
xlabel(sprintf(xlabelbars))
ylabel(sprintf('pix.'))
title(['Distance points-epipolar lines'])
set(gca,'XTick',[1:size(method,2)]);
set(gca,'XLim',[0 size(method,2)+1]);
zoom on

figure;
bar(bartime);
xlabel(sprintf(xlabelbars))
ylabel(sprintf('sec.'))
title(['Time'])
set(gca,'XTick',[1:size(method,2)]);
set(gca,'XLim',[0 size(method,2)+1]);
zoom on
