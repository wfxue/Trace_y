function [model, pfit, residual, zcalc, ssr, exitflag] = ...
    MakeImgModelLSQ(x, y, z, w, p0, lb, ub, model)

%
% DESCRIPTION
% – Part of the semi automatic filament tracer, TraceFilamentMC.m to
% estimate filaments' central axis.
% – Fit a particle geometry to AFM surface topology data in least sqares
% sense
% – Part of TraceFilamentM.m and Trace_y
%
% DEPENDENCIES
% – Optimization Matlab toolbox
% – Part of and used by @AFMimage/TraceFilamentM.m
%
% INPUTS
% x, y, z  –  The image data
% w  –  Weight, has the same dimensions as z or scalar
% p0  –  The initial parameters
% lb, ub  –  Lower/upper bounds
% model  –  Vector of the model geometries to test, i.e. [1 2 ...]. Only 
% supply one number (scalr) to just fit one fixed model. The models are:
%   1: @FibEndLSQ
%   2: @FibSegmentLSQ
%   3: @FibParticleLSQ
%   4: @FibEnd2LSQ
%
% OUTPUTS
% model  –  The model number for the best fit model out of the ones to test
% pfit  –  Fitted model parameters
% residual  –  Residuals
% zcalc  –  Calculated model image
% ssr  –  sum of squared residuals
% exitflag  –  Information about the ptimisation from lsqnonlin
%
% AUTHORS
% Wei-Feng Xue
%
% HISTORY
% Updated with TraceFilamentM.m and see notes in TraceFilamentM.m.
%



%
% Setting up for LSQ
%

% Prep mesh grid
[xg, yg] = meshgrid(x, y);

lsqOpt = optimoptions('lsqnonlin', ...
    'Display', 'off', ...
    'MaxFunEvals', 1000000, 'MaxIter', 500, ...
    'TolX', 1e-14, 'TolFun', 1e-14);

% p: [x, y, sigma, h, theta, l, c2]
% Change for FibSegmentLSQlineConstrained
% p: [x, y, sigma, h, theta, l, c2, line_theta, line_rho]
% not all models use all parameters
modelFcns = {@FibEndLSQ ...
    @FibSegmentLSQ ...
    @FibParticleLSQ ...
    @FibEnd2LSQ ...
    @FibSegmentLineConstrainedLSQ};
modelFitPar = {(1:5) (1:5) (1:6) (1:6) ([3 4 5 9])};

% Incase of angular constraints are bad then alternative way to guess
% initial value
if (ub(5)-lb(5) >= 2*pi)
    fit = 5;
    theta = pi/36:pi/36:2*pi;
    ssrTheta = zeros(length(theta), 1);
    for bbb = 1:length(theta)
        [~, ~, ssrTheta(bbb)] = FibEndLSQ(theta(bbb));
    end
    p0(fit) = theta(find(ssrTheta == min(ssrTheta), 1, 'first'));
end

% Test all the models in [model] input and calc their AIC scores
aicBestFit = Inf;
for aaa = model
    modelFcn = modelFcns{aaa};
    fit = modelFitPar{aaa};
    
    % LSQ fit
    [p, ~, ~, exitflag] = lsqnonlin(modelFcn, p0(fit), lb(fit), ub(fit), lsqOpt);
    [~, ~, ssr, pfit] = modelFcn(p);
    %[residual, zcalc, ssr] = modelFcn(p(fit));
    
    % Check fit quality
    aic = 2*length(fit)+length(z(:))*log(ssr);
    if aic < aicBestFit
        %pBestFit = p;
        pBestFit = pfit;
        %p0(fit) = p(fit);
        aicBestFit = aic;
        bestFitModel = aaa;
        exitflagBestFitModel = exitflag;
    end
    
end

% Best model
model = modelFcns{bestFitModel};
fit = modelFitPar{bestFitModel};
p = pBestFit;
[residual, zcalc, ssr, pfit] = model(p(fit));
residual = residual./w;
model = bestFitModel;
exitflag = exitflagBestFitModel;
%aic = aicBestFit;



%
% LSQ functions
%

    function [residual, zcalc, ssr] = FibEndLSQ(pfit)
        
        p = p0;
        p(fit) = pfit;
        zcalc = FibEnd(xg, yg, p(1), p(2), p(3), p(4), p(5), p(7));
        
        residual = zcalc-z;
        %residual(:) = residual(:).*w(:);
        ssr = sum(residual(:).^2);
    end

    function [residual, zcalc, ssr] = FibSegmentLSQ(pfit)
        
        p = p0;
        p(fit) = pfit;
        zcalc = FibSegment(xg, yg, p(1), p(2), p(3), p(4), p(5), p(7));
        
        residual = zcalc-z;
        %residual(:) = residual(:).*w(:);
        ssr = sum(residual(:).^2);
    end

    function [residual, zcalc, ssr] = FibParticleLSQ(pfit)
        
        p = p0;
        p(fit) = pfit;
        zcalc = FibParticle(xg, yg, p(1), p(2), p(3), p(4), p(5), p(6));
        
        residual = zcalc-z;
        %residual(:) = residual(:).*w(:);
        ssr = sum(residual(:).^2);
    end

    function [residual, zcalc, ssr] = FibEnd2LSQ(pfit)
        
        p = p0;
        p(fit) = pfit;
        zcalc = FibParticle(xg, yg, p(1), p(2), p(3), p(4), p(5), p(6));
        
        residual = zcalc-z;
        %residual(:) = residual(:).*w(:);
        ssr = sum(residual(:).^2);
    end

    function [residual, zcalc, ssr, pout, rmsd] = FibSegmentLineConstrainedLSQ(pfit)
        
        p = p0;
        p(fit) = pfit;

        line_theta = p(8);
        line_rho = p(9);
        x0 = p(1);
        y0 = p(2);

        [xx, yy] = pol2cart(line_theta, line_rho);
        xx = x0+xx;
        yy = y0+yy;

        zcalc = FibSegment(xg, yg, xx, yy, p(3), p(4), p(5), p(7));
        
        residual = zcalc-z;
        %residual(:) = residual(:).*w(:);
        ssr = sum(residual(:).^2);
        pout = p0;
        pout(fit) = pfit;
        pout(1) = xx;
        pout(2) = yy;
        rmsd = sqrt(sum(residual(:).^2)/numel(xg));
    end


%
% Surface geometries
%

    function [z] = FibEnd(xg, yg, x, y, sigma, h, theta, c2)
        
        f = -(xg-x)*sin(theta+pi/2)+(yg-y)*cos(theta+pi/2)+0.5;
        f(f < 0) = 0;
        f(f > 1) = 1;
        a = 1/(2*sigma^2);
        z = h*(f.*exp(-a*((xg-x).^2+(yg-y).^2))+...
            (1-f).*exp(-a*(-(xg-x)*sin(theta)+(yg-y)*cos(theta)+...
            c2*(-(xg-x)*cos(theta)-(yg-y)*sin(theta)).^2).^2));
        
        
    end

    function [z] = FibSegment(xg, yg, x, y, sigma, h, theta, c2)
        
        %[xg, yg] = meshgrid(-ceil(10*sigma):0.5:ceil(10*sigma));
        a = 1/(2*sigma^2);
        z = h*exp(-a*(-(xg-x)*sin(theta)+(yg-y)*cos(theta)+...
            c2*(-(xg-x)*cos(theta)-(yg-y)*sin(theta)).^2).^2);
        %surface(xg, yg, z);
        
        
    end

    function [z] = FibParticle(xg, yg, x, y, sigma, h, theta, l)
        
        xe = x+l*cos(theta);
        ye = y+l*sin(theta);
        f1 = -(xg-x)*sin(theta+pi/2)+(yg-y)*cos(theta+pi/2)+0.5;
        f1( f1 < 0 ) = 0;
        f1( f1 > 1 ) = 1;
        f2 = -(xg-xe)*sin(theta-pi/2)+(yg-ye)*cos(theta-pi/2)+0.5;
        f2( f2 < 0 ) = 0;
        f2( f2 > 1 ) = 1;
        
        a = 1/(2*sigma^2);
        z = h*(f1.*exp(-a*((xg-x).^2+(yg-y).^2))+...
            f2.*exp(-a*((xg-xe).^2+(yg-ye).^2))+...
            (1-f1).*(1-f2).*exp(-a*(-(xg-x)*sin(theta)+(yg-y)*cos(theta)).^2));
        
        
    end



end
