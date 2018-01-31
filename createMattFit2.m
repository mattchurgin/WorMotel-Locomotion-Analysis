function [fitresult, gof] = createMattFit2(ax, fh)
%createMattFit2(AX,FH)
%  Create a fit for activity histogram data captured with WorMotel.
%
%  Data for 'Matt_fit_1' fit:
%      X Input : ax: x values
%      Y Output: fh: y values
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%  NOTE: You may need to adjust fitting initial values to obtain good fits.
%  Initial values are located on line 29


%% Fit
[xData, yData] = prepareCurveData( ax, fh );

% Set up fittype and options.
ft = fittype( 'a*exp(-b*x)+c*exp(-d*(x-e).^2)+f*exp(-g*x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.DiffMaxChange = 0.001;
opts.DiffMinChange = 1e-09;
opts.Display = 'Off';
opts.Lower = [0 0 0 1e-6 200 0 0];
opts.Upper = [Inf .75 Inf 1e-3 1200 Inf 0.75]; 
opts.MaxFunEvals = 1000;
opts.MaxIter = 1000;
%Initial values (adjust these as necessary to obtain good fits
opts.StartPoint = [500 0.05 1000 1e-5 500 100 0.01 ]; 
opts.TolFun = 1e-09;
opts.TolX = 1e-09;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
%figure( 'Name', 'Matt_fit_1' );
h = plot( fitresult, xData, yData );
legend( h, 'fh vs. ax', 'Matt_fit_1', 'Location', 'NorthEast' );
% Label axes
xlabel ax
ylabel fh
grid on


