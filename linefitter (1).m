

% %this is a sample function you can use to see if it works.
% x=1:100;
% y=4*x+5;
% y=y.*(1+rand(size(y)));
% dy=(y.^(2/3)).*rand(size(y));


%x=    ;% put your x data here
%y=    ;put your y data here
%dy=   ;% put your uncertainties here

%w=ones(size(y)); use this if you are fitting without uncertainty
w=1./dy.^2;

%%this is where the meat happens
[fitout,resid,J,cov,mse]=nlinfit( x,  y, @linefunc, [1 1],'Weight',w);
ci = nlparci(fitout,resid,'jacobian',J);

FitSlope=(ci(2,2)+ci(2,1))/2;
DeltaFitSlope=(ci(2,2)-ci(2,1))/2;

FitInt=(ci(1,2)+ci(1,1))/2;
DeltaFitInt=(ci(2,2)-ci(1,1))/2;

%% outputting and plotting
'slope'
[ FitSlope DeltaFitSlope]
'intercept'
[ FitInt DeltaFitInt]


plot(linefunc([FitInt FitSlope],x));
hold on
errorbar(x,y,dy,'o');

hold off