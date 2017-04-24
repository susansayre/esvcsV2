function[h,L,MX,MED,F,U]=myViolin(Y,probDensity,varargin)

% create a violin plot based on a given probability density function. Modified from 

% Hoffmann H, 2015: violin.m - Simple violin plot using matlab default kernel
% density estimation. INRES (University of Bonn), Katzenburgweg 5, 53115 Germany.
% hhoffmann@uni-bonn.de
%
% Modifications designed to handle data specific to this model where the input data is a vector of possible outcomes and
% associated probabilities. We do not need to estimate the kernel density since we can directly calculate the
% probability  based on problem parameters. Challenge is that our input vectors have many values censored at 0
%
% INPUT
%
% Y:     Data to be plotted, being either
%        a) n x m matrix. A 'violin' is plotted for each column m, OR
%        b) 1 x m Cellarry with elements being numerical colums of nx1 length.
%
% varargin:
% xlabel:    xlabel. Set either [] or in the form {'txt1','txt2','txt3',...}
% facecolor: FaceColor. (default [1 0.5 0]); Specify abbrev. or m x 3 matrix (e.g. [1 0 0])
% edgecolor: LineColor. (default 'k'); Specify abbrev. (e.g. 'k' for black); set either [],'' or 'none' if the mean should not be plotted
% facealpha: Alpha value (transparency). default: 0.5
% mc:        Color of the bars indicating the mean. (default 'k'); set either [],'' or 'none' if the mean should not be plotted
% medc:      Color of the bars indicating the median. (default 'r'); set either [],'' or 'none' if the mean should not be plotted
%
%
% OUTPUT
%
% h:     figure handle
% L:     Legend handle
% MX:    Means of groups
% MED:   Medians of groups
% F:	 Final probability values plotted
% U:     Final points plotted
%__________________________________________________________________________

xL=[];
fc=[1 0.5 0];
lc='k';
alp=0.5;
mc='k';
medc='r';
b=[]; %bandwidth
plotlegend=1;
plotmean=1;
plotmedian=1;
x = [];
%_____________________

%convert single columns to cells:
if iscell(Y)==0
    Y = num2cell(Y,1);
end

%get additional input parameters (varargin)
if isempty(find(strcmp(varargin,'xlabel')))==0
    xL = varargin{find(strcmp(varargin,'xlabel'))+1};
end
if isempty(find(strcmp(varargin,'facecolor')))==0
    fc = varargin{find(strcmp(varargin,'facecolor'))+1};
end
if isempty(find(strcmp(varargin,'edgecolor')))==0
    lc = varargin{find(strcmp(varargin,'edgecolor'))+1};
end
if isempty(find(strcmp(varargin,'facealpha')))==0
    alp = varargin{find(strcmp(varargin,'facealpha'))+1};
end
if isempty(find(strcmp(varargin,'mc')))==0
    if isempty(varargin{find(strcmp(varargin,'mc'))+1})==0
        mc = varargin{find(strcmp(varargin,'mc'))+1};
        plotmean = 1;
    else
        plotmean = 0;
    end
end
if isempty(find(strcmp(varargin,'medc')))==0
    if isempty(varargin{find(strcmp(varargin,'medc'))+1})==0
        medc = varargin{find(strcmp(varargin,'medc'))+1};
        plotmedian = 1;
    else
        plotmedian = 0;
    end
end

if isempty(find(strcmp(varargin,'plotlegend')))==0
    plotlegend = varargin{find(strcmp(varargin,'plotlegend'))+1};
end
if isempty(find(strcmp(varargin,'x')))==0
    x = varargin{find(strcmp(varargin,'x'))+1};
end

if any(strcmp(varargin,'ubProb'))
    probAtMax = varargin{find(strcmp(varargin,'ubProb'))+1};
else
	probAtMax = 0;
end

if size(fc,1)==1
    fc=repmat(fc,size(Y,2),1);
end

i=1;
for i=1:size(Y,2)
	%check whether we have a zero
    if Y{i}(1) == 0
		probZero(i) = probDensity{i}(1);
	end
	maxf(i) = max(probDensity{i});
end
maxF = max(maxf);
for i=1:size(Y,2)
    probDensity{i}=probDensity{i}/maxF*0.3;
	cumProb = cumsum(probDensity{i})./sum(probDensity{i});
	belowMedInds = find(cumProb<=0.5);
	[~,lowMedInd] = max(cumProb(belowMedInds));
	lowInd=belowMedInds(lowMedInd);
	aboveMedInds = find(cumProb>=.5);
	[~,highMedInd] = min(cumProb(aboveMedInds));
	highInd=aboveMedInds(highMedInd);
	if any(belowMedInds)
		MED(:,i)=(Y{i}(lowInd)*probDensity{i}(lowInd)+Y{i}(highInd)*probDensity{i}(highInd))./(probDensity{i}(lowInd)+probDensity{i}(highInd));
	else
		MED(:,i)= min(Y{i});
	end
    MX(:,i)= (Y{i}'*probDensity{i})/sum(probDensity{i});
   
end
%%
%-------------------------------------------------------------------------
% Put the figure automatically on a second monitor
% mp = get(0, 'MonitorPositions');
% set(gcf,'Color','w','Position',[mp(end,1)+50 mp(end,2)+50 800 600])
%-------------------------------------------------------------------------
%Check x-value options
if isempty(x)
    x = zeros(size(Y,2));
    setX = 0;
else
    setX = 1;
    if isempty(xL)==0
        disp('_________________________________________________________________')
        warning('Function is not designed for x-axis specification with string label')
        warning('when providing x, xlabel can be set later anyway')
        error('please provide either x or xlabel. not both.')
    end
end

%% Plot the violins
i=1;
for i=i:size(Y,2)
    if isempty(lc) == 1
        if setX == 0
            h(i)=fill([probDensity{i}+i;flipud(i-probDensity{i})],[Y{i};flipud(Y{i})],fc(i,:),'FaceAlpha',alp,'EdgeColor','none');
        else
            h(i)=fill([probDensity{i}+x(i);flipud(x(i)-probDensity{i})],[Y{i};flipud(Y{i})],fc(i,:),'FaceAlpha',alp,'EdgeColor','none');
        end
    else
        if setX == 0
            h(i)=fill([probDensity{i}+i;flipud(i-probDensity{i})],[Y{i};flipud(Y{i})],fc(i,:),'FaceAlpha',alp,'EdgeColor',lc);
        else
            h(i)=fill([probDensity{i}+x(i);flipud(x(i)-probDensity{i})],[Y{i};flipud(Y{i})],fc(i,:),'FaceAlpha',alp,'EdgeColor',lc);
        end
    end
    hold on
    if setX == 0
        if plotmean == 1 && numel(Y{i})>1
            p(1)=plot([interp1(Y{i},probDensity{i}+i,MX(:,i)), interp1(flipud(Y{i}),flipud(i-probDensity{i}),MX(:,i)) ],[MX(:,i) MX(:,i)],mc,'LineWidth',2);
        end
		if any(probZero)
			p(2)=plot([i-probDensity{i}(1) probDensity{i}(1)+i],[0 0],'b','LineWidth',2);
		end
		if any(probAtMax)
			p(3) = plot([i-probAtMax(i) i+probAtMax(i)],[max(Y{i}) max(Y{i})],'b','LineWidth',2);
		end
        if plotmedian == 1 && numel(Y{i})>1
            p(4)=plot([interp1(Y{i},probDensity{i}+i,MED(:,i)), interp1(flipud(Y{i}),flipud(i-probDensity{i}),MED(:,i)) ],[MED(:,i) MED(:,i)],medc,'LineWidth',2,'LineStyle','--');
		end
    elseif setX == 1
        if plotmean == 1
            p(1)=plot([interp1(Y{i},probDensity{i}+i,MX(:,i))+x(i)-i, interp1(flipud(Y{i}),flipud(i-probDensity{i}),MX(:,i))+x(i)-i],[MX(:,i) MX(:,i)],mc,'LineWidth',2);
        end
        if plotmedian == 1
            p(2)=plot([interp1(Y{i},probDensity{i}+i,MED(:,i))+x(i)-i, interp1(flipud(Y{i}),flipud(i-probDensity{i}),MED(:,i))+x(i)-i],[MED(:,i) MED(:,i)],medc,'LineWidth',2);
        end
    end
end

%% Add legend if requested
if plotlegend==1 & plotmean==1 | plotlegend==1 & plotmedian==1
    
    if plotmean==1 & plotmedian==1
        L=legend([p(1) p(4) p(2)],'Mean','Median','Prob at Boundary','location','southOutside');
    elseif plotmean==0 & plotmedian==1
        L=legend([p(2)],'Median');
    elseif plotmean==1 & plotmedian==0
        L=legend([p(1)],'Mean');
    end
    
    set(L,'box','off','FontSize',14)
else
    L=[];
end

%% Set axis
for ii=1:size(Y,2)
	minY(ii) = min(Y{ii});
	maxY(ii) = max(Y{ii});
end
minY = min(minY); maxY = max(maxY);
if setX == 0
    axis([0.5 size(Y,2)+0.5, minY maxY]);
elseif setX == 1
    axis([min(x)-0.05*range(x) max(x)+0.05*range(x), minY maxY]);
end

%% Set x-labels
xL2={''};
i=1;
for i=1:size(xL,2)
    xL2=[xL2,xL{i},{''}];
end
set(gca,'TickLength',[0 0],'FontSize',12)
box on

if isempty(xL)==0
    set(gca,'XtickLabel',xL2)
end
%-------------------------------------------------------------------------
end %of function