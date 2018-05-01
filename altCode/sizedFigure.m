function figHand = sizedFigure(figW,figH,ppi)

%creates a figure that will display at ppi pixels per inch on the screen
if ~exist('ppi','var')
	ppi = 100;
end

figHand = figure();
screenSize = get(0,'screensize');
figHand.Position = [(screenSize(3)-figW*ppi)/2 (screenSize(4)-figH*ppi)/2 figW*ppi figH*ppi];
paperSize = get(figHand,'PaperSize');
figHand.PaperPosition = [(paperSize(1)-figW)/2 (paperSize(2)-figH)/2 figW figH];
