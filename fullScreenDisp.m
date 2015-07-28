function [fh,ah,imh] = fullScreenDisp(screen,calibrate,varargin)

FSOffset = [-8,-8,17,38]; % Windows 7 figure borders - may need adjusting for others
mon = get(0,'Monitor');

if isempty(screen)
    screen = 2;
end

if numel(varargin) == 2
    xOff = varargin{1};
    yOff = varargin{2};
else
    xOff = 0; % The first 'X' pixel on screen two is the correct same as Monitor value
    yOff = 30; % Because the first 'Y' pixel on screen 2 happens to be 30 on this system
end
xpos = mon(screen,1) + xOff;
ypos = mon(screen,2) + yOff;
wpos = mon(screen,3) - mon(screen,1) + 1;
hpos = mon(screen,4) - mon(screen,2) + 1;
pos = [xpos,ypos,wpos,hpos] + FSOffset;

kScreen = zeros(mon(screen,4) - mon(screen,2) + 1,mon(screen,3) - mon(screen,1) + 1,3,'uint8');
wScreen = ones(size(kScreen),'uint8') .* 255;
rScreen = uint8(rand(size(wScreen),'single') .* 255);
ssize = size(kScreen);
if ssize(1) > ssize(2)
    lineWidth = uint16(size(kScreen,2) / 40);
else
    lineWidth = uint16(size(kScreen,1) / 40);
end
center = uint16(size(kScreen,2) / 2);
vLine = kScreen;
vLine(:,center-lineWidth:center+lineWidth,:) = 255;
center = uint16(size(kScreen,1) / 2);
hLine = kScreen;
hLine(center-lineWidth:center+lineWidth,:,:) = 255;

if ~isempty(calibrate) && isa(calibrate,'logical') && calibrate == true
    disp('Assisted calibration not yet implemented.')
end

fh = figure('Visible','off','MenuBar','none','Units','pixels','OuterPosition',pos);
setappdata(fh,'kScreen',kScreen);
setappdata(fh,'wScreen',wScreen);
setappdata(fh,'rScreen',rScreen);
setappdata(fh,'vLine',vLine);
setappdata(fh,'hLine',hLine);

ah = axes('Parent',fh,'Units','normalized','Position',[0,0,1,1]);
imh = image(rScreen,'Parent',ah);
axis('image');
set(ah,'XTick',[],'YTick',[],'XTickMode','manual','YTickMode','manual');
set(fh,'Visible','on')