function [imh,fh,jFrame] = jFullScreen(screen)
mon = get(0,'Monitor');
fh = figure('MenuBar','none','Color','k');
pos = get(fh,'Position');
pos(1:2)= [mon(1,1),mon(1,2)];
pos(3:4) = [mon(screen,3) - mon(screen,1) + 1, mon(screen,4) - mon(screen,2) + 1];
set(fh,'Position',pos);
ah = axes('Parent',fh,'Units','normalized','Position',[0,0,1,1]);
drawnow;
warning off MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame
mjf = get(fh,'JavaFrame');
% mjf.setMaximized(true);

kScreen = zeros(mon(screen,4) - mon(screen,2) + 1,mon(screen,3) - mon(screen,1) + 1,3,'uint8');
wScreen = ones(size(kScreen),'uint8') .* 255;
rScreen = uint8(rand(size(wScreen),'single') .* 255);
ssize = size(rScreen);
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

imh = image(kScreen,'Parent',ah);
setappdata(imh,'kScreen',kScreen);
setappdata(imh,'wScreen',wScreen);
setappdata(imh,'rScreen',rScreen);
setappdata(imh,'vLine',vLine);
setappdata(imh,'hLine',hLine);
set(ah,'XTick',[],'XTickMode','manual','YTick',[],'YTickMode','manual');
axis image;
drawnow
jWindow = mjf.fHG1Client.getWindow;
mjc = jWindow.getContentPane;
% mjr = jWindow.getRootPane;
figTitle = jWindow.getTitle;
jDesktop = com.mathworks.mde.desk.MLDesktop.getInstance;
jFrame = javaObjectEDT(com.mathworks.widgets.desk.DTSingleClientFrame(jDesktop,figTitle));
jFrame.setUndecorated(true);
jFrame.setContentPane(mjc);
jFrame.setLocation(mon(screen,1),mon(screen,2));
jFrame.setMaximized(true);
% mjf.setMinimized(true);
pos(1) = max(mon(:,3)) + 100;
set(fh,'Position',pos);
jFrame.setVisible(true);
set(fh,'CloseRequestFcn',{@jFSClose,jFrame});


function jFSClose(obj,~,jFrame)
dispose(jFrame);
delete(obj);