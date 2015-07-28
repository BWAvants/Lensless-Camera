function captureImage(mouseLocation,delayTimeDisplay,delayTimeCapture,numberOfFrames,nameImg,dirPath,tempPath)
    import java.awt.Robot;
    import java.awt.event.*;
    mouse = Robot;
    set(0,'PointerLocation',mouseLocation);
    pause(delayTimeDisplay)
    
    mouse.mousePress(InputEvent.BUTTON1_MASK);
    mouse.mouseRelease(InputEvent.BUTTON1_MASK);
    pause(delayTimeCapture) 
    l = [];
    while length(l)< numberOfFrames+2
        l = dir(tempPath);
    end
    l = l(3:end);
    try 
        for i = 1:numberOfFrames
            % fprintf('%d',i);
            movefile([tempPath l(i).name],[dirPath nameImg '_' num2str(i) '.bmp']);
        end
    catch err
        pause(delayTimeCapture)

        l = dir(tempPath);
        while length(l) > numberOfFrames + 2
            l = dir(tempPath);
        end
        l = l(3:end);
        for i = 1:numberOfFrames
            % fprintf('%d',i);
            movefile([tempPath l(i).name],[dirPath nameImg '_' num2str(i) '.bmp']);
        end
        % rethrow(err)
    end
    % fprintf('\n');

end