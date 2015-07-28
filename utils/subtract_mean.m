function Iout = subtract_mean(Iin);
% 
% It = It-mean(It,2)*ones(1,size(It,2));
% It = It - ones(size(It,1),1)*mean(It,1);
% 
% It = It - mean(It(:));
 
Iout = 0*Iin; 
for cc = 1:size(Iin,3)
    It = Iin(:,:,cc);
    Iout(:,:,cc) = It-mean(It,2)*ones(1,size(It,2)) - ones(size(It,1),1)*mean(It,1) + mean(It(:));
end