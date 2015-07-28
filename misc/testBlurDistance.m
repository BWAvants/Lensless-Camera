% Experiement to compute 
% distance between sensor and mask (d) 
% and blur size caused by the mask (m)

% monitor distance from the mask is approximately t = 18cm
% size of substrate that is imaged ~ 1.2cm 
% view of substrate on the sensor is approx. 400 pixels
% therefore, size of image is q*1.3/400 in cm

% height of pattern on the monitor = ps or p5 for a square or 5 image
% pixel height on the sensor = qs or q5

% square with shrinking size
qs = [27; 23; 18; 17; 15];
ps = [3.015; 2.23; 1.57; 0.774; 0.37];
As = [ps/18 ones(length(ps),1)];

% 5 with shrinking size
q5 = [90;70;50;33;22;20];
p5 = [13.6; 10.1; 6.7; 3.5; 1.75; 1.35];
A5 = [p5/18 ones(length(p5),1)];

% height of the sensor image is approx..
% ps*d/t+m

dm = [A5; As]\(1.3*[q5; qs]/400) % in cm 

figure(134); plot([(1.3*[q5; qs]/400) [A5;As]*dm])

% system performance is limited by the blur spot size
% in our case spot size is dm(2) = 350 micron (approx.)
% which means that any feature shift on the monitor must be greater than 
% half the spot size ... 
% 
% another explanation is that we can only fit non-overlapping spot-sized
% pixels on the sensor... 
% similarly, visible region on the monitor should map over the entire
% sensor for maximum resolution... 

%% 32R10 mask
% 8, 16, 32... 
pl = [2.97; 1.467; 0.76]; 
ql = [10; 5; 2.5];
[(1.3*ql/400) pl*0.2/18]