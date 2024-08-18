clear all;
addpath('functions');

if ispc
    nasPath = 'Z:/nas1_data/';
    nas3Path = 'Y:/nas3_data/'
    devPath = 'W:/';
elseif isunix
    nasPath = '/home/nas/';
    nas3Path = '/home/nas3/';
    devPath = '/home/';
else
    disp('Unknown operating system.');
end

RIRPATH = [devPath 'data/albert/DB/CONV-TASNet-RIR-v2/target/']
upFs = 3*16000;
N = 200;

RT60 = [0.2 0.3 0.4 0.5 0.6 0.7];

ROOMDIM = {[4.01, 3.65, 2.81]
           [6.02, 3.82, 2.53]
           [3.35, 4.94, 2.83]
           [4.55, 5.76, 2.61]
           [5.33, 3.51, 2.72]
           [4.34, 6.02, 2.55]
           [4.55, 3.54, 2.81]
           [4.76, 3.86, 2.74]
           [4.87, 3.78, 2.57]
           [5.08, 4.02, 2.64]
           [5.94, 7.45, 2.82]
           [7.40, 5.55, 2.87]
           [7.13, 5.78, 2.82]
           [7.23, 5.22, 2.91]
           [8.12, 5.45, 3.01]
           [8.20, 5.55, 2.89]
           [6.63, 4.75, 2.84]
           [6.30, 6.15, 2.92]
           [8.40, 5.45, 2.95]
           [8.10, 6.42, 3.00]};
          

LOCDELTA = {[-0.52 -0.15,-0.11]
            [-0.35 -0.41 0.09]
            [-0.16 0.38 -0.10]
            [-0.45 0.20 0.12]
            [0.45 -0.18 -0.14]
            [0.55 -0.33 0.08]
            [0.31 0.22 0.02]
            [0.27 0.43 -0.03]};

locSensors = [0 0 0];
for rt60 = RT60
    for roomDim_idx = 1:length(ROOMDIM)
        roomDim = cell2mat(ROOMDIM(roomDim_idx));
        for locDelta_idx = 1:length(LOCDELTA)            
            locDelta = cell2mat(LOCDELTA(locDelta_idx));
            centerSensors = [roomDim(1)/2 roomDim(2)/2 1.0] + locDelta;
            rir_path = [RIRPATH 'RT' num2str(rt60) '/ROOM' num2str(roomDim(1)) 'x' num2str(roomDim(2)) 'x' num2str(roomDim(3)) '/LOC' num2str(locDelta(1)) 'x' num2str(locDelta(2)) 'x' num2str(locDelta(3)) '/'];
            mkdir(rir_path)
            disp(['rt60=' num2str(rt60) ', roomDim=' num2str(roomDim_idx) '/' num2str(length(ROOMDIM)) ', locDelta=' num2str(locDelta_idx) '/' num2str(length(LOCDELTA))]);
            tic
            for n = 1:N
                x = rand*(roomDim(1)-0.6)+0.3;
                y = rand*(roomDim(2)-0.6)+0.3;
                z = rand*(roomDim(3)-0.6)+0.3;
                [r,azi,ele] = Cart2Sphe([x-centerSensors(1),y-centerSensors(2),z-centerSensors(3)].');

                generateRIR(upFs, roomDim, rt60, centerSensors, locSensors, r, azi, ele, rir_path);

                if 0
                    figure;
                    plot3(centerSensors(1), centerSensors(2), centerSensors(3), 'ro');
                    hold on;
                    plot3(x, y, z, 'bo');
                    hold off;
                    xlabel('x');
                    ylabel('y');
                    zlabel('z');
                    grid minor;
                    xlim([0 roomDim(1)]);
                    ylim([0 roomDim(2)]);
                    zlim([0 roomDim(3)]);
                    view(3);
                    title(['rt60=' num2str(rt60) ', roomDim=' num2str(roomDim) ', locDelta=' num2str(locDelta)]);
                end                
            end
            toc
        end
    end
end