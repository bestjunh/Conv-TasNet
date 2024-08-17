function [dist, azi, ele] = locDiffuseWall(nx,ny,roomDim,centerSensors,distWall)
    nSource = 2*nx + 2*ny + 4;

    dist = zeros(1,nSource);
    azi = zeros(1,nSource);
    ele = zeros(1,nSource);
    loc = zeros(nSource,3);
    i = 0;
    % Location of coner
    i = i + 1;
    loc(i,:) = [distWall, distWall, centerSensors(3)];
    [dist(i), azi(i), ele(i)] = Cart2Sphe(loc(i,:)-centerSensors);

    i = i + 1;
    loc(i,:) = [distWall, roomDim(2)-distWall, centerSensors(3)];
    [dist(i), azi(i), ele(i)] = Cart2Sphe(loc(i,:)-centerSensors);

    i = i + 1;
    loc(i,:) = [roomDim(1)-distWall, roomDim(2)-distWall, centerSensors(3)];
    [dist(i), azi(i), ele(i)] = Cart2Sphe(loc(i,:)-centerSensors);

    i = i + 1;
    loc(i,:) = [roomDim(1)-distWall, distWall, centerSensors(3)];
    [dist(i), azi(i), ele(i)] = Cart2Sphe(loc(i,:)-centerSensors);

    

    % Location of x-axis
    lenX = roomDim(1) - 2*distWall;
    unitX = lenX/(nx+1);
    for j = 1:nx
        i = i + 1;
        loc(i,:) = [distWall + j*unitX, distWall, centerSensors(3)];
        [dist(i), azi(i), ele(i)] = Cart2Sphe(loc(i,:)-centerSensors);
        i = i + 1;
        loc(i,:) = [distWall + j*unitX, roomDim(2)-distWall, centerSensors(3)];
        [dist(i), azi(i), ele(i)] = Cart2Sphe(loc(i,:)-centerSensors);
    end

    % Location of y-axis
    lenY = roomDim(2) - 2*distWall;
    unitY = lenY/(ny+1);
    for j = 1:ny
        i = i + 1;
        loc(i,:) = [distWall, distWall + j*unitY, centerSensors(3)];
        [dist(i), azi(i), ele(i)] = Cart2Sphe(loc(i,:)-centerSensors);
        i = i + 1;
        loc(i,:) = [roomDim(1)-distWall, distWall + j*unitY, centerSensors(3)];
        [dist(i), azi(i), ele(i)] = Cart2Sphe(loc(i,:)-centerSensors);
    end


    % figure
    if 0
        figure
        hold on
        for i = 1:nSource
            plot3(loc(:,1),loc(:,2),loc(:,3),'o')
        end
        plot3(centerSensors(1),centerSensors(2),centerSensors(3),'x')
        axis equal
        xlim([0 roomDim(1)])
        ylim([0 roomDim(2)])
        zlim([0 roomDim(3)])
        grid on
    end
    
end