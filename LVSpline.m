function [pts, MSE, time] = LVSpline(direc, contourIterations, contourPoints, interpScale, threshLow, threshHigh, initR, resX, resY, resZ)
% The LVSPLINE function uses a provided dataset of MRI cardiac images and
% the Snake function written by Dr. Kroon of University of Twente (available
% on MathWorks File Exchange) to segment the left ventricle. Based on
% these contours, we use a fourth-order polynomial to interpolate the
% x- and y-coordinates for each contour.
% 
% PTS = LVSpline is the interpolation of the acquired contours from MRI
% cardiac images. PTS consists of (x,y,z) tuples, meaning PTS(:,1) is the
% x-points, PTS(:,2) is the y-points, and PTS(:,3) is the z-points.
% 
% PTS = LVSpline(direc)
% Optional input. DIREC defines the folder in which the MRI cardiac
% images are stored. The folder can be listed relative to the current
% folder (a string starting with "./") or the complete folder path
% ("C:/MRIData/Set1", for example). Default value is "./Pics".
% 
% PTS = LVSpline(direc, contourIterations)
% Optional input. CONTOURITERATIONS defines the number of iterations to run
% the active contour algorithm per image. Default value is 300.
% 
% PTS = LVSpline(direc, contourIterations, contourPoints)
% Optional input. CONTOURPOINTS defines the number of points to calculate
% for each contour. Default value is 100.
% 
% PTS = LVSpline(direc, contourIterations, contourPoints, interpScale)
% Optional input. INTERPSCALE defines the number of "slices" to interpolate
% between given images. For example, if interpScale is 4, then there are 3
% additional contours interpolated between the first and second image.
% INTERPSCALE only affects datasets with at least two (2) images. Default
% value is 6.
% 
% PTS = LVSpline(direc, contourIterations, contourPoints, interpScale, threshLow, threshHigh)
% Optional inputs. THRESHLOW and THRESHHIGH define the low and high
% thresholding values, respectively for rescaling the grayscale image.
% Default values are 70 and 100 for low and high thresholds, respectively.
% 
% PTS = LVSpline(direc, contourIterations, contourPoints, interpScale, threshLow, threshHigh, initR)
% Optional input. INITR defines the initial radius for the auto-generated
% contours.  Default value is 5.
% 
% PTS = LVSpline(direc, contourIterations, contourPoints, interpScale, threshLow, threshHigh, initR, resX, resY, resZ)
% Optional inputs. RESX, RESY, RESZ define the image resolution (in mm) in the
% x-, y-, and z- dimensions, respectively. Default values are 1.5256 mm in
% the x- and y-dimensions and 10.02 mm in the z-dimension.
% 
% [PTS, MSE] = LVSpline stores the mean square errors (in mm) for each dimension
% into MSE, in x,y,z order.
% 
% [PTS, MSE, TIME] = LVSpline stores the runtime (in seconds) of LVSPLINE into TIME.
    
    close all;clc; % close all figures and clear Command Window
    t = tic;  % start timer
    time = 0; % total time for program
    
    % define low value for thresholding
    if ~exist('threshLow','var') || isempty(threshLow)
        imlo = 70;
    else
        imlo = threshLow;
    end
    
    % define high value for thresholding
    if ~exist('threshHigh','var') || isempty(threshHigh)
        imhi = 100;
    else
        imhi = threshHigh;
    end
    
    % Define the struct options for the Snake Algorithm
    Options = struct; % create struct object
    
    % define number of iterations for Snake algorithm
    if ~exist('contourIterations','var') || isempty(contourIterations) || contourIterations <= 0
        Options.Iterations = 300; 
    else
        Options.Iterations = contourIterations;
    end
    
    Options.Verbose = true; % define whether the  segmentation figure should be shown (opens in Figure 2)
    
    % define trial-and-error based MRI constants for left ventricle
    Options.Sigma1 = 1;
    Options.Sigma2 = 2;
    
    % number of contour points in final segmentation
    if ~exist('contourPoints','var') || isempty(contourPoints) || contourPoints <= 3
        Options.nPoints = 100; 
    else
        Options.nPoints = contourPoints;
    end
    
    % path to folder where the pictures are stored
    if ~exist('direc','var') || isempty(direc)
    	direc = "./Pics";
    end
    
    folder = dir(direc); % folder information
    nSlices = length(folder)-2; % number of pictures in folder
    
    % number of points to interpolate between pictures
    if ~exist('interpScale','var') || isempty(interpScale) || interpScale <= 1
        interpScale = 6;
    end
    
    nInterp = interpScale*(nSlices-1)+1; % number of total points to interpolate
    O = zeros(Options.nPoints,2,nSlices); % define a matrix of zeros for x- and y-coordinates after image segmentation
    
    
    % initial radius to define a self-generated initial contour
    if ~exist('initR','var') || isempty(initR) || initR <= 0
        r = 5;
    else
        r = initR;
    end
    
    theta = 0:pi/(Options.nPoints-1):2*pi; % define radians for self-generated initial contour
    Z =  zeros(Options.nPoints,nSlices); % define a matrix of zeros for z-coordinates after image segmentation
    
    % For each picture, perform the active contour segmentation based on user choice.
    % First two filenames are "." (for current folder) and ".." for drill-up one
    % folder.
    for i = 3:length(folder)
        file = folder(i); % get an image file from the folder
        name = file.name; % get the filename from the file
        
        disp(name); % show user the image being segmented
        image = imread(direc+"/"+name); % read the image into MATLAB
        image2 = im2double(image); % convert image to grayscale between 0 and 1
        
        % threshold and scale image
        image2 = (image2*255-imlo)/(imhi-imlo);
        image2(image2(:) <= 0) = 0;
        image2(image2(:) >= 1) = 1;
        
        % display image in Figure 1
        figure(1);imagesc(image2);colormap(gray);hold on
        
        choice = 0; % initialize user choice
        disp("Select method of segmenting image:"); % ask user how to segment image
        
        % while the user gives an invalid choice
        while choice == 0
            % options for segmentation
            disp("1. Auto-generate contour with center");
            disp("2. Manually select contour points");
            
            time = time + toc(t); % stop timer for user input
            choice = input("Selection: "); % get user choice
            t = tic; % restart timer
            
            % define different methods of segmenting based on image name
            if choice == 1
                % segment image using a self-generated contour based on
                % user selected center
                disp("select a center point (right-click on image)"); % direct user on how to segment image
                time = time + toc(t); % stop timer for user input
                [yc,xc] = getpts; % get center  of self-generated initial contour based on user input
                t = tic; % restart timer
                x = r.*cos(theta')+xc; % generate x-coordinates for initial contour
                y = r.*sin(theta')+yc; % generate y-coordinates for initial contour
            elseif choice == 2
                % segment image using custom contour defined by user
                % direct user on how to segment image
                disp("Select contour points in a clockwise manner");
                disp("Left-click for all points except for last");
                disp("Right-click for last point)");
                disp("Minimum of 4 points required");
                time = time + toc(t); % stop timer for user input
                [y,x] = getpts; % get x- and y-coordinates from selected points
                t = tic; % restart timer
            else
                % invalid choice. reset choice value and ask user to
                % re-select segmentation option.
                choice = 0;
                disp("Invalid selection. Please try again.")
            end
        end
        
        P = [x, y]; % pair the x- and y-coordinates for initial contour
        Z(:,i-2) = (length(folder)-i)*ones(size(Z,1),1); % define z-coordinates for current image
        [O(:,:,i-2), J] = Snake2D(image2,P,Options); % run Snake algorithm and find segmented contour
        disp(" "); % line breaks to separate information regarding each image
        disp(" ");
    end
    
    % create matrices for interpolated coordinates
    xf = zeros(Options.nPoints,nInterp);
    yf = zeros(Options.nPoints,nInterp);
    zf = zeros(Options.nPoints,nInterp);
    
    % find least-squares solution for each point in z-direction using
    % x(z) = c0*z^4 + c1*z^3 + c2*z^2 + c3*z + c4
    % y(z) = c0*z^4 + c1*z^3 + c2*z^2 + c3*z + c4
    
    % for each paired point in contours, find custom interpolation function (4th order polynomial)
    % and interpolate estimate for x,y locations
    for i=1:size(O,1)
        x = reshape(O(i,1,:),nSlices,1); % get x-coordinates for defined contour point
        y = reshape(O(i,2,:),nSlices,1); % get y-coordinates for defined contour point
        z = Z(i,:);z = z'; % get original z-coordinates for images
        Ax = [z.^4, z.^3, z.^2, z, ones(nSlices,1)]; % define matrix for A*x=z
        cx = Ax\x; % calculate constants for x(z)
        Ay = [z.^4, z.^3, z.^2, z, ones(nSlices,1)]; % define matrix for A*y=z
        cy = Ay\y; % calculate constants for y(z)
        
        newz = 1*(nSlices-1):-1*(nSlices-1)/(nInterp-1):0; newz = newz'; % get linearly interpolated z-coordinates
        zf(i,:) =  newz; % store interpolated z-coordinates
        xf(i,:) = [newz.^4, newz.^3, newz.^2, newz, ones(nInterp,1)]*cx; % get least-squares x-coordinates
        yf(i,:) = [newz.^4, newz.^3, newz.^2, newz, ones(nInterp,1)]*cy; % get least-squares y-coordinates
    end
    
    pts = [xf(:),yf(:),zf(:)]; % store x,y,z tuples for output
    
    % define resolutions in each dimension
    if ~exist('resX','var') || isempty(resX)
        resX = 1.5265;
    end
    if ~exist('resY','var') || isempty(resY)
        resY = 1.5265;
    end
    if ~exist('resZ','var') || isempty(resZ)
        resZ = 10.02;
    end
    
    % display 3D scatter plot of interpolation of segmented left ventricle
    figure(3);
    plot3(pts(:,1).*resX,pts(:,2).*resY,pts(:,3).*resZ,'.','MarkerSize',10);
    xlabel("X Pixel location (mm)");ylabel("Y Pixel location (mm)");zlabel("Ventricle Height (mm)");
    title("Interpolated, Segmented Left Ventricle");
    
    % display 3D scatter of original segmented left ventricle
    figure(4);
    plot3(reshape(O(:,1,:),[],1).*resX,reshape(O(:,2,:),[],1).*resY,Z(:).*resZ,'.','MarkerSize',10);
    xlabel("X Pixel location (mm)");ylabel("Y Pixel location (mm)");zlabel("Ventricle Height (mm)");
    title("Original Segmented Left Ventricle");
    
    % calculate and display mean square errors in pixels values (image slices are understood as z-pixels)
    MSEz = sqrt(sum((reshape(zf(:,1:ceil(nInterp/nSlices):end).*resZ,[],1)-Z(:).*resZ).^2)/length(Z(:)));
    MSEx = sqrt(sum((reshape(xf(:,1:ceil(nInterp/nSlices):end).*resX,[],1)-reshape(O(:,1,:),[],1).*resX).^2)/length(reshape(O(:,1,:),[],1)));
    MSEy = sqrt(sum((reshape(yf(:,1:ceil(nInterp/nSlices):end).*resY,[],1)-reshape(O(:,2,:),[],1).*resY).^2)/length(reshape(O(:,2,:),[],1)));
    MSE = [MSEx, MSEy, MSEz];
    disp("Mean Square Error (MSE) evaluation:")
    disp("MSEz = "+num2str(MSEz)+" mm");
    disp("MSEx = "+num2str(MSEx)+" mm");
    disp("MSEy = "+num2str(MSEy)+" mm");
    
    time = time + toc(t); % stop timer
    
    disp("Elapsed time: "+num2str(time)+" seconds."); % display how long the program ran, excluding user input timings

end