function [saved_tform,needSave] = RegisterImages_FRONTIERS(B1, B2, varargin)
    %{
       REGISTERIMAGES_FRONTIERS Takes two single channel images and
       interactively aligns them. Initial automatic alignment is based on
       minimizing the differences between the two image intensities, and is
       iterated until a minimum is established (either after a threshold of
       overlap or number of cycles)

%B1 is static, B2 is moving. Runs on gray images
    %}

    % Start with a blank transform
    tform = affine2d();

    % Process arguments
    warning on
    doAuto = 1;
    autoAlignOnFullImage = 0;
    doBlur = 0;
    needSave = 1;
    scaleFactor = 1; %1/9

    autoTransformType = 'affine';
    clear startTform;
    if nargin > 2
        for i = 1:(nargin-2)
            %{
                If the value at varargin{i} isn't a string, skipping. If the value at varargin{i}
                is a string but an invalid one from the case list,
                dropping a warning and continuing. Before, it was erroring when
                varargin{i} pointed to the tform element,since it wasn't a string or constant.
                However, it would error after since there was no valid case statement.
            %}
            if ischar(varargin{i})==1
                switch varargin{i}
                    case 'AutomaticAlignment'
                        doAuto = varargin{i+1};
                    case 'StartingTransform'
                        startTform = varargin{i+1};
                    case 'AutoAlignOnFullImage'
                        autoAlignOnFullImage = 1;
                        scaleFactor = 1;
                    case 'ScalingFactor'
                        scaleFactor = varargin{i+1};
                    case 'AllowScaling'
                        autoTransformType = 'similarity';
                    case 'GaussianBlur'
                        doBlur = 1;
                        blurRadius = varargin{i+1};
                    otherwise
                        warning('Name argument of Name-Value pair is invalid. Proceeding with loop.')
                end
            end
        end
    end

    % Create working array
    wrk = {B1, B2};
    movTmp = wrk{2};

    % Resize if applicable
    if ~autoAlignOnFullImage
        for i = 1:2
            wrk{i} = imresize(wrk{i}, scaleFactor);
        end
        movTmp = wrk{2};
    end

    % Do a blur
    if doBlur
        filt = fspecial('gaussian', blurRadius, 0.5);
        for i = 1:2
            wrk{i} = imfilter(wrk{i}, filt);
        end
        movTmp = wrk{2};
    end

    % Start with this transformation
    if exist('startTform', 'var')
        % Convert to affine2d
        tform = affine2d(startTform);

        % Transform and show
        [movTmp, ~] = imwarp(B2, tform, 'OutputView',imref2d(size(B1)));
    
        figure
        imshow(cat(3,imadjust(im2double(movTmp)),imadjust(im2double(movTmp)),imadjust(im2double(B1))));

        % Ask for opinion on quality of alignment
        option = input('This is the existing alignment on file. Does it need to be adjusted? (yes/no) :: ', 's');
        close
        if strcmp(option, 'no')
            saved_tform = tform.T;
            return
        end
        % Otherwise
        option = questdlg('Reset alignment?', 'Reset Alignment', 'yes', 'no', 'no');
        if strcmp(option, 'yes')
            tform = affine2d();
        else
            if doAuto
                option = questdlg('Execute automatic alignment?', 'Autoalign', 'yes', 'no', 'yes');
                if strcmp(option, 'no')
                    doAuto = 0;
                end
            end

            tform.T(3,1:2) = tform.T(3,1:2) * scaleFactor;
            movTmp = imwarp(wrk{2}, tform, 'OutputView',imref2d(size(wrk{1})));
        end
    end

    % Do autoalignment
    if doAuto
        fprintf('Calculating autoalignment transformation... \n');

        [optimizer, metric] = imregconfig('multimodal');
        optimizer.Epsilon = 1.5E-25; %residual - decrease to improve accuracy
        optimizer.GrowthFactor = 1.0125;  %above 1, decreasing makes more accurate
        optimizer.MaximumIterations = 500; %increase to iterate more. 500 originally
        optimizer.InitialRadius=1.25E-3; %decreasing (or increasing) this value actually makes it more accurate, at the expense of time. Originally 6.25E-5

        autoTform = imregtform(adapthisteq(movTmp), adapthisteq(wrk{1}), autoTransformType, optimizer, metric);
        tform.T = tform.T * autoTform.T;
    else
        warning(' --- Automatic alignment has been disabled --- \n' );
    end

    isnt_good='yes';
    %setting up transformation dialogue box


    while strcmp(isnt_good,'yes')==1
        % Transform the original small image with the newest transform
        % FRESH each time. Do not replace the original small image to
        % accumlate transformations, properly calculate the accumulated
        % transformation before transforming each time.
        [movTmp, ~] = imwarp(wrk{2}, tform, 'OutputView',imref2d(size(wrk{1})));

        %equalize image before displaying
        figure
        
        imshow(cat(3,imadjust(im2double(movTmp)),imadjust(im2double(movTmp)),imadjust(im2double(wrk{1}))));
        % Default to yes because the false positive for yes is much better than the false positive for no.
        isnt_good=questdlg('Do the images need to be manually moved?','Manual Tweaking','yes','no','yes');
        if strcmp(isnt_good,'yes')==1

            %% Attempting to create a tform from cpmatch
            %If automatic alignment fails, offer the user the option of
            %manually selecting control point pairs on both images, and
            %generating a transformation matrix based on those
            [m_pts, f_pts]=cpselect(im2double(wrk{2}),im2double(wrk{1}),'wait',true); %[m_pts, f_pts]=cpselect(im2double(wrk{2})*20,im2double(wrk{1})*10,'wait',true);
            tform_cpselect=fitgeotrans(m_pts, f_pts,'nonreflectivesimilarity');
            tform=tform_cpselect;
        end
        close
    end
    
    % Modify the translation portion of the transform to match the
    % non-scaled image's dimentions.
    tform.T(3,1:2) = tform.T(3,1:2) * 1/scaleFactor;

    %output the transformation matrix to be used by
    %image_alignment_FRONTIERS
    saved_tform = tform.T;
    close(gcf);
end