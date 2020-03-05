%Code to pre-process Outfi data 
%Necessary to run prior to running newdemograph.m

function cellList=cellListMod(cellList)
for frame=1:length(cellList)
    for cells=1:length(cellList{frame})
        if ~isempty(cellList{frame}{cells})
            cellList{frame}{cells}=getextradata(cellList{frame}{cells});
        end
    end
end
end

function str=getextradata(str)
% calculates geometrical properties of detected meshes in structure
% "str", corresponding to a single cell in "cellList", recording to the
% same structure. The properties are: "steplength", "steparea",
% "stepvolume" (length, area, and volume of each segment), "length",
% "area", "volume" (foe the whole cell), "lengthvector" - coordinates
% along cell centerline of the centers of each segment.
if isempty(str), return; end
if isfield(str,'mesh') && length(str.mesh)>1
    mesh = str.mesh;
    lng = size(mesh,1)-1;
    if ~isfield(str,'polarity'), str.polarity=0; end
    %length
    str.steplength = edist(mesh(2:end,1)+mesh(2:end,3),mesh(2:end,2)+mesh(2:end,4),...
        mesh(1:end-1,1)+mesh(1:end-1,3),mesh(1:end-1,2)+mesh(1:end-1,4))/2;
    str.length = sum(str.steplength);
    
    str.lengthvector = cumsum(str.steplength)-str.steplength/2;
    
    % area
    steparea = [];
    parfor i=1:lng, steparea=[steparea;polyarea([mesh(i:i+1,1);mesh(i+1:-1:i,3)],[mesh(i:i+1,2);mesh(i+1:-1:i,4)])]; end
    str.area = sum(steparea);
    str.steparea = steparea;
    
    % volume
    d = edist(mesh(:,1),mesh(:,2),mesh(:,3),mesh(:,4));
    str.stepvolume = (d(1:end-1).*d(2:end) + (d(1:end-1)-d(2:end)).^2/3).*str.steplength*pi/4;
    str.volume = sum(str.stepvolume);
elseif isfield(str,'contour') && length(str.contour)>1
    contour = str.contour;
    lng = size(contour,1);
    str.length = sqrt(max(max( (repmat(contour(:,1),1,lng)-repmat(contour(:,1)',lng,1)).^2 + ...
        (repmat(contour(:,2),1,lng)-repmat(contour(:,2)',lng,1)).^2)));
    str.area = polyarea(contour(:,1),contour(:,2));
end
end

function d=edist(x1,y1,x2,y2)
% complementary for "getextradata", computes the length between 2 points
d=sqrt((x2-x1).^2+(y2-y1).^2);
end