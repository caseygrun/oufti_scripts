function [bactmask] = make_cell_mask(meshcoord, img)
    % convert an oufti cell mesh into a binary mask of img, with pixels
    % involving the cell set to 1 and other pixels to 0
    
    % generate a mask identifying this cell, from its mesh
    borders=[flip(meshcoord(:,4)) flip(meshcoord(:,3)); meshcoord(:,2) meshcoord(:,1)];

    % remove points containing +/-infinity
    infrows = find(ismember(abs(borders),[Inf Inf],'rows')==1);

    if isempty(infrows)==0
        borders(infrows,:)=[];
    end

    % create mask identifing the position of this cell
    bactmask=poly2mask(...
        double(borders(:,2)),...
        double(borders(:,1)),...
        size(img,1),...
        size(img,2));
end