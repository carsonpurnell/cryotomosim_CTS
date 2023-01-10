function [out,out2] = ctsutil(op,in,opt)
%out = ctsutil(op,in,opt)
%various utility functions that are used in multiple CTS functions
%'trim' in=vol
%'pad'  in=vol opt=padding, all dims

%'bin'  in=vol opt=binning
%edge blanking in here?
arguments
    op
    in
    opt = 0
end
out2=0;
switch op
    case 'trim' %& iscell(in)
        if iscell(in)
            n=numel(in);
            out = cell(n,1);
            sumvol = sum( cat(4,in{:}) ,4);
            trimr = any(sumvol,[2 3]); trimc = any(sumvol,[1 3]); triml = any(sumvol,[1 2]);
            %cuts out intervening empty space, can screw up individual trim operations potentially
            for i=1:n
                out{i} = in{i}(trimr,trimc,triml);
            end
        else
            %out = in(any(in ~= 0,[2 3]),any(in ~= 0,[1 3]),any(in ~= 0,[1 2]));
            trimr = any(in,[2 3]); trimc = any(in,[1 3]); triml = any(in,[1 2]); 
            out = in(trimr,trimc,triml);
        end
    case 'pad'
        out = padarray(in,[opt,opt,opt]);
    case 'centervol'
        cen = regionprops3(in{1}>0,in{1},'WeightedCentroid').WeightedCentroid; %find center of gravity
        cen = cen([2,1,3]); %fix xy inversion
        difCOM = (cen-size(in{1})/2-0.25)*2;
        out = in;
        for i=1:numel(in)
            out{i} = padarray(out{i},max(round(difCOM*1),0),'post');
            out{i} = padarray(out{i},max(round(difCOM*-1),0),'pre');
        end
    case 'edgeblank'
        in([1:opt,end-opt:end],:,:) = 0;
        in(:,[1:opt,end-opt:end],:) = 0;
        in(:,:,[1:opt,end-opt:end]) = 0;
        out = in;
        
    case 'findloc'
        out = randi(numel(in)); %random start
        while in(out)~=1
            out = randi(numel(in));
        end
        [r,c,l] = ind2sub(size(in),out);
        out2 = out; out = [r,c,l];
        
end