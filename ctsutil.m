function out = ctsutil(op,in,opt)
%out = ctsutil(op,in,opt)
%'trim' in=vol
%'pad'  in=vol opt=padding, all dims
%'bin'  in=vol opt=binning
switch op
    case 'trim'
        out = in(any(in ~= 0,[2 3]),any(in ~= 0,[1 3]),any(in ~= 0,[1 2]));
    case 'pad'
        out = padarray(in,[opt,opt,opt]);
        
end