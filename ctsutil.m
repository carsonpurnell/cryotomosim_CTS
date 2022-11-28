function out = ctsutil(op,in,opt)
switch op
    case 'trim'
        out = in(any(in ~= 0,[2 3]),any(in ~= 0,[1 3]),any(in ~= 0,[1 2]));
    case 'pad'
        out = padarray(in,[opt,opt,opt]);

end