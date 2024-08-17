function H = gen_diffuse_white(mixoutRoot, opts, params)

mkdir(mixoutRoot)
%z = cinf_1D(params.d,params.L,params); 
%z = cinf_3D(params.P,params.L,params); 
%z = sinf_1D(params.d,params.L,params);
z = sinf_3D(params.P,params.L,params);
z = z./max(max(abs(z)));
audiowrite([mixoutRoot, '/x_', num2str(opts.Nd), 'x', num2str(opts.nch), '.wav'], z.', opts.fs);

end

