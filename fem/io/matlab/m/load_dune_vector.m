function A = load_dune_vector(fid)
sz = load_dune_int(fid);
A=zeros(sz,1);
for i=1:sz
  A(i,1)=load_dune_double(fid);
end;