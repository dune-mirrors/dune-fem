function A=store_dune_vector(fid,value)
szvec = size(value);
sz = szvec(1);
store_dune_int(fid,sz);
for i=1:sz
  store_dune_double(fid,value(i,1));
end;
A = 0;