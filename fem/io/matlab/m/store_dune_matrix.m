function A=store_dune_matrix(fid,value)
szvec = size(value);
sz1 = szvec(1);
sz2 = szvec(2);
store_dune_int(fid,sz1);
store_dune_int(fid,sz2);
for j=1:sz1
for i=1:sz2
  store_dune_double(fid,value(j,i));
end;
end;
A = 0;
