function A = load_dune_intmatrix(fid)
rows = load_dune_int(fid);
cols = load_dune_int(fid);
A=zeros(rows,cols);
for i=1:rows
  for j=1:cols
    A(i,j)=load_dune_int(fid);
  end;
end;