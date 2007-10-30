function A = load_dune_grid(fid)
p = load_dune_matrix(fid);
t = load_dune_intmatrix(fid);
disp(p);
disp(t);
A = triagrid( p, t, [] );
