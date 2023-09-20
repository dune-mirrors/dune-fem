import subprocess, datetime, pickle, sys

count = datetime.date.today().toordinal()
res = {}

for np in range(1,5):
    print("==========================")
    print(f"running on {np} procs")
    print("==========================")
    with subprocess.Popen(["mpiexec", "-np",str(np),"python","parallelTest.py",str(count)],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE) as mpi2:
        stdout, stderr = mpi2.communicate()
        res[np] = {"exit":mpi2.returncode}
        if res[np]["exit"] > 0:
            print("---------------------")
            print(stdout.decode('utf-8'))
            print("------")
            print(stderr.decode('utf-8'))
            print("---------------------")

    if res[np]["exit"] == 0:
        with open(f"parallelTest{np}.out","rb") as f:
            failures,iterations = pickle.load(f)
        res[np]["failures"] = failures
        res[np]["iters"] = iterations

        if np > 1: # compare iterations with np-1 (if available)
            for cmpNp in range(np-1,0,-1):
                if res[cmpNp]["exit"] == 0: # we'll compare with this one
                    break
            if res[cmpNp]["exit"] == 0: # none found to compare with - probably everything has broken down
                for rnp, r1 in zip( res[np]["iters"],res[cmpNp]["iters"] ):
                    assert rnp[0] == r1[0] # make sure we are cmp the same simulations
                    # THIS NEEDS LOOKING IN TO
                    if (rnp[1]-r1[1]) > 0.3*r1[1]:  # this is ridiculously high tolerance but lagrange with numpy/istl needs this
                        print(np,f"high iteration count for {rnp[0]} with {rnp[1]}>{r1[1]}")

print([ r["exit"] for r in res.values() ])
sys.exit( sum([ r["exit"] for r in res.values() ]) )
