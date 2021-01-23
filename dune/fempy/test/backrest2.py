#!/usr/bin/env python3

import pickle, numpy

class CP:
    def __init__(self):
        self.items = {}
    def add(self,**kwargs):
        self.items.update(kwargs)

def run(cp=None):
    a = 10
    name = "Test"
    if cp is None:
        a += 10
        # note can't change 'name'
        print("backup")
        cp = CP()
        cp.add( a=a, name=name )
    else:
        print("restore")
        run.__globals__.update(cp.items)
    print(a,name)
    return cp

if __name__ == "__main__":
    cp = run(None)
    run(cp)
