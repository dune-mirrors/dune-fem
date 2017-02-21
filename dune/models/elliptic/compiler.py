from argparse import ArgumentParser

from dune.source.cplusplus import NameSpace, SourceWriter

from .formfiles import compileUFLFile

def main():
    parser = ArgumentParser(description='dune-fem-fluid model compiler')
    parser.add_argument('input', help='name of input .ufl file')
    parser.add_argument('-o', '--output', help='name of output .hh file')
    parser.add_argument('-n', '--namespace', help='C++ namespace for models in')
    parser.add_argument('--no-tempvars', dest='tempVars', action='store_false', help='do not generate temporary variables')
    parser.add_argument('-v', '--version', action='version', version='dune-fem-fluid 2.6')
    try:
        args = parser.parse_args()

        code = compileUFLFile(args.input, tempVars=args.tempVars)
        if args.namespace is not None:
            for namespace in reversed(args.namespace.split('::')):
                code = NameSpace(namespace, code=code)
        if args.output is None:
            if args.input[-4:] != ".ufl":
                args.output = args.input + ".hh"
            else:
                args.output = args.input[:-4] + ".hh"
        SourceWriter(args.output).emit(code)

        return 0
    except Exception as e:
        print(e)
        return 1
