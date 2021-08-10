#!/usr/bin/python
# Generate "wrapper sources" for compiling the same file several times with CMake for UMFPACK
# (original author: Jose Luis Blanco, see
# https://github.com/jlblancoc/suitesparse-metis-for-windows/blob/master/SuiteSparse/makefile2wrappers.py)
# modified by Benjamin DONNOT

import os
import sys
import re
import shutil


# NB it is super important that all the lines of the `makefile2wrappers.txt` files
# have something like "WHATEVER -DDLONG -c WHATEVER", ie that the "-DDLONG" preceed the "-c"
# otherwise it will NOT WORK

#------   MAIN   -------
def main_one_dir(input_dir=".", output_dir="SourceWrappers"):
    with open(os.path.join(input_dir, "makefile2wrappers.txt"), "r") as f:
        lins = f.readlines()

    if os.path.exists(output_dir):
        # remove the content of the directory if it exists
        shutil.rmtree(output_dir)

    # create the output directory if it does not exists
    os.mkdir(output_dir)

    for l in lins:
        l = l.strip()
        if len(l) == 0 or l.startswith('#'):
            continue

        # print('Line: '+l)

        # orig regex
        # defs = re.match(".*\)(.*)-c", l).group(1).strip()  # to retrieve stuff like "-DDINT"
        # If there's no "-o" flag, just compile the file as is:
        if re.search('.*-o.*', l) is not None:
            src = re.match(".*-c(.*)-o", l).group(1).strip()
            out = re.match(".*-o(.*)", l).group(1).strip()
            f = os.path.join(output_dir, out+".c")

            if re.search("cs_convert\..*", src) is not None:
                # special case: this file require a bit more work...
                # courtesy to https://github.com/jlblancoc/suitesparse-metis-for-windows/blob/master/SuiteSparse/CXSparse/SourceWrappers/cs_convert.o.c
                with open(f, "w") as o:
                    o.write('#ifndef NCOMPLEX'+'\n')
                    o.write('#define CS_COMPLEX'+'\n')
                    o.write('#include <'+src+'>'+'\n')
                    o.write('#endif // NCOMPLEX'+'\n')
                continue

            # modified regex:
            my_re = re.search("\s(-D.*\s+)+-c\s", l)  # to retrieve stuff like "-DDINT" if any
            if my_re is not None:
                defs = my_re.group().strip()
                defs = re.sub("\s+-c\s*", "", defs)  # remove the "-c"
            else:
                defs = ""
            # print(' => Creating '+f+'\n')
            with open(f, "w") as o:
                DEFs = defs.strip().split("-D")
                DEFs = [x for x in DEFs if x]  # Remove empty
                # If we have "CS_COMPLEX", wrap the entire thing in a #ifndef NCOMPLEX ... #endif
                has_complex = False
                if 'CS_COMPLEX' in DEFs:
                    has_complex = True
                # TODO: Other COMPLEX flags?

                if has_complex:
                    o.write('#ifndef NCOMPLEX\n')

                # Pre-definitions:
                for d in DEFs:
                    o.write('#define '+d+'\n')

                # The source code itself:
                o.write('#include <'+src+'>'+'\n')
                if has_complex:
                    o.write('#endif // NCOMPLEX\n')
        else:
            src = re.match(".*-c(.*)", l).group(1).strip()
            f = os.path.join(output_dir, os.path.basename(src))
            # print(' => Creating ' + f + '\n')
            if re.search("cs_convert.o.c", src) is not None:
                # special case: this file require a bit more work...
                # courtesy to https://github.com/jlblancoc/suitesparse-metis-for-windows/blob/master/SuiteSparse/CXSparse/SourceWrappers/cs_convert.o.c
                with open(f, "w") as o:
                    o.write('#ifndef NCOMPLEX'+'\n')
                    o.write('#define CS_COMPLEX'+'\n')
                    o.write('#include <'+src+'>'+'\n')
                    o.write('#endif // NCOMPLEX'+'\n')
            else:
                with open(f, "w") as o:
                    o.write('#include <'+src+'>'+'\n')
    return 0


def main():
    this_path = os.path.abspath("./SuiteSparse")
    for dirnm in ("AMD", "BTF", "COLAMD", "CXSparse", "KLU"):
        curr_dir = os.path.join(os.path.join(this_path, dirnm))
        main_one_dir(curr_dir,
                     os.path.join(curr_dir, "SourceWrappers")
                     )


if __name__ == "__main__":
    sys.exit(main())

