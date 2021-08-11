#!/usr/bin/python
# Generate "wrapper sources" for compiling the same file several times with CMake for UMFPACK
# (original author: Jose Luis Blanco, see
# https://github.com/jlblancoc/suitesparse-metis-for-windows/blob/master/SuiteSparse/makefile2wrappers.py)
# modified by Benjamin DONNOT


# In the original "make" commands, some SuiteSparse source files are compiled with different option (such as -DLONG or
# -DCS_COMPLEX etc.) to give different binaries. Following what has been done by Jose Luis Blanco in the
# above mentioned project, we here generate different files, with each a given name, that will:
# 1) define in the file, with the macro #define XXX
# 2) include the original source file with "include <SOURCE_FILE.cc>
# so that it's easy to achieve the same result (as in "make") with cmake.

# the main modifications to the original script are:
# 2) the regex to match the files are different (this is because I generated the "makefile2wrappers.txt" by copy pasting
#    the output of "make" on a linux machine
# 3) this files applies to all subdirectories (eg. AMD, BTF, COLAMD, etc.) of the SuiteSparse,
#    instead of having to modify each one individually.
# 4) PEP-8 convention, especially for this source is "more" respected.

import os
import sys
import re
import shutil


# NB it is super important that all the lines of the `makefile2wrappers.txt` files
# have something like "WHATEVER -DDLONG -c WHATEVER", ie that the "-DDLONG" preceed the "-c"
# otherwise it will NOT WORK


#------   MAIN   -------
def main_one_dir(input_dir=".", output_dir="SourceWrappers", suitsparse_source_dir="../SuiteSparse"):
    with open(os.path.join(input_dir, "makefile2wrappers.txt"), "r") as f:
        lins = f.readlines()
    input_files_path = os.path.abspath(input_dir)
    input_files_path = os.path.split(input_files_path)[0]
    input_files_path = os.path.join(input_files_path, "Source")  # supposes that the original c files are located into
    # a directory named "Source" for each of the SuiteSparse subdirectory (eg. AMD, BTF, etc.)

    # first create the symlink to original suitesparse code
    if not os.path.exists(suitsparse_source_dir) or not os.path.isdir(suitsparse_source_dir):
        raise RuntimeError(f"\"{suitsparse_source_dir}\" should be a suitesparse directory.")
    for dirn in ["Source", "Include", "Lib"]:
        os.symlink(os.path.join(suitsparse_source_dir, dirn),
                   os.path.join(input_dir, dirn))

    # now make the output directory
    if os.path.exists(output_dir):
        # remove the content of the directory if it exists
        shutil.rmtree(output_dir)

    # and now fill this target directory with the modified source files (see description above)
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
                    o.write('#ifndef NCOMPLEX\n')
                    o.write('#define CS_COMPLEX\n')
                    src_wrt = re.sub("../Source/", "", src)
                    o.write(f'#include <{src}>\n')
                    o.write('#endif // NCOMPLEX\n')
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
                    o.write(f'#define {d}\n')

                # The source code itself:
                src_wrt = re.sub("../Source/", "", src)
                o.write(f'#include <{src}>\n')
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
                    o.write('#ifndef NCOMPLEX\n')
                    o.write('#define CS_COMPLEX\n')
                    src_wrt = re.sub("../Source/", "", src)
                    o.write(f'#include <{src}>\n')
                    o.write('#endif // NCOMPLEX\n')
            else:
                with open(f, "w") as o:
                    src_wrt = re.sub("../Source/", "", src)
                    o.write(f'#include <{src}>\n')
    return 0


def main():
    this_path = os.path.abspath("./SuiteSparse")
    suitsparse_source_code = os.path.abspath(".")
    suitsparse_source_code = os.path.split(suitsparse_source_code)[0]
    suitsparse_source_code = os.path.join(suitsparse_source_code, "SuiteSparse")

    # handle the "SuiteSparse_config" directory, where i need to link 2 files
    dirnm = "SuiteSparse_config"
    curr_dir = os.path.join(os.path.join(this_path, dirnm))
    suitsparse_source_dir = os.path.join(suitsparse_source_code, dirnm)
    if not os.path.exists(suitsparse_source_dir) or not os.path.isdir(suitsparse_source_dir):
        raise RuntimeError(f"\"{suitsparse_source_dir}\" should be a suitesparse directory.")
    for fn in ["SuiteSparse_config.c", "SuiteSparse_config.h"]:
        src_file = os.path.join(suitsparse_source_dir, fn)
        if not os.path.exists(src_file) or not\
            os.path.isfile(src_file):
            raise RuntimeError(f"Impossible to find file \"{src_file}\".")
        shutil.copy(src_file, os.path.join(curr_dir, fn))

    # handle the files for all libs
    for dirnm in ("AMD", "BTF", "COLAMD", "CXSparse", "KLU"):
        curr_dir = os.path.join(os.path.join(this_path, dirnm))
        sp_source = os.path.join(suitsparse_source_code, dirnm)
        main_one_dir(curr_dir,
                     os.path.join(curr_dir, "SourceWrappers"),
                     sp_source
                     )


if __name__ == "__main__":
    sys.exit(main())

