"""
This is a file to correct the image locations of the files created by
Weave.jl. Weave.jl makes the image location to be a directory above where the
file gets placed. However, the images are actually located in 'assets/img/'.
At the end of the 'make.jl' script, this script will get called and we will
scrub through the files, find where the images are decalared and fix the
call. The calls originally look something like:
    ![](../assets/img/other-stuff
We will fix this to be:
    ![]({{site.baseurl}}/assets/img/other-stuff
"""

import os
from os import listdir
from os.path import isfile, join
import sys

# Read in the files names from the commandline args and strip off the '.jmd'
# ext and replace with '.md'
file_names = [str(file_name) for file_name in sys.argv[1:]]
for i, file_name in enumerate(file_names):
   file_name, _ = os.path.splitext(file_name)
   file_names[i] = file_name + '.md'

# Go through the '_post' directory and correct the figure calls. The write the
# new files.
for file_name in file_names:
    lines = []
    with open(join("_posts", file_name)) as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if "![](../assets/img/" in line:
                front, back = line.split("..")
                lines[i] = front + "{{site.baseurl}}" + back
    with open(join("_posts", file_name), 'w') as f:
        f.writelines(lines)
