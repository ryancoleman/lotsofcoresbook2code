import os
import shutil
import sys
rootdir = os.path.realpath('../libxc.old')
for root, subFolders, files in os.walk(rootdir):
    for dir in subFolders:
       if dir == '.svn':
           source = os.path.join(root, dir)
           dest = os.path.realpath(source[source.find('libxc.old')+len('libxc.old')+1:])
           shutil.copytree(source, dest)
           print(source, dest)
