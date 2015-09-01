#!/bin/sh

# $Id: mkdocs.sh,v 1.3 2005-08-22 16:09:05 zbigniew Exp $
# Script to build the documentation

# Allow this script to be run from its own directory or its parent directory.

if [ `basename $PWD` = 'docs' ]
then
    docdir=`pwd -P`
    topdir=`echo $docdir | sed -e 's|/docs$||'`
else
    topdir=`pwd -P` 
    docdir=$topdir/docs
fi

# Run doxygen for QIO 

qiodir=$topdir/other_libs/qio/doc
cd $qiodir
doxygen qiodoc
cd -

# The version number
# I assume that the tag is of the form <something>1-2-3<something>

version=`echo '$Name:  $' | sed -e 's/[^1234567890-]//g' -e 'y/-/./'`

tempdoxcfg=/tmp/dox$$

for here in usr ref 
do
	sed -e "s/^PROJECT_NUMBER.*/PROJECT_NUMBER = $version/" $docdir/$here/doxygen.cfg > $tempdoxcfg
	mv $tempdoxcfg $docdir/$here/doxygen.cfg
done

# Do any other doxygen config file editing here.


# Run doxygen  

for here in ref usr 
do
	cd $docdir/$here
	doxygen doxygen.cfg		
done



