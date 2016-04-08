#! /bin/sh

# Test file splitting and recombining

testfile=binary_host_test

#----------------------------------------------------------------------
echo "Testing qio-convert-mesh-singlefs"
#----------------------------------------------------------------------

# Split this file according to the machine dimensions 4 2 1 1
echo "4 4 2 1 1" | qio-convert-mesh-singlefs 0 $testfile
if [ $? -ne 0 ] 
then
  echo "File splitting FAILED"
  exit 1
fi

# Move the file out of the way and recombine the parts
mv $testfile $testfile.bak
echo "4 4 2 1 1" | qio-convert-mesh-singlefs 2 $testfile
if [ $? -ne 0 ] 
then
  echo "File recombination FAILED"
  exit 1
fi

# The resulting file should be identical with the original
diff $testfile $testfile.bak
if [ $? -ne 0 ] 
then
  echo "Test FAILED.  Recombined file differs from original."
  exit 1
fi

echo "PASSED qio-convert-mesh-singlefs"
/bin/rm $testfile.bak

#----------------------------------------------------------------------
echo "Testing qio-convert-mesh-pfs"
#----------------------------------------------------------------------

# Split this file according to the machine dimensions 4 2 2 1
qio-convert-mesh-pfs 0 $testfile < layout_test
if [ $? -ne 0 ] 
then
  echo "File splitting FAILED"
  exit 1
fi

# Move the file out of the way and recombine the parts
mv $testfile $testfile.bak
qio-convert-mesh-pfs 2 $testfile < layout_test
if [ $? -ne 0 ] 
then
  echo "File recombination FAILED"
  exit 1
fi

# The resulting file should be identical with the original
diff $testfile $testfile.bak
if [ $? -ne 0 ] 
then
  echo "Test FAILED.  Recombined file differs from original."
  exit 1
fi

echo "PASSED qio-convert-mesh-pfs"
/bin/rm $testfile.bak
/bin/rm -r path??

#----------------------------------------------------------------------
echo "Testing qio-convert-mesh-ppfs"
#----------------------------------------------------------------------

# Split this file according to the machine dimensions 4 2 2 1
# and I/O machine dimensions 2 2 1 1
qio-convert-mesh-pfs 0 $testfile < layout_test_ppfs
if [ $? -ne 0 ] 
then
  echo "File splitting FAILED"
  exit 1
fi

# Move the file out of the way and recombine the parts
mv $testfile $testfile.bak
qio-convert-mesh-pfs 2 $testfile < layout_test_ppfs
if [ $? -ne 0 ] 
then
  echo "File recombination FAILED"
  exit 1
fi

# The resulting file should be identical with the original
diff $testfile $testfile.bak
if [ $? -ne 0 ] 
then
  echo "Test FAILED.  Recombined file differs from original."
  exit 1
fi

echo "PASSED qio-convert-mesh-pfs"
/bin/rm $testfile.bak
/bin/rm -r path??

