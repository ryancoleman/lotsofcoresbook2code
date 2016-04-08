"""Functions for in-place manipulation of bundletrajectories.

This module defines a number of functions that can be used to 
extract and delete data from BundleTrajectories directly on 
disk.  The functions are intended for large-scale MD output, 
so they avoid copying the potentially large amounts of data.
In stead, data is either directly deleted in-place; or copies
are made by creating a new directory structure, but hardlinking
the data files.  Hard links makes it possible to delete the
original data without invalidating the copy.
"""

from ase.io.bundletrajectory import BundleTrajectory
import os
import cPickle as pickle

def copy_frames(inbundle, outbundle, start=0, end=None, step=1,
                verbose=False):
    """Copies selected frame from one bundle to the next."""
    if not (isinstance(start, int) and
            (isinstance(end, int) or end is None) and
            isinstance(step, int)):
        raise TypeError("copy_frames: start, end and step must be integers.")
    metadata, nframes = read_bundle_info(inbundle)
    if start < 0:
        start += nframes
    if end is None:
        end = nframes
    if end < 0:
        end += nframes
    if start < 0 or (start > nframes-1 and end > 0):
        raise ValueError("copy_frames: Invalid start value.")
    if end < 0 or (end > nframes-1 and end < 0):
        raise ValueError("copy_frames: Invalid end value.")
    if step == 0:
        raise ValueError("copy_frames: Invalid step value (zero)")
    frames = range(start, end, step)
    if verbose:
        print "Copying the frames", frames
    
    # Make the new bundle directory
    os.mkdir(outbundle)
    f = open(os.path.join(outbundle, 'metadata'), 'w')
    pickle.dump(metadata, f, -1)
    f.close()
    
    for nout, nin in enumerate(frames):
        if verbose:
            print "F%i -> F%i" % (nin, nout)
        indir = os.path.join(inbundle, "F"+str(nin))
        outdir = os.path.join(outbundle, "F"+str(nout))
        os.mkdir(outdir)
        names = os.listdir(indir)
        for name in names:
            fromfile = os.path.join(indir, name)
            tofile = os.path.join(outdir, name)
            os.link(fromfile, tofile)
        if nout == 0 and nin != 0:
            if verbose:
                print "F0 -> F0 (supplemental)"
            # Data for first frame must be supplemented with
            # data from the first frame of the source bundle.
            firstnames = os.listdir(os.path.join(inbundle, "F0"))
            n_from_first = 0
            for name in firstnames:
                if name not in names:
                    if verbose:
                        print "   ", name
                    fromfile = os.path.join(inbundle, "F0", name)
                    tofile = os.path.join(outdir, name)
                    os.link(fromfile, tofile)
                    n_from_first += 1
            # Also, the smalldata.pickle stuff must be updated.
            # At the same time, check that the number of fragments
            # has not changed, if the data is written in a fragmented
            # way AND it looks like we got such data from F0
            assert metadata['backend'] == "pickle"
            f = open(os.path.join(inbundle, "F0", "smalldata.pickle"))
            data0 = pickle.load(f)
            f = open(os.path.join(indir, "smalldata.pickle"))
            data1 = pickle.load(f)
            if metadata['subtype'] == 'split' and n_from_first >= data0['fragments']:
                if data0['fragments'] != data1['fragments']:
                    raise RuntimeError("Cannot combine data from F0 and F%i since the number of fragments has changed"
                                       % (nin,))
            data0.update(data1)  # Data in frame overrides data from frame 0.
            smallname = os.path.join(outdir, "smalldata.pickle")
            os.unlink(smallname) 
            f = open(smallname, "w")
            pickle.dump(data0, f, -1)
            f.close()
    # Finally, write the number of frames
    f = open(os.path.join(outbundle, 'frames'), 'w')
    f.write(str(len(frames))+'\n')
    f.close()
    
# Helper functions

def read_bundle_info(name):
    """Read global info about a bundle.
    
    Returns (metadata, nframes)
    """
    if not os.path.isdir(name):
        raise IOError("No directory (bundle) named '%' found." % (name,))
    metaname = os.path.join(name, 'metadata')
    if not os.path.isfile(metaname):
        raise IOError("'%s' does not appear to be a BundleTrajectory (no %s)" 
                      % (name, metaname))
    f = open(metaname)
    mdata = pickle.load(f)
    f.close()
    if 'format' not in mdata or mdata['format'] != 'BundleTrajectory':
        raise IOError("'%s' does not appear to be a BundleTrajectory" % (name,))
    if mdata['version'] != 1:
        raise IOError("Cannot manipulate BundleTrajectories with version number %s"
                      % (mdata['version'],))
    f = open(os.path.join(name, "frames"))
    nframes = int(f.read())
    if nframes == 0:
        raise IOError("'%s' is an empty BundleTrajectory" % (name,))
    return mdata, nframes



if __name__ == '__main__':
    import sys
    inname, outname = sys.argv[1:3]
    if len(sys.argv) > 3:
        start = int(sys.argv[3])
    else:
        start = 0
    if len(sys.argv) > 4:
        end = int(sys.argv[4])
    else:
        end = -1
    if len(sys.argv) > 5:
        step = int(sys.argv[5])
    else:
        step = 1
    copy_frames(inname, outname, start, end, step, verbose=1)
    
