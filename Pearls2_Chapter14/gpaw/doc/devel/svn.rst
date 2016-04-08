.. _svn:

===
SVN
===

Browse svn online_.

.. _online: http://svn.fysik.dtu.dk/projects/gpaw/


Take a look at this `SVN cheat sheet`_

.. _SVN cheat sheet: ../_static/svn-refcard.pdf



Working with branches
=====================

Creating a new branch::

  $ svn copy https://svn.fysik.dtu.dk/projects/gpaw/trunk https://svn.fysik.dtu.dk/projects/gpaw/branches/mixing -m "Experimental density mixing branch"

Merge changes from trunk into branch::

  $ svn merge https://svn.fysik.dtu.dk/projects/gpaw/trunk
  $ svn ci -m "Merged changes from trunk into branch."

Merge branch to trunk::

  $ cd <root directory of trunk>
  $ svn up
  At revision 957.
  $ svn merge --reintegrate https://svn.fysik.dtu.dk/projects/gpaw/branches/new-interface
  $ svn ci -m "Merged branch to trunk."


Reverting a bad commit
======================

Go back to revision 2748::

  $ svn merge -r BASE:2748 .

And commit the change::

  $ svn ci -m "Reverting repository to r2748."
