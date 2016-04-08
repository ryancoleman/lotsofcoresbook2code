.. _newrelease:

===========
New release
===========

When it is time for a new release of the code, here is what you have to do:

* **Warning:** use only three digits release numbers, e.g. *0.7.2*,

* Checkout the :ref:`latest_development_release`,

* then :ref:`running_tests`.

* If a new ase release is required to pass the tests
  modify ``required_ase_version`` and ``required_ase_svnversion``
  in :trac:`gpaw/version.py`, and checkin the changes.

* ``svn up`` and :ref:`running_tests` again.

* Make a tag in svn, using the current version number
  (to make sure **not** to include changes done by other developers
  in the meantime!)::

    svn copy -r 6972 https://svn.fysik.dtu.dk/projects/gpaw/trunk https://svn.fysik.dtu.dk/projects/gpaw/tags/0.7.2 -m "Version 0.7.2"

  **Note** the resulting tag's revision ``tags_revision``.

* **Checkout** the source, specifying the version number in the directory name::

   svn co -r tags_revision https://svn.fysik.dtu.dk/projects/gpaw/tags/0.7.2 gpaw-0.7.2

* Create the tar file::

   cd gpaw-0.7.2
   rm -f MANIFEST gpaw/svnversion.py*
   python setup.py sdist

  Note that the ``tags_revision`` is put into the name of the
  tar file automatically. Make sure that you are getting only
  ``tags_revision`` in the tar file name! Any changes to the source
  will be reflected as a mixed or modified revision tag!

* Put the tar file on web2 (set it read-able for all)::

   scp dist/gpaw-0.7.2."tags_revision".tar.gz root@web2:/var/www/wiki/gpaw-files

* Add a link to the new GPAW on :ref:`news` and update the information
  on the :ref:`download` page.

* Update the :ref:`releasenotes` including the compatible ASE release!.

* Increase the version number in gpaw/version.py, and commit the change::

    cd ~/gpaw
    svn ci -m "Version 0.8.0"

  Now the trunk is ready for work on the new version.

* Send announcement email to the ``gpaw-developers`` mailing list (see :ref:`mailing_lists`).
