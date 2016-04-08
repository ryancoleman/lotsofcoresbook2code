.. _license_info:

=======
License
=======

.. contents::

Human-readable version
======================

ASE is released under the `GNU Lesser General Public License`_ (LGPL).
In short, this means

* **NO WARRANTY:** We provide the code free of charge. You cannot sue
  us if it does not work, gives wrong results, or even if it damages
  your computer and/or scientific reputation.

* You may use the code for whatever you like.  You may study and modify it.

* You may distribute distribute modified or unmodified versions of ASE
  as long as you do it under the LGPL_ (or GPL_) licence. You may
  distribute unmodified versions of ASE together with software under
  other licences (even commercial) as long as ASE itself is clearly
  identified as being under the LGPL_, but if you modify ASE for such
  purposes you are *required* to make the modifications available under
  the LGPL_.

Note that we appreciate that you send modifications, bug fixes and
improvements back to us instead of just distributing them, but the
license has no such requirements.

You can read more about the `LGPL on Wikipedia`_.


Legal version of the license
============================

Please read the full text of the `GNU Lesser General Public License`_
(as published by the Free Software Foundation).


What happens when ASE Calculators are under another license?
============================================================

We are sometimes asked if it is problematic to use ASE together with
calculators under other licenses, for example GPL_. It is clear that a
program under the GPL_ can use a library under the LGPL_, whereas a
program under the LGPL_ cannot be derived from (and link) a library
under the GPL_. Does this cause a problem if someone uses ASE with a
calculator such as GPAW_ licensed under the GPL_? We do not think so,
for the following reasons:

1. The LGPL_ and GPL_ do not limit how you *use* the codes, only how
   you *distribute* them.

2. ASE does not require any specific calculator to function, but many
   calculators require ASE to function, supporting the interpretation
   that ASE is a library for the calculator.

3. Although ASE includes a few cases where it imports calculators such
   as GPAW_ and Asap_, these can be regarded as "hooks" helping ASE to
   support these calculators, ASE does not depend on these calculators
   for its functionality.

4. The LGPL_ / GPL_ concept of "derived work" relies on the concept of
   "linking" which only makes sense in compiled languages. It is
   generally agreed that it is unproblematic when an interpreted
   language uses different modules under different licenses. See e.g.
   this `statement by Fedora`_: *Mere use of independent modules in a
   true interpreted language environment (like Perl or Python) is not a
   situation where Fedora is generally concerned about license
   compatibility, as long as those multiply licensed modules are not
   compiled together into a single binary and there is no code copying
   between the two.*

5. The actual executable doing the linkage is not ASE, but Python.
   However, nobody doubts that it is OK for Python (which has a very
   permissible license) to load modules licensed under the GPL_. Probably
   because of point 1 above.

6. Point 5 is not valid when running parallel GPAW_ or Asap_
   calculations. In these cases GPAW_ and Asap_ provide specially built
   Python executables with the GPAW_ or Asap_ code built-in, i.e. derived
   work based on Python but licensed under the GPL_ (or LGPL_ for Asap).
   In these cases it is absolutely clear that it is GPAW_ or Asap_
   loading ASE, not the other way around; so there are no problems.


.. _`GNU Lesser General Public License`: http://www.gnu.org/licenses/lgpl.html
.. _LGPL: http://www.gnu.org/licenses/lgpl.html
.. _GPL: http://www.gnu.org/licenses/gpl.html
.. _`LGPL on Wikipedia`: http://en.wikipedia.org/wiki/GNU_Lesser_General_Public_License
.. _GPAW: https://wiki.fysik.dtu.dk/gpaw
.. _Asap: https://wiki.fysik.dtu.dk/asap
.. _`statement by Fedora`: https://fedoraproject.org/wiki/Licensing:FAQ?rd=Licensing/FAQ#Linking_and_multiple_licenses

