.. _devel:

===========
Development
===========

GPAW development can be done by anyone! Just take a look at the
:ref:`todolist` list and find something that suits your talents!

The primary source of information is still the :ref:`manual` and
:ref:`documentation`, but as a developer you might need additional
information which can be found here. For example the :ref:`code_overview`.

As a developer, you should subscribe to all GPAW related :ref:`mailing_lists`.
We would also like to encourage you to join our channel for :ref:`irc`.

Now you are ready to to perfom a :ref:`developer_installation` and
start development!

.. toctree::
   :maxdepth: 2

   projects/projects
   developer_installation
   workshop

.. note --- below toctrees are defined in separate files to make sure that the line spacing doesn't get very large (which is of course a bad hack)

Development topics
==================

When committing significant changes to the code, remember to add a
note in the :ref:`releasenotes` at the top (current svn) - the version
to become the next release.

.. toctree::
   :maxdepth: 1

   codingstandard
   c_extension
   writing_documentation
   formulas
   debugging
   profiling
   testing
   ase_optimize/ase_optimize
   bugs
   newrelease
   svn
   technology
   todo

* The latest report_ from PyLint_ about the GPAW coding standard.

.. spacer

* Details about supported :ref:`platforms_and_architectures`.

.. _report: http://dcwww.camd.dtu.dk/~s052580/pylint/gpaw/
.. _PyLint: http://www.logilab.org/857


.. _the_big_picture:
.. _code_overview:

Code Overview
=============

Keep this picture under your pillow:


.. image:: bigpicture.png
   :target: ../bigpicture.svg

The developer guide provides an overview of the PAW quantities and how
the corresponding objects are defined in the code:

.. toctree::
   :maxdepth: 2

   overview
   developersguide
   proposals/proposals
   paw
   symmetry
   wavefunctions
   setups
   density_and_hamiltonian
   planewaves
   communicators
   others


The GPAW logo
=============

The GPAW-logo is available in the odg_ and svg_ file formats:
gpaw-logo.odg_, gpaw-logo.svg_

.. _odg: http://www.openoffice.org/product/draw.html
.. _svg: http://en.wikipedia.org/wiki/Scalable_Vector_Graphics
.. _gpaw-logo.odg: ../_static/gpaw-logo.odg
.. _gpaw-logo.svg: ../_static/gpaw-logo.svg


Statistics
==========

The image below shows the development in the volume of the code as per
March 18 2014.

.. image:: stat.png

*Documentation* refers solely the contents of this homepage. Inline
documentation is included in the other line counts.


Contributing to GPAW
====================

Getting commit access to our SVN repository works the same way as for
the `ASE project`_.  Here is the list of current committers:

=========  =========================  ========================================
id         real name
=========  =========================  ========================================
anpet      Andrew Peterson            andrew_peterson:brown,edu
aoj        Ari Ojanperä               ari,ojanpera:tkk,fi
askhl      Ask Hjorth Larsen          asklarsen:gmail,com
dlandis    David Landis               dlandis:fysik,dtu,dk
dmowbray   Duncan Mowbray             duncan,mowbray:gmail,com
dulak      Marcin Dulak               dulak:fysik,dtu,dk
eojons     Elvar Örn Jónsson          elvar,jonsson:fysik,dtu,dk
falh       Falco Huser                falh:fysik,dtu,dk
fras       Filip Anselm Rasmussen     fras:fysik,dtu,dk
georg      Poul Georg Moses           gtritsaris:seas,harvard,edu
getri      George Tritsaris           getri:fysik,dtu,dk
hahansen   Heine Anton Hansen         hahansen:fysik,dtu,dk
haiping    Haiping Lin                H,Lin1:liverpool,ac,uk
hhk05      Henrik Kristoffersen       hhk05:inano,dk
hildur     Hildur Guõmundsdóttir      hildur,gud:gmail,com
jensj      Jens Jørgen Mortensen      jensj:fysik,dtu,dk
jesswe     Jess Wellendorff Pedersen  jesswe:fysik,dtu,dk
jingzhe    Jingzhe Chen               jingzhe:fysik,dtu,dk
jsm        Jess Stausholm-Møller      jsm:phys,au,dk
jussie     Jussi Enkovaara            jussi,enkovaara:csc,fi
juya       Jun Yan                    juya:fysik,dtu,dk
kelkkane   Andre Kelkkanen            kelkkane:fysik,dtu,dk
kelu       Keld Troen Lundgaard       kelu:fysik,dtu,dk
kiran      Kirsten Andersen           kiran:fysik,dtu,dk
kkaa       Kristen Kaasbjerg          kkaa:fysik,dtu,dk
ksaha      Kamal Saha                 ?
kuismam    Mikael Kuisma              mikael,kuisma:tut,fi
lara       Lara Ferrighi              laraf:phys,au,dk
linli      Lin Li                     lli2010:stanford,edu
lauri      Lauri Lethovaara           lauri,lehtovaara:iki,fi
lopeza     Olga Lopez                 lopez:cc,jyu,fi
marsalek   Ondrej Marsalek            ondrej,marsalek:gmail,com
mathiasl   Mathias Ljungberg          mathiasl:physto,se
miwalter   Michael Walter             Michael,Walter:fmf,uni-freiburg,de
mogje      Morten Gjerding            mortengjerding:gmail,com
mohpa      Mohnish Pandey             mohpa:fysik,dtu,dk
moses      Poul Georg Moses           poulgeorgmoses:gmail,com
mvanin     Marco Vanin                mvanin:fysik,dtu,dk
naromero   Nichols Romero             naromero:alcf,anl,gov
peterklue  Peter Kluepfel             peter:theochem,org
s032082    Christian Glinsvad         s032082:fysik,dtu,dk
s093017    Per Simmendefeldt Schmidt  s093017:student,dtu,dk
sakko      Arto Sakko                 arto,sakko:aalto,fi
schiotz    Jakob Schiotz              schiotz:fysik,dtu,dk
shakala    Samuli Hakala              samuli,hakala:tkk,fi
simola     Simone Latini              simola:fysik,dtu,dk
strange    Mikkel Strange             mikkel,strange:gmail,com
thahn      Torsten Hahn               der,hahn,torsten:googlemail,com
tjiang     Tao Jiang                  tjiang:fysik,dtu,dk
tolsen     Thomas Olsen               tolsen:fysik,dtu,dk
tprossi    Tuomas Rossi               tuomas,rossi:aalto,fi
=========  =========================  ========================================


Former committers:

=========  =========================  ========================================
id         real name
=========  =========================  ========================================
carstenr   Carsten Rostgaard          carstenr:fysik,dtu,dk 
jstenlund  Jonathan Stenlund          jonathan,stenlund:abo,fi 
madsbk     Mads Burgdorff Kristensen  madsbk:diku,dk 
s042606    Janosch Michael Rauba      s042606:fysik,dtu,dk 
s052580    Troels Kofoed Jacobsen     s052580:fysik,dtu,dk 
=========  =========================  ========================================


.. _ASE project: https://wiki.fysik.dtu.dk/ase/development/contribute.html
