.. _literature:

----------
Literature
----------


Links to guides and manual pages
--------------------------------

* The GPAW calculator :ref:`manual`

* The :ref:`devel` pages

* The :ref:`guide for developers <developersguide>`

* The code :ref:`overview`

* The :ref:`algorithms` used in the code


.. _literature_reports_presentations_and_theses:

Reports, presentations, and theses using gpaw
---------------------------------------------

* Summer-school 2014 talk about `PAW and GPAW`_

* A short note on the basics of PAW: `paw note`_

* A master thesis on the inclusion of non-local exact exchange in the
  PAW formalism, and the implementation in gpaw: `exact exchange`_

* A master thesis on the inclusion of a localized basis in the PAW
  formalism, plus implementation and test results in GPAW: `lcao`_

* A master thesis on the inclusion of localized basis sets in the PAW
  formalism, focusing on basis set generation and force calculations:
  `localized basis sets`_

* A course report on a project involving the optimization of the
  setups (equivalent of pseudopotentials) in gpaw: `setup
  optimization`_

* Slides from a talk about PAW: `introduction to PAW slides`_

* Slides from a talk about GPAW development: `gpaw for developers`_

* Slides from a mini symposium during early development stage: `early gpaw`_

.. _paw note: ../paw_note.pdf
.. _exact exchange: ../_static/rostgaard_master.pdf
.. _lcao: ../_static/marco_master.pdf
.. _localized basis sets: ../_static/askhl_master.pdf
.. _setup optimization: ../_static/askhl_10302_report.pdf
.. _introduction to PAW slides: ../_static/mortensen_paw.pdf
.. _gpaw for developers: ../_static/mortensen_gpaw-dev.pdf
.. _early gpaw: ../_static/mortensen_mini2003talk.pdf
.. _PAW and GPAW: ../_static/ss14.pdf


.. _paw_papers:

Articles on the PAW formalism
-----------------------------

The original article introducing the PAW formalism:
   | P. E. Blöchl
   | `Projector augmented-wave method`__
   | Physical Review B, Vol. **50**, 17953, 1994

   __ http://dx.doi.org/10.1103/PhysRevB.50.17953

A different formulation of PAW by Kresse and Joubert designed to make the transistion from USPP to PAW easy.
  | G. Kresse and D. Joubert
  | `From ultrasoft pseudopotentials to the projector augmented-wave method`__
  | Physical Review B, Vol. **59**, 1758, 1999

  __ http://dx.doi.org/10.1103/PhysRevB.59.1758

A second, more pedagogical, article on PAW by Blöchl and co-workers.
  | P. E. Blöchl, C. J. Först, and J. Schimpl
  | `Projector Augmented Wave Method: ab-initio molecular dynamics with full wave functions`__
  | Bulletin of Materials Science, Vol. **26**, 33, 2003

  __ http://www.ias.ac.in/matersci/


.. _gpaw_publications:

Citations of the GPAW method papers
-----------------------------------

.. image:: citations.png
   :width: 750

(updated on May 18, 2013)

The five method papers are:

gpaw1:
    \J. J. Mortensen, L. B. Hansen, and K. W. Jacobsen

    `Real-space grid implementation of the projector augmented
    wave method`__

    Physical Review B, Vol. **71**, 035109 (2005)

    __ http://dx.doi.org/10.1103/PhysRevB.71.035109

tddft:
    \M. Walter, H. Häkkinen, L. Lehtovaara, M. Puska, J. Enkovaara,
    C. Rostgaard, and J. J. Mortensen

    `Time-dependent density-functional theory in the projector
    augmented-wave method`__

    Journal of Chemical Physics, Vol. **128**, 244101 (2008)

    __ http://dx.doi.org/10.1063/1.2943138

lcao:
    \A. H. Larsen, M. Vanin, J. J. Mortensen, K. S. Thygesen, and
    K. W. Jacobsen

    `Localized atomic basis set in the projector augmented wave method`__

    Physical Review B, Vol. **80**,  195112 (2009)

    __ http://dx.doi.org/10.1103/PhysRevB.80.195112

gpaw2:
    \J. Enkovaara, C. Rostgaard, J. J. Mortensen, J. Chen, M. Dulak,
    L. Ferrighi, J. Gavnholt, C. Glinsvad, V. Haikola, H. A. Hansen,
    H. H. Kristoffersen, M. Kuisma, A. H. Larsen, L. Lehtovaara,
    M. Ljungberg, O. Lopez-Acevedo, P. G. Moses, J. Ojanen, T. Olsen,
    V. Petzold, N. A. Romero, J. Stausholm, M. Strange, G. A. Tritsaris,
    M. Vanin, M. Walter, B. Hammer, H. Häkkinen, G. K. H. Madsen,
    R. M. Nieminen, J. K. Nørskov, M. Puska, T. T. Rantala,
    J. Schiøtz, K. S. Thygesen, and K. W. Jacobsen   

    `Electronic structure calculations with GPAW: a real-space
    implementation of the projector augmented-wave method`__ 

    \J. Phys.: Condens. Matter **22**, 253202 (2010)

    __ http://stacks.iop.org/0953-8984/22/253202

response:
    Jun Yan, Jens. J. Mortensen, Karsten W. Jacobsen, and Kristian S. Thygesen

    `Linear density response function in the projector augmented wave method:
    Applications to solids, surfaces, and interfaces`__

    Phys. Rev. B **83**, 245122 (2011)

    __ http://prb.aps.org/abstract/PRB/v83/i24/e245122


All citing articles:

.. csv-table::
   :file: citations.csv
   :header: #, title
   :widths: 1, 15

