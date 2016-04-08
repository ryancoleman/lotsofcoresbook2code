.. _translate:

Translate ASE
=============

You can contribute by translating the ASE GUI, :file:`ase-gui`, into your
language.


How to translate
----------------

If any of the below steps prove difficult, be sure to ask on the
:ref:`developer mailing list <mailing_lists>`.  These steps should work on GNU/Linux.

* :ref:`Download <download_and_install>` ASE.
* Go to :file:`ase/gui/po`.  There is a directory of the form :file:`{ll}` or :file:`{ll}_{LL}` for each language, where :file:`{ll}` is the `language code`_ and :file:`{LL}` the `country code`_.  The latter is necessary only if variants of the same language are spoken in multiple countries.
* If your language is not there already, run :file:`LANG={ll} make init`, substituting the desired language code.  If necessary append the country code as well: :file:`LANG={ll_LL} ...`.
* There should now be a template file for your language, :file:`{ll}/LC_MESSAGES/ag.po`, which can be filled out.

You can edit the po-file with any text editor.  It is easiest with a dedicated po-editor such as gtranslator, poedit or the gettext mode in EMACS (the package ":file:`gettext-el`").  Fill out the missing :file:`msgstr` entries like this::

  #: ../energyforces.py:61
  msgid "Calculate potential energy and the force on all atoms"
  msgstr "Beregn potentiel energi og kræfter på alle atomer"

Your editor may wholly or partially hide some of the difficult
formatting syntax.  This next example shows a more syntactically
complex case in case you need it::

  #: ../calculator.py:107
  msgid ""
  "GPAW implements Density Functional Theory using a\n"
  "<b>G</b>rid-based real-space representation of the wave\n"
  "functions, and the <b>P</b>rojector <b>A</b>ugmented <b>W</b>ave\n"
  "method for handling the core regions.  \n"
  msgstr ""
  "GPAW implementerer tæthedsfunktionalteori med en <b>G</b>itterbaseret\n"
  "repræsentation af bølgefunktioner i det reelle rum, samt\n"
  "<b>P</b>rojector <b>A</b>ugmented <b>W</b>ave-metoden til behandling\n"
  "af regionen omkring atomkerner.  \n"

If you are maintaining an existing translation, there may be some
"fuzzy" messages.  These are translations that were written
previously but now need to be reviewed, maybe because the original
string has been slightly modified.  Edit them as appropriate and remove the
"fuzzy" flag.

There will be a few special constructs such as string substitution
codes :file:`%(number)d` or :file:`%s`.  These should remain unchanged
in the translation as they are replaced by numbers or text at runtime.
An underscore like in :file:`msgid "_File"` indicates that `F` is a
shortcut key.  Conflicting shortcut keys are not a big problem, but
avoid them if you see them.  Finally, some messages may have a lot of
whitespace in them.  This is due to bad programming style; just try to
get approximately the same spacing in your translation.

Already after writing a few translations, you can check that the
translation works as expected by following the instructions in the
next section.

Check and commit your translation
---------------------------------

* You can check the syntax by running :file:`msgfmt -cv ag.po`.  This will
  report any syntax errors.

* You can test your translation in :file:`ase-gui` directly.  First issue
  the command :file:`make` in :file:`ase/gui/po`, then reinstall ASE
  using the usual procedure.  The translations will then be in the
  newly installed ASE.  If you translate into the same language as
  your computer's locale, you should see the translations when you
  start :file:`ase-gui` normally.  If you translate ASE into another
  language, then run :file:`LANG={ll}_{LL}.UTF-8 ase-gui`.  On some
  operating systems you may need to run
  :file:`LANGUAGE={ll}_{LL}.UTF-8 ase-gui` instead.

Depending on your operating system, you may need to install
:file:`gettext` or :file:`locales`.

Send the partially or completely translated po-file to the developers
mailing list and ask to have it committed.  In fact, we will be quite thrilled
if you send an e-mail even before you start, and be sure to send one
whenever you have questions.

.. note::

  Certain uncommon languages such as Lojban, Anglo-Saxon or Klingon
  may not be compatible with our current build system.  Please let us
  know if you want to translate ASE into such languages.

Maintaining translations
------------------------

Messages will once in a while be added or changed in the ASE.  Running
:file:`make` in :file:`ase/gui/po` automatically synchronizes all templates with
the messages in the current source tree while maximally reusing the
existing translations.  Some strings may be marked "fuzzy", indicating
that they need review by translators (this happens e.g. if an English
message is changed only slightly).  One can then update the few fuzzy
or untranslated messages.  The obvious time to do this is shortly
before a new stable release.

If you are a committer, please run :file:`make` before committing and
briefly check by running the translated ase-gui that nothing is obviously horrible.

.. _language code: http://www.gnu.org/software/gettext/manual/gettext.html#Language-Codes
.. _country code: http://www.gnu.org/software/gettext/manual/gettext.html#Country-Codes
