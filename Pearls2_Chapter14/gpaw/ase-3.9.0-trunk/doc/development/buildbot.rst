.. _buildbot: http://trac.buildbot.net/

========
Buildbot
========

buildbot_ is a continuous integration system.
The server called ``buildmaster`` defines and schedules processes of building
and testing of software to be performed by so called builders on machines
called ``buildslaves`` (or buildbots).
One buildslave can have several builders associated (for example a given
operating system version ``buildslave`` performs testing of several software
branches ``builders``).
One builder can also be associated with several ``buildslaves`` (several,
"identical" machines running the same operating system),
in which case ``buildmaster`` will schedule the given builder to switch
periodically between the ``buildslaves``.
At the time of writing (2012) buildbot_ does not offer satisfactory security
mechanisms, and due to clear text passwords stored/transferred both on
``buildmaster`` and ``buildslaves`` external security measures like firewall,
proxy, ssh/stunnel should be used.

The sections below describe the configuration of ``buildmaster`` and ``buildslaves``.

Configuration of buildmaster
============================

The ``buildmaster`` should be running as unprivileged user
(usually called ``buildmaster``). The user and group need to be created.

In order to install buildbot on RHEL6 system you need::

  yum install python-twisted python-jinja2 python-tempita

Additional dependencies + buildbot itself need to be installed manually::

  wget https://downloads.sourceforge.net/project/sqlalchemy/sqlalchemy/0.7.9/SQLAlchemy-0.7.9.tar.gz
  wget http://sqlalchemy-migrate.googlecode.com/files/sqlalchemy-migrate-0.7.2.tar.gz

It is sufficient to unpack the files and set :envvar:`PYTHONPATH`,
and :envvar:`PATH` accordingly.
Buildbot may need to be patched for twisted compatibility
http://www.npcglib.org/~stathis/blog/2012/10/20/bug-buildbot-dies-with-exceptions-importerror-cannot-import-name-noresource/

Proceed with buildbot configuration:

* create the master::

    buildbot create-master --relocatable python-ase

* reconfigure the master with :svn:`doc/development/master.cfg`::

    cp master.cfg python-ase
    buildbot checkconfig python-ase/master.cfg
    buildbot start python-ase

* consider setting up a crontab which starts ``buildmaster``
  after the server reboot (this solution is not reliable,
  deploy rather init scripts - see below)::

    # run every 5 minutes
    */5 * * * * sh ~/python-ase-restart.sh

  with the following ``python-ase-restart.sh``::

    bot=python-ase
    # no pid file - assume buildbot is not running
    if ! test -f $bot/twistd.pid;
    then
        buildbot restart $bot || buildbot start $bot
    else
        # pid file exists but buildbot is not running
        pid=`cat $bot/twistd.pid`
        if test -z `ps -p $pid -o pid=`;
        then
            buildbot restart $bot
        fi
    fi

* or create a system V init script under ``/etc/init.d``

  .. literalinclude:: python-ase-buildmaster

  The service is added with::

    chkconfig --add python-ase-buildmaster

  started with::

    service python-ase-buildmaster start

  end enabled for boot time with::

    chkconfig python-ase-buildmaster on

  See also an example of ``systemd`` script in the section below.

* consider protecting the ``buildmaster`` web interface by, e.g.
  apache reverse proxy (http://httpd.apache.org/docs/2.4/mod/mod_proxy.html).
  The basic configuration file may look like

  .. literalinclude:: buildbot.conf

Installation and configuration of buildslaves
=============================================

Some examples are given below. The configuration should be performed
under an unprivileged ``buildslave``, or buildbot user account.

Installation
------------

Fedora
++++++

Install with::

  yum install buildbot-slave

You can configure ``systemd`` service by creating
:file:`python-ase-fedora-18-x86_64-gcc-2.7.service` file
under ``/usr/lib/systemd/system``.

.. literalinclude:: python-ase-fedora-18-x86_64-gcc-2.7.service

``python-ase-fedora+18+x86_64+gcc+2.7-start.sh`` simply exports
the necessary environment variables (if needed) and starts ``buildslave``
(use the full path), e.g.::

   #!/bin/sh
   bot=/home/buildslave-username/python-ase-fedora+18+x86_64+gcc+2.7
   buildslave start $bot

Choose ``User`` and ``Group`` under which ``buildslave`` will be running.
The service is started with::

  systemctl start python-ase-fedora-18-x86_64-gcc-2.7.service

In order to force the service to be started at boot time create a link::

  cd /etc/systemd/system/multi-user.target.wants
  ln -s /usr/lib/systemd/system/python-ase-fedora-18-x86_64-gcc-2.7.service .

OS X
++++

Configure `Homebrew <https://wiki.fysik.dtu.dk/gpaw/install/MacOSX/homebrew.html>`_, and::

  brew install subversion

Configure a virtualenv, and then::

  pip install numpy
  pip install buildbot-slave

RHEL5
+++++

Install recent ``python-setuptools`` in order to get ``easy_install``::

  mkdir ~/buildbot-slave-el5
  export PATH=$HOME/buildbot-slave-el5:${PATH}
  export PYTHONPATH=$HOME/buildbot-slave-el5:${PYTHONPATH}
  wget http://pypi.python.org/packages/2.4/s/setuptools/setuptools-0.6c11-py2.4.egg
  sh setuptools-0.6c11-py2.4.egg --install-dir=$HOME/buildbot-slave-el5

then::

  easy_install --install-dir=$HOME/buildbot-slave-el5 zope.interface==3.6.7
  easy_install --install-dir=$HOME/buildbot-slave-el5 twisted==9.0.0  # ignore errors
  easy_install --install-dir=$HOME/buildbot-slave-el5 buildbot-slave

RHEL6
+++++

Install ``build-slave`` and dependencies::

  mkdir ~/buildbot-slave-el6
  export PATH=$HOME/buildbot-slave-el6:${PATH}
  export PYTHONPATH=$HOME/buildbot-slave-el6:${PYTHONPATH}
  easy_install --install-dir=$HOME/buildbot-slave-el6 buildbot-slave

RHEL7
+++++

Install ``build-slave`` and dependencies::

  mkdir ~/buildbot-slave-el7
  export PATH=$HOME/buildbot-slave-el7:${PATH}
  export PYTHONPATH=$HOME/buildbot-slave-el7:${PYTHONPATH}
  easy_install --install-dir=$HOME/buildbot-slave-el7 buildbot-slave

Windows
+++++++

``build-slave`` can be installed and configured to start as a service on
Windows `<http://trac.buildbot.net/wiki/RunningBuildbotOnWindows>`_.
This involves several steps:

1. install Python(x,y) from https://code.google.com/p/pythonxy/wiki/Downloads

2. configure distutils to use mingw. First enable showing file extensions::

    Open a folder with IE -> Folder and search options-> View -> Folder Options:
    Check: Show hidden files, ...; uncheck: Hide extensions for known file types.

   Then, in notepad, create C:\\python27\\lib\\distutils\\distutils.cfg,
   containing::

    [build]
    compiler=mingw32

3. install ``build-slave`` on the command line::

    C:\python27\scripts\easy_install.exe buildbot-slave

4. create a local (domain ``computer-name``) user that will run the
   ``buildbot`` service::

     control panel->administrative tools->computer management->local users and groups->users->new user: buildslave-username.
     Click the created user: member of: administrators->check names

5. grant ``buildslave-username`` the ability to run the services.
   Login as ``computer-name\buildslave-username``:
   Run ``secpol.msc`` on the command line as administrator
   (task bar->cmd->right click: run as administrator):

   - Select the "Local Policies" folder
   - Select the "User Rights Assignment" folder
   - Double click "Log on as a service"
   - Use "Add User or Group..." to add your user here.

   Select the correct "from this location"
   (may require login as the current domain administrator) and Enter the
   object names: ``computer-name\buildslave-username``.

6. on the command line install the ``buildbot`` service::

    buildbot_service.py --user computer-name\buildslave-username --password thepassword --startup auto install

7. start the service (for the moment it does not start any ``buildslave``,
   because they are not configured yet)::

     Start->Control Panel> Administrative Tools->Services->Buildbot (Start)

   There are additional steps mentioned in the buildbot wiki, but it seems
   just to work on Windows 7.

8. run regedit as administrator (type "regedit" on the command line) and add
   "directories" parameter of the String Value type, containing paths to all
   your ``buildslave`` instances (they will be configured in the Configuration
   section below)::

     HKEY_LOCAL_MACHINE\System\CurrentControlSet\services\Buildbot->paramaters->new (String Value): directories C:\python-ase-windows+7+AMD64+msc+2.7;C:\proj2

9. configure ``buildslave`` instance as described in the Configuration section
   below and start the service again (point 7.). Test that ``buildslave``
   comes online, and verify that the service starts after reboot. 

Configuration
-------------

After having installed the buildbot create a name which will identify
your ``buildslave``.
Obtain the first part of the name for your ``buildslave`` by
running :svn:`doc/development/master.cfg`::

  python master.cfg

This will return a string like ``redhat+6+x86_64+gcc+2.6``
(OS+OSversion+Bitness+Ccompiler+PythonVersion).
Append it with something that identifies your ``buildslave``,
for example ``gpu``.
For a very special system you can use a name like
``mycluster+x86_64+open64+2.5 gpu``.
Note that ASE buildbot will perform verification of python version based on
the ``buildslave`` name so stay consistent!

Generate a strong ``buildslave`` password with :ref:`devel_passwd_create`.
Don't reuse any valuable passwords. You don't need to remember it,
buildbot stores it plain text!

Create the ``buildslave`` with::

  buildslave create-slave python-ase-redhat+6+x86_64+gcc+2.6 buildbot.fysik.dtu.dk:ASE_BB_PORT "redhat+6+x86_64+gcc+2.6 gpu" "password"

ASE_BB_PORT is the port ASE buildbot uses for communication with the ``buildslave``.
You will have to tell us the name and password of your ``buildslave``.
Please contact ase-developers list :ref:`mailing_lists`, but
don't send the name and password there!

Edit the ``python-ase-redhat+6+x86_64+gcc+2.6/info/{admin,info}`` files:
describe your ``buildslave`` configuration relevant for the builder process
in the ``info`` file.

Note that before starting the slave you need to perform an temporary
svn checkout of ASE in order to accept the certificate permanently.

Start the ``buildslave`` with::

  buildslave start python-ase-redhat+6+x86_64+gcc+2.6

Don't forget to configure a crontab job or a service as described in the
previous sections.

By default all slaves run the continuous integration for the trunk.
If you prefer your ``buildslave`` works also on one of the branches, write
this in the email to ase-developers :ref:`mailing_lists`.
