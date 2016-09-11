Installing Connor
==================
Connor requires python 2.7 or later and has been tested with:

* Python 2.7 and 3.4
* pysam 0.8.4 and 0.9.0
* OSX and \*nix RHEL6/7

Connor does not work in Windows OS because it depends on the python library
pysam, which is not supported on Windows.

Prerequisites
-------------
* Python (2.7 or later) and pip. See https://www.python.org/downloads/ for more details on
  downloading a recent version of Python 2/3.
.. note:: Pip installs all required libraries; see [Installing] below.

* cython
* pysam
* sortedcontainters


Installing
----------

* The simplest way to install Connor is through PyPI:

  ``$ pip install connor``

* If you don't have root permissions, you can install connor locally:

  ``$ pip install connor --user``

  Following a --user install, you may need to adjust your path settings to
  include $HOME/.local/bin. 

* You can also install directly from source from github:

  ``$ pip install git+https://github.com/umich-brcf-bioinf/Connor``

Uninstalling
------------
``$ pip uninstall connor``


Advanced / Connor developers
----------------------------
Connor has been tested with virtualenv and Conda. (Note Conda's pysam version
lags behind so if you're using conda and installing modules manually, we
recommend you install pysam via 'pip install' instead of 'conda install').

If you already have prerequisite modules installed, you can also clone from
github and run directly from the source like so:

   ``$ git clone https://github.com/umich-brcf-bioinf/Connor``

   ``$ connor/connor-runner.py input.bam output.bam``

For running automated tests, Connor also requires:
 * nosetests
 * testfixtures


