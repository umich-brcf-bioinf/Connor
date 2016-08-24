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
.. note:: Pip installs all required libraries; see [Installing] below.

* cython
* pysam


Installing
----------

* The simplest way to install Connor is through PyPI:

  *(This option only available after general release.)*

  ``$ pip install connor``

* You can also install directly from source from github:

  ``$ pip install git+https://github.com/umich-brcf-bioinf/Connor``

* If you don't have root permissions, you can install connor locally:

  ``$ pip install git+https://github.com/umich-brcf-bioinf/Connor --user``

  Following a --user install, you may need to adjust your path settings to
  include $HOME/.local/bin. 


Uninstalling
------------
``$ pip uninstall connor``


Advanced / Connor developers
----------------------------
If you already have prerequisite modules installed, you can also clone from
github and run directly from the source like so:

   ``$ git clone https://github.com/umich-brcf-bioinf/Connor``

   ``$ connor/connor-runner.py input.bam output.bam``

For running automated tests, Connor also requires:
 * nosetests
 * testfixtures


