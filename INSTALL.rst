Installing Connor
==================
TBD

Prerequisites
-------------
TBD



Installing
----------
The easiest way to install Connor is through PyPI. Get pip if it's
not available in your system:

You can install from source from github:

``$ pip install git+https://github.com/umich-brcf-bioinf/Connor``

If you don't have root permissions, you can install locally:

``$ pip install git+https://github.com/umich-brcf-bioinf/Connor --user``

Following the pip install, you may need to adjust your path settings to include home/.local/bin. 
Then you can execute so:

``$ connor input.bam output.bam``

If you already have prerequisite modules installed, you can also clone from github and run directly from the source like so:

``$ git clone https://github.com/umich-brcf-bioinf/Connor``

``$ connor/connor-runner.py input.bam output.bam``

