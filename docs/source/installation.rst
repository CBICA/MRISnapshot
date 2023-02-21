************
Installation
************

MRISnapshot can be installed using `pip`.
To start the installation, start a terminal and check the installed version
of Python.
The output should show a version `3.8` or newer.
If the command is not recognized, Python is not installed or the `python` command
is not found by the system.
If Python is installed--for instance in Windows, it is possible to call the command
with the full path.

.. tabs::

   .. code-tab:: shell Linux

         /usr/bin/python --version

.. code-block:: shell

    python -m venv .env
    source .env/bin/activate
    python -m pip install https://github.com/CBICA/NiBAx
