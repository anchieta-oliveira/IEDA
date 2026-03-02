Getting Started
===============

Installation
------------

Clone the repository:

.. code-block:: bash

    git clone https://github.com/anchieta-oliveira/IEDA

Enter the project folder:

.. code-block:: bash

    cd IEDA/

Create a conda environment:

.. code-block:: bash

    conda create -n IEDA python=3.12
    conda activate IEDA

Install the package:

.. code-block:: bash

    pip install .

----

Quickstart
----------

Display usage options:

.. code-block:: bash

    IEDA

Compute IED matrix:

.. code-block:: bash

    IEDA matrix --pdb file.pdb --qm file.out --qm_sof orca

Generate 3D map:

.. code-block:: bash

    map_3D --pdb file.pdb --ied IED.npy --intrachain False
