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
