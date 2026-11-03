:html_theme.sidebar_secondary.remove:

prepmd
======

.. ccpbiosim-hero::
   :logo-light: _static/logos/prepmd-logo-black.svg
   :logo-dark: _static/logos/prepmd-logo-white.svg
   :logo-alt: prepmd logo
   :tagline: prepmd
   :button: Installation | Get started
   :github: CCPBioSim/prepmd

   Automatically prepare PDB and mmCIF structures for molecular dynamics simulations

A utility to automatically prepare structures from the PDB for molecular dynamics
simulation and perform minimisations and simple MD simulations.

Quick start
-----------

.. tab-set::

   .. tab-item:: Install

      .. code-block:: bash

         conda install -c ccpbiosim prepmd

   .. tab-item:: Prepare

      .. code-block:: bash

         prepmd 6xov 6xov_processed.cif

   .. tab-item:: Run

      .. code-block:: bash

         runmd 6xov_processed.cif --traj_out traj.xtc --md_steps 5000

Explore the documentation
-------------------------

.. grid:: 1 2 2 3
   :gutter: 3

   .. grid-item-card:: 📦 Installation
      :link: Installation
      :link-type: doc

      Install prepmd with Conda and add a MODELLER licence.

   .. grid-item-card:: 🚀 Quickstart
      :link: Quickstart
      :link-type: doc

      Prepare a structure and run a short simulation.

   .. grid-item-card:: 🧬 Prepare Structures
      :link: Prep-Structure
      :link-type: doc

      Prepare structure files for simulation with prepmd.

   .. grid-item-card:: ⚙️ Run Simulations
      :link: Run-with-runmd
      :link-type: doc

      Run MD simulations with runmd.

   .. grid-item-card:: 🐍 Python API
      :link: Python-API
      :link-type: doc

      Use prepmd from Python.

   .. grid-item-card:: 📖 API Reference
      :link: reference/prepmd/index
      :link-type: doc

      Reference for the prepmd package.

.. toctree::
   :maxdepth: 2
   :hidden:

   Installation <Installation>
   Quickstart <Quickstart>
   Prepare <Prep-Structure>
   Run <Run-with-runmd>
   What Next <Whats-next>
   Python API <Python-API>
   Developer Notes <Developer-Notes>
   API Reference <reference/prepmd/index>
