.. _stereo_combination:

Combining the stereo events: ``sst1mpipe-stereo``
=================================================

The SST-1M telescopes have no array trigger: each telescope writes its own events, one raw data
(zfits) file per run and per telescope. Each run is processed into a DL1 file by ``sst1mpipe-process``
(see :ref:`data_processing`), the DL1 files of the runs of a night are merged for each telescope
(``ctapipe-merge``, see :ref:`merging_step`), and the coincident events of the merged files of the
telescopes are combined into stereo events by ``sst1mpipe-stereo``, a `ctapipe.core.Tool`:

.. code-block:: text

    SST1M1_20260120_merged.dl1.h5 ─┐
                                   ├──► sst1mpipe-stereo ──► SST1M_20260120_merged.stereo.dl1.h5
    SST1M2_20260120_merged.dl1.h5 ─┘

The coincident events are matched and combined by a `~sst1mpipe.stereo.StereoEventCombiner`. The
stereo reconstruction (stereo DL2) is a later step.


Usage
-----

.. code-block:: bash

    sst1mpipe-stereo \
        --input SST1M1_20260120_merged.dl1.h5 \
        --input SST1M2_20260120_merged.dl1.h5 \
        --output SST1M_20260120_merged.stereo.dl1.h5 \
        --config sst1mpipe/resources/config/sst1mpipe_data_config.json

``sst1mpipe-stereo --help`` lists the options (``--help-all`` all the settings). The default settings
are in ``sst1mpipe_data_config.json``.

.. warning::

    The combination does not work properly with more than two telescopes, see
    `~sst1mpipe.stereo.StereoEventCombiner`.


Combiners
---------

The matching method is chosen with ``--combiner``: each combiner implements the matching of the events
of two telescopes, and the combination of each part of the event is a method of the base class, which a
child combiner can override.

.. autoclass:: sst1mpipe.stereo.StereoEventCombiner
   :members: match, match_pair, combine, combine_index, combine_trigger, combine_r0, combine_r1,
             combine_dl0, combine_dl1, combine_dl2, combine_muon, combine_monitoring,
             combine_array_pointing, combine_simulation

.. autoclass:: sst1mpipe.stereo.WhiteRabbitClosestCombiner
   :members: match_pair

.. autoclass:: sst1mpipe.stereo.SlidingWindowCombiner
   :members: match_pair, find_time_offset, n_coincidences

.. autoclass:: sst1mpipe.stereo.SWATEventIdCombiner
   :members: match_pair


Examples
--------

Matching of the events of the merged files of a night, without the tool:

.. code-block:: python

    from sst1mpipe.scripts.sst1mpipe_stereo_tool import read_telescope_triggers
    from sst1mpipe.stereo import StereoEventCombiner

    triggers = {
        21: read_telescope_triggers(["SST1M1_20260120_merged.dl1.h5"]),
        22: read_telescope_triggers(["SST1M2_20260120_merged.dl1.h5"]),
    }
    combiner = StereoEventCombiner.from_name("SWATEventIdCombiner")
    groups = combiner.match(triggers)
    print(triggers[22][groups[22]][["obs_id", "event_id"]])

A child combiner overriding a part of the combination (it must be imported before the tool is
configured, to be found by ``--combiner``):

.. code-block:: python

    from sst1mpipe.stereo import WhiteRabbitClosestCombiner


    class LowestTelescopeIds(WhiteRabbitClosestCombiner):
        """stereo events with the ids of the event of the lowest telescope id"""

        def combine_index(self, stereo, events):
            stereo.index = events[min(events)].index


Tool
----

.. autoclass:: sst1mpipe.scripts.sst1mpipe_stereo_tool.StereoTool

.. autofunction:: sst1mpipe.scripts.sst1mpipe_stereo_tool.read_telescope_triggers

.. autoclass:: sst1mpipe.scripts.sst1mpipe_stereo_tool.TelescopeEventReader
   :members: get

.. autoclass:: sst1mpipe.scripts.sst1mpipe_stereo_tool.CombinedRuns

.. autofunction:: sst1mpipe.io.append_monitoring_tables

.. autofunction:: sst1mpipe.io.containers.to_sst1m_event

.. autofunction:: sst1mpipe.stereo.nearest_index
