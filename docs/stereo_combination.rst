.. _stereo_combination:

Combining the stereo events: ``sst1mpipe-stereo``
=================================================

The SST-1M telescopes have no array trigger: each telescope triggers and writes its own events,
one DL1 file per run and per telescope (``sst1mpipe-process``, see :ref:`data_processing`). The
events of a shower seen by both telescopes are found afterwards, from the trigger times or the
array event ids of the events, and combined into stereo events by the ctapipe tool
``sst1mpipe-stereo``:

.. code-block:: text

    SST1M1_20260120_1179_W1.dl1.h5 ─┐
    SST1M1_20260120_1180_W1.dl1.h5 ─┤
                                    ├──► sst1mpipe-stereo ──► SST1M_20260120_1179.stereo.dl1.h5
    SST1M2_20260120_1102_W1.dl1.h5 ─┤        (StereoEventCombiner: match, combine)
    SST1M2_20260120_1103_W1.dl1.h5 ─┘

The matching and the combination of the events are done by a ``StereoEventCombiner``
(``sst1mpipe.stereo``), a configurable ctapipe component chosen with ``--combiner``. The stereo
reconstruction (stereo DL2, e.g. the geometry of the showers) is a later step.

.. contents::
   :local:
   :depth: 2


Usage
-----

.. code-block:: bash

    sst1mpipe-stereo \
        --input SST1M1_20260120_1179_W1.dl1.h5 \
        --input SST1M2_20260120_1102_W1.dl1.h5 \
        --output SST1M_20260120_1179.stereo.dl1.h5 \
        --combiner WhiteRabbitClosestCombiner

The input files are DL1 files of ``sst1mpipe-process``, of one telescope each, any number of files
per telescope (e.g. all the runs of a night), of at least two telescopes including the reference
telescope. The telescope of each file is read from its subarray.

.. list-table:: Main options
   :header-rows: 1
   :widths: 35 65

   * - Option
     - Description
   * - ``-i``, ``--input``
     - Input DL1 file, repeated for each file (``StereoTool.input_files``).
   * - ``-o``, ``--output``
     - Output file (``DataWriter.output_path``).
   * - ``--combiner``
     - Matching method: ``WhiteRabbitClosestCombiner`` (default), ``SlidingWindowCombiner`` or
       ``SWATEventIdCombiner`` (``StereoTool.combiner_type``), see :ref:`stereo_matching`.
   * - ``--reference-tel``
     - Reference telescope (``StereoEventCombiner.reference_tel_id``, 21), see
       :ref:`stereo_reference_telescope`.
   * - ``--keep-mono``
     - Also write the events which are not combined (mono events), see :ref:`stereo_mono_events`.
   * - ``--overwrite``
     - Overwrite the output file (``DataWriter.overwrite``).
   * - ``--progress``
     - Show a progress bar.
   * - ``--StereoTool.copy_monitoring``
     - Copy the monitoring data of the input files (true by default), see :ref:`stereo_output`.
   * - ``--DataWriter.write_dl1_images``, ``--DataWriter.write_dl1_parameters``, ...
     - Data levels written, see :ref:`stereo_output`.

``--help`` lists the main options and ``--help-all`` all the settings, including those of each
combiner (e.g. ``--WhiteRabbitClosestCombiner.max_time_diff``). As for the other tools, the
settings can be given in a configuration file (``--config``, JSON or YAML, see
:ref:`processing_config`):

.. code-block:: json

    {
      "StereoTool": {"combiner_type": "SlidingWindowCombiner", "keep_mono": true},
      "StereoEventCombiner": {"reference_tel_id": 21},
      "SlidingWindowCombiner": {"window_half_width": "500 us", "offset_step": "10 us"},
      "DataWriter": {"write_dl1_images": false}
    }

The settings of ``StereoEventCombiner`` apply to all the combiners.

The simulations are not combined: their stereo events are already combined by sim_telarray (array
trigger, same ``event_id`` in all the telescopes), and ``sst1mpipe-stereo`` refuses simulated files.


.. _stereo_matching:

Matching of the coincident events
---------------------------------

The coincident events are matched from the trigger tables of the DL1 files
(``/dl1/event/telescope/trigger``: ``obs_id``, ``event_id``, trigger ``time`` with the ns
precision of White Rabbit, ``event_type``), by the method of the combiner. For all the
combiners:

* Only the events triggered by the cameras (``event_type`` SUBARRAY) are matched. The pedestal
  events (SKY_PEDESTAL, ...) are never combined; ``sst1mpipe-process`` does not write them anyway.
* The events of each telescope are matched to those of the reference telescope (see
  :ref:`stereo_reference_telescope`).
* Each event is in at most one stereo event.

The number of events matched for each telescope is given in the log.

``WhiteRabbitClosestCombiner`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For data with White Rabbit timestamps (taken after 25 September 2023). Each event of the reference
telescope is matched to the event of the other telescope with the closest trigger time, if:

* their trigger times differ by less than ``max_time_diff`` (10 µs), and
* it is also the closest event of the reference telescope to that event (mutual closest events).

The trigger times of coincident events differ by a few tens of ns (32 and 48 ns for the two
coincident events of the test data), while the events of a telescope are hundreds of µs apart at the
usual trigger rates: 10 µs does not give random coincidences.

``SlidingWindowCombiner``
~~~~~~~~~~~~~~~~~~~~~~~~~

For data without White Rabbit timestamps (taken before 25 September 2023), whose clocks have an
offset. An event of the other telescope is coincident with an event of the reference telescope if
its trigger time is within ``window_half_width`` (1 ms, strictly) of the trigger time of the event of
the reference telescope plus the offset of the clocks (time of the other telescope - time of the
reference telescope).

* The offset is ``time_offset`` if given. Otherwise it is searched from ``offset_start`` (-10 ms) to
  ``offset_stop`` (10 ms, excluded) in steps of ``offset_step`` (20 µs): the offset giving the most
  coincidences is used, and given in the log.
* An event with several events of the other telescope in its window is ambiguous and not matched,
  in both directions. The number of ambiguous events is given in the log.

.. warning::

    With its default 1 ms window, the method gives random coincidences at the usual trigger rates.
    On the test data (with White Rabbit timestamps), it finds a "clock offset" of 3.36 ms and 7 events,
    while the two telescopes have 2 coincident events. Use ``WhiteRabbitClosestCombiner`` for data with
    White Rabbit timestamps, or a fixed ``time_offset`` with a narrow window
    (``--SlidingWindowCombiner.time_offset="0 s" --SlidingWindowCombiner.window_half_width="10 us"``).

``SWATEventIdCombiner``
~~~~~~~~~~~~~~~~~~~~~~~

For data taken after 30 January 2024: SWAT gives the coincident events of the telescopes the same
array event id (``arrayEvtNum`` of the raw data), which is the ``event_id`` of the events in the DL1
files. Only the event ids are used:

* the events without the same id are not combined, whatever their trigger times, and nothing else is
  tried if the files have no SWAT event ids (the events then have their camera event numbers, which
  do not match);
* an id found several times in the events of a telescope is not matched (warning in the log).


.. _stereo_reference_telescope:

Reference telescope
-------------------

A stereo event is an event of the reference telescope (``--reference-tel``, 21 by default) with the
coincident events of the other telescopes. Swapping the reference telescope:

* changes the ``obs_id`` and ``event_id`` of the stereo events, which are those of the event of the
  reference telescope. With SWAT event ids, the ``event_id`` of the telescopes are the same, and only
  the ``obs_id`` (run of the reference telescope) changes;
* gives the same stereo events with ``WhiteRabbitClosestCombiner`` (mutual closest events) and
  ``SWATEventIdCombiner`` (same ids), whose matching of two telescopes is symmetric, and with
  ``SlidingWindowCombiner`` with a fixed ``time_offset`` of opposite sign;
* can give other stereo events with ``SlidingWindowCombiner`` when the offset is searched: the
  search keeps the first offset with the most coincidences, from ``offset_start``. On the test data
  (random coincidences of the 1 ms window), the offset found is 3.36 ms with 7 events with tel 21 as
  reference, and -3.8 ms with 5 events with tel 22.

More than two telescopes
~~~~~~~~~~~~~~~~~~~~~~~~

The events of each other telescope are matched to those of the reference telescope, separately. A
stereo event has the events of all the telescopes matched to the same event of the reference
telescope (e.g. 21, 22 and 23). A coincidence of other telescopes without an event of the reference
telescope (e.g. 22 and 23 without 21) is not a stereo event: its events are mono events.

``StereoEventCombiner.match`` returns the stereo events as groups of events (one row per stereo
event, the row of the event of each telescope or -1): a matching of all the telescopes at once can
be implemented by overriding ``match`` in a new combiner, without changing the tool.


.. _stereo_combination_of_events:

Combination of the events
-------------------------

The events of the telescopes are read from their DL1 files with the ctapipe ``HDF5EventSource``, and
each stereo event is a new ``SST1MArrayEventContainer`` (``sst1mpipe.io.containers``) made by
``StereoEventCombiner.combine``, with one method per part of the event:

.. list-table::
   :header-rows: 1
   :widths: 25 30 45

   * - Method
     - Part of the stereo event
     - Content
   * - ``combine_index``
     - ``index``
     - ``obs_id`` and ``event_id`` of the event of the reference telescope.
   * - ``combine_trigger``
     - ``trigger``
     - ``tels_with_trigger``: the telescopes of the events; ``time``: the earliest trigger time of the
       events; ``event_type``: SUBARRAY; ``tel[tel_id]``: the trigger of each telescope (time, event
       type, number of trigger pixels).
   * - ``combine_r0``, ``combine_r1``, ``combine_dl0``
     - ``r0.tel``, ``r1.tel``, ``dl0.tel``
     - Data of each telescope, from its event (if the DL1 files have them).
   * - ``combine_dl1``
     - ``dl1.tel``
     - DL1 of each telescope (image, peak time, image mask, parameters), from its event.
   * - ``combine_dl2``
     - ``dl2.tel``
     - DL2 of each telescope, from its event. The stereo DL2 (``dl2.stereo``) is not filled.
   * - ``combine_muon``
     - ``muon.tel``
     - Muon parameters of each telescope, from its event.
   * - ``combine_monitoring``
     - ``monitoring.tel``
     - Monitoring of each telescope (pointing, camera monitoring), from its event.
   * - ``combine_array_pointing``
     - ``monitoring.pointing``
     - **Not called**: the array pointing of the stereo events is not filled. The method keeps the
       array pointing of the events if they all have the same, and raises ``NotImplementedError``
       otherwise: the combination of different array pointings is not defined yet.
   * - ``combine_simulation``
     - ``simulation``
     - Raises ``NotImplementedError`` (simulated events are combined by sim_telarray).

The data of each telescope are those of its mono event, unchanged. A new combiner can change the
combination of any part of the event by overriding its method, see :ref:`stereo_new_combiner`.


.. _stereo_mono_events:

Mono events
-----------

The events which are not combined (mono events: events without a coincident event, events of a
type other than SUBARRAY, events of other telescopes without an event of the reference telescope)
are not written by default. With ``--keep-mono``, they are written as they are (converted to
``SST1MArrayEventContainer``, with their own ``obs_id``, ``event_id`` and event type), in time order
with the stereo events. The numbers of stereo events and of mono events of each telescope are
given in the log.

Each event of the output file must have its own ``obs_id`` and ``event_id``: if a mono event of
another telescope has the ids of another output event (e.g. the same file given twice), the tool
stops with an error.


.. _stereo_output:

Output file
-----------

The events are written in time order by the ctapipe ``DataWriter``: the data levels written are its
settings (``--DataWriter.write_dl1_images``, ``write_dl1_parameters``, ``write_dl2``,
``write_r0_waveforms``, ...), and only the data levels of the input files can be written. The output
file has:

* the subarray trigger (``/dl1/event/subarray/trigger``: ``obs_id``, ``event_id``, time,
  ``tels_with_trigger``, ``event_type``) and the trigger of each telescope
  (``/dl1/event/telescope/trigger``) of the events;
* the DL1 (and the data levels asked for) of each telescope, in the tables of each telescope
  (``/dl1/event/telescope/parameters/tel_021``, ...);
* the subarray of all the telescopes (positions of the DL1 files), and the observation and
  scheduling blocks of the runs of all the input files;
* with ``--StereoTool.copy_monitoring`` (default), the monitoring data of the input files with events
  in the output (pointing of the telescopes, statistics of the pedestal events, camera calibration
  coefficients: groups ``/r0/monitoring``, ``/dl0/monitoring`` and ``/dl1/monitoring``), the tables of
  each telescope appended run after run (``sst1mpipe.io.append_monitoring_tables``).

The output file has no stereo reconstruction (stereo DL2), and no summary of the processing of the
runs (``/dl1/info`` of the DL1 files).

.. note::

    When a DL1 file with monitoring data is read, the ``HDF5EventSource`` interpolates the pointing of
    the telescope at the time of each event, and stops at an event out of the time range of the
    pointing table. ``sst1mpipe-stereo`` does not write this pointing (the monitoring tables are
    copied instead): it sets ``PointingInterpolator.bounds_error`` to false, so that such an event gets a
    NaN pointing, unless ``--PointingInterpolator.bounds_error`` is given.

Reading the stereo files
~~~~~~~~~~~~~~~~~~~~~~~~

The stereo files are read like the DL1 files of ``sst1mpipe-process``: by the ctapipe
``HDF5EventSource`` (stereo events with the DL1 of both telescopes), ``ctapipe.io.TableLoader``, or
``load_dl1_sst1m`` for each telescope (``tel="tel_021"``, with the pointing of the telescope from the
copied monitoring data).

.. warning::

    ``sst1mpipe_dl1_dl2 --stereo`` (and ``load_dl1_sst1m(..., stereo=True)``) do not read the files of
    ``sst1mpipe-stereo`` yet. They need the stereo geometry reconstruction
    (``/dl2/event/telescope/impact/HillasReconstructor/tel_XXX`` and
    ``/dl2/event/subarray/geometry/HillasReconstructor``: impact distances, ``h_max`` and core, features
    of the stereo models), done by a later step, and the ``/dl1/info`` table.


.. _stereo_new_combiner:

Writing a new combiner
----------------------

A combiner is a subclass of ``StereoEventCombiner`` which implements ``match_pair``: the coincident
events of the reference telescope and of another telescope, from their events of type SUBARRAY sorted
by trigger time (tables with the columns ``obs_id``, ``event_id``, ``tel_id``, ``event_type`` and
``time_ns``, the trigger time in integer ns). It returns the rows of the coincident events in each
table. Any part of the combination can be changed by overriding its ``combine_<part>`` method, and the
matching of all the telescopes by overriding ``match``. Its configuration options are traits, as for
any ctapipe component:

.. code-block:: python

    import astropy.units as u
    import numpy as np
    from ctapipe.core.traits import AstroQuantity

    from sst1mpipe.stereo import StereoEventCombiner


    class FixedOffsetCombiner(StereoEventCombiner):
        """Events of the other telescope within ``max_time_diff`` of the reference time + ``offset``"""

        offset = AstroQuantity(default_value=0 * u.s, physical_type=u.physical.time).tag(config=True)
        max_time_diff = AstroQuantity(default_value=10 * u.us, physical_type=u.physical.time).tag(config=True)

        def match_pair(self, reference, other):
            time_1 = np.asarray(reference["time_ns"]) + int(self.offset.to_value(u.ns))
            time_2 = np.asarray(other["time_ns"])
            other_index = np.clip(np.searchsorted(time_2, time_1), 0, len(time_2) - 1)
            close = np.abs(time_2[other_index] - time_1) < self.max_time_diff.to_value(u.ns)
            reference_index = np.flatnonzero(close)
            # each event of the other telescope at most once
            _, first = np.unique(other_index[close], return_index=True)
            return reference_index[first], other_index[close][first]

        def combine_index(self, stereo, events):
            """ids of the event of the lowest telescope id instead of the reference telescope"""
            stereo.index = events[min(events)].index

The class must be imported before the tool is configured to be found by name, e.g. with the
Python API:

.. code-block:: python

    from ctapipe.core import run_tool
    from sst1mpipe.scripts.sst1mpipe_stereo_tool import StereoTool

    run_tool(StereoTool(), argv=[
        "--input=SST1M1_20260120_1179_W1.dl1.h5",
        "--input=SST1M2_20260120_1102_W1.dl1.h5",
        "--output=stereo.dl1.h5",
        "--combiner=FixedOffsetCombiner",
        "--FixedOffsetCombiner.offset=3 ms",
    ], raises=True)


Using the combiners in Python
-----------------------------

The combiners can be used without the tool, on trigger tables (e.g. read by
``read_telescope_triggers`` of ``sst1mpipe.scripts.sst1mpipe_stereo_tool``) and events:

.. code-block:: python

    from sst1mpipe.scripts.sst1mpipe_stereo_tool import read_telescope_triggers
    from sst1mpipe.stereo import StereoEventCombiner

    triggers = {
        21: read_telescope_triggers(["SST1M1_20260120_1179_W1.dl1.h5"]),
        22: read_telescope_triggers(["SST1M2_20260120_1102_W1.dl1.h5"]),
    }
    combiner = StereoEventCombiner.from_name("SWATEventIdCombiner")
    groups = combiner.match(triggers)
    # obs_id and event_id of the events of tel 22 in the stereo events
    print(triggers[22][groups[22]][["obs_id", "event_id"]])


API
---

Combiners
~~~~~~~~~

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

.. autofunction:: sst1mpipe.stereo.to_sst1m_event

.. autofunction:: sst1mpipe.stereo.nearest_index

Tool
~~~~

.. autoclass:: sst1mpipe.scripts.sst1mpipe_stereo_tool.StereoTool

.. autofunction:: sst1mpipe.scripts.sst1mpipe_stereo_tool.read_telescope_triggers

.. autoclass:: sst1mpipe.scripts.sst1mpipe_stereo_tool.TelescopeEventReader
   :members: get

.. autoclass:: sst1mpipe.scripts.sst1mpipe_stereo_tool.CombinedRuns

.. autofunction:: sst1mpipe.io.append_monitoring_tables
