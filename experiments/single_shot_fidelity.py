# Copyright 2016-2017 Raytheon BBN Technologies
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0

from auspex.exp_factory import QubitExpFactory
from auspex.experiment import Experiment
import auspex.config as auspex_config

from QGL import *
from QGL import config as QGLconfig
from QGL.BasicSequences.helpers import create_cal_seqs, time_descriptor, cal_descriptor

import os
import json
import numpy as np
import networkx as nx

class QubitFidelityExperimentFactory(object):
   
    @staticmethod
    def create(qubit_name, notebook=False, expname=None):
        
        # We need to create a qubit and the metafile in order to properly 
        # setup the filter graph

        qubit    = QubitFactory(qubit_name)
        filename = "SingleShot/SingleShot"

        sequence = [[Id(qubit), MEAS(qubit)], [X(qubit), MEAS(qubit)]]
        axis_descriptor = [{
            'name': 'state',
            'unit': 'state',
            'points': ["0", "1"],
            'partition': 1
        }]

        seq_files = compile_to_hardware(sequence, fileName=filename, axis_descriptor=axis_descriptor)
        meta_file = os.path.join(QGLconfig.AWGDir, filename + '-meta.json') # This is all we need for now

        with open(auspex_config.instrumentLibFile, 'r') as FID:
            instrument_settings = json.load(FID)

        with open(auspex_config.measurementLibFile, 'r') as FID:
            measurement_settings = json.load(FID)

        with open(auspex_config.sweepLibFile, 'r') as FID:
            sweep_settings = json.load(FID)

        # Create a mapping from qubits to data writers
        qubit_to_writer = {}

        # Use the meta info to modify the other JSON
        with open(meta_file, 'r') as FID:
            meta_info = json.load(FID)

        # Construct a graph of all instruments in order to properly enabled those
        # associated with the meta_file. We only need to use string representations
        # here, not actual filter and instrument objects.

        # Strip any colons
        def strip_conn_name(text):
            if ':' in text:
                return text.split(":")[0]
            return text

        # Graph edges for the measurement filters
        edges = [(strip_conn_name(pars['data_source']), name) for name, pars in measurement_settings["filterDict"].items()]
        dag = nx.DiGraph()
        dag.add_edges_from(edges)

        inst_to_enable = []
        filt_to_enable = []

        # Find any writer endpoints of the receiver channels
        for receiver_text, num_segments in meta_info['receivers'].items():
            dig_name, chan_name = receiver_text.split("-")

            # Enable this digitizer
            inst_to_enable.append(dig_name)

            # Set number of segments in the digitizer
            instrument_settings['instrDict'][dig_name]['nbr_segments'] = num_segments
            
            # Find descendants of the channel selector
            chan_descendants = nx.descendants(dag, chan_name)
            
            # Find endpoints within the descendants
            endpoints = [n for n in chan_descendants if dag.in_degree(n) == 1 and dag.out_degree(n) == 0]
            
            # Find endpoints of various types
            writers = [e for e in endpoints if measurement_settings["filterDict"][e]["x__class__"] == "WriteToHDF5" and
                                               measurement_settings["filterDict"][e]["enabled"]]
            plotters = [e for e in endpoints if measurement_settings["filterDict"][e]["x__class__"] == "Plotter" and
                                                measurement_settings["filterDict"][e]["enabled"]]
            averagers = [e for e in endpoints if measurement_settings["filterDict"][e]["x__class__"] == "Averager" and
                                                 measurement_settings["filterDict"][e]["enabled"]]
            integrators = [e for e in endpoints if measurement_settings["filterDict"][e]["x__class__"] == "KernelIntegrator" and
                                                   measurement_settings["filterDict"][e]["enabled"]]                                
            # The user should only have one writer enabled, otherwise we will be confused.
            if len(writers) > 1:
                raise Exception("More than one viable data writer was found for a receiver channel {}. Please enabled only one!".format(receiver_text))
            if len(writers) == 0:
                raise Exception("No viable data writer was found for receiver channel {}. Please enabled only one!".format(receiver_text))

            # Replace file writers with buffers
            buffers = []
            for w in writers:
                label = measurement_settings["filterDict"][w]["label"]
                buff = {
                        "data_source": measurement_settings["filterDict"][w]["data_source"],
                        "enabled": True,
                        "label": label,
                        "x__class__": "DataBuffer",
                        "x__module__": "MeasFilters"
                        }
                # Remove the writer
                measurement_settings["filterDict"].pop(measurement_settings["filterDict"][w]["label"])
                # Substitute the buffer
                measurement_settings["filterDict"][label] = buff
                # Store buffer name for local use
                buffers.append(label)
            writers = buffers

            # For now we assume a single qubit
            # TODO: have meta info give the relationships of qubits to receivers so we don't need to dig in the channel lib
            with open(auspex_config.channelLibFile, 'r') as FID:
                chan_settings = json.load(FID)
            for chan in chan_settings['channelDict']:
                if 'receiverChan' in chan_settings['channelDict'][chan] and  chan_settings['channelDict'][chan]['receiverChan'] == receiver_text:
                    qubit_to_writer[chan.strip('M-')] = writers[0]

            # Trace back our ancestors
            writer_ancestors = nx.ancestors(dag, writers[0])

            # We will have gotten the digitizer, which should be removed since we're already taking care of it
            writer_ancestors.remove(dig_name)

            # We also need to remove any integrators and averagers, since we want single shots
            # Step back in the graph and remove any such filters
            for filt in averagers + integrators:
                if filt in writer_ancestors:
                    # Disable the filter
                    measurement_settings['filterDict'][filt]['enabled'] = False

                    # Fix the graph: find the incoming connector
                    in_filt   = dag.in_edges(filt)[0][0]
                    out_filts = [e[1] for e in dag.out_edges(filt)]
                    # Remove it from the graph
                    dag.remove_node(filt)
                    # Bridge over the gap
                    dag.add_edges_from([(in_filt, out_filt) for out_filt in out_filts])

                    # Update the connected filters to match
                    for out_filt in out_filts:
                        measurement_settings['filterDict'][out_filt]['data_source'] = measurement_settings['filterDict'][filt]['data_source']

            instrument_settings['instrDict'][dig_name]['nbr_segments'] = num_segments

            plotters = [e for e in endpoints if measurement_settings["filterDict"][e]["x__class__"] == "Plotter" and
                                               measurement_settings["filterDict"][e]["enabled"]]
            if plotters:
                plotter_ancestors = set().union(*[nx.ancestors(dag, pl) for pl in plotters])
                plotter_ancestors.remove(dig_name)
            else:
                plotter_ancestors = []

            filt_to_enable.extend(set().union(writers, writer_ancestors, plotters, plotter_ancestors))

            writer_to_qubit = {v: k for k, v in qubit_to_writer.items()}

            # Disable digitizers and APSs and then build ourself back up with the relevant nodes
            for instr_name in instrument_settings['instrDict'].keys():
                if instrument_settings['instrDict'][instr_name]["x__module__"] in ['instruments.Digitizers', 'instruments.APS', 'instruments.APS2']:
                    instrument_settings['instrDict'][instr_name]['enabled'] = False
            for instr_name in inst_to_enable:
                instrument_settings['instrDict'][instr_name]['enabled'] = True

            for meas_name in measurement_settings['filterDict'].keys():
                measurement_settings['filterDict'][meas_name]['enabled'] = False
            for meas_name in filt_to_enable:
                measurement_settings['filterDict'][meas_name]['enabled'] = True
                #label measurement with qubit name (assuming the convention "M-"+qubit_name)
                # if not calibration and measurement_settings['filterDict'][meas_name]["x__class__"] == "WriteToHDF5":
                #     measurement_settings['filterDict'][meas_name]['groupname'] = writer_to_qubit[meas_name].strip('M-')

            # First enable any instruments and set the sequence files
            for instr_name, seq_file in meta_info['instruments'].items():
                instrument_settings['instrDict'][instr_name]['enabled']  = True
                instrument_settings['instrDict'][instr_name]['seq_file'] = seq_file

            # Set the appropriate sweep
            desc = meta_info["axis_descriptor"]
            sweep_settings["sweepDict"] = {"SegmentSweep": {
                                            "axisLabel": "{} ({})".format(desc[0]["name"], desc[0]["unit"]),
                                            "enabled": True,
                                            "label": "SegmentSweep",
                                            "meta_file": meta_file,
                                            "meta_info": meta_info,
                                            "x__class__": "SegmentNum",
                                            "x__module__": "Sweeps"
                                            }
                                          }

            # Replace the sweep order with just the metafile sweep
            sweep_settings["sweepOrder"] = ["SegmentSweep"]

        class QubitExperiment(Experiment):
            """Experiment with a specialized run method for qubit experiments run via factory below."""
            def init_instruments(self):
                for name, instr in self._instruments.items():
                    instr_par = self.instrument_settings['instrDict'][name]
                    logger.debug("Setting instr %s with params %s.", name, instr_par)
                    instr.set_all(instr_par)

                self.digitizers = [v for _, v in self._instruments.items() if v.instrument_type == "Digitizer"]
                self.awgs       = [v for _, v in self._instruments.items() if v.instrument_type == "AWG"]

                # Swap the master AWG so it is last in the list
                master_awg_idx = next(ct for ct,awg in enumerate(self.awgs) if self.instrument_settings['instrDict'][awg.name]['is_master'])
                self.awgs[-1], self.awgs[master_awg_idx] = self.awgs[master_awg_idx], self.awgs[-1]

                # attach digitizer stream sockets to output connectors
                for chan, dig in self.chan_to_dig.items():
                    socket = dig.get_socket(chan)
                    oc = self.chan_to_oc[chan]
                    self.loop.add_reader(socket, dig.receive_data, chan, oc)

            def shutdown_instruments(self):
                # remove socket readers
                for chan, dig in self.chan_to_dig.items():
                    socket = dig.get_socket(chan)
                    self.loop.remove_reader(socket)
                for name, instr in self._instruments.items():
                    instr.disconnect()

            async def run(self):
                """This is run for each step in a sweep."""

                # Recompile the sequence
                seq_files = compile_to_hardware(self.sequence, fileName=self.filename, axis_descriptor=self.axis_descriptor)
                for dig_name, seq_file in meta_info['instruments']:
                    self._instruments[dig_name].seq_file = seq_file

                # Arm the instruments
                for dig in self.digitizers:
                    dig.acquire()
                for awg in self.awgs:
                    awg.run()

                # Wait for all of the acquisitions to complete
                timeout = 10
                await asyncio.wait([dig.wait_for_acquisition(timeout)
                    for dig in self.digitizers])

                for dig in self.digitizers:
                    dig.stop()
                for awg in self.awgs:
                    awg.stop()

                # hack to try to get plots to finish updating before we exit
                await asyncio.sleep(2)

        experiment = QubitExperiment()
        experiment.instrument_settings  = instrument_settings
        experiment.measurement_settings = measurement_settings
        experiment.sweep_settings       = sweep_settings
        experiment.run_in_notebook = notebook
        experiment.name = expname

        experiment.qubit_to_writer = qubit_to_writer
        
        # Specific to single shot measurements
        experiment.qubit           = qubit
        experiment.filename        = filename
        experiment.axis_descriptor = axis_descriptor
        experiment.sequence        = sequence
        experiment.meta_info       = meta_info

        QubitExpFactory.load_instruments(experiment)
        QubitExpFactory.load_segment_sweeps(experiment)
        QubitExpFactory.load_filters(experiment)
        QubitExpFactory.load_parameter_sweeps(experiment)

        return experiment
