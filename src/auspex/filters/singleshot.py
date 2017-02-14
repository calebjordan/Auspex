# Copyright 2016 Raytheon BBN Technologies
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0

import numpy as np

from .filter import Filter
from auspex.parameter import Parameter, FloatParameter, IntParameter, BoolParameter
from auspex.stream import DataStreamDescriptor, InputConnector, OutputConnector
from auspex.log import logger

class SingleShot(Filter):
    
    sink      = InputConnector()
    fidelity  = OutputConnector()
    histogram = OutputConnector()
    kernel    = OutputConnector()

    log_regression       = BoolParameter(default=False)
    zero_mean            = BoolParameter(default=False)
    opt_integration_time = BoolParameter(default=True)

    """Find the optimally matched filter for integration the results."""
    def __init__(self, **kwargs):
        super(SingleShot, self).__init__(**kwargs)
        if len(kwargs) > 0:
            self.log_regression.value = kwargs['log_regression']
            self.zero_mean.value = kwargs['zero_mean']

        self.quince_parameters = [self.log_regression, self.opt_integration_time, self.zero_mean]

    def update_descriptors(self):
        logger.debug('Updating Single "%s" descriptors based on input descriptor: %s.', self.name, self.sink.descriptor)
        
        # We will always be averaging at the round-robin level, and integrating
        # along the time axis. Both of these axes will be destroyed.

        descriptor_in = self.sink.descriptor
        names = [a.name for a in descriptor_in.axes]

        self.axis_num  = names.index("round_robins")
        self.data_dims = descriptor_in.data_dims()
        self.avg_dims  = self.data_dims[self.axis_num+1:]

        # If we get multiple final average simultaneously
        self.reshape_dims = [-1] + self.data_dims[self.axis_num:]
        self.mean_axis = self.axis_num - len(self.data_dims)

        # This should collapse along all of the internal axes
        # Including time, round_robins, and segments
        fidelity_descriptor  = descriptor_in.copy()
        self.num_time_points = fidelity_descriptor.pop_axis("time").num_points()
        self.num_averages    = fidelity_descriptor.pop_axis("round_robins").num_points()

        # The kernel will have a time axis
        kernel_descriptor = descriptor_in.copy()
        kernel_descriptor.pop_axis("round_robins")

        # The histograms have a relatively different set of axes
        kernel_descriptor = descriptor_in.copy()zz

        # Set output descriptors
        self.fidelity.set_descriptor(fidelity_descriptor)
        self.histogram.set_descriptor(histogram_descriptor)
        self.kernel.set_descriptor(kernel_descriptor)

    def final_init(self):
        self.kernel_sum = np.zeros()

    async def process_data(self, data):
        # Assume we get entire shots at a time
        pass


        # Wait until we have all of the data before proceeding
