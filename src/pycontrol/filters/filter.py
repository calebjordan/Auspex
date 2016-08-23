import asyncio
import zlib
import pickle
from concurrent.futures import FIRST_COMPLETED

from pycontrol.stream import DataStream, InputConnector, OutputConnector
from pycontrol.logging import logger

class MetaFilter(type):
    """Meta class to bake the input/output connectors into a Filter class description
    """
    def __init__(self, name, bases, dct):
        type.__init__(self, name, bases, dct)
        logger.debug("Adding connectors to %s", name)
        self._input_connectors  = []
        self._output_connectors = []
        for k,v in dct.items():
            if isinstance(v, InputConnector):
                logger.debug("Found '%s' input connector.", k)
                self._input_connectors.append(k)
            elif isinstance(v, OutputConnector):
                logger.debug("Found '%s' output connector.", k)
                self._output_connectors.append(k)

class Filter(metaclass=MetaFilter):
    """Any node on the graph that takes input streams with optional output streams"""
    def __init__(self, name=None):
        self.name = name
        self.input_connectors = {}
        self.output_connectors = {}

        for ic in self._input_connectors:
            a = InputConnector(name=ic, parent=self)
            a.parent = self
            self.input_connectors[ic] = a
            setattr(self, ic, a)
        for oc in self._output_connectors:
            a = OutputConnector(name=oc, parent=self)
            a.parent = self
            self.output_connectors[oc] = a
            setattr(self, oc, a)

    def __repr__(self):
        return "<Filter(name={})>".format(self.name)

    def update_descriptors(self):
        self.descriptor = list(self.input_connectors.values())[0].descriptor
        logger.debug("Starting descriptor update in filter %s, where the descriptor is %s",
                self.name, self.descriptor)
        for oc in self.output_connectors.values():
            oc.descriptor = self.descriptor
            oc.update_descriptors()

    async def run(self):
        """
        Generic run method which waits on a single stream and calls `process_data` on any new_data
        """
        logger.debug('Running "%s" run loop', self.name)

        input_stream = getattr(self, self._input_connectors[0]).input_streams[0]

        while True:

            message = await input_stream.queue.get()
            message_type = message['type']
            message_data = message['data']
            message_comp = message['compression']
            
            if message_comp == 'zlib':
                message_data = pickle.loads(zlib.decompress(message_data))
            # If we receive a message
            if message['type'] == 'event':
                logger.debug('%s "%s" received event "%s"', self.__class__.__name__, self.name, message_data)
                
                # Propagate along the graph
                for oc in self.output_connectors.values():
                    for os in oc.output_streams:
                        logger.debug('%s "%s" pushed event "%s" to %s, %s', self.__class__.__name__, self.name, message_data, oc, os)
                        await os.queue.put(message)
                
                # Check to see if we're done
                if message['data'] == 'done':
                    break

            elif message['type'] == 'data':
                logger.debug('%s "%s" received %d points.', self.__class__.__name__, self.name, message_data.size)
                logger.debug("Now has %d of %d points.", input_stream.points_taken, input_stream.num_points())
                await self.process_data(message_data)

    async def process_data(self, data):
        """Generic pass through.  """
        return data
