from argparse import Namespace, ArgumentParser


class Subapplication(object):
    '''
    Super class for different parts of application
    '''

    def __init__(self, key):
        self.key = key

    def modifyParser(self, parser: ArgumentParser):
        '''
        Adds the subapplication to the default parser

        Parameters
        ----------
        parser: ArgumentParser
            the parser to add the subapplication to
        '''
        parser.set_defaults(application=lambda x: self.run(x))

    def run(self, args: Namespace):
        '''
        Run the subapplication

        Parameters
        ----------
        args: Namespace
            arguments from parser
        '''
        pass