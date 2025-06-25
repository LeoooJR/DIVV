from argparse import Action

class ValidateProcessAction(Action):

    def __call__(self, parser, namespace, values, option_string = None):
        
        if values < 0:

            raise parser.error("Number of processes alocated must be an positive unsigned integer.")
        
        setattr(namespace, self.dest, values)