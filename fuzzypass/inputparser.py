
class InputParser:
    
    def __init__(self, filepath, filetype, program):
        self.filetype = filetype
        self.program = program
        self.filepath = filepath
    
    def read(self):
        
        print("Read file")
