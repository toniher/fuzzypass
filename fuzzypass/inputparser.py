
class InputParser:
    
    list_seqs = []
    
    def __init__(self, filepath, filetype="json", program="blastp"):
        self.filepath = filepath
        self.filetype = filetype
        self.program = program
    
    def read(self):
        
        print("Read file")
        print( self.filepath )
        print( len( self.list_seqs ) )
