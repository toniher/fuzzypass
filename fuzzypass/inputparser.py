import json

class InputParser:
    
    list_seqs = []
    
    def __init__(self, filepath, filetype="json", program="blastp"):
        self.filepath = filepath
        self.filetype = filetype
        self.program = program
    
    def read(self):
        
        print("Read file")
        print( self.filepath )
        
        if self.filetype == "json" :

            json_struct = {}

            # TODO: Handle if file exists
            with open(self.filepath, 'r') as f:
                # TODO: Handle if not JSON
                json_struct = json.load(f)
            
            if json_struct :
                
                if "BlastOutput2" in json_struct :
                    
                    if "report" in json_struct["BlastOutput2"] :
                        
                        self.list_seqs.append( self.handle_report( json_struct["BlastOutput2"]["report"] ) )

                    else :
                        
                         for i in range(0, len( json_struct["BlastOutput2"] ) ):
                    
                            if "report" in  json_struct["BlastOutput2"][i] :
                                self.list_seqs.append( self.handle_report( json_struct["BlastOutput2"][i]["report"] ) )
                    
                else :
                    
                    print("No BLAST")
        
        print( self.list_seqs )
        
    def handle_report( self, report ):

        iterations = []
        
        if "results" in report :
                
            if "iterations" in report["results"] :
                                
                for i in range(0, len( report["results"]["iterations"] ) ):
                
                    print( "iter "+str( i ) )
                    iterations.append( self.handle_results( report["results"]["iterations"][i] ) )
                
            else :
                iterations.append( self.handle_results( report["results"] ) )
                    
        return iterations
    
    
    def handle_results( self, results ) :
        
        hits = []
        
        if "search" in results :
            
            qtitle = None
            qid = None
            qlen = None
            
            if "query_title" in results["search"] :
                
                qtitle = results["search"]["query_title"]
                
            if "query_id" in results["search"] :

                qid = results["search"]["query_id"]

            if "query_len" in results["search"] :
                
                qlen = int( results["search"]["query_len"] )

            if "hits" in results["search"] :
                
                for h in range(0, len( results["search"]["hits"] ) ) :
                    
                    hit = results["search"]["hits"][h]
                    
                    desc_title = None
                    desc_acc = None
                    desc_id = None
                    desc_taxid = None
                    
                    hit_len = 0
                    
                    if "len" in hit :
                        hit_len = hit["len"]
                    
                    if len( hit["description"] ) > 0 :
                    
                        # We take only the first one
                        desc = hit["description"][0]
                        
                        if "title" in desc :
                            desc_title = desc["title"]
                        
                        if "accession" in desc :
                            desc_acc = desc["accession"]
                        
                        if "id" in desc :
                            desc_id = desc["id"]
                        
                        if "taxid" in desc :
                            desc_taxid = int( desc["taxid"] )


                    if len( hit["hsps"] ) > 0 :
                        
                        # We take only the first one
                        hsp = hit["hsps"][0]
                        
                        hsp_len = 0
                        
                        if "align_len" in hsp :
                            hsp_len = int( hsp["align_len"] )
                            
                        if "qseq" in hsp :
                            
                            hsp_qseq = hsp["qseq"]

                            
                        if "hseq" in hsp :
                            
                            hsp_hseq = hsp["hseq"]

                            
                        if "qstart" in hsp :
                                                       
                            hsp_qstart = int( hsp["query_from"] )

                            
                        if "qend" in hsp :
                            
                            hsp_qend = int( hsp["query_to"] )

                            
                        if "hstart" in hsp :
                            
                            hsp_hstart = int( hsp["hit_from"] )

                        if "hend" in hsp :
                            
                            hsp_hend = int( hsp["hit_to"] )
                            
                        if "bits" in hsp :
                            
                            hsp_bits = float( hsp["bit_score"] )
                            
                        if "score" in hsp :
                            
                            hsp_score = float( hsp["score"] )

                            
                        if "evalue" in hsp :
                            
                            hsp_evalue = float( hsp["evalue"] )

                            
                        if "identical" in hsp :
                            
                            hsp_identical = int( hsp["identity"] )

                            
                        if "conserved" in hsp :
                            
                            hsp_conserved = int( hsp["positive"] )
  
                            
                        if "mseq" in hsp :
                            
                            hsp_mseq = hsp["midline"]

                            
                        if "gaps" in hsp :
                            
                            hsp_gaps = int( hsp["gaps"] )
    
                        if hsp_len == 0 :
                            
                            hsp_len = hsp_hend - hsp_hstart + 1

                    print( h )
                    hits.append( h )
        return hits
