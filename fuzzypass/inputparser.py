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
                
                    iterations.append( self.handle_results( report["results"]["iterations"][i] ) )
                
            else :
                iterations.append( self.handle_results( report["results"] ) )
                    
        return iterations
    
    
    def handle_results( self, results ) :
        
        result = {}
        
        if "iter_num" in results :
            result["iter_num"] = results["iter_num"]
        
        if "search" in results :
            
            if "query_title" in results["search"] :
                
                result["query_title"] = results["search"]["query_title"]
                
            if "query_id" in results["search"] :

                result["query_id"] = results["search"]["query_id"]

            if "query_len" in results["search"] :
                
                result["query_len"] = int( results["search"]["query_len"] )


            if "hits" in results["search"] :

                result["hits"] = []
                
                for h in range(0, len( results["search"]["hits"] ) ) :
                    
                    hstruct = {}

                    hstruct["num"] = h + 1
                    
                    hit = results["search"]["hits"][h]
                    
                    result["hit_len"] = 0
                    
                    if "len" in hit :
                        result["hit_len"] = int( hit["len"] )
                    
                    hstruct["description"] = []
                    
                    if len( hit["description"] ) > 0 :
                                        
                        # We take only the first one
                        desc = hit["description"][0]
                        
                        dstruct = {}

                        if "title" in desc :
                            dstruct["title"] = desc["title"]
                        
                        if "accession" in desc :
                            dstruct["accession"] = desc["accession"]

                        if "id" in desc :
                            dstruct["id"] = desc["id"]
                        
                        if "taxid" in desc :
                            dstruct["taxid"] = int( desc["taxid"] )
                            
                        hstruct["description"].append( dstruct )


                    if len( hit["hsps"] ) > 0 :
                        
                        # We take only the first one
                        hsp = hit["hsps"][0]
                        
                        hstruct["length"] = 0
                        
                        if "align_len" in hsp :
                            hstruct["length"] = int( hsp["align_len"] )
                            
                        if "qseq" in hsp :
                            
                            hstruct["qseq"] = hsp["qseq"]

                            
                        if "hseq" in hsp :
                            
                           hstruct["hseq"] = hsp["hseq"]

                            
                        if "qstart" in hsp :
                                                       
                            hstruct["qstart"] = int( hsp["query_from"] )

                            
                        if "qend" in hsp :
                            
                            hstruct["qend"] = int( hsp["query_to"] )

                            
                        if "hstart" in hsp :
                            
                            hstruct["hstart"] = int( hsp["hit_from"] )

                        if "hend" in hsp :
                            
                            hstruct["hend"] = int( hsp["hit_to"] )
                            
                        if "bits" in hsp :
                            
                            hstruct["bits"] = float( hsp["bit_score"] )
                            
                        if "score" in hsp :
                            
                            hstruct["score"] = float( hsp["score"] )

                            
                        if "evalue" in hsp :
                            
                            hstruct["evalue"] = float( hsp["evalue"] )

                            
                        if "identical" in hsp :
                            
                            hstruct["identical"] = int( hsp["identity"] )

                            
                        if "conserved" in hsp :
                            
                            hstruct["conserved"] = int( hsp["positive"] )
  
                            
                        if "mseq" in hsp :
                            
                            hstruct["mseq"] = hsp["midline"]

                            
                        if "gaps" in hsp :
                            
                            hstruct["gaps"] = int( hsp["gaps"] )
    
                        if hstruct["length"] == 0 :
                            
                            hstruct["length"] = hstruct["hend"] - hstruct["hstart"] + 1

                    result["hits"].append( hstruct )
                    
        return result
