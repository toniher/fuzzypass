from . import metrics
import math

class FuzzyFier:
    
    list_fuzz = [  ]
    discarded_val = -1
    min_aln_size = 5
    
    def __init__(self, config):
        print("Fuzzifier")

    def calculate(self, list_seqs):
        
        for r in range(0, len( list_seqs ) ):
            
            self.list_fuzz.append( [] )
            
            for i in range(0, len( list_seqs[r] ) ):
                    
                MaxVal=0;
                DistMax=0;
                KDMax=0;
                FlxMax=0;
                
                self.list_fuzz[r].append( [] )

                if "hits" in list_seqs[r][i] :
                
                    for h in range(0, len( list_seqs[r][i]["hits"] ) ):
    
                        self.list_fuzz[r][i].append( {} )
    
                        qseq = list_seqs[r][i]["hits"][h]["qseq"]
                        hseq = list_seqs[r][i]["hits"][h]["hseq"]
    
                        self.list_fuzz[r][i][h]["length"] = list_seqs[r][i]["hits"][h]["hit_len"]
    
                        qseq = qseq.replace( ".", "-" )
                        qseq = qseq.replace( " ", "-" )
                        hseq = hseq.replace( ".", "-" )
                        hseq = hseq.replace( " ", "-" )
                        
                        if len( hseq ) < self.min_aln_size :
                            self.list_fuzz[r][i][h]["Dis"] = self.discarded_val
                        else :
                            pearson = metrics.DistPearson( qseq, hseq )
                            pearson.calculate()
                            self.list_fuzz[r][i][h]["Dis"] = pearson.distance
                            
                            KDF = metrics.KDF( qseq, hseq )
                            KDF.calculate()
                            self.list_fuzz[r][i][h]["KD"] = KDF.KD
                            self.list_fuzz[r][i][h]["Flx"] = KDF.Flx
                         
                            if self.list_fuzz[r][i][h]["length"] > MaxVal :
                                MaxVal = self.list_fuzz[r][i][h]["length"]
                                
                            if self.list_fuzz[r][i][h]["Dis"] > DistMax :
                                DistMax = self.list_fuzz[r][i][h]["Dis"]

                            if self.list_fuzz[r][i][h]["KD"] > KDMax :
                                KDMax = self.list_fuzz[r][i][h]["KD"]
                                
                            if self.list_fuzz[r][i][h]["Flx"] > FlxMax :
                                FlxMax = self.list_fuzz[r][i][h]["Flx"]            
        
        for r in range(0, len( self.list_fuzz ) ):
            
            for i in range(0, len( self.list_fuzz[r] ) ):

                for h in range(0, len( self.list_fuzz[r][i] ) ):

                    self.list_fuzz[r][i][h]["Fuz"] = 0

                    if self.list_fuzz[r][i][h]["Dis"] == self.discarded_val :
                        self.list_fuzz[r][i][h]["Fuz"] = 0
                        self.list_fuzz[r][i][h]["Flx"] = self.discarded_val
                        self.list_fuzz[r][i][h]["KD"] = self.discarded_val

                    else :
                    
                        Len = 100*self.list_fuzz[r][i][h]["length"]/MaxVal
                        EucDist = 100*self.list_fuzz[r][i][h]["Dis"]/DistMax
                        Transm = 100-(100*self.list_fuzz[r][i][h]["Flx"]/FlxMax)
                        KD = 100*self.list_fuzz[r][i][h]["KD"]/KDMax
            
                        if math.isnan(EucDist) or math.isnan(Transm) or math.isnan(KD) :
                            self.list_fuzz[r][i][h]["Fuz"] = 0
                        else :
                            self.go_fuzzy( Len, EucDist, Transm, KD )

    def go_fuzzy( self, Len, Eucdist, Transm, KD ) :
        
        print( "go_fuzzy" )

        
