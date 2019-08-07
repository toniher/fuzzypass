from . import metrics

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
                            
                            if self.list_fuzz[r][i][h]["Dis"] > DistMax :
                                DistMax = self.list_fuzz[r][i][h]["Dis"]

                        
            print( self.list_fuzz )
        
	# for (int iter = 0; iter < numiter; ++iter ) {
	# 	
	# 	fuzzvalues.push_back(vector<pair<int, int> >());
	# 
	# 	//Creates Fuzzy Structure
	# 	FC = new TFunctionFC();
	# 
	# 	double MaxVal, DistMax, KDMax, FlxMax;
	# 	MaxVal=0;
	# 	DistMax=0;
	# 	KDMax=0;
	# 	FlxMax=0;
	# 
	# 	int numprot = Prote[iter].size();
	# 	
	# 	for (int i=0; i<numprot; i++) {
	# 
	# 		// Turn to capital letters seqs
	# 		string qseq = Prote[iter][i].qseq;
	# 		string hseq = Prote[iter][i].hseq;
	# 		transform(qseq.begin(), qseq.end(),qseq.begin(), ::toupper);
	# 		transform(hseq.begin(), hseq.end(),hseq.begin(), ::toupper);
	# 		replace( qseq.begin(), qseq.end(), '.', '-');
	# 		replace( hseq.begin(), hseq.end(), '.', '-');
	# 		replace( qseq.begin(), qseq.end(), ' ', '-');
	# 		replace( hseq.begin(), hseq.end(), ' ', '-');
	# 
	# 		if ( hseq.length() < minAlnSize ) {
	# 			Prote[iter][i].Dis = discardedVal;
	# 		} else {
	# 		
	# 			DistPearson( qseq, hseq, Prote[iter][i].Dis, Prote[iter][i].r);
	# 
	# 			KDF((qseq).c_str(), (hseq).c_str(), Prote[iter][i].KD, Prote[iter][i].Flx);
	# 
	# 			if (Prote[iter][i].length > MaxVal) {MaxVal= Prote[iter][i].length;}
	# 			if (Prote[iter][i].Dis > DistMax) {DistMax=Prote[iter][i].Dis;}
	# 			if (Prote[iter][i].KD > KDMax) {KDMax= Prote[iter][i].KD;}
	# 			if (Prote[iter][i].Flx > FlxMax) {FlxMax= Prote[iter][i].Flx;}
	# 		}
	# 
	# 	}
	# 	
	# 	for (int i=0; i<numprot; i++) {
	# 	
	# 		if ( Prote[iter][i].Dis == discardedVal ) {
	# 			Prote[iter][i].Fuz = 0;
	# 			Prote[iter][i].Flx = discardedVal;
	# 			Prote[iter][i].KD = discardedVal;
	# 		} else {
	# 		
	# 			double Len = 100*Prote[iter][i].length/MaxVal;
	# 			double EucDist = 100*Prote[iter][i].Dis/DistMax;
	# 			double Transm = 100-(100*Prote[iter][i].Flx/FlxMax);
	# 			double KD = 100*Prote[iter][i].KD/KDMax;
	# 
	# 			if ( isnan(EucDist) || isnan(Transm) || isnan(KD) ) {
	# 				Prote[iter][i].Fuz = 0; // NO FUZZY
	# 			} else {
	# 				Prote[iter][i].Fuz = OutPut( Len, EucDist, Transm, KD, JSONConf);
	# 			}
	# 		}
	# 	
	# 		fuzzvalues[iter].push_back( make_pair( i + 1, Prote[iter][i].Fuz ) );
	# 
	# 	}
	# 			
	# 	delete FC;
	# 
	# }