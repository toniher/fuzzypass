import numpy as np

class DistPearson:
    
    distance = -1
    rcoeff = -1
    
    def __init__(self, ori=[], other=[], seqtype="prot"):

        self.ori = ori
        self.other = other
        self.seqtype = seqtype
        
        self.numdiv = 22;
        self.orddiv = np.asarray( list( "-ACDEFGHIKLMNPQRSTVWYX") )

        self.X = 0
        self.Y = 0
        self.X2 = 0
        self.Y2 = 0
        self.XY = 0
        
        self.ori_acu = np.zeros( (self.numdiv,), dtype=int )
        self.other_acu = np.zeros( (self.numdiv,), dtype=int )


    def calculate( self ):
        
        ori = np.array(["A", "B", "C"])
        other = np.array(["Z", "Y", "X"])
          
        it = np.nditer(ori, flags=['f_index'])
        while not it.finished:
            
            loc = np.where(self.orddiv == it[0])
            print( loc )
            if len( loc ) > 0 :
                self.ori_acu[ loc[0] ] = self.ori_acu[ loc[0] ] + 1
            it.iternext()
        
        it = np.nditer(other, flags=['f_index'])
        while not it.finished:
            print("%s <%d>" % (it[0], it.index), end=' ')

            loc = np.where(self.orddiv == it[0])
            print( loc )
            if len( loc ) > 0 :
                self.other_acu[ loc[0] ] = self.other_acu[ loc[0] ] + 1
                
            it.iternext()
        
        if [ 11 > 10 ] :
            
            for x in range(1, 22):
                
                print( self.orddiv[x] )
                
                self.distance = self.distance + pow( self.other_acu[x]-self.ori_acu[x], 2 )
                self.X = self.X + self.other_acu[x];
                self.Y = self.Y + self.ori_acu[x];
                self.X2 = self.X2 + pow(self.other_acu[x], 2);
                self.Y2 = self.Y2 + pow(self.ori_acu[x], 2);
                self.XY = self.XY + ( self.other_acu[x] * self.ori_acu[x] )
        
        
            self.rcoeff = ((20*self.X2)-pow(self.X,2))*((20*self.Y2)-pow(self.Y,2))
            
            if self.rcoeff > 0 :
                self.rcoeff = ((20*self.XY)-(self.X*self.Y))/pow(self.rcoeff,0.5)

        
        print( self.ori_acu )
        print( self.other_acu )

class KDF:
    
    kd = -1
    flex = -1
    
    def __init__(self, ori=[], other=[], seqtype="prot"):

        self.ori = ori
        self.other = other
        self.seqtype = seqtype
        
        self.AAKD = np.asarray( list( "ARNDCQEGHILKMFPSTWYV" ) );
        self.Rigid = np.asarray( list( "ALHVYIFCWM" ) );
        self.KDp = np.array( [ 0,1.8,-4.5,-3.5,-3.5,2.5,-3.5,-3.5,-0.4,-3.2,4.5,3.8,-3.9,1.9,2.8,-1.6,-0.8,-0.7,-0.9,-1.3,4.2 ] );
        self.AAFlex = np.asarray( list( "ARNDCQEGHILKMFPSTWYV" ) );
        self.Flx1R = np.array( [ 0,0.946,1.028,1.006,1.089,0.878,1.028,1.036,1.042,0.952,0.892,0.961,1.082,0.862,0.912,1.085,1.048,1.051,0.917,0.93,0.927 ] );
        self.Flx0R = np.array( [ 0,1.041,1.038,1.117,1.033,0.96,1.165,1.094,1.142,0.982,1.002,0.967,1.093,0.947,0.93,1.055,1.169,1.073,0.925,0.961,0.982 ] );
        self.Flx2R = np.array( [ 0,0.892,0.901,0.93,0.932,0.925,0.885,0.933,0.923,0.894,0.872,0.921,1.057,0.804,0.914,0.932,0.923,0.934,0.803,0.837,0.913 ] );


    def calculate( self ):

        ori = np.array(["A", "B", "C"])
        other = np.array(["Z", "Y", "X"])
        LQ = length( ori )
        
        
        if LQ > 0 :
        
            FlQ = np.empty( LQ + 1 )
            FlS = np.empty( LQ + 1 )
        
            if np.where(self.Rigid == ori[0]) > 0 :
                
                AAFlexLoc = np.where(self.AAFlex == ori[1])
                FlQ[0] = Flx1R[AAFlexLoc]
                
            else :
                
                AAFlexLoc = np.where(self.AAFlex == ori[1])
                FlQ[0] = Flx0R[AAFlexLoc]
                
        
            if np.where(self.Rigid == other[0]) > 0 :
                
                AAFlexLoc = np.where(self.AAFlex == ori[1])
                FlS[0] = Flx1R[AAFlexLoc]
                
            else :
                
                AAFlexLoc = np.where(self.AAFlex == ori[1])
                FlS[0] = Flx0R[AAFlexLoc]
        
            if np.where(self.Rigid == ori[ LQ - 2 ]) > 0 :
                
                AAFlexLoc = np.where(self.AAFlex == ori[ LQ - 1 ])
                FlQ[ LQ - 1 ] = Flx1R[AAFlexLoc]
                
            else :
                
                AAFlexLoc = np.where(self.AAFlex == ori[ LQ - 1 ])
                FlQ[ LQ - 1 ] = Flx0R[AAFlexLoc]
                
            
            if np.where(self.Rigid == other[ LQ - 2 ]) > 0 :
                
                AAFlexLoc = np.where(self.AAFlex == ori[ LQ - 1 ])
                FlS[ LQ - 1 ] = Flx1R[AAFlexLoc]
                
            else :
                
                AAFlexLoc = np.where(self.AAFlex == other[ LQ - 1 ])
                FlS[ LQ - 1 ] = Flx0R[AAFlexLoc]  
                
        
        print("Calculate")
        
        

# 
# //Kye-Dolittle Function
# //---------------------------------------------------------------------------
# void KDF(string OriProt, string OtherProt, double &KD, double &Flx){
# 
# 	int i, LQ, Rig;
# 	double Dif, *FlQ, *FlS, Query, Subje;
# 	LQ   = string(OriProt).length();
# 	Dif  = 0;
# 
# 	if(LQ>0){
# 	
# 			FlQ = (double *) calloc(LQ+1, sizeof(double));
# 			FlS = (double *) calloc(LQ+1, sizeof(double));
# 	
# 			if(Rigid.find(OriProt.substr(0,1))>0){
# 		
# 			   FlQ[0] = Flx1R[AAFlex.find(OriProt.substr(1,1))];
# 			}
# 	
# 	else{
# 			   FlQ[0] = Flx0R[AAFlex.find(OriProt.substr(1,1))];
# 			}
# 	
# 	if(Rigid.find(OtherProt.substr(0,1))>0){
# 			   FlS[0] = Flx1R[AAFlex.find(OtherProt.substr(1,1))];
# 			}
# 	
# 	else{
# 			   FlS[0] = Flx0R[AAFlex.find(OtherProt.substr(1,1))];
# 		}
# 	
# 			//-----------------------------------------------------------
# 	//
# 			if(Rigid.find(OriProt.substr(LQ-2,1))>0){
# 			   FlQ[LQ-1] = Flx1R[AAFlex.find(OriProt.substr(LQ-1,1))];
# 			}
# 	
# 	else{
# 			   FlQ[LQ-1] = Flx0R[AAFlex.find(OriProt.substr(LQ-1,1))];
# 			}
# 	
# 			if(Rigid.find(OtherProt.substr(LQ-2,1))>0){
# 			   FlS[LQ-1] = Flx1R[AAFlex.find(OtherProt.substr(LQ-1,1))];
# 			}
# 
# 	else{
# 			   FlS[LQ-1] = Flx0R[AAFlex.find(OtherProt.substr(LQ-1,1))];
# 			}
# 
# 			//-----------------------------------------------------------
# 	//
# 			for(i=1;i<(LQ-1);i++){
# 		
# 				//-------------------------------------------------------
# 		//
# 				if(Rigid.find(OtherProt.substr(i,1))>0){
# 					Rig = 1;
# 				}
# 		
# 		else{
# 					Rig = 0;
# 				}
# 		
# 				if(Rigid.find(OtherProt.substr(i+2,1))>0){
# 					Rig++;
# 				}
# 		
# 				switch (Rig) {
# 			
# 					case 0 :
# 						 FlS[i] = Flx0R[AAFlex.find(OtherProt.substr(i+1,1))];
# 						 break;
# 					case 1 :
# 						 FlS[i] = Flx1R[AAFlex.find(OtherProt.substr(i+1,1))];
# 						 break;
# 					case 2 :
# 						 FlS[i] = Flx2R[AAFlex.find(OtherProt.substr(i+1,1))];
# 						 break;
# 				}
# 		
# 				//-------------------------------------------------------
# 		//
# 				if(Rigid.find(OriProt.substr(i,1))>0){
# 					Rig = 1;
# 				}
# 		
# 		else{
# 					Rig = 0;
# 				}
# 				
# 		if(Rigid.find(OriProt.substr(i+2,1))>0){
# 					Rig++;
# 				}
# 		
# 				switch (Rig) {
# 					case 0 :
# 						 FlQ[i] = Flx0R[AAFlex.find(OriProt.substr(i+1,1))];
# 						 break;
# 					case 1 :
# 						 FlQ[i] = Flx1R[AAFlex.find(OriProt.substr(i+1,1))];
# 						 break;
# 					case 2 :
# 						 FlQ[i] = Flx2R[AAFlex.find(OriProt.substr(i+1,1))];
# 						 break;
# 				}
# 		
# 				//-------------------------------------------------------
# 				
# 		Dif    = Dif+pow(KDp[AAKD.find(OriProt.substr(i+1,1))]
# 								-KDp[AAKD.find(OtherProt.substr(i+1,1))],2);
# 		}
# 	
# 			KD   = Dif /(double)LQ;
# 			Dif  = 0;
# 	
# 			for(i=3;i<(LQ-3);i++){
# 		
# 				Query = 0.25*FlQ[i-3]+0.5*FlQ[i-2]+0.75*FlQ[i-1]+FlQ[i]+
# 						0.25*FlQ[i+3]+0.5*FlQ[i+2]+0.75*FlQ[i+1];
# 				Subje = 0.25*FlS[i-3]+0.5*FlS[i-2]+0.75*FlS[i-1]+FlS[i]+
# 						0.25*FlS[i+3]+0.5*FlS[i+2]+0.75*FlS[i+1];
# 				Dif   = Dif+pow((Query/4)-(Subje/4),2);
# 			}
# 	
# 			Flx     = Dif /(double)LQ;
# 		free(FlS);
# 		free(FlQ);
# 	}
# }

