import numpy as np

class DistPearson:
    
    distance = 0
    
    
    def __init__(self, ori=[], other=[], seqtype="prot"):

        self.ori = ori
        self.other = other
        self.seqtype = seqtype
        
        self.numdiv = 22;
        self.orddiv = "-ACDEFGHIKLMNPQRSTVWYX";

        self.X = 0
        self.Y = 0
        self.X2 = 0
        self.Y2 = 0
        self.XY = 0
        
        self.ori_acu = np.zeros( (self.numdiv,), dtype=int )
        self.other_acu = np.zeros( (self.numdiv,), dtype=int )


    def calculate( self ):
        
        '''
        it = np.nditer(a, flags=['f_index'])
        while not it.finished:
            print("%d <%d>" % (it[0], it.index), end=' ')
            it.iternext()
        '''
        return self.distance


'''
void DistPearson(string OriProt, string OtherProt, double &Dis, double &r){

	long OriAcu[22], OtherAcu[22], i, Lon;
	double X, Y, X2, Y2, XY;
	for(i=0;i<NumAA;i++){
		OriAcu[i]   = 0;
		OtherAcu[i] = 0;
	}

	Lon=OriProt.length();

	for(i=0;i<Lon+1;i++){
		OriAcu[AA.find(OriProt.substr(i,1))]++;
	}

	Lon=OtherProt.length();

	for(i=0;i<Lon+1;i++){
		OtherAcu[AA.find(OtherProt.substr(i,1))]++;
	}

	Dis=0;
	X=Y=X2=Y2=XY=0;

	if(Lon>11){
	
		for(i=1;i<NumAA;i++){
			Dis  += pow(OtherAcu[i]-OriAcu[i],2);
			X    += OtherAcu[i];
			Y    += OriAcu[i];
			X2   += pow(OtherAcu[i],2);
			Y2   += pow(OriAcu[i],2);
			XY   += OtherAcu[i]*OriAcu[i];
		}
	
		r = ((20*X2)-pow(X,2))*((20*Y2)-pow(Y,2));
	
		if(r>0){
			r = ((20*XY)-(X*Y))/pow(r,0.5);
		}
	
		else{
			cerr << "Error: overflow calculating R pearson coef.";
			//r = pow(10,99);
		}
			
		Dis=pow(Dis,0.5);
	}
}
'''