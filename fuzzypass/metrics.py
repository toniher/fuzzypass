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

