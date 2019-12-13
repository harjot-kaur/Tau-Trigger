from array import array
from math import sqrt

class getbinning2017:

	def __init__(self):
		pass

	def getBinning(self):
		edgesdiTau = []
		for i in range( 20, 50, 2) :
			edgesdiTau.append( float(i) )
		for i in range( 50, 75, 5 ) :
			edgesdiTau.append( float(i) )
		for i in range( 75, 100, 25) :
   			edgesdiTau.append( float(i) )
		for i in range(100, 150, 50 ) :
			edgesdiTau.append( float(i) )
		for i in range(150, 500, 300 ) :
 	  	 	edgesdiTau.append( float(i) )
 	  	 
 	  	edgesMuTau =[]
		for i in range( 20, 50, 2) :
			edgesMuTau.append( float(i) )
		for i in range( 50, 75, 5 ) :
			edgesMuTau.append( float(i) )
		for i in range( 75, 100, 25) :
   			edgesMuTau.append( float(i) )
		for i in range(100, 150, 50 ) :
			edgesMuTau.append( float(i) )
		for i in range(150, 500, 300 ) :
 	  	 	edgesMuTau.append( float(i) )
 	  	 	
 	  	edgesETau =[]
		for i in range( 20, 24, 1) :
			edgesETau.append( float(i) )
		for i in range( 24, 30, 3) :
			edgesETau.append( float(i) )
		for i in range( 30, 55, 5 ) :
			edgesETau.append( float(i) )
		for i in range( 55, 70, 15 ) :
			edgesETau.append( float(i) )
		for i in range( 70, 90, 20) :
			edgesETau.append( float(i) )
		for i in range(90, 160, 50 ) :
			edgesETau.append( float(i) )
		for i in range(160, 600, 350 ) :
			edgesETau.append( float(i) )
			
		edgesdict ={"ditau":edgesdiTau, "mutau":edgesMuTau, "etau":edgesETau}
 	 	return edgesdict
 	 	
 	def getBinningDM(self):
 	 	edgesDiTauDM0 = []
 	 	for i in range( 20, 30, 2) :
			edgesDiTauDM0.append( float(i) )
		for i in range( 30, 55, 5) :
			edgesDiTauDM0.append( float(i) )
		for i in range( 55, 75, 20 ) :
			edgesDiTauDM0.append( float(i) )
		for i in range( 75, 112, 37) :
			edgesDiTauDM0.append( float(i) )
		for i in range(112, 500, 300 ) :
			edgesDiTauDM0.append( float(i) )

                edgesDiTauDM1 = []	
		for i in range( 20, 30, 2) :
			edgesDiTauDM1.append( float(i) )
		for i in range( 30, 55, 5) :
			edgesDiTauDM1.append( float(i) )
		for i in range( 55, 75, 10 ) :
			edgesDiTauDM1.append( float(i) )
		for i in range( 75, 90, 25) :
			edgesDiTauDM1.append( float(i) )
		for i in range(90, 112, 22 ) :
			edgesDiTauDM1.append( float(i) )
		for i in range(112, 500, 300 ) :
			edgesDiTauDM1.append( float(i) )
				
		edgesDiTauDM10Data = []	
		for i in range(20, 30, 2) :
			edgesDiTauDM10Data.append( float(i) ) 
		for i in range(30, 50, 5) :
			edgesDiTauDM10Data.append( float(i) )
		for i in range(50, 70, 10 ) :
			edgesDiTauDM10Data.append( float(i) )
		for i in range(70, 85, 15 ) :
			edgesDiTauDM10Data.append( float(i) )
		for i in range(85, 115, 30 ) :
			edgesDiTauDM10Data.append( float(i) )
		for i in range(115, 500, 300 ) :
			edgesDiTauDM10Data.append( float(i) )

                edgesDiTauDM10MC = []				
		for i in range(20, 30, 2) :
			edgesDiTauDM10MC.append( float(i) )
		for i in range(30, 50, 5) :
			edgesDiTauDM10MC.append( float(i) )
		for i in range(50, 58, 8 ) :
			edgesDiTauDM10MC.append( float(i) )
		for i in range(58, 78, 20 ) :
			edgesDiTauDM10MC.append( float(i) )
		for i in range(78, 120, 42 ) :
			edgesDiTauDM10MC.append( float(i) )
		for i in range(120, 500, 300 ) :
			edgesDiTauDM10MC.append( float(i) )

                edgesDiTauDM0 = {"DATA":edgesDiTauDM0, "MC": edgesDiTauDM0}
                edgesDiTauDM1 = {"DATA":edgesDiTauDM1, "MC": edgesDiTauDM1}
                edgesDiTauDM10 = {"DATA":edgesDiTauDM10Data, "MC": edgesDiTauDM10MC}

                edgesDiTauDM = {"dm0": edgesDiTauDM0, "dm1": edgesDiTauDM1, "dm10": edgesDiTauDM10}

                edgesMuTauDM0 = []
 	 	for i in range( 20, 30, 2) :
			edgesMuTauDM0.append( float(i) )
		for i in range( 30, 55, 5) :
			edgesMuTauDM0.append( float(i) )
		for i in range( 55, 65, 10 ) :
			edgesMuTauDM0.append( float(i) )
		for i in range( 65, 80, 15) :
			edgesMuTauDM0.append( float(i) )
		for i in range(80, 110, 30 ) :
			edgesMuTauDM0.append( float(i) )
		for i in range(110, 500, 300 ) :
			edgesMuTauDM0.append( float(i) )

                edgesMuTauDM1 = []
		for i in range( 20, 30, 2) :
			edgesMuTauDM1.append( float(i) )
		for i in range( 30, 55, 5) :
			edgesMuTauDM1.append( float(i) )
		for i in range( 55, 65, 10 ) :
			edgesMuTauDM1.append( float(i) )
		for i in range( 65, 80, 15) :
			edgesMuTauDM1.append( float(i) )
		for i in range(80, 110, 30 ) :
			edgesMuTauDM1.append( float(i) )
		for i in range(110, 500, 300 ) :
			edgesMuTauDM1.append( float(i) )

                edgesMuTauDM10 = []
		for i in range(20, 24, 1) :
			edgesMuTauDM10.append( float(i) )
		for i in range(24, 30, 2) :
			edgesMuTauDM10.append( float(i) )
		for i in range(30, 40, 5 ) :
			edgesMuTauDM10.append( float(i) )
		for i in range(40, 58, 9 ) :
			edgesMuTauDM10.append( float(i) )
		for i in range(58, 80, 22 ) :
			edgesMuTauDM10.append( float(i) )
		for i in range(80, 115, 35) :
			edgesMuTauDM10.append( float(i) )
		for i in range(115,  500, 300 ) :
			edgesMuTauDM10.append( float(i) )

                edgesMuTauDM0 = {"DATA":edgesMuTauDM0, "MC": edgesMuTauDM0}
                edgesMuTauDM1 = {"DATA":edgesMuTauDM1, "MC": edgesMuTauDM1}
                edgesMuTauDM10 = {"DATA":edgesMuTauDM10, "MC": edgesMuTauDM10}

                edgesMuTauDM = {"dm0": edgesMuTauDM0, "dm1": edgesMuTauDM1, "dm10": edgesMuTauDM10}

                edgesETauDM0 = []
 	 	for i in range( 20, 24, 2) :
			edgesETauDM0.append( float(i) )
		for i in range( 24, 30, 3) :
			edgesETauDM0.append( float(i) )
		for i in range( 30, 40, 5) :
			edgesETauDM0.append( float(i) )
		for i in range( 40, 60, 10 ) :
			edgesETauDM0.append( float(i) )
		for i in range( 60, 75, 15) :
			edgesETauDM0.append( float(i) )
		for i in range(75, 95, 20 ) :
			edgesETauDM0.append( float(i) )
		for i in range(95, 500, 300 ) :
			edgesETauDM0.append( float(i) ) 

                edgesETauDM1 = []
		for i in range( 20, 30, 2) :
			edgesETauDM1.append( float(i) )
		for i in range( 30, 40, 5) :
			edgesETauDM1.append( float(i) )
		for i in range( 40, 55, 15 ) :
			edgesETauDM1.append( float(i) )
		for i in range( 55, 75, 20) :
			edgesETauDM1.append( float(i) )
		for i in range(75, 125, 50 ) :
			edgesETauDM1.append( float(i) )
		for i in range(125, 500, 350 ) :
			edgesETauDM1.append( float(i) )

                edgesETauDM10 = []
		for i in range( 20, 22, 1) :
			edgesETauDM10.append( float(i) )
		for i in range( 22, 30, 2) :
			edgesETauDM10.append( float(i) )
		for i in range( 30, 55, 5 ) :
			edgesETauDM10.append( float(i) )
		for i in range( 55, 75, 20 ) :
			edgesETauDM10.append( float(i) )
		for i in range(75, 100, 25 ) :
			edgesETauDM10.append( float(i) )
		for i in range(100, 600, 300 ) :
			edgesETauDM10.append( float(i) )

                edgesETauDM0 = {"DATA":edgesETauDM0, "MC": edgesETauDM0}
		edgesETauDM1 = {"DATA":edgesETauDM1, "MC": edgesETauDM1}
		edgesETauDM10 = {"DATA":edgesETauDM10, "MC": edgesETauDM10}

		edgesETauDM = {"dm0": edgesETauDM0, "dm1": edgesETauDM1, "dm10": edgesETauDM10}
 	 	edgesDM =  {"ditau": edgesDiTauDM, "mutau": edgesMuTauDM, "etau": edgesETauDM}
 	 	return edgesDM

