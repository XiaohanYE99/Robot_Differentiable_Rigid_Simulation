#Local-Global
if [ ! -d "Ripple" ] 
then
	./mainDeformedEnvironment range=5 resolution=10,10,2 expression='sin(x)*sin(y)' name=Ripple wACAP=0.1 > Ripple.dat
fi
if [ ! -d "Cone" ] 
then
	./mainDeformedEnvironment range=5 resolution=10,10,2 expression='(x*x+y*y)/10' name=Cone wACAP=0.1 > Cone.dat
fi
if [ ! -d "Stair" ] 
then
	./mainDeformedEnvironment range=5 resolution=10,10,2 expression='if(x<0)0;else 2;' name=Stair wACAP=0.1 > Stair.dat
fi
if [ ! -d "Hill" ] 
then
	./mainDeformedEnvironment range=5 resolution=10,10,2 expression='if(x<-2 or x>2 or y<-2 or y>2)0;else 1;' name=Hill wACAP=0.1 > Hill.dat
fi
if [ ! -d "Hill2" ] 
then
	./mainDeformedEnvironment range=5 resolution=10,10,2 expression='if(x<-2 or x>2)0;else 2;' name=Hill2 wACAP=0.1 > Hill2.dat
fi
if [ ! -d "Pit" ] 
then
	./mainDeformedEnvironment range=5 resolution=10,10,2 expression='if(x<-2 or x>2 or y<-2 or y>2)0;else -2;' name=Pit wACAP=0.1 > Pit.dat
fi

#Newton
if [ ! -d "RippleNewton" ] 
then
	./mainDeformedEnvironment range=5 resolution=10,10,2 expression='sin(x)*sin(y)' name=RippleNewton wACAP=0.1 useLocalGlobal=0 > RippleNewton.dat
fi
if [ ! -d "ConeNewton" ] 
then
	./mainDeformedEnvironment range=5 resolution=10,10,2 expression='(x*x+y*y)/10' name=ConeNewton wACAP=0.1 useLocalGlobal=0 > ConeNewton.dat
fi
if [ ! -d "StairNewton" ] 
then
	./mainDeformedEnvironment range=5 resolution=10,10,2 expression='if(x<0)0;else 2;' name=StairNewton wACAP=0.1 useLocalGlobal=0 > StairNewton.dat
fi
if [ ! -d "HillNewton" ] 
then
	./mainDeformedEnvironment range=5 resolution=10,10,2 expression='if(x<-2 or x>2 or y<-2 or y>2)0;else 1;' name=HillNewton wACAP=0.1 useLocalGlobal=0 > HillNewton.dat
fi
if [ ! -d "Hill2Newton" ] 
then
	./mainDeformedEnvironment range=5 resolution=10,10,2 expression='if(x<-2 or x>2)0;else 2;' name=Hill2Newton wACAP=0.1 useLocalGlobal=0 > Hill2Newton.dat
fi
if [ ! -d "PitNewton" ] 
then
	./mainDeformedEnvironment range=5 resolution=10,10,2 expression='if(x<-2 or x>2 or y<-2 or y>2)0;else -2;' name=PitNewton wACAP=0.1 useLocalGlobal=0 > PitNewton.dat
fi
