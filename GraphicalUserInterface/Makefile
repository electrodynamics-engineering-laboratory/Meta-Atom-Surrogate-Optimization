#Author: Joe Haun
#Assignment: EEL Metamodeling Suite
#File: Makefile
#Date: June 30, 2019
#Update: June 30, 2019

FLAGS = `wx-config --cxxflags --libs all`

all: EELMetamodelGUI.x

EELMetamodelGUI.x: ElectrodynamicsEngineeringLaboratoryGUI.cpp
	g++ ElectrodynamicsEngineeringLaboratoryGUI.cpp $(FLAGS) -o EELMetamodelGUI.x

test: TestFile.cpp
	g++ TestFile.cpp $(FLAGS) -o Test.x
	./Test.x

clean:
	rm *~ *.o *.x *#

run: EELMetamodelGUI.x
	./EELMetamodelGUI.x
