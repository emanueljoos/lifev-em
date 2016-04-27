/*
 * RossiModel14.cpp
 *
 *  Created on: 13/mag/2014
 *      Author: srossi
 */

#include <lifev/em/solver/activation/activeStressModels/ActiveStressRossiModel14.hpp>
#include<iostream>

namespace LifeV
{

ActiveStressRossiModel14::ActiveStressRossiModel14 (Real beta, Real mu, Real Tmax ) :
    M_coefficientBeta (beta),
    M_coefficientMu (mu),
    M_maximumActiveTenstion (Tmax)
{
}



void
ActiveStressRossiModel14::setParameters( EMData& data )
{
	M_maximumActiveTenstion = data.activationParameter<Real>("MaxActiveTension");
	M_coefficientBeta = data.activationParameter<Real>("ActiveStress_Beta");
	M_coefficientMu = data.activationParameter<Real>("ActiveStress_Mu");
	M_caseArea = data.activationParameter<int>("case_area");
        M_weakessFactor = data.activationParameter<double>("weakness_factor");
	M_xMax = data.activationParameter<double>("x_max");
	M_yMax = data.activationParameter<double>("y_max");
	M_zMax = data.activationParameter<double>("z_max");
	M_xMin = data.activationParameter<double>("x_min");
	M_yMin = data.activationParameter<double>("y_min");
	M_zMin = data.activationParameter<double>("z_min");
	M_aMajor = data.activationParameter<double>("a_major");
	M_bMinor = data.activationParameter<double>("b_minor");
	M_xMid = data.activationParameter<double>("x_midpoint");
	M_yMid = data.activationParameter<double>("y_midpoint");
	
}

void
ActiveStressRossiModel14::setup(EMData& data, const MapEpetra& map)
{
	setParameters(data);
    this->M_fiberActivationPtr.reset ( new vector_Type ( map ) );
    super::setup(data, map);

}

void
ActiveStressRossiModel14::setupActivationPtrs(	vectorPtr_Type& fiberActivationPtr,
												vectorPtr_Type& sheetActivationPtr,
												vectorPtr_Type& normalActivationPtr )
{
	fiberActivationPtr = 	this->M_fiberActivationPtr;
	sheetActivationPtr.reset();
	normalActivationPtr.reset();
}

void
ActiveStressRossiModel14::updateActivation(	vectorPtr_Type& fiberActivationPtr,
											vectorPtr_Type& sheetActivationPtr,
											vectorPtr_Type& normalActivationPtr )
{
	setupActivationPtrs(fiberActivationPtr, sheetActivationPtr, normalActivationPtr);
}

void ActiveStressRossiModel14::solveModel2 (Real& timeStep,boost::shared_ptr<RegionMesh<LinearTetra> > fullMeshPtr)
{

//determine the number of meshPoints:
UInt num_MeshPoints = fullMeshPtr -> numPoints();
std::cout<<"number of Mesh Points  :"<<num_MeshPoints;
//Shows important things about the Mesh:
//fullMeshPtr -> showMe();

//determine the indices in the mesh which lie in the region:
std::vector<UInt> indices;
if (M_caseArea  == 1)
{
for (UInt i = 0;i<num_MeshPoints;i++)
	{
		double x_value = fullMeshPtr -> point(i).x();

		double y_value = fullMeshPtr -> point(i).y();

		double z_value = fullMeshPtr -> point(i).z();
		
		if(( M_xMin<x_value) && (x_value<M_xMax) && (M_yMin<y_value) && (y_value<M_yMax) && (M_zMin<z_value) && (z_value<M_zMax)   )
			{
				 indices.push_back(i);	
								
			}	
	 } 

}

else if (M_caseArea == 2)
{
        for (UInt i=0;i<num_MeshPoints;i++)
        {
                double x_value = fullMeshPtr -> point(i).x();

                double y_value = fullMeshPtr -> point(i).y();

                double z_value = fullMeshPtr -> point(i).z();

                double elliptic_value = (x_value-M_xMid)*(x_value-M_xMid)/(M_aMajor*M_aMajor)+(y_value-M_yMid)*(y_value-M_yMid)/(M_bMinor*M_bMinor);

                        if( (elliptic_value<=1)  && (M_zMin<z_value) && (z_value<M_zMax)   )
                        {
                                 indices.push_back(i);

                        }
        }
}
else{}

int size_of_indicesvector=indices.size();
//solving de Activation equation with the help of explicit Euler
    VectorEpetra rhs ( *( this->M_electroSolution.at(0) ) ); 
    rhs *= rhs.operator > (0.0);
    rhs *= M_coefficientBeta;
    rhs -= 2.0 * *super::M_fiberActivationPtr;
    rhs *= (timeStep * M_maximumActiveTenstion / M_coefficientMu);
  //  //rhs -= 2.0 * M_maximumActiveTenstion * *super::M_fiberActivationPtr;
  //  //    //rhs *= (timeStep / M_coefficientMu); 
  //  //        *super::M_fiberActivationPtr += rhs;
   
//weak the found region by multiplying the corresponding rhs with a factor:
int size_of_rhs=rhs.blockMap().NumMyElements(); 

for (int i =0; i<size_of_rhs;i++)
{
	int index_of_rhs = rhs.blockMap().GID (i);
	
	for (int j = 0;j<size_of_indicesvector;j++)
		{
			int dummy = indices[j];
	
				if(index_of_rhs == dummy)
					{	
						rhs [index_of_rhs] *=  M_weakessFactor;					
					}
		}
}

//update rhs:
*super::M_fiberActivationPtr += rhs;
}
void ActiveStressRossiModel14::solveModel (Real& timeStep)
{
    VectorEpetra rhs ( *( this->M_electroSolution.at(0) ) );
    rhs *= rhs.operator > (0.0);
    rhs *= M_coefficientBeta;
    rhs -= 2.0 * *super::M_fiberActivationPtr;
    rhs *= (timeStep * M_maximumActiveTenstion / M_coefficientMu);
    //rhs -= 2.0 * M_maximumActiveTenstion * *super::M_fiberActivationPtr;
    //rhs *= (timeStep / M_coefficientMu);

   *super::M_fiberActivationPtr += rhs;
   
}

} /* namespace LifeV */
