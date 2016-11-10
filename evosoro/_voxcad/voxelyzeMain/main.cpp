#include <iostream>
#include "VX_Object.h"
#include "VX_Environment.h"
#include "VX_Sim.h"
#include "VX_SimGA.h"
#include "VX_MeshUtil.h"
#include "VXS_SimGLView.h"

int main(int argc, char *argv[])
{
	char* InputFile;
	//create the main objects
	CVXC_Structure structure;
	CVX_Object Object;
	CVX_Environment Environment;
	CVX_SimGA Simulator;
	CVX_MeshUtil DeformableMesh;
	long int Step = 0;
	vfloat Time = 0.0; //in seconds
	bool print_scrn = false;
	bool computeShapeDescriptors = false;

	// TODO: add --genVolumeInfo option?
	//first, parse inputs. Use as: -f followed by the filename of the .vxa file that describes the simulation. Can also follow this with -p to cause console output to occur
	if (argc < 3) 
	{ // Check the value of argc. If not enough parameters have been passed, inform user and exit.
		std::cout << "\nInput file required. Quitting.\n";
		return(0);	//return, indicating via code (0) that we did not complete the simulation
	} 
	else 
	{ // if we got enough parameters...
		for (int i = 1; i < argc; i++) 
		{ 
			if (strcmp(argv[i],"-f") == 0) 
			{
				InputFile = argv[i + 1];	// We know the next argument *should* be the filename:
			} 
			else if (strcmp(argv[i],"-p") == 0) 
			{
				print_scrn=true;	//decide if output to the console is desired
			}
			else if (strcmp(argv[i],"--computeShapeDescriptors") == 0)
			{
				computeShapeDescriptors = true;
			}
		}
	} 

	//setup main object
	Simulator.pEnv = &Environment;	//connect Simulation to environment
	Environment.pObj = &Object;		//connect environment to object
	Simulator.setInternalMesh(&DeformableMesh);

	//import the configuration file
	if (!Simulator.LoadVXAFile(InputFile)){
		if (print_scrn) std::cout << "\nProblem importing VXA file. Quitting\n";
		return(0);	//return, indicating via code (0) that we did not complete the simulation
		}
	std::string ReturnMessage;
	if (print_scrn) std::cout << "\nImporting Environment into simulator...\n";

	Simulator.Import(&Environment, 0, &ReturnMessage);
	if (print_scrn) std::cout << "Simulation import return message:\n" << ReturnMessage << "\n";
	
	Simulator.pEnv->UpdateCurTemp(Time);	//set the starting temperature (nac: pointer removed for debugging)


	double hulVolumeStart, hulVolumeEnd, robotVolumeStart, robotVolumeEnd, sComplexityStart, sComplexityEnd;

	DeformableMesh.initializeDeformableMesh(&Simulator); // Initialize internal mesh and link it to the simulation

	if(computeShapeDescriptors)
	{		
		hulVolumeStart = DeformableMesh.computeAndStoreQHullStart();
		robotVolumeStart = DeformableMesh.computeAndStoreRobotVolumeStart();
		sComplexityStart = DeformableMesh.computeInitialShapeComplexity();

		if (print_scrn)
		{
			std::cout << "Robot mesh has " << DeformableMesh.getMeshVertNum() << " vertices and " << DeformableMesh.getMeshFacetsNum() << " facets" << std::endl
					  << "Init robot volume: " << robotVolumeStart << std::endl 
					  << "Init convex hull volume: " << hulVolumeStart << std::endl << std::endl
					  << "Init shape complexity: " << sComplexityStart << std::endl;

			DeformableMesh.printAllMeshInfo();
		}
	}

	while (not Simulator.StopConditionMet())
	{
		// do some reporting via the stdoutput if required:
		if (Step%100 == 0.0 && print_scrn) //Only output every n time steps
		{
			std::cout << "Time: " << Time << std::endl;
			std::cout << "CM: " << Simulator.GetCM().Length() << std::endl << std::endl;
			
			// std::cout << " \tVox 0 X: " << Vox0Pos.x << "mm" << "\tVox 0 Y: " << Vox0Pos.y << "mm" << "\tVox 0 Z: " << Vox0Pos.z << "mm\n";	//just display the position of the first voxel in the voxelarray
			std::cout << "Vox[0]  Scale: " << Simulator.VoxArray[0].GetCurScale() << std::endl;
			std::cout << "Vox[0]  TempAmp: " << Simulator.VoxArray[0].TempAmplitude << std::endl;
			std::cout << "Vox[0]  TempPer: " << Simulator.VoxArray[0].TempPeriod << std::endl;
			std::cout << "Vox[0]  phaseOffset: " << Simulator.VoxArray[0].phaseOffset << std::endl;
			// std::cout << "Vox[5]  Scale: " << Simulator.VoxArray[5].GetCurScale() << std::endl;
			// std::cout << "Vox[10] Scale: " << Simulator.VoxArray[10].GetCurScale() << std::endl;
		}

		//do the actual simulation step
		Simulator.TimeStep(&ReturnMessage);
		Step += 1;	//increment the step counter
		Time += Simulator.dt;	//update the sim tim after the step
		Simulator.pEnv->UpdateCurTemp(Time);	//pass in the global time, and a pointer to the local object so its material temps can be modified (nac: pointer removed for debugging)	
	}

	if(computeShapeDescriptors)
	{
		hulVolumeEnd = DeformableMesh.computeAndStoreQHullEnd();
		robotVolumeEnd = DeformableMesh.computeAndStoreRobotVolumeEnd();
		sComplexityEnd = DeformableMesh.computeFinalShapeComplexity();

		if (print_scrn)
		{
			std::cout << "Final robot volume: " << robotVolumeEnd << std::endl 
					  << "Final convex hull volume: " << hulVolumeEnd << std::endl << std::endl
					  << "Final shape complexity: " << sComplexityEnd << std::endl;
  			DeformableMesh.printAllMeshInfo();
  		}
	}

	if (print_scrn) std::cout << "Ended at: " << Time << std::endl;

	Simulator.SaveResultFile(Simulator.FitnessFileName);

	return 1;	//code for successful completion  // could return fitness value if greater efficiency is desired
}
