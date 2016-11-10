/*******************************************************************************
Copyright (c) 2010, Jonathan Hiller (Cornell University)
If used in publication cite "J. Hiller and H. Lipson "Dynamic Simulation of Soft Heterogeneous Objects" In press. (2011)"

This file is part of Voxelyze.
Voxelyze is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Voxelyze is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
See <http://www.opensource.org/licenses/lgpl-3.0.html> for license details.
*******************************************************************************/

#include "VX_SimGA.h"
#include <iostream>

CVX_SimGA::CVX_SimGA()
{
	Fitness = 0.0f;
	TrackVoxel = 0;
	FitnessFileName = "";
//	print_scrn = false;
	WriteFitnessFile = false;
	FitnessType = FT_NONE;	//no reporting is default

}

void CVX_SimGA::SaveResultFile(std::string filename)
{
	CXML_Rip XML;
	WriteResultFile(&XML);
	XML.SaveFile(filename);
}


void CVX_SimGA::WriteResultFile(CXML_Rip* pXML)
{


	float normFinalDist; // , normRegimeDist;
	normFinalDist = pow(pow(SS.CurCM.x-IniCM.x,2)+pow(SS.CurCM.y-IniCM.y,2),0.5)/LocalVXC.GetLatticeDim();

//	if(pEnv->pObj->GetUsingFinalVoxelSize())
//	{
//		normFinalDist = pow(pow(SS.EndOfLifetimeCM.x-IniCM.x,2)+pow(SS.EndOfLifetimeCM.y-IniCM.y,2),0.5)/LocalVXC.GetLatticeDim();
////		normRegimeDist = pow(pow(SS.CurCM.x-SS.EndOfLifetimeCM.x,2)+pow(SS.CurCM.y-SS.EndOfLifetimeCM.y,2),0.5)/LocalVXC.GetLatticeDim();
//	}
//	else
//	{
//		normFinalDist = pow(pow(SS.CurCM.x-IniCM.x,2)+pow(SS.CurCM.y-IniCM.y,2),0.5)/LocalVXC.GetLatticeDim();
////		normRegimeDist = pow(pow(SS.CurCM.x-SS.RegimeStartCM.x,2)+pow(SS.CurCM.y-SS.RegimeStartCM.y,2),0.5)/LocalVXC.GetLatticeDim();
//	}
//
//	pXML->DownLevel("Voxelyze_Sim_Result");
//	pXML->SetElAttribute("Version", "1.0");
//	pXML->DownLevel("Fitness");
//	pXML->Element("NormFinalDist", normFinalDist);
//	pXML->Element("NormRegimeDist", normRegimeDist);
//
//	pXML->UpLevel();
//	pXML->UpLevel();


	pXML->DownLevel("Voxelyze_Sim_Result");

		pXML->SetElAttribute("Version", "1.0");
		pXML->DownLevel("Fitness");
			pXML->Element("NormFinalDist", normFinalDist);
//			pXML->Element("NormRegimeDist", normRegimeDist);
		pXML->UpLevel();

		if (pEnv->getTimeBetweenTraces() > 0)
		{
			pXML->DownLevel("CMTrace");
				for(std::vector<vfloat>::size_type i = 0; i != SS.CMTraceTime.size(); ++i)
				{
	   				pXML->DownLevel("TraceStep");
	   					pXML->Element("Time",SS.CMTraceTime[i]);
						pXML->Element("TraceX",SS.CMTrace[i].x);
						pXML->Element("TraceY",SS.CMTrace[i].y);
						pXML->Element("TraceZ",SS.CMTrace[i].z);
					pXML->UpLevel();
				}
			pXML->UpLevel();
		}

	pXML->UpLevel();


}

void CVX_SimGA::WriteAdditionalSimXML(CXML_Rip* pXML)
{
	pXML->DownLevel("GA");
		pXML->Element("Fitness", Fitness);
		pXML->Element("FitnessType", (int)FitnessType);
		pXML->Element("TrackVoxel", TrackVoxel);
		pXML->Element("FitnessFileName", FitnessFileName);
		pXML->Element("WriteFitnessFile", WriteFitnessFile);
	pXML->UpLevel();
}

bool CVX_SimGA::ReadAdditionalSimXML(CXML_Rip* pXML, std::string* RetMessage)
{
	if (pXML->FindElement("GA")){
		int TmpInt;
		if (!pXML->FindLoadElement("Fitness", &Fitness)) Fitness = 0;
		if (pXML->FindLoadElement("FitnessType", &TmpInt)) FitnessType=(FitnessTypes)TmpInt; else Fitness = 0;
		if (!pXML->FindLoadElement("TrackVoxel", &TrackVoxel)) TrackVoxel = 0;
		if (!pXML->FindLoadElement("FitnessFileName", &FitnessFileName)) FitnessFileName = "";
		if (!pXML->FindLoadElement("QhullTmpFile", &QhullTmpFile)) QhullTmpFile = "";		
		if (!pXML->FindLoadElement("CurvaturesTmpFile", &CurvaturesTmpFile)) CurvaturesTmpFile = "";		
		if (!pXML->FindLoadElement("WriteFitnessFile", &WriteFitnessFile)) WriteFitnessFile = true;
		pXML->UpLevel();
	}

	return true;
}