/*******************************************************************************
Copyright (c) 2010, Jonathan Hiller (Cornell University)
If used in publication cite "J. Hiller and H. Lipson "Dynamic Simulation of Soft Heterogeneous Objects" In press. (2011)"

This file is part of Voxelyze.
Voxelyze is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Voxelyze is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
See <http://www.opensource.org/licenses/lgpl-3.0.html> for license details.
*******************************************************************************/

#include "VX_MeshUtil.h"
#include "VX_Object.h"
#include "VX_FEA.h"
#include "VX_Sim.h"
#include "VXS_SimGLView.h"
#include "Utils/MarchCube.h"
#include <iostream>


CVX_MeshUtil::CVX_MeshUtil(void)
{
	pSim = NULL;
	pSimView = NULL;
	pFEA = NULL;
	qhullVolumeStart = -1;
	qhullVolumeEnd = -1;
	robotVolumeStart = -1;
	robotVolumeEnd = -1;
	shapeComplexityStart = -1;
	shapeComplexityEnd = -1;
	usingGUI = false;
	tipVoxel = 0;
}

CVX_MeshUtil::~CVX_MeshUtil(void)
{
	Clear();
}

void CVX_MeshUtil::initializeDeformableMesh(CVX_Sim* pSimIn)
{
	LinkSimVoxels(pSimIn, NULL);
	DefMesh.DrawSmooth=false;

}

void CVX_MeshUtil::Clear(void)
{
	DefMesh.Clear();
	CalcVerts.clear();
	FacetToSIndex.clear();
	pSim = NULL;
	pFEA = NULL;
}


void CVX_MeshUtil::LinkSimExisting(CVX_Sim* pSimIn, CVXS_SimGLView* pSimViewIn, CMesh* pMeshIn) //links the simulation to an existing surface mesh so it can be deformed correctly
{
	if (pSimIn != NULL) pSim = pSimIn;
	if (pSimViewIn != NULL) pSimView = pSimViewIn;
	if (pMeshIn != NULL) DefMesh = *pMeshIn;
	DefMesh.tipVoxel = tipVoxel;
	if (!DefMesh.Exists()) return; //cna't inport if no mesh...

	vfloat FuzzDist = 1.5; //multiple of Lattice dimension to include voxels (radius)

	CVX_Object* pObj = pSimIn->pEnv->pObj; //temp, for convenience
	Vec3D<> WS = pObj->GetWorkSpace();
	Vec3D<> tmpVert, VoxCtr, Offset;
	int xi, yi, zi, Vi;
	vfloat dist;

	CalcVerts.clear();
	for (int i=0; i<(int)(DefMesh.Vertices.size()); i++){
		CVertexCalc tmpVertCalc;
		tmpVert = DefMesh.Vertices[i].v;

		//check for all voxel centerpoints within one voxel of the voxel this one would be contained in?
		//see which voxe this point is containted in...
		xi = (int)(tmpVert.x/WS.x*pObj->GetVXDim());
		yi = (int)(tmpVert.y/WS.y*pObj->GetVYDim());
		zi = (int)(tmpVert.z/WS.z*pObj->GetVZDim());
		int D2S = int(FuzzDist*(pObj->GetLatticeDim()/pObj->GetLatDimEnv().Min())); //distance to search... make sure we look far enough if some dimensions are smaller than lattice dimension
		for (int x=xi-D2S; x<=xi+D2S; x++){
			for (int y=yi-D2S; y<=yi+D2S; y++){
				for (int z=zi-D2S; z<=zi+D2S; z++){
					Vi = pObj->GetIndex(x, y, z);
					if (Vi != -1 && pObj->GetMat(Vi) != 0){ //if valid voxel location
						VoxCtr = pObj->GetXYZ(Vi);
						Offset = tmpVert-VoxCtr;
						dist = Offset.Length()/pObj->GetLatticeDim();
						if (dist<FuzzDist){ //FINALLY! make a connection!
							tmpVertCalc.ConVoxels.push_back(CVertexComp(Vi, Offset, (FuzzDist-dist)/FuzzDist));
						}
					}
				}
			}
		}
		CalcVerts.push_back(tmpVertCalc);
	}

	DefMesh.DrawSmooth = true;

}

void CVX_MeshUtil::LinkSimSmooth(CVX_Sim* pSimIn, CVXS_SimGLView* pSimViewIn) //creates a voxel mesh and links to simulation
{
	CMesh GeneratedSmoothMesh;

	CArray3Df OccupancyArray(pSimIn->pEnv->pObj->GetVXDim(), pSimIn->pEnv->pObj->GetVYDim(), pSimIn->pEnv->pObj->GetVZDim()); 
	int NumPossibleVox = pSimIn->pEnv->pObj->GetStArraySize();
	for (int g=0; g<NumPossibleVox; g++){
		if (pSimIn->pEnv->pObj->Structure.GetData(g)>0) OccupancyArray[g] = 1.0;
	}
	CMarchCube::SingleMaterial(&GeneratedSmoothMesh, &OccupancyArray, 0.5, pSimIn->pEnv->pObj->GetLatticeDim());
	LinkSimExisting(pSimIn, pSimViewIn, &GeneratedSmoothMesh);
}

void CVX_MeshUtil::LinkSimVoxels(CVX_Sim* pSimIn, CVXS_SimGLView* pSimViewIn)
{
	if (pSimIn != NULL) pSim = pSimIn;
	if (pSimViewIn != NULL) pSimView = pSimViewIn;

	CMesh tmpMesh;

	CVX_Object* pDM = pSim->pEnv->pObj; //for convenience of access...

	//make all the potential lattice points:
//	std::vector<CVertexCalc> tmpVerts; //3D array [x][y][z]
	CalcVertsAll.clear();
	int tVX = pDM->GetVXDim()+1;
	int tVY = pDM->GetVYDim()+1;
	int tVZ = pDM->GetVZDim()+1;
	Vec3D<> Offset, tmp;

	FacetToSIndex.clear();
	CalcVertsAll.resize(tVX*tVY*tVZ); //each row in x, rolling around to consecutive Y's, then z.

	int NumVox = pSim->NumVox();
	for (int i=0; i<NumVox; i++){ //for every voxel in the simulation...
		CVXS_Voxel& CurVox = pSim->VoxArray[i];
		int CurXIndex = pSim->StoXIndexMap[i];
//		int CurXIndex = CurVox.GetVxcIndex();
		int cX, cY, cZ, tX, tY, tZ;
		pDM->GetXYZNom(&cX, &cY, &cZ, CurXIndex); //get the nominal coords if this voxel...
//		bool IsSurfaceVox = false;

		for (int i=0; i<8; i++){
			int tmpVertIndex = D3IndCorner(cX, cY, cZ, tVX-1, tVY-1, i);
//			GetCornerDir(i, &tmp);
//			Offset = (pDM->GetLatDimEnv()/2).Scale(tmp); //vector from center of voxel to this corner point
			CalcVertsAll[tmpVertIndex].ConVoxels.push_back(CVertexComp(CurXIndex, i));
		}

//		//look at all 8 corners, add link back to this voxel to all exterior corners
//		for (int x=-1; x<=1; x+=2){ //x = -1, then x = +1
//			for (int y=-1; y<=1; y+=2){ //y = -1, then y = +1
//				for (int z=-1; z<=1; z+=2){ //z = -1, then z = +1
////					if (pDM->IsVoxRelativeIndex(CurXIndex, x, y, z) && pDM->IsVoxRelativeIndex(CurXIndex, 0, y, z) && pDM->IsVoxRelativeIndex(CurXIndex, x, 0, z) && pDM->IsVoxRelativeIndex(CurXIndex, x, y, 0) && pDM->IsVoxRelativeIndex(CurXIndex, 0, 0, z) && pDM->IsVoxRelativeIndex(CurXIndex, x, 0, 0) && pDM->IsVoxRelativeIndex(CurXIndex, 0, y, 0)) continue; //if interior, then continue to next point
//
//					tX = (x == -1) ? 0 : 1; //tX to 0 or 1
//					tY = (y == -1) ? 0 : 1; //tY to 0 or 1
//					tZ = (z == -1) ? 0 : 1; //tZ to 0 or 1
//					int tmpVertIndex = D3Ind(cX+tX, cY+tY, cZ+tZ, tVX, tVY); // (cZ+tZ)*tVX*tVY + (cY+tY)*tVX + (cX+tX); //which vertex in out n+1 array are we dealing with?
//
//					Vec3D Offset = (pDM->GetLatDimEnv()/2).Scale(Vec3D(x, y, z)); //vector from center of voxel to this corner point
//					CalcVertsAll[tmpVertIndex].ConVoxels.push_back(CVertexComp(CurXIndex, Offset, 1.0f));
//				}
//			}
//		}

//		if (IsSurfaceVox){
			//look at all faces, add facets (no non-existant vertex indices.. but we'll shift down
			if (!pDM->IsVoxRelativeIndex(CurXIndex, 1, 0, 0)){ //if we have an open face +X
				tmpMesh.Facets.push_back(CFacet(Vec3D<>(1, 0, 0), D3Ind(cX+1, cY, cZ, tVX, tVY), D3Ind(cX+1, cY+1, cZ, tVX, tVY), D3Ind(cX+1, cY+1, cZ+1, tVX, tVY), CurXIndex));
				tmpMesh.Facets.push_back(CFacet(Vec3D<>(1, 0, 0), D3Ind(cX+1, cY, cZ, tVX, tVY), D3Ind(cX+1, cY+1, cZ+1, tVX, tVY), D3Ind(cX+1, cY, cZ+1, tVX, tVY), CurXIndex));
			}
			if (!pDM->IsVoxRelativeIndex(CurXIndex, -1, 0, 0)){ //if we have an open face -X
				tmpMesh.Facets.push_back(CFacet(Vec3D<>(-1, 0, 0), D3Ind(cX, cY, cZ, tVX, tVY), D3Ind(cX, cY+1, cZ+1, tVX, tVY), D3Ind(cX, cY+1, cZ, tVX, tVY), CurXIndex));
				tmpMesh.Facets.push_back(CFacet(Vec3D<>(-1, 0, 0), D3Ind(cX, cY, cZ, tVX, tVY), D3Ind(cX, cY, cZ+1, tVX, tVY), D3Ind(cX, cY+1, cZ+1, tVX, tVY), CurXIndex));

			}
			if (!pDM->IsVoxRelativeIndex(CurXIndex, 0, 1, 0)){ //if we have an open face + Y
				tmpMesh.Facets.push_back(CFacet(Vec3D<>(0, 1, 0), D3Ind(cX, cY+1, cZ, tVX, tVY), D3Ind(cX, cY+1, cZ+1, tVX, tVY), D3Ind(cX+1, cY+1, cZ+1, tVX, tVY), CurXIndex));
				tmpMesh.Facets.push_back(CFacet(Vec3D<>(0, 1, 0), D3Ind(cX, cY+1, cZ, tVX, tVY), D3Ind(cX+1, cY+1, cZ+1, tVX, tVY), D3Ind(cX+1, cY+1, cZ, tVX, tVY), CurXIndex));
			}
			if (!pDM->IsVoxRelativeIndex(CurXIndex, 0, -1, 0)){ //if we have an open face + Y
				tmpMesh.Facets.push_back(CFacet(Vec3D<>(0, -1, 0), D3Ind(cX, cY, cZ, tVX, tVY), D3Ind(cX+1, cY, cZ+1, tVX, tVY), D3Ind(cX, cY, cZ+1, tVX, tVY), CurXIndex));
				tmpMesh.Facets.push_back(CFacet(Vec3D<>(0, -1, 0), D3Ind(cX, cY, cZ, tVX, tVY), D3Ind(cX+1, cY, cZ, tVX, tVY), D3Ind(cX+1, cY, cZ+1, tVX, tVY), CurXIndex));
			}
			if (!pDM->IsVoxRelativeIndex(CurXIndex, 0, 0, 1)){ //if we have an open face + Y
				tmpMesh.Facets.push_back(CFacet(Vec3D<>(0, 0, 1), D3Ind(cX, cY, cZ+1, tVX, tVY), D3Ind(cX+1, cY, cZ+1, tVX, tVY), D3Ind(cX+1, cY+1, cZ+1, tVX, tVY), CurXIndex));
				tmpMesh.Facets.push_back(CFacet(Vec3D<>(0, 0, 1), D3Ind(cX, cY, cZ+1, tVX, tVY), D3Ind(cX+1, cY+1, cZ+1, tVX, tVY), D3Ind(cX, cY+1, cZ+1, tVX, tVY), CurXIndex));
			}
			if (!pDM->IsVoxRelativeIndex(CurXIndex, 0, 0, -1)){ //if we have an open face + Y
				tmpMesh.Facets.push_back(CFacet(Vec3D<>(0, 0, -1), D3Ind(cX, cY, cZ, tVX, tVY), D3Ind(cX+1, cY+1, cZ, tVX, tVY), D3Ind(cX+1, cY, cZ, tVX, tVY), CurXIndex));
				tmpMesh.Facets.push_back(CFacet(Vec3D<>(0, 0, -1), D3Ind(cX, cY, cZ, tVX, tVY), D3Ind(cX, cY+1, cZ, tVX, tVY), D3Ind(cX+1, cY+1, cZ, tVX, tVY), CurXIndex));
			}


			while (FacetToSIndex.size() < tmpMesh.Facets.size()) //update our linking list with all the facets we just added
				FacetToSIndex.push_back(pSim->XtoSIndexMap[CurXIndex]);


			//look at all 12 edges... (shared lines added twice, but oh well... we will delete them later
			//vertical lines
			if (!pDM->IsVoxRelativeIndex(CurXIndex, 1, 0, 0) || !pDM->IsVoxRelativeIndex(CurXIndex, 0, 1, 0))
				tmpMesh.Lines.push_back(CLine(D3Ind(cX+1, cY+1, cZ, tVX, tVY), D3Ind(cX+1, cY+1, cZ+1, tVX, tVY)));
			if (!pDM->IsVoxRelativeIndex(CurXIndex, -1, 0, 0) || !pDM->IsVoxRelativeIndex(CurXIndex, 0, 1, 0))
				tmpMesh.Lines.push_back(CLine(D3Ind(cX, cY+1, cZ, tVX, tVY), D3Ind(cX, cY+1, cZ+1, tVX, tVY)));
			if (!pDM->IsVoxRelativeIndex(CurXIndex, -1, 0, 0) || !pDM->IsVoxRelativeIndex(CurXIndex, 0, -1, 0))
				tmpMesh.Lines.push_back(CLine(D3Ind(cX, cY, cZ, tVX, tVY), D3Ind(cX, cY, cZ+1, tVX, tVY)));
			if (!pDM->IsVoxRelativeIndex(CurXIndex, 1, 0, 0) || !pDM->IsVoxRelativeIndex(CurXIndex, 0, -1, 0))
				tmpMesh.Lines.push_back(CLine(D3Ind(cX+1, cY, cZ, tVX, tVY), D3Ind(cX+1, cY, cZ+1, tVX, tVY)));

			//top lines
			if (!pDM->IsVoxRelativeIndex(CurXIndex, 1, 0, 0) || !pDM->IsVoxRelativeIndex(CurXIndex, 0, 0, 1))
				tmpMesh.Lines.push_back(CLine(D3Ind(cX+1, cY, cZ+1, tVX, tVY), D3Ind(cX+1, cY+1, cZ+1, tVX, tVY)));
			if (!pDM->IsVoxRelativeIndex(CurXIndex, 0, 1, 0) || !pDM->IsVoxRelativeIndex(CurXIndex, 0, 0, 1))
				tmpMesh.Lines.push_back(CLine(D3Ind(cX, cY+1, cZ+1, tVX, tVY), D3Ind(cX+1, cY+1, cZ+1, tVX, tVY)));
			if (!pDM->IsVoxRelativeIndex(CurXIndex, -1, 0, 0) || !pDM->IsVoxRelativeIndex(CurXIndex, 0, 0, 1))
				tmpMesh.Lines.push_back(CLine(D3Ind(cX, cY, cZ+1, tVX, tVY), D3Ind(cX, cY+1, cZ+1, tVX, tVY)));
			if (!pDM->IsVoxRelativeIndex(CurXIndex, 0, -1, 0) || !pDM->IsVoxRelativeIndex(CurXIndex, 0, 0, 1))
				tmpMesh.Lines.push_back(CLine(D3Ind(cX, cY, cZ+1, tVX, tVY), D3Ind(cX+1, cY, cZ+1, tVX, tVY)));

			//bottom lines
			if (!pDM->IsVoxRelativeIndex(CurXIndex, 1, 0, 0) || !pDM->IsVoxRelativeIndex(CurXIndex, 0, 0, -1))
				tmpMesh.Lines.push_back(CLine(D3Ind(cX+1, cY, cZ, tVX, tVY), D3Ind(cX+1, cY+1, cZ, tVX, tVY)));
			if (!pDM->IsVoxRelativeIndex(CurXIndex, 0, 1, 0) || !pDM->IsVoxRelativeIndex(CurXIndex, 0, 0, -1))
				tmpMesh.Lines.push_back(CLine(D3Ind(cX, cY+1, cZ, tVX, tVY), D3Ind(cX+1, cY+1, cZ, tVX, tVY)));
			if (!pDM->IsVoxRelativeIndex(CurXIndex, -1, 0, 0) || !pDM->IsVoxRelativeIndex(CurXIndex, 0, 0, -1))
				tmpMesh.Lines.push_back(CLine(D3Ind(cX, cY, cZ, tVX, tVY), D3Ind(cX, cY+1, cZ, tVX, tVY)));
			if (!pDM->IsVoxRelativeIndex(CurXIndex, 0, -1, 0) || !pDM->IsVoxRelativeIndex(CurXIndex, 0, 0, -1))
				tmpMesh.Lines.push_back(CLine(D3Ind(cX, cY, cZ, tVX, tVY), D3Ind(cX+1, cY, cZ, tVX, tVY)));
//		}
	}


	Vec3D<> tmpOffset = pDM->GetLatDimEnv()/2;
	CVX_Object tmpDM = *pDM;
	tmpDM.Resize(tVX, tVY, tVZ);
	//add all vertices..
	for (int k=0; k<tVZ; k++){
		for (int j=0; j<tVY; j++){
			for (int i=0; i<tVX; i++){
				Vec3D<> pos = tmpDM.GetXYZ(i, j, k)-tmpOffset;
				tmpMesh.Vertices.push_back(CVertex(pos));
			}
		}
	}

	tmpMesh.RemoveDupLines();

	//tmpMesh is now complete, but with way too many vertices, most unused.

	//commit to DefMesh, minus all the extra vertices!!!
	DefMesh.Clear();
	CalcVerts.clear();

	int* Map = new int[tVX*tVY*tVZ];
	int NewInd = 0;
	for (int i=0; i<tVX*tVY*tVZ; i++){
		if (CalcVertsAll[i].ConVoxels.size() != 0 && CalcVertsAll[i].ConVoxels.size() != 8){
			DefMesh.Vertices.push_back(tmpMesh.Vertices[i]);
			CalcVerts.push_back(CalcVertsAll[i]);
			Map[i] = NewInd++;

		}
		else Map[i] = -1;
	}

	DefMesh.Facets = tmpMesh.Facets;
	for (int i=0; i<(int)DefMesh.Facets.size(); i++){
		for (int j=0; j<3; j++){
			DefMesh.Facets[i].vi[j] = Map[tmpMesh.Facets[i].vi[j]];
		}
	}

	DefMesh.Lines = tmpMesh.Lines;
	for (int i=0; i<(int)DefMesh.Lines.size(); i++){
		for (int j=0; j<2; j++){
			DefMesh.Lines[i].vi[j] = Map[tmpMesh.Lines[i].vi[j]];
		}
	}


	/*
	updateDeformableMesh(); 
	std::cout << "Tip voxel is " << pSim->getTipVoxel() << std::endl;

	for(int i = 0; i < DefMesh.Facets.size(); i++)
	{	
		if(FacetToSIndex[i] == pSim->getTipVoxel())
		{
			DefMesh.Facets[i].drawPointingDirection = true;
			//std::cout << "Facet no. " << i << " belongs to voxel no. " << pSim->getTipVoxel() << ": will draw normals for debug"<< std::endl;

			// store a pointer to the correct normal so that sim (via pSim) and simGA (calling methods of internalMesh) can get the correct pointing dir
		}else
			DefMesh.Facets[i].drawPointingDirection = false;
	}
*/
}

void CVX_MeshUtil::UpdateMesh(int CurSel) //updates mesh based on linked FEA/Relaxation
{

	if (pSim){/*
		//update vertex positions
		Vec3D avgPos, ThisPos;
		vfloat TotWeight, Ra, Ga, Ba;
		for (int i=0; i<(int)CalcVerts.size(); i++){
			avgPos = Vec3D(0,0,0);
			TotWeight = 0.0;
#ifdef USE_OPEN_GL
			CColor Tmp;
			Ra = 0.0;
			Ga = 0.0;
			Ba = 0.0;
#endif

			for (int j=0; j<(int)CalcVerts[i].ConVoxels.size(); j++){
				CVXS_Voxel& CurVox = pSim->VoxArray[pSim->XtoSIndexMap[CalcVerts[i].ConVoxels[j].XIndex]];

				//temp, until scale is updated by temp!
				vfloat TempFact = 1.0;
				if(pSim->pEnv->TempEnabled) TempFact = (1+(pSim->pEnv->CurTemp-pSim->pEnv->TempBase)*CurVox.GetCTE()); //pSim->LocalVXC.GetBaseMat(CurVox.MatIndex)->GetCTE());

				ThisPos = CurVox.GetCurPos() + (CurVox.Angle*CQuat(TempFact*CalcVerts[i].ConVoxels[j].Off)*CurVox.Angle.Conjugate()).ToVec(); //todo:scale!
				avgPos += CalcVerts[i].ConVoxels[j].Weight*ThisPos;

#ifdef USE_OPEN_GL
				Tmp = pSim->GetCurColor(pSim->XtoSIndexMap[CalcVerts[i].ConVoxels[j].XIndex], CurSel);
				Ra += CalcVerts[i].ConVoxels[j].Weight*Tmp.r;
				Ga += CalcVerts[i].ConVoxels[j].Weight*Tmp.g;
				Ba += CalcVerts[i].ConVoxels[j].Weight*Tmp.b;
#endif

				TotWeight += CalcVerts[i].ConVoxels[j].Weight;
			}

			DefMesh.Vertices[i].DrawOffset = avgPos/TotWeight - DefMesh.Vertices[i].v; //yes, we're subtracting this out just to add it back.

#ifdef USE_OPEN_GL
			DefMesh.Vertices[i].VColor = CColor(Ra/TotWeight, Ga/TotWeight, Ba/TotWeight, 1.0);
#endif
		}*/

		Vec3D<> NewPos;
		for (int i=0; i<(int)CalcVerts.size(); i++){
			GetCurVLoc(CalcVerts[i], &NewPos);
			DefMesh.Vertices[i].DrawOffset = NewPos - DefMesh.Vertices[i].v; //yes, we're subtracting this out just to add it back.

#ifdef USE_OPEN_GL
			GetCurVCol(CalcVerts[i], &DefMesh.Vertices[i].VColor, CurSel);
#endif
		}


		//update colors that aren't by vertex!
#ifdef USE_OPEN_GL
		if(pSimView)
			if (FacetToSIndex.size() != 0)
			{ //if pulling color for full fact
				for (int i=0; i<(int)DefMesh.Facets.size(); i++){
					DefMesh.Facets[i].FColor = pSimView->GetCurVoxColor(FacetToSIndex[i], CurSel);
					//DefMesh.Facets[i].FColor = pSim->GetCurColor(i, CurSel);
				}
			}
		//else { //color by vertex!!!...

		//}
#endif

		//update normals!
		if (!DefMesh.DrawSmooth) //if drawing faces...
			DefMesh.CalcFaceNormals();
		else
			DefMesh.CalcVertNormals();
		

	}
	else if (pFEA){
		//todo...
	}
}

void CVX_MeshUtil::GetCurVLoc(CVertexCalc& VertCalc, Vec3D<>* pLocOut)
{
	Vec3D<> avgPos(0,0,0);
	Vec3D<> ThisPos;
	vfloat TotWeight = 0;

	for (int j=0; j<(int)VertCalc.ConVoxels.size(); j++){
		CVertexComp CurComp = VertCalc.ConVoxels[j];
		CVXS_Voxel& CurVox = pSim->VoxArray[pSim->XtoSIndexMap[CurComp.XIndex]];

		Vec3D<> ThisOffset; //offset to this point
		if (CurComp.Corner != -1){
			switch (CurComp.Corner){
			case NNN: ThisOffset = CurVox.GetCornerNeg(); break;
			case NNP: ThisOffset = Vec3D<>(CurVox.GetCornerNeg().x, CurVox.GetCornerNeg().y, CurVox.GetCornerPos().z); break;
			case NPN: ThisOffset = Vec3D<>(CurVox.GetCornerNeg().x, CurVox.GetCornerPos().y, CurVox.GetCornerNeg().z); break;
			case NPP: ThisOffset = Vec3D<>(CurVox.GetCornerNeg().x, CurVox.GetCornerPos().y, CurVox.GetCornerPos().z); break;
			case PNN: ThisOffset = Vec3D<>(CurVox.GetCornerPos().x, CurVox.GetCornerNeg().y, CurVox.GetCornerNeg().z); break;
			case PNP: ThisOffset = Vec3D<>(CurVox.GetCornerPos().x, CurVox.GetCornerNeg().y, CurVox.GetCornerPos().z); break;
			case PPN: ThisOffset = Vec3D<>(CurVox.GetCornerPos().x, CurVox.GetCornerPos().y, CurVox.GetCornerNeg().z); break;
			case PPP: ThisOffset = CurVox.GetCornerPos(); break;
			}


		}
		else {
			vfloat ScaleFact = CurVox.GetCurScale() / CurVox.GetNominalSize(); //Assumes square, isotropic expansion
			ThisOffset = ScaleFact*VertCalc.ConVoxels[j].Off;

		//ThisPos = CurVox.GetCurPos() + (CurVox.GetCurAngle()*CQuat<>(ScaleFact*VertCalc.ConVoxels[j].Off)*CurVox.GetCurAngle().Conjugate()).ToVec(); //todo:scale!
		}
		
//		ThisPos = CurVox.GetCurPos() + (CurVox.GetCurAngle()*CQuat<>(ThisOffset)*CurVox.GetCurAngle().Conjugate()).ToVec(); //todo:scale!
		ThisPos = CurVox.GetCurPos() + CurVox.GetCurAngle().RotateVec3D(ThisOffset);

		avgPos += VertCalc.ConVoxels[j].Weight*ThisPos;
		TotWeight += VertCalc.ConVoxels[j].Weight;
	}

	*pLocOut = avgPos/TotWeight; 
}

void CVX_MeshUtil::GetCurVCol(CVertexCalc& VertCalc, CColor* pColOut, int CurSel)
{
#ifdef USE_OPEN_GL
	vfloat TotWeight = 0.0;
	vfloat Ra = 0.0f;
	vfloat Ga = 0.0f;
	vfloat Ba = 0.0f;
	CColor Tmp;

	for (int j=0; j<(int)VertCalc.ConVoxels.size(); j++){
		Tmp = pSimView->GetCurVoxColor(pSim->XtoSIndexMap[VertCalc.ConVoxels[j].XIndex], CurSel);
		Ra += VertCalc.ConVoxels[j].Weight*Tmp.r;
		Ga += VertCalc.ConVoxels[j].Weight*Tmp.g;
		Ba += VertCalc.ConVoxels[j].Weight*Tmp.b;

		TotWeight += VertCalc.ConVoxels[j].Weight;
	}


	*pColOut = CColor(Ra/TotWeight, Ga/TotWeight, Ba/TotWeight, 1.0);
#endif
}

int CVX_MeshUtil::D3IndCorner(int XInd, int YInd, int ZInd, const int XSize, const int YSize, int const Corner) //returns the vertex index in the [x+1, y+1, z+1] CalcVertsAll array corresponding to the specified corner of this voxel.
{
	switch (Corner){
//	case NNN: *pOut = Vec3D(-1,-1,-1); break;
	case NNP: ZInd++; break;
	case NPN: YInd++; break;
	case NPP: YInd++; ZInd++; break;
	case PNN: XInd++; break;
	case PNP: XInd++; ZInd++; break;
	case PPN: XInd++; YInd++; break;
	case PPP: XInd++; YInd++; ZInd++; break;
	}
	return D3Ind(XInd, YInd, ZInd, XSize+1, YSize+1);
}

void CVX_MeshUtil::GetCornerDir(int Corner, Vec3D<>* pOut) //returns vector with component +/- 1 dependeind on which corner
{
	switch (Corner){
	case NNN: *pOut = Vec3D<>(-1,-1,-1); break;
	case NNP: *pOut = Vec3D<>(-1,-1, 1); break;
	case NPN: *pOut = Vec3D<>(-1, 1,-1); break;
	case NPP: *pOut = Vec3D<>(-1, 1, 1); break;
	case PNN: *pOut = Vec3D<>( 1,-1,-1); break;
	case PNP: *pOut = Vec3D<>( 1,-1, 1); break;
	case PPN: *pOut = Vec3D<>( 1, 1,-1); break;
	case PPP: *pOut = Vec3D<>( 1, 1, 1); break;
	default: *pOut = Vec3D<>(0,0,0);
	}
}

void CVX_MeshUtil::Draw(void)
{
#ifdef USE_OPEN_GL
	DefMesh.Draw();
#endif
}


//Misc
bool CVX_MeshUtil::ToStl(std::string BasePath, CVX_Object* pObj, bool WantDefMes)
{
	if (BasePath == "") return false;

	CMesh tmpMesh;
	int XSize = pObj->GetVXDim();
	int YSize = pObj->GetVYDim();

	switch (pObj->Voxel.GetVoxName()){
	case VS_BOX:
		for (int m=1; m<=(int)pObj->Palette.size(); m++){ //for each material (skip null material)
			if (pObj->GetNumVox(m) == 0) continue; //see if there's material here, continue if not

			tmpMesh.Clear();
			int CheckIndex = 0;
			Vec3D<> Center;
			int NomX, NomY, NomZ;
			Vec3D<> LatDims = pObj->GetLatDimEnv();
			//vertices of the cube in order of (---, --+, -+-, -++, +--, +-+, ++-, +++)
			Vec3D<> V[8] = {Vec3D<>(0,0,0)};

			for (int i=0; i<pObj->GetStArraySize(); i++){ //for every element of the master array
				if (pObj->GetMat(i) != m) continue;  //if there's not a voxel here of the correct material, keep going

				pObj->GetXYZ(&Center, i); //loads the center variable with coordinates of center of cube
				pObj->GetXYZNom(&NomX, &NomY, &NomZ, i);

				if (WantDefMes){
					GetCurVLoc(CalcVertsAll[D3IndCorner(NomX, NomY, NomZ, XSize, YSize, NNN)], &V[NNN]);
					GetCurVLoc(CalcVertsAll[D3IndCorner(NomX, NomY, NomZ, XSize, YSize, NNP)], &V[NNP]);
					GetCurVLoc(CalcVertsAll[D3IndCorner(NomX, NomY, NomZ, XSize, YSize, NPN)], &V[NPN]);
					GetCurVLoc(CalcVertsAll[D3IndCorner(NomX, NomY, NomZ, XSize, YSize, NPP)], &V[NPP]);
					GetCurVLoc(CalcVertsAll[D3IndCorner(NomX, NomY, NomZ, XSize, YSize, PNN)], &V[PNN]);
					GetCurVLoc(CalcVertsAll[D3IndCorner(NomX, NomY, NomZ, XSize, YSize, PNP)], &V[PNP]);
					GetCurVLoc(CalcVertsAll[D3IndCorner(NomX, NomY, NomZ, XSize, YSize, PPN)], &V[PPN]);
					GetCurVLoc(CalcVertsAll[D3IndCorner(NomX, NomY, NomZ, XSize, YSize, PPP)], &V[PPP]);

					//V[NNN] = 
				}
				else{ //not deformed mesh
				V[NNN] = Center+Vec3D<>(-LatDims.x/2,-LatDims.y/2,-LatDims.z/2);
				V[NNP] = Center+Vec3D<>(-LatDims.x/2,-LatDims.y/2, LatDims.z/2);
				V[NPN] = Center+Vec3D<>(-LatDims.x/2, LatDims.y/2,-LatDims.z/2);
				V[NPP] = Center+Vec3D<>(-LatDims.x/2, LatDims.y/2, LatDims.z/2);
				V[PNN] = Center+Vec3D<>( LatDims.x/2,-LatDims.y/2,-LatDims.z/2);
				V[PNP] = Center+Vec3D<>( LatDims.x/2,-LatDims.y/2, LatDims.z/2);
				V[PPN] = Center+Vec3D<>( LatDims.x/2, LatDims.y/2,-LatDims.z/2);
				V[PPP] = Center+Vec3D<>( LatDims.x/2, LatDims.y/2, LatDims.z/2);
				}

				for (int j=0; j<6; j++){ //for each face of the cube
					switch (j){
						case 0: //top! (+Z)
							CheckIndex = pObj->GetIndex(NomX, NomY, NomZ+1);
							if (CheckIndex == -1 || pObj->GetMat(CheckIndex) != m) //if its outside range or there's no cube there...
								tmpMesh.AddQuadFacet(V[PPP], V[NPP], V[NNP], V[PNP]);
						break;
						case 1: //bottom!
							CheckIndex = pObj->GetIndex(NomX, NomY, NomZ-1);
							if (CheckIndex == -1 || pObj->GetMat(CheckIndex) != m)
								tmpMesh.AddQuadFacet(V[PNN], V[NNN], V[NPN], V[PPN]);
						break;
						case 2: //right!
							CheckIndex = pObj->GetIndex(NomX+1, NomY, NomZ);
							if (CheckIndex == -1 || pObj->GetMat(CheckIndex) != m)
								tmpMesh.AddQuadFacet(V[PNN], V[PPN], V[PPP], V[PNP]);
						break;

						case 3: //left!
							CheckIndex = pObj->GetIndex(NomX-1, NomY, NomZ);
							if (CheckIndex == -1 || pObj->GetMat(CheckIndex) != m)
								tmpMesh.AddQuadFacet(V[NNN], V[NNP], V[NPP], V[NPN]);
						break;

						case 4: //front!
							CheckIndex = pObj->GetIndex(NomX, NomY+1, NomZ);
							if (CheckIndex == -1 || pObj->GetMat(CheckIndex) != m)
								tmpMesh.AddQuadFacet(V[NPN], V[NPP], V[PPP], V[PPN]);
						break;

						case 5: //back!
							CheckIndex = pObj->GetIndex(NomX, NomY-1, NomZ);
							if (CheckIndex == -1 || pObj->GetMat(CheckIndex) != m)
								tmpMesh.AddQuadFacet(V[PNN], V[PNP], V[NNP], V[NNN]);
						break;
					}
				}
			}
			std::string tmp = "_" + pObj->Palette[m].GetName();
			std::string ThisPath = BasePath;
			ThisPath.insert(ThisPath.size()-4, tmp);
			tmpMesh.SaveSTL(ThisPath);
		}
		
		break;
	case VS_SPHERE:
		for (int m=1; m<=(int)pObj->Palette.size(); m++){ //for each material (skip null material)
			if (pObj->GetNumVox(m) == 0) continue; //see if there's material here, continue if not

			tmpMesh.Clear();
			Vec3D<> Center, p1, p2, p3, p4;
			vfloat a1, a2, x1, x2, Rad1, Rad2;

			for (int i=0; i<pObj->GetStArraySize(); i++){ //for every element of the master array
				if (pObj->GetMat(i) != m) continue;  //if there's not a voxel here of the correct material, keep going
				pObj->GetXYZ(&Center, i); //loads the center variable with coordinates of center of cube
				
				//adds a sphere to the mesh:
				vfloat pi = atan(1.0)*4.0;
				vfloat da = pi/12;
				vfloat Scale = pObj->GetLatticeDim()/2.01;

				for (a2=0; a2<pi; a2+=da) {
					x1 = cos(a2);
					Rad1 = sin(a2);
					x2 = cos(a2+da);
					Rad2 = sin(a2+da);
					for (a1=0; a1<pi*2.01; a1+=da) {
						p1 = Vec3D<>(x1, Rad1*cos(a1), Rad1*sin(a1));
						p2 = Vec3D<>(x2, Rad2*cos(a1), Rad2*sin(a1));
						p3 = Vec3D<>(x1, Rad1*cos(a1+da), Rad1*sin(a1+da));
						p4 = Vec3D<>(x2, Rad2*cos(a1+da), Rad2*sin(a1+da));

						p1 = p1*Scale + Center;
						p2 = p2*Scale + Center;
						p3 = p3*Scale + Center;
						p4 = p4*Scale + Center;


						if (Rad1<1e-6)
							tmpMesh.AddFacet(p1, p2, p4);
						else if (Rad2<1e-6)
							tmpMesh.AddFacet(p1, p2, p3);
						else {
							tmpMesh.AddFacet(p1, p2, p3);
							tmpMesh.AddFacet(p2, p4, p3);
						}

					}
				}

			}


			std::string tmp = "_" + pObj->Palette[m].GetName();
			std::string ThisPath = BasePath;
			ThisPath.insert(ThisPath.size()-4, tmp);
			tmpMesh.SaveSTL(ThisPath);
		}
		break;

	default:
		return false;
	}
	return true;
}



bool CVX_MeshUtil::FromStl(CMesh* pMeshIn, CVX_Object* pObj, int MatIndex)
{
	Vec3D<> Loc;
	for (int i=0; i<pObj->GetStArraySize(); i++){
		Loc = pObj->GetXYZ(i);
		if(pMeshIn->IsInside(&Loc)) pObj->SetMat(i, MatIndex);
	}

	return true;

	/*
//	if (pObj->Structure.X_Voxels == 0)
//		pObj->InitializeMatter(); //initialize defaults if haven't loaded yet...

	Vec3D Point;
	//pObj->Cl
	int TotNumFacets = (int)pMeshIn->Facets.size();
	int* pZList = new int[TotNumFacets]; //array of facet indices that touch current Z plane...
	for (int i=0; i<TotNumFacets; i++) pZList[i] = -1; //initialize all to -1
	
	vfloat* pIntrscts = new vfloat[1000]; //the most intersections we expect to see...

	int Xv, Yv, Zv;
	int CurNumInt; //current number of intersections...

	pObj->GetVDim(&Xv, &Yv, &Zv);
	for (int k=0; k<Zv; k++){ //iterate through z
		pObj->GetXYZ(&Point, 0, 0, k); //this height
		vfloat ZHeight = Point.z;

		//clear the previous list:
		int iter = 0;
		while (pZList[iter] != -1 && iter < TotNumFacets)
			pZList[iter++] = -1;

		//add any Facets whose Z coordinates are not all above or all below this Z plane
		int NumZFacets = 0;
		bool IsAbove, IsBelow;

		//TRACE("ZHeight: %g\n", ZHeight);
		for (int m=0; m<TotNumFacets; m++){
			IsAbove = true; IsBelow = true;

			for (int n=0; n<3; n++){
				if(pMeshIn->Vertices[pMeshIn->Facets[m].vi[n]].v.z > ZHeight) IsBelow = false;
				if(pMeshIn->Vertices[pMeshIn->Facets[m].vi[n]].v.z < ZHeight) IsAbove = false;
			}
			if (!IsAbove && !IsBelow){ //if this facet is not fully above or fully below our ZPlane
				pZList[NumZFacets++] = m; //add to out current ZList...
				//TRACE("V1.z = %g, V2.z = %g, V3.z = %g\n", pModel->Vertices[pModel->Facets[m].vi[0]].v.z, pModel->Vertices[pModel->Facets[m].vi[1]].v.z, pModel->Vertices[pModel->Facets[m].vi[2]].v.z);
			}
		}

		
		for (int j=0; j<Yv; j++){ //iterate through y
			pObj->GetXYZ(&Point, 0, j, k); //this location
			CurNumInt = pMeshIn->GetXIntersections(Point.z, Point.y, pIntrscts, NumZFacets, pZList); //fill in pIntrscts

			int IntrInd = 0;
			for (int i=0; i<Xv; i++){ //iterate through x
				pObj->GetXYZ(&Point, i, j, k); //this location

				while (pIntrscts[IntrInd] < Point.x && IntrInd < CurNumInt) //step through our array of intersections until
					IntrInd++;

				if (IntrInd%2 == 1) pObj->SetVoxel(i, j, k, MatIndex); //if we've gone through an odd number of intersections, put material here
				//else pObj->SetVoxel(i, j, k, 0); //otherwise empty...

			}
		}
	}

	delete [] pZList;
	pZList = NULL;
	delete [] pIntrscts;
	pIntrscts = NULL;

	return true;
*/
}


void CVX_MeshUtil::printMeshNormals()
{
	std::cout << " -----------------------------------------------------------"  << std::endl;
	std::cout << "| 		PRINTING DEFORMABLE MESH NORMALS 				  |" << std::endl;				  
	std::cout << " -----------------------------------------------------------"  << std::endl;
	std::vector<CFacet>::iterator VIt;
	for(VIt=DefMesh.Facets.begin(); VIt != DefMesh.Facets.end(); VIt++)
		std::cout << VIt->n.x  << " " << VIt->n.y  << " " << VIt->n.z << std::endl;
	std::cout << " -----------------------------------------------------------"  << std::endl;

}

void CVX_MeshUtil::printMeshVertices()
{
	std::cout << " -----------------------------------------------------------"  << std::endl;
	std::cout << "| 		PRINTING DEFORMABLE MESH VERTICES 				  |" << std::endl;				  
	std::cout << " -----------------------------------------------------------"  << std::endl;
	std::vector<CVertex>::iterator VIt;
	for(VIt=DefMesh.Vertices.begin(); VIt != DefMesh.Vertices.end(); VIt++)
		std::cout << VIt->DrawOffset.x + VIt->v.x  << " " << VIt->DrawOffset.y + VIt->v.y  << " " << VIt->DrawOffset.z + VIt->v.z << std::endl;
	std::cout << " -----------------------------------------------------------"  << std::endl;
}

void CVX_MeshUtil::printMeshFacets()
{
	std::cout << " ---------------------------------------------------------------------"  << std::endl;
	std::cout << "| 	PRINTING DEFORMABLE MESH FACETS (TERNS OF VERTICES INDICES) 	|" << std::endl;				  
	std::cout << " ---------------------------------------------------------------------"  << std::endl;
	std::vector<CFacet>::iterator FIt;
	for(FIt=DefMesh.Facets.begin(); FIt != DefMesh.Facets.end(); FIt++)
		std::cout << FIt->vi[0]  << " " << FIt->vi[1]  << " " << FIt->vi[2] << std::endl;
	std::cout << " -----------------------------------------------------------"  << std::endl;
}

void CVX_MeshUtil::printAllMeshInfo()
{
	printMeshVertices();
	printMeshFacets();
	printMeshNormals();
}


bool CVX_MeshUtil::writeQhullInputFile()
{
	// Let's dump current position of all vertices to a Qhull input file.
	// can be used to invoke qhull and compute the convex hull of the robot at this moment.
	if(pSim->getQhullTmpFile() == "")
	{
		if(usingGUI)
			pSim->setQhullTmpFile(std::string("/tmp/VoxCad_qhull_tmp.txt"));
		else
		{
			std::cout << "[**ERROR CVX_MeshUtil::writeQhullInputFile **] qhull temp file is not defined and we're not using the GUI! Quitting."<< std::endl;
			return 1;
		}
	}

	std::ofstream ofs;
	ofs.open (pSim->getQhullTmpFile().c_str(), std::ofstream::out);// | std::ofstream::app);

	if(ofs.fail())
	{
		std::cout << "ERROR: CVX_MeshUtil::writeQhullInputFile() failed creating a file for qhull!" << std::endl;
		return 1; 
	}

	// Now writing a qhull input file. Later qhull will be invocated on this file to compute the convex hull of the softbot
	int numVertices = getMeshVertNum();

	ofs << "3" << std::endl; 			// Points dimensionality 
	ofs << numVertices << std::endl; 	// Number of points
	// Let's iterate over the mesh and print to file the x,y,z position of each vertex in the Qhull input file format
	std::vector<CVertex>::iterator VIt;
	for(VIt=DefMesh.Vertices.begin(); VIt != DefMesh.Vertices.end(); VIt++)
		ofs << VIt->v.x + VIt->DrawOffset.x << " " << VIt->v.y + VIt->DrawOffset.y << " " << VIt->v.z + VIt->DrawOffset.z << std::endl;


// FC: OLD VERSION USING VOXELS' CENTERS. NOW USING DEFORMED MESH's VERTICES FOR BETTER ACCURACY	
//	int nVx = NumVox();
//	ofs << "3" << std::endl; 	 // Point dimensionality
//	ofs << nVx << std::endl;  // Points number

//	// Now iterating on all voxels and append the x,y,z coordinates to the file

//	for (int i = 0; i < nVx; i++) // TODO: INSTEAD OF WRITING CENTERS, APPEND ALL VERTICES OF THIS VOXEL, IN THEIR DEFORMED STATE
//		ofs << VoxArray[i].GetCurPos().x << " " << VoxArray[i].GetCurPos().y << " " << VoxArray[i].GetCurPos().z << std::endl;

	ofs.close();
	return 0;

}

double CVX_MeshUtil::invokeQhull()
{
	if(pSim->getQhullTmpFile() == "")
	{	
		if(usingGUI)
			pSim->setQhullTmpFile(std::string("/tmp/VoxCad_qhull_tmp.txt"));
		else
		{
			std::cout << "WARNING: CVX_MeshUtil::invokeQhull() called but pSim->getQhullTmpFile() not defined. Check your vxa file for the relevant tag." << std::endl;
			return -1;
		}

	}

	if(writeQhullInputFile() != 0)
	{
		std::cout << "CVX_MeshUtil ERROR: writeQhullInputFile failed in invokeQhull function. Returning -1 volume." << std::endl;
		return -1;
	}

	FILE *fp;
	char* r;
	char buf[1024];
	
	int maxAttempts = 60;
	int attempts = 0;
	double result = -1;

	while(attempts++ < maxAttempts)
	{
		//cout << "Invoking qhull to compute convex hull volume - attempt " << attempts << endl;
		system("sleep 1"); // wait a bit for the file to be created and populated by writeQhull, just in case

		std::string qhullCmd1, qhullCmd2;
		qhullCmd1 = "./qhull FS TI "+pSim->getQhullTmpFile()+" | sed 1d | tr -s \\ | cut -f3 -d \\ "; // invoke qhull with input file QhullTmpFile and get the volume of the convex hull
		qhullCmd2 = "qhull FS TI "+pSim->getQhullTmpFile()+" | sed 1d | tr -s \\ | cut -f3 -d \\ "; // invoke qhull with input file QhullTmpFile and get the volume of the convex hull

		fp = popen(qhullCmd1.c_str(), "r");
		r = fgets(buf, 1024, fp);  
		if (!fp || !r)
		{
			fp = popen(qhullCmd2.c_str(), "r");
			r = fgets(buf, 1024, fp);  
			if (!fp || !r)
				continue;
		}

		result = atof(buf); // we were able to read. Let's convert to double

/*		if(result <= 0) // if atof failed
		{
			std::cout << "CVX_MeshUtil ERROR: Qhull returned a volume <= 0. Probably an error" << std::endl;
			return -1;		
		}
		else
			break;
*/
		fclose(fp);
		break;
		
	}

	if(result <= 0)
	{
		std::cout << "CVX_MeshUtil ERROR: Cannot get a valid result from qhull in 20 attempts. Returning -1." << std::endl;
		result = -1;
	}
	
//	cout << "The value returned by qhull is: " << result << endl;	
	// we can now delete the QhullTmpFile
	std::string removeCmd;
	removeCmd = "rm \""+pSim->getQhullTmpFile()+"\"";
	system(removeCmd.c_str());
	return result;

}

double CVX_MeshUtil::computeCurrentRobotVolume()
{

	// FC
	// Volume computation based on the following paper:
	// Zhang, Cha, and Tsuhan Chen. "Efficient feature extraction for 2D/3D objects in mesh representation."
	// Image Processing, 2001. Proceedings. 2001 International Conference on. Vol. 3. IEEE, 2001.
	// 
	// See also http://n-e-r-v-o-u-s.com/blog/?p=4415
	// for more concise formulas
	// 
	updateDeformableMesh(); 

	double volume = 0.0;

	// Defining useful iterators, just for convenience
	std::vector<CFacet>::iterator FIt;
	std::vector<CVertex>::iterator V = DefMesh.Vertices.begin();

	for(FIt=DefMesh.Facets.begin(); FIt != DefMesh.Facets.end(); FIt++)
	{
		// Vertices indices for this facet
		int idxV1, idxV2, idxV3;
		idxV1 = FIt->vi[0];
		idxV2 = FIt->vi[1];
		idxV3 = FIt->vi[2];

		// Let's get the vectors associated with the vertices of this facet
		Vec3D<> V1 = V[idxV1].v + V[idxV1].DrawOffset;
		Vec3D<> V2 = V[idxV2].v + V[idxV2].DrawOffset;
		Vec3D<> V3 = V[idxV3].v + V[idxV3].DrawOffset;		

		volume += (1.0/6.0)*((V1.Cross(V2)).Dot(V3));

		//
		// If we want to be sure that the sign is correct, we can try using explicitly the normal, this way (pseudocode):
		// Vec3D<> N = FIt->n; // facet normal
		// s = sign(dot(V1, N))
		// volume += (1.0/6.0)*s*abs(((V1.Cross(V2)).Dot(V3)));
		//
		// Looks like it's not necessary here, as vertices appear to be ordered in a conventional manner
		//
	}

	return volume;
}

double CVX_MeshUtil::computeShapeComplexity(std::string tmpFileName)
{
	int vertexCount = DefMesh.Vertices.size();
	double angleExcesses[vertexCount]; // This will contain the approximation of curvatures at each vertex.


	// For each vertex we need to find all triangles it is part of, 
	// and compute angle at that vertex
	for(int V = 0; V < vertexCount; V++)
	{
		//double area = 0;
		double angleSum = 0;

		// Now we iterate over facets to find triangles enclosing this vertex
		std::vector<CFacet>::iterator F;
		for(F=DefMesh.Facets.begin(); F != DefMesh.Facets.end(); F++)
		{			
			// F is a facet, and the following are vertices indices composing the facet
			int a = F->vi[0];
			int b = F->vi[1];
			int c = F->vi[2];
			// we need to find if 'v' is one of the three, and which one, 
			// so that we can measure the angle at 'v'
				
			// Vectors associated to each vertex, in the global reference frame
			Vec3D<> va, vb, vc;
			va = DefMesh.Vertices[a].v + DefMesh.Vertices[a].DrawOffset;
			vb = DefMesh.Vertices[b].v + DefMesh.Vertices[b].DrawOffset;
			vc = DefMesh.Vertices[c].v + DefMesh.Vertices[c].DrawOffset;

			// The angle at V will be the angle between v1 and v2, once they are correctly defined
			// bsed on va,vb,vc
			Vec3D<> v1, v2;

			if(V == a)
			{
				v1 = vb-va;
				v2 = vc-va;			
			}
			else if (V == b)
			{
			
				v1 = va-vb;
				v2 = vc-vb;
			}
			else if (V == c)
			{
				v1 = va-vc;
				v2 = vb-vc;
			}
			else
				continue;

			// Let's normalize v1 and v2
			v1.Normalize();
			v2.Normalize();

			// Now we just need to compute the angle between v1 and v2 (a12)...
			double a12 = acos(v1.Dot(v2));

			// ... and accumulate it
			angleSum += a12;
		}

		// Finally, angle excess for V is
		angleExcesses[V] = (2.0*PI) - angleSum;
	}

	// Time to print angleExcesses to file 'tmpFileName'

	std::ofstream curvatureFile;
	curvatureFile.open(tmpFileName.c_str());

	if(!curvatureFile.is_open())
	{
		std::cout << "[V_MeshUtil.cpp] Cannot open file "<< tmpFileName << " to log mesh curvatures!" << std::endl;
		return -1;
	}	


	for(int V = 0; V < vertexCount; V++)
		curvatureFile << angleExcesses[V] << "\t";		

	curvatureFile.close();

	// Call python program that will compute the entropy of the values in tmpFileName and override the file with the compyted entropy
	char command[200];
	sprintf(command,"python curvatureEntropy.py %s", tmpFileName.c_str());
	system(command);
	
	double curvatureEntropy = -1.0;
	

	int maxAttempts = 60;
	int attempts = 0;

	while(attempts++ < maxAttempts)
	{
		system("sleep 1");

		std::ifstream entropyIn(tmpFileName.c_str());

		if (!entropyIn.is_open())
			continue;		

		entropyIn >> curvatureEntropy;
		break;
	}


	if(curvatureEntropy <= 0)
		std::cout << "[V_MeshUtil.cpp] [WARNING:] Could not get an admissible entropy value from python script curvatureEntropy.py in 60 attempts. CVX_MeshUtil::computeShapeComplexity" << std::endl;

	// Tmp file can now be removed.
	char command2[200];
	sprintf(command2,"rm %s", tmpFileName.c_str());
	system(command2);

	return curvatureEntropy;
}

double CVX_MeshUtil::computeInitialShapeComplexity()
{
	if(pSim)
	{
		if(pSim->getCurvaturesTmpFile() == "")
		{			
			if(usingGUI)
				pSim->setCurvaturesTmpFile(std::string("/tmp/VoxCad_curvatures_tmp.txt"));
			else
			{
				std::cout << "[VX_MeshUtil.cpp] ERROR: curvatures tmp filename is empty, and not using a GUI. Returning -1!" << std::endl;
				return -1;
			}
		}
		shapeComplexityStart = computeShapeComplexity(pSim->getCurvaturesTmpFile());
		return shapeComplexityStart;
	}

	std::cout << "[VX_MeshUtil.cpp] ERROR: Pointer to pSim null! Returning -1!" << std::endl;

	return -1;
}

double CVX_MeshUtil::computeFinalShapeComplexity()
{
	if(pSim)
	{
		if(pSim->getCurvaturesTmpFile() == "")
			if(usingGUI)
				pSim->setCurvaturesTmpFile(std::string("/tmp/VoxCad_curvatures_tmp.txt"));
			else
			{
				std::cout << "[VX_MeshUtil.cpp] ERROR: curvatures tmp filename is empty, and not using a GUI. Returning -1!" << std::endl;
				return -1;
			}

		shapeComplexityEnd = computeShapeComplexity(pSim->getCurvaturesTmpFile());
		return shapeComplexityEnd;
	}

	std::cout << "[VX_MeshUtil.cpp] ERROR: Pointer to pSim null! Returning -1!" << std::endl;
	return -1;
}