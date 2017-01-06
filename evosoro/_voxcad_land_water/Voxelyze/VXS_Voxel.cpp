/*******************************************************************************
Copyright (c) 2010, Jonathan Hiller (Cornell University)
If used in publication cite "J. Hiller and H. Lipson "Dynamic Simulation of Soft Heterogeneous Objects" In press. (2011)"

This file is part of Voxelyze.
Voxelyze is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Voxelyze is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
See <http://www.opensource.org/licenses/lgpl-3.0.html> for license details.
*******************************************************************************/

#include "VXS_Voxel.h"
#include "VXS_Bond.h"
#include "VX_Sim.h"
#include <iostream>

CVXS_Voxel::CVXS_Voxel(CVX_Sim* pSimIn, int SIndexIn, int XIndexIn, int MatIndexIn, Vec3D<>& NominalPositionIn, vfloat OriginalScaleIn) : CVX_Voxel(pSimIn, SIndexIn, XIndexIn, MatIndexIn, NominalPositionIn, OriginalScaleIn)
{
	ResetVoxel(); //sets state variables to zero


	ExternalInputScale = 1.0;

//	CornerPosCur = Vec3D<>(OriginalScaleIn/2, OriginalScaleIn/2, OriginalScaleIn/2);
//	CornerNegCur = Vec3D<>(-OriginalScaleIn/2, -OriginalScaleIn/2, -OriginalScaleIn/2);

	SetColor(0,0,0,1);

	// Here we should be able to override Vox_E (elastic modulus)
	// if(pSimIn->pEnv->pObj->GetEvolvingStiffness())
	// {

	// 	std::cout << "Material stiffness should be " << Vox_E << std::endl;
	// 	Vox_E = evolvedStiffness; // already assigned to the voxel in VX_Sim.cpp
	// 	std::cout << "Evolving stiffness, CPPN-dictated one is "<< Vox_E << std::endl;		
	// }	

}

CVXS_Voxel::~CVXS_Voxel(void)
{
}

CVXS_Voxel& CVXS_Voxel::operator=(const CVXS_Voxel& VIn)
{
	CVX_Voxel::operator=(VIn);
	
	ExternalInputScale=VIn.ExternalInputScale;
	
	InputForce = Vec3D<>(0,0,0);

	Pos = VIn.Pos;
	LinMom = VIn.LinMom;
	Angle = VIn.Angle;
	AngMom = VIn.AngMom;
	Vel = VIn.Vel;
	KineticEnergy = VIn.KineticEnergy;
	AngVel = VIn.AngVel;
	Pressure = VIn.Pressure;
	Scale=VIn.Scale;

	StaticFricFlag = VIn.StaticFricFlag;
	VYielded = VIn.VYielded;
	VBroken = VIn.VBroken;

	ColBondInds = VIn.ColBondInds;
	UpdateColBondPointers();

	m_Red = VIn.m_Red;
	m_Green = VIn.m_Green;
	m_Blue = VIn.m_Blue;
	m_Trans = VIn.m_Trans;

//	SizeCurrent = VIn.SizeCurrent;
	CornerPosCur = VIn.CornerPosCur;
	CornerNegCur = VIn.CornerNegCur;

	ForceCurrent = VIn.ForceCurrent;

	return *this;
}

void CVXS_Voxel::ResetVoxel(void) //resets this voxel to its defualt (imported) state.
{
	LinMom = Vec3D<double>(0,0,0);
	Angle = CQuat<double>(1.0, 0, 0, 0);
	AngMom = Vec3D<double>(0,0,0);
	Scale = 0;
	Vel = Vec3D<>(0,0,0);
	KineticEnergy = 0;
	AngVel = Vec3D<>(0,0,0);
	Pressure=0;


	Pos = GetNominalPosition(); //only position and size need to be set
	Scale = GetNominalSize();


	InputForce = Vec3D<>(0,0,0); //?


	StaticFricFlag = false;
	VYielded = false;
	VBroken = false;

//	SizeCurrent = Vec3D<>(Scale, Scale, Scale);
	CornerPosCur = Vec3D<>(Scale/2, Scale/2, Scale/2);
	CornerNegCur = Vec3D<>(-Scale/2, -Scale/2, -Scale/2);
//	CornerPosCur = Vec3D<>(0,0,0);
//	CornerNegCur = Vec3D<>(0,0,0);
	ForceCurrent = Vec3D<>(0,0,0);
	StrainPosDirsCur = Vec3D<>(0,0,0);
	StrainNegDirsCur = Vec3D<>(0,0,0);
}


bool CVXS_Voxel::LinkColBond(int CBondIndex) //simulation bond index...
{
	if (!pSim || CBondIndex >= pSim->BondArrayCollision.size()) return false;

	ColBondInds.push_back(CBondIndex);
	ColBondPointers.push_back(&(pSim->BondArrayCollision[CBondIndex]));

	return true;
}

void CVXS_Voxel::UnlinkColBonds(void)
{
	ColBondInds.clear();
	ColBondPointers.clear();
}


void CVXS_Voxel::UpdateColBondPointers() //updates all links (pointers) to bonds according top current p_Sim
{
	int NumColBonds = ColBondInds.size();
	if (NumColBonds == 0) return;
	ColBondPointers.resize(NumColBonds);
	for (int i=0; i<NumColBonds; i++){
		ColBondPointers[i] = &(pSim->BondArrayCollision[ColBondInds[i]]);
	}
}





//http://klas-physics.googlecode.com/svn/trunk/src/general/Integrator.cpp (reference)
void CVXS_Voxel::EulerStep()
{
	double dt = pSim->dt;
	//bool EqMode = p_Sim->IsEquilibriumEnabled();
	if (IS_ALL_FIXED(DofFixed) & !pSim->IsFeatureEnabled(VXSFEAT_VOLUME_EFFECTS)){ //if fixed, just update the position and forces acting on it (for correct simulation-wide summing
		LinMom = Vec3D<double>(0,0,0);
		Pos = NominalPosition + ExternalInputScale*ExternalDisp;
		AngMom = Vec3D<double>(0,0,0);
		Angle.FromRotationVector(Vec3D<double>(ExternalInputScale*ExternalTDisp));
	}
	else {
		Vec3D<> ForceTot = CalcTotalForce(); //TotVoxForce;

		//DISPLACEMENT
		LinMom = LinMom + ForceTot*dt;
		Vec3D<double> Disp(LinMom*(dt*_massInv)); //vector of what the voxel moves

//		if(pSim->IsMaxVelLimitEnabled()){ //check to make sure we're not going over the speed limit!
		if(pSim->IsFeatureEnabled(VXSFEAT_MAX_VELOCITY)){ //check to make sure we're not going over the speed limit!
			vfloat DispMag = Disp.Length();
			vfloat MaxDisp = pSim->GetMaxVoxVelLimit()*NominalSize; // p_Sim->pEnv->pObj->GetLatticeDim();
			if (DispMag>MaxDisp) Disp *= (MaxDisp/DispMag);
		}
		Pos += Disp; //update position (source of noise in float mode???

		if (IS_FIXED(DOF_X, DofFixed)){Pos.x = NominalPosition.x + ExternalInputScale*ExternalDisp.x; LinMom.x = 0;}
		if (IS_FIXED(DOF_Y, DofFixed)){Pos.y = NominalPosition.y + ExternalInputScale*ExternalDisp.y; LinMom.y = 0;}
		if (IS_FIXED(DOF_Z, DofFixed)){Pos.z = NominalPosition.z + ExternalInputScale*ExternalDisp.z; LinMom.z = 0;}

		//ANGLE
		Vec3D<> TotVoxMoment = CalcTotalMoment(); //debug
		



		AngMom = AngMom + TotVoxMoment*dt;

		if (pSim->IsFeatureEnabled(VXSFEAT_VOLUME_EFFECTS)) AngMom /= 1.01; //TODO: remove angmom altogehter???
		else {
			vfloat AngMomFact = (1 - 10*pSim->GetSlowDampZ() * _inertiaInv *_2xSqIxExSxSxS*dt);
			AngMom *= AngMomFact; 
		}
//		if ((AngMomXPos && AngMom.x < 0) || (!AngMomXPos && AngMom.x > 0)) AngMom.x = 0;
//		if ((AngMomYPos && AngMom.y < 0) || (!AngMomYPos && AngMom.y > 0)) AngMom.y = 0;
//		if ((AngMomZPos && AngMom.z < 0) || (AngMomZNeg && AngMom.z > 0)) AngMom.z = 0;


		//convert Angular velocity to quaternion form ("Spin")
		Vec3D<double> dSAngVel(AngMom * _inertiaInv);
		CQuat<double> Spin = 0.5 * CQuat<double>(0, dSAngVel.x, dSAngVel.y, dSAngVel.z) * Angle; //current "angular velocity"

		Angle += CQuat<double>(Spin*dt); //see above
		Angle.NormalizeFast(); //Through profiling, quicker to normalize every time than check to see if needed then do it...

	//	TODO: Only constrain fixed angles if one is non-zero! (support symmetry boundary conditions while still only doing this calculation) (only works if all angles are constrained for now...)
		if (IS_FIXED(DOF_TX, DofFixed) && IS_FIXED(DOF_TY, DofFixed) && IS_FIXED(DOF_TZ, DofFixed)){
			Angle.FromRotationVector(Vec3D<double>(ExternalInputScale*ExternalTDisp));
			AngMom = Vec3D<>(0,0,0);
		}
	}

	//SCALE
	//	ScaleMom = ScaleMom + CalcTotalScaleForce()*p_Sim->dt;
	vfloat TempFact = 1.0;
	if(pSim->IsFeatureEnabled(VXSFEAT_TEMPERATURE) and pSim->CurTime >= pSim->GetInitCmTime()){ 
		//TempFact = (1+(p_Sim->pEnv->CurTemp-p_Sim->pEnv->TempBase)*GetCTE()); //LocalVXC.GetBaseMat(VoxArray[i].MatIndex)->GetCTE());
		
//		vfloat ThisTemp = p_Sim->pEnv->pObj->GetBaseMat(GetMaterial())->GetCurMatTemp();
//		vfloat ThisTemp = pSim->pEnv->pObj->GetBaseMat(GetMaterialIndex())->GetCurMatTemp();
		vfloat ThisTemp = _pMat->GetCurMatTemp();

		vfloat ThisCTE = GetCTE();
		vfloat TempBase =  pSim->pEnv->GetTempBase();

		// std::cout << "Scale: " << Simulator.VoxArray[0].GetCurScale() << std::endl;
		// std::cout << "TempAmp: " << TempAmplitude << std::endl;
		// std::cout << "TempPer: " << TempPeriod << std::endl;
		// std::cout << "TempAmp: " << pSim->pEnv->GetTempAmplitude() << std::endl;
		// std::cout << "TempPer: " << pSim->pEnv->GetTempPeriod() << std::endl;
		// std::cout << "phase offset: " << phaseOffset << std::endl;

		// TempFact = (1+(ThisTemp - TempBase)*ThisCTE);	//To allow selective temperature actuation for each different material
		TempFact = (1+(TempAmplitude*sin(2*3.1415926f * (pSim->CurTime/TempPeriod + phaseOffset)))*ThisCTE); 
		// if (pEnv->TempPeriod == 0) return; //avoid NaNs.
		// 	CurTemp = TempBase + TempAmplitude*sin(2*3.1415926/TempPeriod*time);
	}
	if (TempFact < MIN_TEMP_FACTOR) TempFact = MIN_TEMP_FACTOR;
	Scale = TempFact*NominalSize;
	

	//Recalculate secondary:
	AngVel = AngMom * _inertiaInv;
	Vel = LinMom * _massInv;
	if(pSim->StatToCalc & CALCSTAT_KINE) KineticEnergy = 0.5*Mass*Vel.Length2() + 0.5*Inertia*AngVel.Length2(); //1/2 m v^2
	if(pSim->StatToCalc & CALCSTAT_PRESSURE) {
		//vfloat VolumetricStrain = StrainPosDirsCur.x/2 + StrainPosDirsCur.y/2 + StrainPosDirsCur.z/2 + StrainNegDirsCur.x/2 + StrainNegDirsCur.y/2 + StrainNegDirsCur.z/2;
		vfloat VolumetricStrain = GetVoxelStrain(AXIS_X) + GetVoxelStrain(AXIS_Y) + GetVoxelStrain(AXIS_Z);
		Pressure = - Vox_E*VolumetricStrain/(3*(1-2*_pMat->GetPoissonsRatio())); //http://www.colorado.edu/engineering/CAS/courses.d/Structures.d/IAST.Lect05.d/IAST.Lect05.pdf
	}
}


void CVXS_Voxel::SetColor(float r, float g, float b, float a)
{
	m_Red = r;
	m_Green = g;
	m_Blue = b;
	m_Trans = a;
}	

void CVXS_Voxel::SetStrainDir(BondDir Bond, vfloat StrainIn)
{
	switch (Bond){
	case BD_PX: StrainPosDirsCur.x = StrainIn; break;
	case BD_PY: StrainPosDirsCur.y = StrainIn; break;
	case BD_PZ: StrainPosDirsCur.z = StrainIn; break;
	case BD_NX: StrainNegDirsCur.x = StrainIn; break;
	case BD_NY: StrainNegDirsCur.y = StrainIn; break;
	case BD_NZ: StrainNegDirsCur.z = StrainIn; break;
	}
}

vfloat CVXS_Voxel::GetVoxelStrain(Axis DesiredAxis)
{
	bool pd, nd; //positive and negative directions
	switch (DesiredAxis){
		case AXIS_X:
			pd = InternalBondPointers[BD_PX]!=NULL, nd = InternalBondPointers[BD_NX]!=NULL;
			if (!pd && !nd) return 0;
			else if (pd && !nd) return StrainPosDirsCur.x;
			else if (!pd && nd) return StrainNegDirsCur.x;
			else return 0.5*(StrainPosDirsCur.x + StrainNegDirsCur.x);
			break;
		case AXIS_Y:
			pd = InternalBondPointers[BD_PY]!=NULL, nd = InternalBondPointers[BD_NY]!=NULL;
			if (!pd && !nd) return 0;
			else if (pd && !nd) return StrainPosDirsCur.y;
			else if (!pd && nd) return StrainNegDirsCur.y;
			else return 0.5*(StrainPosDirsCur.y + StrainNegDirsCur.y);
			break;
		case AXIS_Z:
			pd = InternalBondPointers[BD_PZ]!=NULL, nd = InternalBondPointers[BD_NZ]!=NULL;
			if (!pd && !nd) return 0;
			else if (pd && !nd) return StrainPosDirsCur.z;
			else if (!pd && nd) return StrainNegDirsCur.z;
			else return 0.5*(StrainPosDirsCur.z + StrainNegDirsCur.z);
			break;
		default: return 0;

	}


}

Vec3D<> CVXS_Voxel::CalcTotalForce()
{
//	THE NEXT optimization target
	//INTERNAL forces
	Vec3D<> TotalForce = Vec3D<>(0,0,0);
	TotalForce += -pSim->GetSlowDampZ() * Vel*_2xSqMxExS; //(2*sqrt(Mass*GetEMod()*Scale.x));

	//POSSIONS!
//	Vec3D<> pStrain = Vec3D<>(0,0,0), nStrain = Vec3D<>(0,0,0),
	Vec3D<> CurLocStrain = Vec3D<>(0,0,0);
	//End Poissons


	//Forces from permanent bonds:
	for (int i=0; i<6; i++){
		CVXS_Bond* pThisBond = InternalBondPointers[i];
		if (!pThisBond) continue;

		if (IAmInternalVox2(i)) TotalForce += pThisBond->GetForce2();
		else TotalForce += pThisBond->GetForce1();


		//if (pSim->IsFeatureEnabled(VXSFEAT_VOLUME_EFFECTS)){
		//	bool IAmVox1 = !IAmInternalVox2(i);
		//	switch (pThisBond->GetBondAxis()){
		//		case AXIS_X: IAmVox1 ? pStrain.x = pThisBond->GetStrainV1() : nStrain.x = pThisBond->GetStrainV2(); break; 
		//		case AXIS_Y: IAmVox1 ? pStrain.y = pThisBond->GetStrainV1() : nStrain.y = pThisBond->GetStrainV2(); break;
		//		case AXIS_Z: IAmVox1 ? pStrain.z = pThisBond->GetStrainV1() : nStrain.z = pThisBond->GetStrainV2(); break;
		//	}
		//}

	}

	//Forces from collision bonds: To optimize!
	if (pSim->IsFeatureEnabled(VXSFEAT_COLLISIONS)){
		int NumColBond = ColBondPointers.size();
		for (int i=0; i<NumColBond; i++){
			if (IAmVox2Col(i)) TotalForce += ColBondPointers[i]->GetForce2();
			else TotalForce += ColBondPointers[i]->GetForce1();

//			CVXS_Bond* pThisBond = ColBondPointers[i];
//			bool IAmVox1 = !IAmVox2Col(i); //IsMe(pThisBond->GetpV1()); //otherwise vox 2 of the bond

//			if (IAmVox1) TotalForce -= pThisBond->GetForce1(); //Force on Vox 1 from this bond
//			else TotalForce -= pThisBond->GetForce2(); //Force on Vox 2 from this bond
		}
	}

	//Forced from input bond
	TotalForce -= InputForce;



	//From gravity
	// if (pSim->IsFeatureEnabled(VXSFEAT_GRAVITY))
	if (pSim->IsFeatureEnabled(VXSFEAT_GRAVITY) and not pSim->fluidEnvironment) // FC: assuming neutral buoyancy, simply disabling gravity
		TotalForce.z += Mass*pSim->pEnv->GetGravityAccel();

	//EXTERNAL forces
	TotalForce += ExternalInputScale*ExternalForce; //add in any external forces....

	// // // nac: water:
	//std::cout << "TotalForce: " << TotalForce.x << ", " << TotalForce.y << ", " << TotalForce.z << " --> " << TotalForce.Length() << std::endl;
	//std::cout << "DragForce: " << DragForce.x << ", " << DragForce.y << ", " << DragForce.z << " --> " << DragForce.Length() << std::endl;
	//std::cout << std::endl;

//	std::cout << "TotalForce:" << TotalForce.Length() << std::endl;
//	std::cout << "Drag Force:" << DragForce.Length() << std::endl;
	if (pSim->fluidEnvironment)
		TotalForce += DragForce;


	if (pSim->IsFeatureEnabled(VXSFEAT_VOLUME_EFFECTS)){
		//http://www.colorado.edu/engineering/CAS/courses.d/Structures.d/IAST.Lect05.d/IAST.Lect05.pdf

		vfloat mu = GetPoisson();

		bool px = (InternalBondPointers[BD_PX] != NULL);
		bool nx = (InternalBondPointers[BD_NX] != NULL);
		bool py = (InternalBondPointers[BD_PY] != NULL);
		bool ny = (InternalBondPointers[BD_NY] != NULL);
		bool pz = (InternalBondPointers[BD_PZ] != NULL);
		bool nz = (InternalBondPointers[BD_NZ] != NULL);

		bool Tx = px && nx || ((px || nx) && (IS_FIXED(DOF_X, DofFixed) || ExternalForce.x != 0)); //if bond on both sides or pulling against a fixed or forced constraint
		bool Ty = py && ny || ((py || ny) && (IS_FIXED(DOF_Y, DofFixed) || ExternalForce.y != 0)); //if bond on both sides or pulling against a fixed or forced constraint
		bool Tz = pz && nz || ((pz || nz) && (IS_FIXED(DOF_Z, DofFixed) || ExternalForce.z != 0)); //if bond on both sides or pulling against a fixed or forced constraint
		if (Tx){
			CurLocStrain.x = GetVoxelStrain(AXIS_X);
//			if (px && !nx) CurLocStrain.x = StrainPosDirsCur.x;
//			else if (!px && nx) CurLocStrain.x = StrainNegDirsCur.x;
//			else CurLocStrain.x = 0.5*(StrainPosDirsCur.x + StrainNegDirsCur.x);
		}
		if (Ty){
			CurLocStrain.y = GetVoxelStrain(AXIS_Y);

//			if (py && !ny) CurLocStrain.y = StrainPosDirsCur.y;
//			else if (!py && ny) CurLocStrain.y = StrainNegDirsCur.y;
//			else CurLocStrain.y = 0.5*(StrainPosDirsCur.y + StrainNegDirsCur.y);
		}
		if (Tz){
			CurLocStrain.z = GetVoxelStrain(AXIS_Z);

//			if (pz && !nz) CurLocStrain.z = StrainPosDirsCur.z;
//			else if (!pz && nz) CurLocStrain.z = StrainNegDirsCur.z;
//			else CurLocStrain.z = 0.5*(StrainPosDirsCur.z + StrainNegDirsCur.z);
		}


		if (!Tx && !Ty && !Tz) CurLocStrain = Vec3D<>(0,0,0); //if nothing pushing or pulling, no strain on this bond!
		//else if (!Tx && Ty && Tz) CurLocStrain.x = pow(1+CurLocStrain.y, -mu)-1 + pow(1+CurLocStrain.z, -mu)-1;
		//else if (Tx && !Ty && Tz) CurLocStrain.y = pow(1+CurLocStrain.x, -mu)-1 + pow(1+CurLocStrain.z, -mu)-1;
		//else if (Tx && Ty && !Tz) CurLocStrain.z = pow(1+CurLocStrain.x, -mu)-1 + pow(1+CurLocStrain.y, -mu)-1;
		else if (!Tx && Ty && Tz) CurLocStrain.x = pow(1+CurLocStrain.y + CurLocStrain.z, -mu)-1;
		else if (Tx && !Ty && Tz) CurLocStrain.y = pow(1+CurLocStrain.x + CurLocStrain.z, -mu)-1; //??
		else if (Tx && Ty && !Tz) CurLocStrain.z = pow(1+CurLocStrain.x + CurLocStrain.y, -mu)-1;
		else if (!Tx && !Ty && Tz) CurLocStrain.x = CurLocStrain.y = pow(1+CurLocStrain.z, -mu)-1;
		else if (!Tx && Ty && !Tz) CurLocStrain.x = CurLocStrain.z = pow(1+CurLocStrain.y, -mu)-1;
		else if (Tx && !Ty && !Tz) CurLocStrain.y = CurLocStrain.z = pow(1+CurLocStrain.x, -mu)-1;
		//else if (Tx && Ty && Tz) //we already have everything!

//		SizeCurrent.x = (1+CurLocStrain.x)*NominalSize;
//		SizeCurrent.y = (1+CurLocStrain.y)*NominalSize;
//		SizeCurrent.z = (1+CurLocStrain.z)*NominalSize;


		//TODO: get down to a force for this bond, then dump it (and current stiffness) to the bond
//		for (int i=0; i<NumLocBond; i++){
		for (int i=0; i<6; i++){
//			CVXS_Bond* pThisBond = GetBond(i);
			CVXS_Bond* pThisBond = InternalBondPointers[i];
			if (!pThisBond) continue;

			bool IAmVox1 = !IAmInternalVox2(i); //IsMe(pThisBond->GetpV1()); //otherwise vox 2 of the bond
			switch (pThisBond->GetBondAxis()){
			case AXIS_X:
				if (IAmVox1) {pThisBond->TStrainSum1 = CurLocStrain.y + CurLocStrain.z; pThisBond->CSArea1 = (1+CurLocStrain.y)*(1+CurLocStrain.z)*NominalSize*NominalSize;}
				else {pThisBond->TStrainSum2 = CurLocStrain.y + CurLocStrain.z; pThisBond->CSArea2 = (1+CurLocStrain.y)*(1+CurLocStrain.z)*NominalSize*NominalSize;}
				break;
			case AXIS_Y:
				if (IAmVox1) {pThisBond->TStrainSum1 = CurLocStrain.x + CurLocStrain.z; pThisBond->CSArea1 = (1+CurLocStrain.x)*(1+CurLocStrain.z)*NominalSize*NominalSize;}
				else {pThisBond->TStrainSum2 = CurLocStrain.x + CurLocStrain.z; pThisBond->CSArea2 = (1+CurLocStrain.x)*(1+CurLocStrain.z)*NominalSize*NominalSize;}
				break;
			case AXIS_Z:
				if (IAmVox1) {pThisBond->TStrainSum1 = CurLocStrain.y + CurLocStrain.x;  pThisBond->CSArea1 = (1+CurLocStrain.y)*(1+CurLocStrain.x)*NominalSize*NominalSize;}
				else {pThisBond->TStrainSum2 = CurLocStrain.y + CurLocStrain.x;  pThisBond->CSArea2 = (1+CurLocStrain.y)*(1+CurLocStrain.x)*NominalSize*NominalSize;}
				break;
			}
		}
	}
	else { //volume effects off
//		SizeCurrent = Vec3D<>(NominalSize, NominalSize, NominalSize);
		//for (int i=0; i<NumLocBond; i++){
		//	CVXS_Bond* pThisBond = GetBond(i);
		for (int i=0; i<6; i++){ //update for collision bonds?
			CVXS_Bond* pThisBond = InternalBondPointers[i];
			if (pThisBond){
				pThisBond->CSArea1 = pThisBond->CSArea2 = NominalSize*NominalSize;
				//pThisBond->CSArea2 = NominalSize*NominalSize;
			}
		}
	}

//	if(pSim->pEnv->IsFloorEnabled()){
	// if(pSim->IsFeatureEnabled(VXSFEAT_FLOOR))
	if(pSim->IsFeatureEnabled(VXSFEAT_FLOOR) and not pSim->fluidEnvironment)
	{
		TotalForce += CalcFloorEffect(Vec3D<vfloat>(TotalForce));
		//else StaticFricFlag = false;
		if (StaticFricFlag) {TotalForce.x = 0; TotalForce.y = 0;} //no lateral movement if static friction in effect
	}

	CornerPosCur = (Vec3D<>(1,1,1)+StrainPosDirsCur)*NominalSize/2;
	CornerNegCur = -(Vec3D<>(1,1,1)+StrainNegDirsCur)*NominalSize/2;

	//Enforce fixed degrees of freedom (put no force on them so they don't move)
//	if (IS_FIXED(DOF_X, DofFixed) && WithRestraint) TotalForce.x=0;
//	if (IS_FIXED(DOF_Y, DofFixed) && WithRestraint) TotalForce.y=0;
//	if (IS_FIXED(DOF_Z, DofFixed) && WithRestraint) TotalForce.z=0;

	ForceCurrent=TotalForce;
	return ForceCurrent;
}

Vec3D<> CVXS_Voxel::CalcTotalMoment(void)
{
	Vec3D<> TotalMoment(0,0,0);
//	for (int i=0; i<GetNumLocalBonds(); i++) {
	//permanent bonds
	for (int i=0; i<6; i++){ //update for collision bonds?
		CVXS_Bond* pThisBond = InternalBondPointers[i];
		if (pThisBond){
			if (IAmInternalVox2(i)){ TotalMoment -= InternalBondPointers[i]->GetMoment2(); } //if this is voxel 2		//add moments from bond
			else { TotalMoment -= InternalBondPointers[i]->GetMoment1(); } //if this is voxel 1
		}
	}

	//EXTERNAL moments
//	TotalMoment += -pSim->GetSlowDampZ() * AngVel *_2xSqIxExSxSxS; 
	TotalMoment += ExternalInputScale*ExternalTorque; //add in any external forces....

	if (IS_FIXED(DOF_TX, DofFixed)) TotalMoment.x=0;
	if (IS_FIXED(DOF_TY, DofFixed)) TotalMoment.y=0;
	if (IS_FIXED(DOF_TZ, DofFixed)) TotalMoment.z=0;

	return TotalMoment;
}

vfloat CVXS_Voxel::GetCurGroundPenetration() //how far into the ground penetrating (penetration is positive, no penetration is zero)
{
	vfloat Penetration = 0.5*Scale - Pos.z;
	return Penetration <= 0 ? 0 : Penetration;
}

inline bool CVXS_Voxel::IAmVox2Col(const int BondDirColIndex) const 
{
	return (this == ColBondPointers[BondDirColIndex]->GetpV2());
} //returns true if this voxel is Vox2 of the specified bond


Vec3D<> CVXS_Voxel::CalcFloorEffect(Vec3D<> TotalVoxForce) //calculates the object's interaction with a floor. should be calculated AFTER all other forces for static friction to work right...
{
	Vec3D<> FloorForce(0,0,0); //the force added by floor interactions...

	StaticFricFlag = false; //assume not under static friction unless we decide otherwise
	vfloat CurPenetration = GetCurGroundPenetration();


	if (CurPenetration>0){ 
		vfloat LocA1 = GetLinearStiffness(); //p_Sim->LocalVXC.GetBaseMat(MatIndex)->GetElasticMod()*2*NominalSize; 
		vfloat LocUDynamic = _pMat->GetuDynamic();
		vfloat LocUStatic = _pMat->GetuStatic();

		vfloat NormalForce = LocA1 * CurPenetration; //positive for penetration...
		FloorForce.z += NormalForce; //force resisting penetration
	
		//do vertical damping here...
		FloorForce.z -= pSim->GetCollisionDampZ()*_2xSqMxExS*Vel.z;  //critically damp force for this bond to ground
//		FloorForce.z -= p_Sim->GetCollisionDampZ()*2*Mass*sqrt(LocA1/Mass)*Vel.z;  //critically damp force for this bond to ground
		
		//lateral friction
		vfloat SurfaceVel = sqrt(Vel.x*Vel.x + Vel.y*Vel.y); //velocity along the floor...
		vfloat SurfaceVelAngle = atan2(Vel.y, Vel.x); //angle of sliding along floor...
		vfloat SurfaceForce = sqrt(TotalVoxForce.x*TotalVoxForce.x + TotalVoxForce.y*TotalVoxForce.y);
		vfloat dFrictionForce = LocUDynamic*NormalForce; 
		Vec3D<> FricForceToAdd = -Vec3D<>(cos(SurfaceVelAngle)*dFrictionForce, sin(SurfaceVelAngle)*dFrictionForce, 0); //always acts in direction opposed to velocity in DYNAMIC friction mode
		//alwyas acts in direction opposite to force in STATIC friction mode

		if (Vel.x == 0 && Vel.y == 0){ //STATIC FRICTION: if this point is stopped and in the static friction mode...
			if (SurfaceForce < LocUStatic*NormalForce) StaticFricFlag = true; //if we don't have enough to break static friction
		}
		else { //DYNAMIC FRICTION
			if (dFrictionForce*pSim->dt < Mass*SurfaceVel){ //check momentum. if we are not possibly coming to a stop this timestep, add in the friction force
				FloorForce += FricForceToAdd;
			}
			else { //if we are coming to a stop, don't overshoot the stop. Set to zero and zero the momentum to get static friction to kick in.
				StaticFricFlag = true;
				LinMom.x = 0; //fully stop the voxel here! (caution...)
				LinMom.y = 0;
			}
		}
		
	}
	return FloorForce;

}

Vec3D<> CVXS_Voxel::CalcGndDampEffect() //damps everything to ground as qucik as possible...
{
	vfloat tmp = 0;
//	for (int i=0; i<GetNumLocalBonds(); i++) tmp+=sqrt(GetBond(i)->GetLinearStiffness()*Mass);
	for (int i=0; i<6; i++) tmp+=sqrt(InternalBondPointers[i]->GetLinearStiffness()*Mass);

	return -pSim->GetSlowDampZ()*2*tmp*Vel;
}


vfloat CVXS_Voxel::GetMaxBondStrain(void) const
{
	vfloat MxSt = 0;
	for (int i=0; i<6; i++){
		if (InternalBondPointers[i] != NULL){
			vfloat TSt = InternalBondPointers[i]->GetEngStrain();
			if (TSt>MxSt) MxSt = TSt; 
		}
	}
	return MxSt;
}

vfloat CVXS_Voxel::GetMaxBondStrainE(void) const
{
	vfloat MxSt = 0;
	for (int i=0; i<6; i++){
		if (InternalBondPointers[i] != NULL){
			vfloat TSt = InternalBondPointers[i]->GetStrainEnergy();
			if (TSt>MxSt) MxSt = TSt; 
		}
	}
	return MxSt;
}

vfloat CVXS_Voxel::GetMaxBondStress(void) const
{
	vfloat MxSt = 0;
//	for (int i=0; i<NumLocalBonds; i++){
	for (int i=0; i<6; i++){
		if (InternalBondPointers[i] != NULL){
			vfloat TSt = InternalBondPointers[i]->GetEngStress();
			if (TSt>MxSt) MxSt = TSt; 
		}
	}
	return MxSt;

}

vfloat CVXS_Voxel::CalcVoxMatStress(const vfloat StrainIn, bool* const IsPastYielded, bool* const IsPastFail) const
{
		return _pMat->GetModelStress(StrainIn, IsPastYielded, IsPastFail, Vox_E);
}