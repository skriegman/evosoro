/*******************************************************************************
Copyright (c) 2010, Jonathan Hiller (Cornell University)
If used in publication cite "J. Hiller and H. Lipson "Dynamic Simulation of Soft Heterogeneous Objects" In press. (2011)"

This file is part of Voxelyze.
Voxelyze is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Voxelyze is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
See <http://www.opensource.org/licenses/lgpl-3.0.html> for license details.
*******************************************************************************/

#ifndef CVX_OBJECT_H
#define CVX_OBJECT_H

#include <vector>
#include <iostream>
#include "Utils/Vec3D.h" //use this for portability, instead of Vec3D()
#include "Utils/XML_Rip.h" 

/*Written by: Jonathan Hiller (jdh74) */
class CVXC_Lattice;
class CVXC_Voxel;
class CVXC_Structure;
class CVXC_Material;

//!Structure compression modes. Defines the compression type that will be used to store the raw voxel data
enum CompType {
	CP_RAWDATA, /*!< @deprecated binary data.*/
	CP_ASCIIREADABLE, /*!< Ascii readable data. Index 0 begins at ascii '0', each subsequent index is the next ascii value.*/
	CP_QTZLIB, /*!< @deprecated QT's internal ZLIB format data. Use CP_ZLIB.*/
	CP_BASE64, /*!< Base64 data. Converts raw data to a base 64 ascii-readable character. Extends length slightly.*/
	CP_ZLIB  /*!< ZLib compress data. Zlib must be installed and USE_ZLIB_COMPRESSION must be manually defined pre-processor.*/
};

//The maximum number of layers we'll ever have
#define LAYERMAX 99999 

//!Voxel lattice information class.	
/*!Defines voxel tiling scheme. By default a rectangular lattice is used, but incorporating line and layer offsets allows more complex lattices such as hexagonal close packed (HCP) and Face-centered cubic (FCC).*/
class CVXC_Lattice //container for information about the lattice of possible voxel locations
{
public:
	CVXC_Lattice(void); //!< Constructor
	~CVXC_Lattice(void); //!< Destructor
	CVXC_Lattice(const CVXC_Lattice& Lattice) {*this = Lattice;} //!< Copy constructor
	CVXC_Lattice& operator=(const CVXC_Lattice& RefLat); //!< Overload "="

	//I/O functions for load/save
	void ReadXML(CXML_Rip* pXML); //!< Reads all data from XML tree.
	void WriteXML(CXML_Rip* pXML); //!< Stores all data to XML tree.

	//functions to get info about lattice
	vfloat GetLatticeDim(void) const {return Lattice_Dim;} //!< Returns the base lattice dimension in meters.
	Vec3D<> GetDimAdj(void) const {return Vec3D<>(X_Dim_Adj, Y_Dim_Adj, Z_Dim_Adj);} //!< Returns the scaling factors of the lattice spacing in X, Y, and Z.
	vfloat GetXDimAdj(void) const {return X_Dim_Adj;} //!< Returns the scaling factor if the lattice in the X direction.
	vfloat GetYDimAdj(void) const {return Y_Dim_Adj;} //!< Returns the scaling factor if the lattice in the Y direction.
	vfloat GetZDimAdj(void) const {return Z_Dim_Adj;} //!< Returns the scaling factor if the lattice in the Z direction.
	vfloat GetXLiO(void) const {return X_Line_Offset;} //!< Returns the line X direction offset percentage.
	vfloat GetYLiO(void) const {return Y_Line_Offset;} //!< Returns the line Y direction offset percentage.
	vfloat GetXLaO(void) const {return X_Layer_Offset;} //!< Returns the layer X direction offset percentage.
	vfloat GetYLaO(void) const {return Y_Layer_Offset;} //!< Returns the layer Y direction offset percentage.
	vfloat GetMaxOffsetX(int yV=5, int zV=5) const; //!<Returns the maximum X offset in the lattice. 
	vfloat GetMaxOffsetY(int xV=5, int zV=5) const; //!<Returns the maximum Y offset in the lattice.

	//Funtions to modify the lattice
	void SetLattice(vfloat LatDimIn, vfloat XDAdjIn = 1.0, vfloat YDAdjIn = 1.0, vfloat ZDAdjIn = 1.0, vfloat XLiOIn = 0.0, vfloat YLiOIn = 0.0, vfloat XLaOIn = 0.0, vfloat YLaOIn = 0.0); //!< Sets the lattice to the specified parameters.
	void ClearLattice(void); //!< Returns lattice to default state.

	void SetLatticeDim(vfloat Lattice_Dim_In) {Lattice_Dim = Lattice_Dim_In;} //!< Sets the base lattice dimension in meters. @param[in] Lattice_Dim_In The base lattice dimension in meters.
	void SetXDimAdj(vfloat X_Dim_Adj_In) {X_Dim_Adj = X_Dim_Adj_In;} //!< Sets the scaling factor if the lattice in the X direction. @param[in] X_Dim_Adj_In The lattice dimension X adjustment percentage. Range [0.0, 1.0].
	void SetYDimAdj(vfloat Y_Dim_Adj_In) {Y_Dim_Adj = Y_Dim_Adj_In;} //!< Sets the scaling factor if the lattice in the Y direction. @param[in] Y_Dim_Adj_In The lattice dimension Y adjustment percentage. Range [0.0, 1.0].
	void SetZDimAdj(vfloat Z_Dim_Adj_In) {Z_Dim_Adj = Z_Dim_Adj_In;} //!< Sets the scaling factor if the lattice in the Z direction. @param[in] Z_Dim_Adj_In The lattice dimension Z adjustment percentage. Range [0.0, 1.0].
	void SetXLiO(vfloat X_Line_Offset_In) {X_Line_Offset = X_Line_Offset_In;} //!< Sets the line X direction offset percentage. @param[in] X_Line_Offset_In The lattice line X adjustment percentage. Range [0.0, 1.0].
	void SetYLiO(vfloat Y_Line_Offset_In) {Y_Line_Offset = Y_Line_Offset_In;} //!< Sets the line Y direction offset percentage. @param[in] Y_Line_Offset_In The lattice line Y adjustment percentage. Range [0.0, 1.0].
	void SetXLaO(vfloat X_Layer_Offset_In) {X_Layer_Offset = X_Layer_Offset_In;} //!< Sets the layer X direction offset percentage. @param[in] X_Layer_Offset_In The lattice layer X adjustment percentage. Range [0.0, 1.0].
	void SetYLaO(vfloat Y_Layer_Offset_In) {Y_Layer_Offset = Y_Layer_Offset_In;} //!< Sets the layer Y direction offset percentage. @param[in] Y_Layer_Offset_In The lattice layer Y adjustment percentage. Range [0.0, 1.0].

protected:
	//variables from file:
	vfloat Lattice_Dim; //!< The base lattice dimension. The lattice dimension defines the distance between voxels in meters. This distance can be modified separately in the X, Y, and Z dimensions by X_Dim_Adj, Y_Dim_Adj and Z_Dim_Adj respectively.
	vfloat X_Dim_Adj; //!< The percentage to adjust the lattice dimension in the X direction. [0.0 to 1.0]
	vfloat Y_Dim_Adj; //!< The percentage to adjust the lattice dimension in the Y direction. [0.0 to 1.0]
	vfloat Z_Dim_Adj; //!< The percentage to adjust the lattice dimension in the Z direction. [0.0 to 1.0]
	vfloat X_Line_Offset; //!< The amount to offset each subsequent line in X (as Y increases) as a percentage of Lattice_Dim. [0.0 to 1.0]
	vfloat Y_Line_Offset; //!< The amount to offset each subsequent line in Y (as Y increases) as a percentage of Lattice_Dim. [0.0 to 1.0]
	vfloat X_Layer_Offset; //!< The amount to offset each subsequent layer in X (as Z increases) as a percentage of Lattice_Dim. [0.0 to 1.0]
	vfloat Y_Layer_Offset; //!< The amount to offset each subsequent layer in Y (as Z increases) as a percentage of Lattice_Dim. [0.0 to 1.0]
};

enum VoxShapes {
	VS_SPHERE, 
	VS_BOX, 
	VS_CYLINDER, 
	VS_VOXFILE
};

class CVXC_Voxel
{
public:
	CVXC_Voxel(void) {ClearVoxel();};
	~CVXC_Voxel(void) {};
	CVXC_Voxel(const CVXC_Voxel& Voxel) {*this = Voxel;} //copy constructor
	CVXC_Voxel& operator=(const CVXC_Voxel& RefVoxel); //overload "=" 

	//I/O function for save/loading
	void ReadXML(CXML_Rip* pXML);
	void WriteXML(CXML_Rip* pXML);

	//functions to modify voxel
	void SetVoxName(int VoxNameIn){Vox_Name = VoxNameIn;}; //sets the voxel type
	void SetVoxName(int VoxNameIn, vfloat XSQ, vfloat YSQ, vfloat ZSQ){Vox_Name = VoxNameIn; X_Squeeze=XSQ; Y_Squeeze=YSQ; Z_Squeeze=ZSQ;};
	void SetXSqueeze(vfloat XSqueezeIn){X_Squeeze = XSqueezeIn;};
	void SetYSqueeze(vfloat YSqueezeIn){Y_Squeeze = YSqueezeIn;};
	void SetZSqueeze(vfloat ZSqueezeIn){Z_Squeeze = ZSqueezeIn;};
	void ClearVoxel(void);

	//functions to get information about voxel
	inline int GetVoxName() const {return Vox_Name;};
	inline vfloat GetXSqueeze(void) const {return X_Squeeze;};
	inline vfloat GetYSqueeze(void) const {return Y_Squeeze;};
	inline vfloat GetZSqueeze(void) const {return Z_Squeeze;};

#ifdef USE_OPEN_GL
	void DrawVoxel(Vec3D<>* Center, vfloat Lat_Dim); //draw it (OpenGL)
	void DrawVoxel2D(Vec3D<>* Center, vfloat Lat_Dim, Vec3D<>* Normal, bool Fill = false); //draw 2D view as seen from viewpoint pointed to by normal
#endif

protected:
	//variable from file:
	int Vox_Name;
	std::string File;
	vfloat X_Squeeze;
	vfloat Y_Squeeze;
	vfloat Z_Squeeze;
};

static const std::string base64_chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
static const char defaultReturn = -1;

class CVXC_Structure //contains voxel location information in vast 1D array
{
public:
	CVXC_Structure() {pData = NULL; ClearStructure();};
	CVXC_Structure(int xV, int yV, int zV) {pData = NULL; CreateStructure(xV, yV, zV);};
	~CVXC_Structure() {DeleteData();};
	CVXC_Structure(const CVXC_Structure& RefStruct) {*this = RefStruct;}; //copy constructor
	CVXC_Structure& operator=(const CVXC_Structure& RefStruct); //overload "=" 
	char& operator [](int i) const {return GetData(i);}

	//I/O function for save/loading
	bool WriteXML(CXML_Rip* pXML, int Compression = CP_ASCIIREADABLE, std::string* RetMessage = NULL);
	bool ReadXML(CXML_Rip* pXML, std::string Version = "", std::string* RetMessage = NULL); //version -1 means latest.
	static inline bool is_base64(unsigned char c) {return (isalnum(c) || (c == '+') || (c == '/'));}
	std::string ToBase64(unsigned char const* , unsigned int len);
	std::string FromBase64(std::string const& s);

	//Get information about the structure:
	inline char& GetData(int Index) const {if (DataInit) return pData[Index]; else return (char&)defaultReturn;} //Gets the material index here (this should be the only place we access pData)
	inline int GetArraySize(void) const {return m_SizeOfArray;}; //gets the size of the master array of voxels
	inline int GetVXDim(void) const {return X_Voxels;}; //get number of voxels in each dimension
	inline int GetVYDim(void) const {return Y_Voxels;};
	inline int GetVZDim(void) const {return Z_Voxels;};
	int GetIndex(int x, int y, int z) const; //returns the index of the array from xyz indices
	bool GetXYZNom(int* x, int* y, int* z, int index) const; //gets the physical position of the voxels from index
	bool ContainsMatIndex(int MatIndex) const {for (int i=0; i<GetArraySize(); i++)if (GetData(i) == MatIndex) return true; return false;}; //returns true if the structure contains this material index

	//Functions to modify the lattice:
	inline void SetData(int Index, char Data){if (DataInit) pData[Index] = Data;} //sets the material index here
	void CreateStructure(int xV, int yV, int zV); //creates empty structure with these dimensions
	void Resize(int xS, int yS, int zS); //resizes the structure preserving data
	void ClearStructure(); //completly erases, frees, and destroys the voxel array
	void ResetStructure(); //erases all voxel imformation within voxel array
	void ReplaceMaterial(int Matindex, int ReplaceWith = 0, bool ShiftDown = false); //converts a material to ReplaceWith, option to decrement all higher indices

	// nac: neuralNet:
	// inline void SetWeight(int Index1, int Index2, double Weight){pSynapseWeigths[Index1][Index2] = Weight;} // nac: sets the material index here
	// inline double GetSynapseWeightArray(int Index1, int Index2) const {return pSynapseWeigths[Index1][Index2];}
	// inline void SetSynapseWeight(int Index, double Weight) {pSynapseWeigths[Index] = Weight;} // nac: sets the material index here
	// inline double GetSynapseWeight(int Index) const {return pSynapseWeigths[Index];}
	// inline void InitSynpaseArray(int Size) {pSynapseWeigths = new double[Size];}
	// int getNumNuerons() {return NumNuerons;}

	// nac: phase offsets
	inline void SetPhaseOffset(int Index, double Value) {pPhaseOffsets[Index] = Value;} // nac: sets the material index here
	inline double GetPhaseOffset(int Index) const {return pPhaseOffsets[Index];}
	inline void InitPhaseOffsetArray(int Size) {pPhaseOffsets = new double[Size];}
	inline bool GetUsingPhaseOffset(void) {return usingPhaseOffset;}

	inline void SetFinalPhaseOffset(int Index, double Value) {pFinalPhaseOffsets[Index] = Value;} // nac: sets the material index here
	inline double GetFinalPhaseOffset(int Index) const {return pFinalPhaseOffsets[Index];}
	inline void InitFinalPhaseOffsetArray(int Size) {pFinalPhaseOffsets = new double[Size];}
	inline bool GetUsingFinalPhaseOffset(void) {return usingFinalPhaseOffset;}

	// damping phase functions
	inline void SetTempAmpDamp(int Index, double Value) {pTempAmpDamp[Index] = Value;}
	inline double GetTempAmpDamp(int Index) const {return pTempAmpDamp[Index];}
	inline void InitTempAmpDampArray(int Size) {pTempAmpDamp = new double[Size];}
	inline bool GetUsingTempAmpDamp(void) {return usingTempAmpDamp;}

	inline void SetFinalTempAmpDamp(int Index, double Value) {pFinalTempAmpDamp[Index] = Value;}
	inline double GetFinalTempAmpDamp(int Index) const {return pFinalTempAmpDamp[Index];}
	inline void InitFinalTempAmpDampArray(int Size) {pFinalTempAmpDamp = new double[Size];}
	inline bool GetUsingFinalTempAmpDamp(void) {return usingFinalTempAmpDamp;}

	// EvolvedVoxelSize functions
	inline void SetInitialVoxelSize(int Index, double Value) {pInitialVoxelSize[Index] = Value;}
	inline double GetInitialVoxelSize(int Index) const {return pInitialVoxelSize[Index];}
	inline void InitInitialVoxelSizeArray(int Size) {pInitialVoxelSize = new double[Size];}
	inline bool GetUsingInitialVoxelSize(void) {return usingInitialVoxelSize;}

	// FinalVoxelSize functions
	inline void SetFinalVoxelSize(int Index, double Value) {pFinalVoxelSize[Index] = Value;}
	inline double GetFinalVoxelSize(int Index) const {return pFinalVoxelSize[Index];}
	inline void InitFinalVoxelSizeArray(int Size) {pFinalVoxelSize = new double[Size];}
	inline bool GetUsingFinalVoxelSize(void) {return usingFinalVoxelSize;}


     // GrowthTime functions
    inline void SetGrowthTime(int Index, double Value) {pGrowthTime[Index] = Value;}
    inline double GetGrowthTime(int Index) const {return pGrowthTime[Index];}
    inline void InitGrowthTimeArray(int Size) {pGrowthTime = new double[Size];}
    inline bool GetUsingGrowthTime(void) {return usingGrowthTime;}

    // StartGrowthTime functions
    inline void SetStartGrowthTime(int Index, double Value) {pStartGrowthTime[Index] = Value;}
    inline double GetStartGrowthTime(int Index) const {return pStartGrowthTime[Index];}
    inline void InitStartGrowthTimeArray(int Size) {pStartGrowthTime = new double[Size];}
    inline bool GetUsingStartGrowthTime(void) {return usingStartGrowthTime;}


	// stiffness functions
	inline void SetStiffness(int Index, double Value) {pStiffness[Index] = Value;} 
	inline double GetStiffness(int Index) const {return pStiffness[Index];}
	inline void InitStiffnessArray(int Size) {pStiffness = new double[Size];}
	inline bool GetEvolvingStiffness(void) {return evolvingStiffness;}

	// stiffness plasticity functions
	inline void SetStiffnessPlasticityRate(int Index, double Value) {pStiffnessPlasticityRate[Index] = Value;} 
	inline double GetStiffnessPlasticityRate(int Index) const {return pStiffnessPlasticityRate[Index];}
	inline void InitStiffnessPlasticityRateArray(int Size) {pStiffnessPlasticityRate = new double[Size];}
	inline bool GetUsingStiffnessPlasticity(void) {return stiffnessPlasticity;}

	// evoDevo parameters
	inline void SetKP(int Index, double Value) {pKP[Index] = Value;} 
	inline double GetKP(int Index) const {return pKP[Index];}
	inline void InitKPArray(int Size) {pKP = new double[Size];}
	inline void SetKI(int Index, double Value) {pKI[Index] = Value;} 
	inline double GetKI(int Index) const {return pKI[Index];}
	inline void InitKIArray(int Size) {pKI = new double[Size];}
	inline void SetANTIWINDUP(int Index, double Value) {pANTIWINDUP[Index] = Value;} 
	inline double GetANTIWINDUP(int Index) const {return pANTIWINDUP[Index];}
	inline void InitANTIWINDUPArray(int Size) {pANTIWINDUP = new double[Size];}	
	inline bool GetEvoDevoParameters(void) {return evoDevoParameters;}


	// int getNumNuerons() {return NumNuerons;}

protected:
	//variable from file
	std::string Compression;
	char* pData; //the main voxel array. This is an array of chars; the entries correspond with the material IDs, and the position in the array corresponds with the position in the 3D structure, where the array is ordered: starting at (x0,x0,z0), proceeding to (xn,y0,z0), next to (xn,yn,z0), and on to (xn,yn,zn)

	// std::vector< std::vector<double> > pSynapseWeigths;
	// double* pSynapseWeigths;
	// int NumNuerons;

	bool usingPhaseOffset; 
	double* pPhaseOffsets;

    bool usingFinalPhaseOffset;
	double* pFinalPhaseOffsets;

	bool usingTempAmpDamp;
    double* pTempAmpDamp;

    bool usingFinalTempAmpDamp;
    double* pFinalTempAmpDamp;

	bool usingInitialVoxelSize;
	double* pInitialVoxelSize;

	bool usingFinalVoxelSize;
	double* pFinalVoxelSize;

	bool usingGrowthTime;
    double* pGrowthTime;

    bool usingStartGrowthTime;
    double* pStartGrowthTime;

	bool evolvingStiffness;
	double* pStiffness;

	bool stiffnessPlasticity;
	double* pStiffnessPlasticityRate;


	bool evoDevoParameters;
	double* pKP;
	double* pKI;
	double* pANTIWINDUP;


	int X_Voxels;
	int Y_Voxels;
	int Z_Voxels;

	//other variables:
	int m_SizeOfArray; //keep track of the number of elements in the voxel array

	//non-stored variables
	bool DataInit; //flag for if the dynamically allocated structure 

	//functions that should only be used locally
	void DeleteData(void); //sandbox the creation and destruction...
	void IniData(int Size);
};

//!Voxel object class.	
/*!Contains all information to describe a multi-material voxel object.*/
class CVX_Object
{
public:
	//Constructors and operators
	CVX_Object(void); //!< Constructor
	~CVX_Object(void); //!< Destructor
	CVX_Object(const CVX_Object& RefObj) {*this = RefObj;} //!< Copy constructor
	CVX_Object& operator=(const CVX_Object& RefObj); //!< Overload "=" 

	//Member classes (see individual descriptiona )
	CVXC_Lattice Lattice; //!< The one and only CVXC_Lattice object for this voxel object. Describes the lattice to be used.
	CVXC_Voxel Voxel; //!< The one and only CVXC_Voxel object for this voxel object. Describes the way voxels of this object should be handled in a visualization or simulation.
	CVXC_Structure Structure; //!< The one and only CVXC_Structure object for this voxel object. Contains the voxel bitmap describing where voxels are present and their material indices for this object.
	std::vector<CVXC_Material> Palette; //!< The palette of materials available for use in this object.

	//I/O function for save/loading
	void SaveVXCFile(std::string filename, int Compression = CP_ASCIIREADABLE); //!< Saves all voxel object information from this voxel object into an XML-based format.
	bool LoadVXCFile(std::string filename); //!< Attempts to load all voxel object information from the specified file.
	void SaveVXPFile(std::string filename); //!< Saves all material palette information from this voxel object into an XML-based format.
	bool LoadVXPFile(std::string filename); //!< Attempts to load all material palette information from the specified file.
	bool ExportKV6File(std::string filename); //!< Exports kV6 file compatible with slab6

	void WriteXML(CXML_Rip* pXML, int Compression = CP_ASCIIREADABLE, std::string* RetMessage = NULL); //!< Writes all voxel object information to an XML ripping stream.
	void WritePaletteXML(CXML_Rip* pXML, int Compression = CP_ASCIIREADABLE); //!< Writes just material palette information to an XML ripping stream.
	void ReadXML(CXML_Rip* pXML, bool OnlyPal = false, std::string* RetMessage = NULL); //!< Attempts to extract voxel object from XML tree.
	void ReadPaletteXML(CXML_Rip* pXML, std::string Version = ""); //!< Attempts to extract material palette from XML tree.

	//Creating and editing the whole VXC object
	void InitializeMatter(void); //!< Initializes voxel object with default values.
	void InitializeMatter(vfloat iLattice_Dim, int xV, int yV, int zV); //!< Initializes voxel object with the specifies sizes and default lattice.
	void InitializeMatter(CVXC_Lattice* pLattice, int xV, int yV, int zV); //!< initializes voxel object with the specified lattice and size.
	void ClearMatter(void); //!< Clears the entire voxel object. Must be re-initialized before subsequent use.
	void Transform(Vec3D<> Trans); //!< Moves all voxels within the workspace the specified displacement.
	void Resize(CVXC_Structure* Structure) {Resize(Structure->GetVXDim(), Structure->GetVYDim(), Structure->GetVZDim());} //!< Resizes the voxel object structure to size of the provided structure. Voxel data remaining within the new structure area is preserved. @param[in] Structure CVXC_Structure object to extract the desired size from.
	void Resize(int xS, int yS, int zS) {Structure.Resize(xS, yS, zS);} //!< Resizes the voxel object structure to specified number of voxels in each dimension. Voxel data remaining within the new structure area is preserved. @param[in] xS Number of desired voxels in X direction. @param[in] yS Number of desired voxels in Y direction. @param[in] zS Number of desired voxels in Z direction.

	//Basic editing of VXC object
	bool SetMat(int x, int y, int z, int MatIndex) {return SetMat(GetIndex(x, y, z), MatIndex);} //!< Sets a single voxel to the specified material. Returns true if succesful. Returns false if indices are outside of workspace or material index is not contained within current palette. @param[in] x Integer X index of voxel location to set. @param[in] y Integer Y index of voxel location to set. @param[in] z Integer Z index of voxel location to set.  @param[in] MatIndex Specifies the index within the material palette to set this voxel to.
	bool SetMat(int StructIndex, int MatIndex); //!< Sets a single voxel to the specified material.
	bool SetMatFill(int MatIndex); //!< Fills the entire workspace to the specifies material.
	inline int GetMat(int x, int y, int z) const {return GetMat(GetIndex(x, y, z));} //!< Returns the base material index at this location.  @param[in] x Integer X index of voxel location to get. @param[in] y Integer Y index of voxel location to get. @param[in] z Integer Z index of voxel location to get. 
	int GetMat(int StructIndex) const; //!< Returns the base material index at this location.
	inline CVXC_Material* GetBaseMatAt(int StructIndex) {return &Palette[GetMat(StructIndex)];} //!< Returns a temporary pointer to material at the specified structure location. @param[in] StructIndex The global voxel structure index to query.
	inline CVXC_Material* GetBaseMat(int MatIndex) {return &Palette[MatIndex];} //!< Returns a temporary pointer to the specified palette material index. @param[in] MatIndex The desired material palette index.
	inline CVXC_Material* GetLeafMat(int StructIndex, bool* pVisible = NULL) {return &Palette[GetLeafMatIndex(StructIndex, pVisible)];} //!< Returns a temporary pointer to the leaf material at this location after evaluating all sub materials. @param[in] StructIndex The global voxel structure index to query. @param[out] pVisible Set to false if this material is currently to be hidden for visualization. Otherwise true.

	//dep
	int GetSubMatIndex(int StructIndex, bool* pVisible = NULL) {return GetSubMatIndex(StructIndex, Structure.GetData(StructIndex), pVisible);} //!< Returns the material index of the next sub-material at this location. Entry point into sub-material recursion process. @param[in] StructIndex The global voxel structure index to query. @param[out] pVisible Set to false if this material is currently to be hidden for visualization. Otherwise true.
	int GetSubMatIndex(int StructIndex, int MatIndex, bool* pVisible = NULL); //!< Returns the material index of the next sub-material at this location.
	//end dep
	
	int GetSubMatIndex(int* pXIndex, int* pYIndex, int* pZIndex, int MatIndex, bool* pVisible = NULL);
	int GetLeafMatIndex(int StructIndex, bool* pVisible = NULL); //!< Returns material index of the leaf material at this location after evaluating all sub materials.


	//Get information about the workspace:
	inline vfloat GetLatticeDim(void) const {return Lattice.GetLatticeDim();} //!< Returns the current lattice dimension (base voxel size) in meters.
	inline Vec3D<> GetLatDimEnv(void) const {return GetLatticeDim()*Lattice.GetDimAdj();} //!< Returns the voxel lattice spacing in meters. Includes the X, Y, and Z dimension adjustments.
	inline Vec3D<> GetLatDimEnvPerc(void) const {return GetLatDimEnv().ScaleInv(GetWorkSpace());} //!< Returns the voxel lattice spacing as a percentage of the entire workspace dimension in each respective dimension. Includes the X, Y, and Z dimension adjustments. All components have a range [0.0 to 1.0]
	inline void GetVDim(int* X, int* Y, int* Z) const {*X = GetVXDim(); *Y = GetVYDim(); *Z = GetVZDim();} //!< Provides the number of voxels of the workspace in all directions. @param[out] X Number of voxels in the X direction. @param[out] Y Number of voxels in the Y direction. @param[out] Z Number of voxels in the Z direction.
	inline int GetVXDim(void) const {return Structure.GetVXDim();} //!< Returns the number of voxels in the X direction.
	inline int GetVYDim(void) const {return Structure.GetVYDim();} //!< Returns the number of voxels in the Y direction.
	inline int GetVZDim(void) const {return Structure.GetVZDim();} //!< Returns the number of voxels in the Z direction.
	inline int GetStArraySize(void) const {return Structure.GetArraySize();} //!< Returns the number of elements in the global voxel structure array.
	void GetWorkSpace(Vec3D<>* Dim) const; //!< Provides the calculated workspace dimensions in meters.
	Vec3D<> GetWorkSpace() const; //!< Returns the calculated workspace dimensions in meters.
	
	//Get information about a voxel (location) within the workspace
	bool GetXYZ(Vec3D<>* Point, int index, bool WithOff = true) const; //!< Calculates the XYZ coordinates of a voxel.
	inline bool GetXYZ(Vec3D<>* Point, int x, int y, int z, bool WithOff = true) const {return GetXYZ(Point, GetIndex(x, y, z), WithOff);} //!< Calculates the XYZ coordinates of a voxel. Returns true if a valid location was specified. Otherwise false. The resulting position is given in meters from the origin. @param[out] Point The X, Y, and Z components of the voxel location is stored in this structure. @param[in] x Integer X index of voxel location to get. @param[in] y Integer Y index of voxel location to get. @param[in] z Integer Z index of voxel location to get. 
	inline Vec3D<> GetXYZ(int index, bool WithOff = true) const {Vec3D<> toReturn; if (GetXYZ(&toReturn, index, WithOff)) return toReturn; else return Vec3D<>(-1, -1, -1);} //!< Returns the XYZ coordinates of a voxel. The resulting position if the voxel at the specified global index is given in meters from the origin. @param[in] index The global structural index to query.
	inline Vec3D<> GetXYZ(int x, int y, int z, bool WithOff = true) const {return GetXYZ(GetIndex(x, y, z), WithOff);} //!< Returns the XYZ coordinates of a voxel. The resulting position if the voxel at the specified X, Y, and Z indices is given in meters from the origin. @param[in] x Integer X index of voxel location to get. @param[in] y Integer Y index of voxel location to get. @param[in] z Integer Z index of voxel location to get. 
	inline Vec3D<> GetXYZPerc(int index) const {return GetXYZ(index).ScaleInv(GetWorkSpace());} //!< Returns the location of a voxel normalized by the size of workspace. All components are in the range [0.0 to 1.0]. @param[in] index The global structural index to query.
	inline int GetIndex(int x, int y, int z) const {return Structure.GetIndex(x, y, z);} //!< Returns the global structure index from XYZ indices. @param[in] x Integer X index of voxel location to get. @param[in] y Integer Y index of voxel location to get. @param[in] z Integer Z index of voxel location to get. 
	inline int GetRelativeIndex(int BaseIndex, int rX, int rY, int rZ) {int bX, bY, bZ; GetXYZNom(&bX, &bY, &bZ, BaseIndex); return GetIndex(bX+rX, bY+rY, bZ+rZ);} //!< Returns the global index of the voxel at the specified relative offset from another voxel. Returns -1 if the resulting location is outside of the workspace. The index is returned regardless of whether a voxel is present there or not. @param[in] BaseIndex The global structural index to add offsets to. @param[in] rX Integer X offset from the base voxel to query. @param[in] rY Integer Y offset from the base voxel to query. @param[in] rZ Integer Z offset from the base voxel to query.
	inline int GetRelativeMatIndex(int BaseIndex, int rX, int rY, int rZ) {return GetMat(GetRelativeIndex(BaseIndex, rX, rY, rZ));} //!< Returns the material index of the voxel at the specified relative offset from another voxel. Returns -1 if invalid material/location. @param[in] BaseIndex The global structural index to add offsets to. @param[in] rX Integer X offset from the base voxel to query. @param[in] rY Integer Y offset from the base voxel to query. @param[in] rZ Integer Z offset from the base voxel to query.
	inline bool IsVoxRelativeIndex(int BaseIndex, int rX, int rY, int rZ) {int tmpMat = GetMat(GetRelativeIndex(BaseIndex, rX, rY, rZ)); if (tmpMat == 0 || tmpMat == -1) return false; else return true;} //!< Returns true if a voxel is present at the specified relative offset from another voxel. Returns false if empty location, invalid material, or out-of-bounds location. @param[in] BaseIndex The global structural index to add offsets to. @param[in] rX Integer X offset from the base voxel to query. @param[in] rY Integer Y offset from the base voxel to query. @param[in] rZ Integer Z offset from the base voxel to query.
	inline bool GetXYZNom(int* x, int* y, int* z, int index) const {return Structure.GetXYZNom(x, y, z, index);} //!< Provides the XYZ indices of the specified global structure index. @param[out] x X index of specified voxel. @param[out] y Y index of specified voxel. @param[out] z Z index of specified voxel. @param[in] index The global structural index to query.
	bool IsAdjacent(int Index1, int Index2, bool IncludeEmpty = false, Vec3D<>* RelativeLoc = NULL); //!< Returns true if the two specified voxels are touching, Otherwise false.

	//Creating and editting Palette
	void LoadDefaultPalette(void); //!< Loads a default palette.
	void ClearPalette(); //!< Erases the material palette, leaving only the erase (void) material.
	int AddMat(std::string Name, std::string* RetMessage = NULL) {return AddMat(Name, 1000000.0, 0.35, RetMessage);} //!< Appends a material to the palette. Returns the index within the palette that this material was created at. This index may change if materials are deleted from the palette. All colors and physical properties are set to defaults. @param[in] Name The name of the material to add. If the name is already in use a variant will be generated. @param[out] RetMessage Pointer to an initialized string. Messages generated in this function will be appended to the string.
	int AddMat(CVXC_Material& MatToAdd, bool ForceBasic = false, std::string* RetMessage = NULL); //!< Appends a material to the palette.
	int AddMat(std::string Name, vfloat EMod, vfloat PRatio, std::string* RetMessage = NULL); //!< Appends a material to the palette.
	bool MatNameExists(std::string Name);
	bool ReplaceMaterial(int IndexToReplace, int NewIndex, bool DeleteReplaced = false, std::string* RetMessage = NULL);
	bool FlattenMaterial(int IndexToFlatten, std::string* RetMessage = NULL); //converts a compound material into is subsidiary materials
	bool DeleteMat(int Index, bool FlattenFirst = false, std::string* RetMessage = NULL); //deletes a material and shifts the rest down
	int GetNumMaterials(void) {return (int)Palette.size();};
	int GetNumLeafMatInUse(void);
	bool IsInRecursivePath(int MatIndexStart, int MatIndexRef); //returns true if MatIndexRef is part of the recusive path of MatIndexStart

	//Statistics of structure
	int GetNumVox(void);
	int GetNumVox(int MatIndex, bool OnlyLeafMat = false);
	vfloat GetSurfaceArea(); //in m2
	vfloat GetVolume(); //in m3
	vfloat GetWeight(void); //based on densities

	//strings for outputting info
	void GetVXCInfoStr(std::string* pString = NULL); //fills pString with overall info about the VXC
	void GetVoxInfoStr(int VoxIndex, std::string* pString = NULL); //fills pString with info about a specific voxel

	//DEPRECATED function names to access voxel data
	bool SetVoxel(int x, int y, int z, int MatIndex) {return SetMat(x, y, z, MatIndex);} //!< @deprecated Sets a voxel to a material.
	bool SetVoxel(int Index, int MatIndex) {return SetMat(Index, MatIndex);} //!< @deprecated Sets a voxel to a material.
	int GetVoxel(int x, int y, int z) {return GetMat(x, y, z);}//!< @deprecated Returns material index at location.
	int GetVoxel(int Index) {return GetMat(Index);} //!< @deprecated Returns material index at location.

	// nac: neural net:
	// inline void SetSynapseWeight(int Index, double Weight) {Structure.SetSynapseWeight(Index,Weight);} // nac: sets the material index here
	// inline double GetSynapseWeight(int Index) const {return Structure.GetSynapseWeight(Index);}
	// inline void InitSynpaseArray(int Size) {Structure.InitSynpaseArray(Size);}
	// int getNumNuerons() {return Structure.getNumNuerons();}
	inline void SetPhaseOffset(int Index, double Value) {Structure.SetPhaseOffset(Index,Value);}
	inline double GetPhaseOffset(int Index) const {return Structure.GetPhaseOffset(Index);}
	inline void InitPhaseOffsetArray(int Size) {Structure.InitPhaseOffsetArray(Size);}
	inline bool GetUsingPhaseOffset(void) {return Structure.GetUsingPhaseOffset();}

	inline void SetFinalPhaseOffset(int Index, double Value) {Structure.SetFinalPhaseOffset(Index,Value);}
	inline double GetFinalPhaseOffset(int Index) const {return Structure.GetFinalPhaseOffset(Index);}
	inline void InitFinalPhaseOffsetArray(int Size) {Structure.InitFinalPhaseOffsetArray(Size);}
	inline bool GetUsingFinalPhaseOffset(void) {return Structure.GetUsingFinalPhaseOffset();}

	// damping actuation functions
    inline void SetTempAmpDamp(int Index, double Value) { Structure.SetTempAmpDamp(Index, Value);}
	inline double GetTempAmpDamp(int Index) const {return Structure.GetTempAmpDamp(Index);}
	inline void InitTempAmpDampArray(int Size) { Structure.InitTempAmpDampArray(Size); }
	inline bool GetUsingTempAmpDamp(void) {return Structure.GetUsingTempAmpDamp(); }

	inline void SetFinalTempAmpDamp(int Index, double Value) { Structure.SetFinalTempAmpDamp(Index, Value);}
	inline double GetFinalTempAmpDamp(int Index) const {return Structure.GetFinalTempAmpDamp(Index);}
	inline void InitFinalTempAmpDampArray(int Size) { Structure.InitFinalTempAmpDampArray(Size); }
	inline bool GetUsingFinalTempAmpDamp(void) {return Structure.GetUsingFinalTempAmpDamp(); }

	// EvolvedVoxelSize functions
	inline void SetInitialVoxelSize(int Index, double Value) { Structure.SetInitialVoxelSize(Index, Value);}
	inline double GetInitialVoxelSize(int Index) const {return Structure.GetInitialVoxelSize(Index);}
	inline void InitInitialVoxelSizeArray(int Size) { Structure.InitInitialVoxelSizeArray(Size); }
	inline bool GetUsingInitialVoxelSize(void) {return Structure.GetUsingInitialVoxelSize(); }

	// FinalVoxelSize functions
	inline void SetFinalVoxelSize(int Index, double Value) { Structure.SetFinalVoxelSize(Index, Value);}
	inline double GetFinalVoxelSize(int Index) const {return Structure.GetFinalVoxelSize(Index);}
	inline void InitFinalVoxelSizeArray(int Size) { Structure.InitFinalVoxelSizeArray(Size); }
	inline bool GetUsingFinalVoxelSize(void) {return Structure.GetUsingFinalVoxelSize(); }

	// GrowthTime functions
    inline void SetGrowthTime(int Index, double Value) { Structure.SetGrowthTime(Index, Value);}
    inline double GetGrowthTime(int Index) const {return Structure.GetGrowthTime(Index);}
    inline void InitGrowthTimeArray(int Size) { Structure.InitGrowthTimeArray(Size); }
    inline bool GetUsingGrowthTime(void) {return Structure.GetUsingGrowthTime(); }

    // StartGrowthTime functions
    inline void SetStartGrowthTime(int Index, double Value) { Structure.SetStartGrowthTime(Index, Value);}
    inline double GetStartGrowthTime(int Index) const {return Structure.GetStartGrowthTime(Index);}
    inline void InitStartGrowthTimeArray(int Size) { Structure.InitStartGrowthTimeArray(Size); }
    inline bool GetUsingStartGrowthTime(void) {return Structure.GetUsingStartGrowthTime(); }

	// stiffness functions
	inline void SetStiffness(int Index, double Value) {Structure.SetStiffness(Index,Value);} 
	inline double GetStiffness(int Index) const {return Structure.GetStiffness(Index);}
	inline void InitStiffnessArray(int Size) {Structure.InitStiffnessArray(Size);}
	inline bool GetEvolvingStiffness(void) {return Structure.GetEvolvingStiffness();}

	// stiffness plasticity functions
	inline void SetStiffnessPlasticityRate(int Index, double Value) {Structure.SetStiffnessPlasticityRate(Index,Value);}
	inline double GetStiffnessPlasticityRate(int Index) const {return Structure.GetStiffnessPlasticityRate(Index);}
	inline void InitStiffnessPlasticityRateArray(int Size) {Structure.InitStiffnessPlasticityRateArray(Size);}
	inline bool GetUsingStiffnessPlasticity(void) {return Structure.GetUsingStiffnessPlasticity();}

	// evoDevoParameters function
	inline void SetKP(int Index, double Value) { Structure.SetKP(Index,Value);  } 
	inline double GetKP(int Index) const {return Structure.GetKP(Index); }
	inline void InitKPArray(int Size) { Structure.InitKPArray(Size); }
	inline void SetKI(int Index, double Value) { Structure.SetKI(Index,Value);  } 
	inline double GetKI(int Index) const {return Structure.GetKI(Index); }
	inline void InitKIArray(int Size) { Structure.InitKIArray(Size); }
	inline void SetANTIWINDUP(int Index, double Value) { Structure.SetANTIWINDUP(Index,Value);  } 
	inline double GetANTIWINDUP(int Index) const {return Structure.GetANTIWINDUP(Index); }
	inline void InitANTIWINDUPArray(int Size) { Structure.InitANTIWINDUPArray(Size); }
	inline bool GetEvoDevoParameters(void) {return Structure.GetEvoDevoParameters();}

	// int getNumNuerons() {return Structure.getNumNuerons();}

	//OpenGL drawing code
#ifdef USE_OPEN_GL
	//call from within OpenGL drawing class to draw matter
	void Draw(int Selected = -1, int XMin=-1, int XMax=LAYERMAX, int YMin=-1, int YMax=LAYERMAX, int ZMin=-1, int ZMax=LAYERMAX); //draws voxels
	void DrawVoxel(Vec3D<>* Center, vfloat Lat_Dim); //Low level (assumes correct GL color, etc have been called
	void DrawSingleVoxel(int StructIndex); //draws voxel with the correct material color, etc.
	void DrawVoxel2D(Vec3D<>* Center, vfloat Lat_Dim, Vec3D<>* Normal, bool Fill = false) {Voxel.DrawVoxel2D(Center, Lat_Dim, Normal, Fill);};
	

#endif
};

enum MatMode {SINGLE, INTERNAL, EXTERNAL, DITHER}; //don't change these orders!
enum MatModelMode {MDL_LINEAR, MDL_LINEAR_FAIL, MDL_BILINEAR, MDL_DATA};
enum FailMode {FM_MAXSTRESS, FM_MAXSTRAIN};
enum RotateAxis {RAX_X, RAX_Y, RAX_Z};
enum RotateAmount {RAM_0, RAM_90, RAM_180, RAM_270};

class CVXC_Material
{
public:
	CVXC_Material(void);
	CVXC_Material(std::string Name, float r, float g, float b, float a = 1.0, vfloat EMod = 0.0, vfloat Poiss = 0.0); //allow us to create easily
	CVXC_Material(const CVXC_Material& RefMat) {*this = RefMat;}; //copy constructor
	~CVXC_Material(void) {DeleteIntStruct();};
	CVXC_Material& operator=(const CVXC_Material& RefMat); //overload "=" 

	//I/O function for save/loading
	void WriteXML(CXML_Rip* pXML, int Compression = CP_ASCIIREADABLE);
	void ReadXML(CXML_Rip* pXML, std::string Version = "", std::string* RetMessage = NULL);

	//Functions to modify material
	void ClearMaterial(void); //reset everything to defaults (grey single material)
	void SetName(std::string NameIn){Name = NameIn;}
	void SetMatType(int NewMatType);
	void SetMatModel(int NewMatModel);
	void SetFailModel(int NewFailModel);
	void SetVisible(bool VisIn) {Visible=VisIn;}
	//Color:
	void SetColor(float r, float g, float b, float a) {Red = r; Green = g; Blue = b; Alpha = a;}
	void SetColor(float r, float g, float b) {SetColor(r, g, b, 1.0f);}
	void SetRed(float r){Red=r;}
	void SetGreen(float g){Green=g;}
	void SetBlue(float b){Blue=b;}
	void SetAlpha(float a){Alpha=a;}
	//Physical properties:
	void SetElasticMod(vfloat EMIn){Elastic_Mod = EMIn;}
	void SetPlasticMod(vfloat PMIn){Plastic_Mod = PMIn;}
	void SetYieldStress(vfloat YSIn){Yield_Stress = YSIn;}
	void SetFailStress(vfloat FSIn){Fail_Stress = FSIn;}
	void SetFailStrain(vfloat FSIn){Fail_Strain = FSIn;}
	bool SetSSData(std::vector<vfloat>* pStrain, std::vector<vfloat>* pStress, std::string* RetMsg); //imports point data to use for tension sterss/strain curve...
	bool ValidateSSData(std::string* RetMsg = NULL); //Does all the checks in the imported stress/strain data to make sure its suitable for simulating

	void SetPoissonsRatio(vfloat PIn){Poissons_Ratio = PIn;}
	void SetDensity(vfloat DensIn){Density = DensIn;}
	void SetCTE(vfloat CTEIn){CTE = CTEIn;}
	void SetuStatic(vfloat uStaticIn){uStatic = uStaticIn;}
	void SetuDynamic(vfloat uDynamicIn){uDynamic = uDynamicIn;}
	void SetMatTempPhase(vfloat MatTmpPhIn) {MaterialTempPhase = MatTmpPhIn;} 
	void SetCurMatTemp(vfloat CurMatTempIn) {CurMaterialTemp = CurMatTempIn;}	//set the local material temperature

	//Dither:
	void SetRandInd1(int RandInd1) {RandIndex1 = RandInd1;}
	void SetRandInd2(int RandInd2) {RandIndex2 = RandInd2;}
	void SetRandPerc1(vfloat RandPerc1) {PercIndex1 = RandPerc1;}
	//Structure:
	void SetSubXOffset(int XOff) {X_Offset = XOff;}
	void SetSubYOffset(int YOff) {Y_Offset = YOff;}
	void SetSubZOffset(int ZOff) {Z_Offset = ZOff;}
	void SetSubXSize(int XDim){if (MatType == INTERNAL) pStructure->Resize(XDim, pStructure->GetVYDim(), pStructure->GetVZDim());}
	void SetSubYSize(int YDim){if (MatType == INTERNAL) pStructure->Resize(pStructure->GetVXDim(), YDim, pStructure->GetVZDim());}
	void SetSubZSize(int ZDim){if (MatType == INTERNAL) pStructure->Resize(pStructure->GetVXDim(), pStructure->GetVYDim(), ZDim);}
	void SetSubRotateAxis(int RotateAxis) {StructRotateAxis = RotateAxis;}
	void SetSubRotateAmount(int RotateAmount) {StructRotateAmount = RotateAmount;}

	//functions to get info about material
	inline const std::string& GetName(void) const {return Name;}
	inline int GetMatType(void) const {return MatType;}
	inline int GetMatModel(void) const {return MatModel;}
	inline int GetFailModel(void) const {return FailModel;}
	inline bool IsVisible(void) const {return Visible;}
	//Color:
	void GetColorf(float* R, float* G, float* B, float* A = NULL) {if(R) *R=Red; if(G) *G=Green; if(B) *B=Blue; if(A) *A=Alpha;}
	void GetColori(int* R, int* G, int* B, int* A = NULL) {if(R) *R=(int)(Red*255); if(G) *G=(int)(Green*255); if(B) *B=(int)(Blue*255); if(A) *A=(int)(Alpha*255);}
	inline float GetRedf(void) const {return Red;}
	inline int GetRedi(void) const {return (int)(Red*255);}
	inline float GetGreenf(void) const {return Green;}
	inline int GetGreeni(void) const {return (int)(Green*255);}
	inline float GetBluef(void) const {return Blue;}
	inline int GetBluei(void) const {return (int)(Blue*255);}
	inline float GetAlphaf(void) const {return Alpha;};
	inline int GetAlphai(void) const {return (int)(Alpha*255);}

	//Model:
	inline vfloat GetElasticMod(void) const {return Elastic_Mod;}
	inline vfloat GetPlasticMod(void) const {return Plastic_Mod;}
	inline vfloat GetYieldStress(void) const {return Yield_Stress;}
	inline vfloat GetFailStress(void) const {return Fail_Stress;}
	inline vfloat GetFailStrain(void) const {return Fail_Strain;}
	vfloat GetModelStiffness(vfloat StrainIn); //returns the stiffness of this material at a given strain
	vfloat GetModelStress(const vfloat StrainIn, bool* const pIsPastYielded, bool* const pIsPastFail) const ; //returns the stress of the current material model at the given strain. Also returns flags for if this is yielded and/or failed.

	//physical:
	inline vfloat GetPoissonsRatio(void) const {return Poissons_Ratio;}
	inline vfloat GetDensity(void) const {return Density;}
	inline vfloat GetCTE(void) const {return CTE;}
	inline vfloat GetuStatic(void) const {return uStatic;}
	inline vfloat GetuDynamic(void) const {return uDynamic;}
	inline double GetCurMatTemp(void) {return CurMaterialTemp;}	//return the local material temperature
	inline double GetMatTempPhase(void) {return MaterialTempPhase;}	//return the phase response of the material to temperature changes (param allows certain materials to warm/cool before/after others) 
	
	//Dither:
	inline int GetRandInd1(void) const {return RandIndex1;};
	inline int GetRandInd2(void) const {return RandIndex2;};
	inline vfloat GetRandPerc1(void) const {return PercIndex1;};
	//Structure:
	inline bool HasLocalStructure(void) const {if (pStructure) return true; else return false;};
	inline CVXC_Structure* GetPStructure(void) const {return pStructure;};
	inline CVXC_Structure& GetStructure(void) const {return *pStructure;};
	inline int GetSubXOffset(void) const {return X_Offset;};
	inline int GetSubYOffset(void) const {return Y_Offset;};
	inline int GetSubZOffset(void) const {return Z_Offset;};
	inline int GetSubRotateAxis(void) const {return StructRotateAxis;};
	inline int GetSubRotateAmount(void) const {return StructRotateAmount;};
	inline int GetSubXSize(void) const {if (MatType == INTERNAL) return pStructure->GetVXDim(); else return 1;};
	inline int GetSubYSize(void) const {if (MatType == INTERNAL) return pStructure->GetVYDim(); else return 1;};
	inline int GetSubZSize(void) const {if (MatType == INTERNAL) return pStructure->GetVZDim(); else return 1;};

#ifdef USE_OPEN_GL
	void SetGLColor(void);
#endif

protected:
	std::string Name; //material name
	int MatType; //keep track of which material this is!
	int MatModel; //keep track of whcih material model we want to use...
	int FailModel; //keep track of criteria for material breaking

	//Display
	float Red, Green, Blue, Alpha;
	bool Visible;
	/*Physical properties
	* MaterialTemp - local storage of the temperature of each material (units in degrees C)
	* MaterialTempPhase - the relative response phase of the temperature change of each material - allows certain materials to heat up before other to get more complicated movements (units in radians)
	*/
	vfloat Density, Elastic_Mod, Yield_Stress, Plastic_Mod, Fail_Stress, Fail_Strain, Poissons_Ratio, CTE, CurMaterialTemp, MaterialTempPhase, uStatic, uDynamic;
	std::vector<vfloat> DStress, DStrain; //arbitrary stress/strain values for material model (in tension only for now)

	//for internal material definition
	CVXC_Structure* pStructure; //create dynamically to keep things small...
	void AllocateIntStruct(void) {DeleteIntStruct(); pStructure = new CVXC_Structure();};
	void AllocateIntStruct(int xD, int yD, int zD) {DeleteIntStruct(); pStructure = new CVXC_Structure(xD, yD, zD);};
	void DeleteIntStruct(void){if (pStructure != NULL){delete pStructure; pStructure = NULL;}};
	int X_Offset, Y_Offset, Z_Offset;
	int StructRotateAxis;
	int StructRotateAmount; //This is currently an enumeration!

	//for random dither material definition
	int RandIndex1, RandIndex2;
	vfloat PercIndex1; //ranges from 0 to 1;
};

//PRSM:
//Pseudo-Random Spatial Map: prsm.h: Spatial random number generator based on the maximally equiditributedcombined Tausworthe-88 generator. By Daniel Ly and Hod Lipson (2010)
//#include <climits>

typedef struct {unsigned long int s1, s2, s3;} taus_state;
unsigned long int rand_seed(unsigned long int x);
unsigned long int taus_get(taus_state* state);
vfloat prsm(vfloat x, vfloat y, vfloat z=0, int k=0);

#endif //CVX_OBJECT_H





