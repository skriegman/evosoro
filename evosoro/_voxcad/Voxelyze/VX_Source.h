/*#ifndef VX_SOURCE_H
#define VX_SOURCE_H

#include <vector>
#include <string>
#include "Utils/Vec3D.h" //use this for portability, instead of Vec3D()
*/
/* Author of this file: Francesco Corucci - f.corucci@sssup.it */
/*
// Known environmental sources
enum SourceType {
	HEAT
};

class VX_Source //container for information about environmental sources of information
{
public:

	// Constructor
	VX_Source(std::string source_name, SourceType source_type, bool source_enabled, float x, float y, float z) // Constructor
	{
		this->name = source_name;
		this->type = source_type;
		this->enabled = source_enabled;
		this->position.setX(x);
		this->position.setY(y);
		this->position.setZ(z);
	}

	// Get methods
	std::string getSourceName(){ return this->name; }
	SourceType getSourceType() { return this->type; }
	bool isSourceEnabled(){ return this->enabled; }
	Vec3D<float> getSourcePosition(){return this->position;}

	// Set methods
	bool setSourceName(std::string source_name){ this->name = source_name; }
	bool setSourceType(SourceType source_type){ this->type = source_type; }
	bool setSourceEnabled(bool source_enabled){ this->enabled = source_enabled; }
	bool setSourcePosition(Vec3D<float> source_position){ 	this->position = source_position; }	

protected:
	std::string name;
	SourceType type;
	bool enabled;
	Vec3D<float> position; // The position of the source in the 3D space
	
};

#endif //VX_SOURCE_H




*/
