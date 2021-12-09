#pragma once
#include "icVector.H"
#include <vector>

class LineSegment
{
public:

	// fields
	icVector3 start, end;
	double len;
	double R = 0;
	double G = 0;
	double B = 0;
	double thickness = 1.0;

	// constructors

	LineSegment(icVector3 start_in, icVector3 end_in)
	{
		start = start_in;
		end = end_in;
		len = length(end - start);
	}

	LineSegment(double sx, double sy, double sz, double ex, double ey, double ez)
	{
		start = icVector3(sx, sy, sz);
		end = icVector3(ex, ey, ez);
		len = length(end - start);
	}

	// methods

	icVector3 midpoint()
	{
		icVector3 diff = end - start;
		return start + (0.5 * diff);
	}
};

// PolyLine is a list of connected line segments
typedef std::vector<LineSegment> PolyLine;