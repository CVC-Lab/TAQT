#pragma once

class Contour2D;

class Contour2DExtractor
{
public:

	Contour2DExtractor(void)
	{
	}

	virtual ~Contour2DExtractor(void)
	{
	}

	virtual Contour2D* extractContour2D(float isoval) = 0;
};
