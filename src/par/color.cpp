#include "color.h"

/*
   Calculate RGB from HSV, reverse of RGB2HSV()
   Hue is in degrees
   Lightness is between 0 and 1
   Saturation is between 0 and 1
*/
RGB HSV2RGB(HSV c1)
{
   RGB c2,sat;

   while (c1.h < 0)
      c1.h += 360;
   while (c1.h > 360)
      c1.h -= 360;

   if (c1.h < 120) {
      sat.r = (120 - c1.h) / 60.0f;
      sat.g = c1.h / 60.0f;
      sat.b = 0;
   } else if (c1.h < 240) {
      sat.r = 0;
      sat.g = (240 - c1.h) / 60.0f;
      sat.b = (c1.h - 120) / 60.0f;
   } else {
      sat.r = (c1.h - 240) / 60.0f;
      sat.g = 0;
      sat.b = (360 - c1.h) / 60.0f;
   }
   sat.r = MIN(sat.r,1);
   sat.g = MIN(sat.g,1);
   sat.b = MIN(sat.b,1);

   c2.r = (1 - c1.s + c1.s * sat.r) * c1.v;
   c2.g = (1 - c1.s + c1.s * sat.g) * c1.v;
   c2.b = (1 - c1.s + c1.s * sat.b) * c1.v;

   return(c2);
}

/*
   Calculate HSV from RGB
   Hue is in degrees
   Lightness is betweeen 0 and 1
   Saturation is between 0 and 1
*/
HSV RGB2HSV(RGB c1)
{
   float themin,themax,delta;
   HSV c2;

   themin = MIN(c1.r,MIN(c1.g,c1.b));
   themax = MAX(c1.r,MAX(c1.g,c1.b));
   delta = themax - themin;
   c2.v = themax;
   c2.s = 0;
   if (themax > 0)
      c2.s = delta / themax;
   c2.h = 0;
   if (delta > 0) {
      if (themax == c1.r && themax != c1.g)
         c2.h += (c1.g - c1.b) / delta;
      if (themax == c1.g && themax != c1.b)
         c2.h += (2 + (c1.b - c1.r) / delta);
      if (themax == c1.b && themax != c1.r)
         c2.h += (4 + (c1.r - c1.g) / delta);
      c2.h *= 60;
   }
   return(c2);
}
