using System;
using System.Collections.Generic;
using Rhino;
using Rhino.Geometry;

namespace RoomSurveyor
{
    /// <summary>
    /// Algorithms based on the rotating calliper proposed by Shamos, Houle and Toussaint
    /// </summary>
    public static class RotatingCallipers
    {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="poly">A simple polygon, i.e. non-self intersecting</param>
        /// <returns></returns>
        public static double Diameter2D(Polyline poly)
        {
            double diameter;

            //1 - We make sure the polygon is convex by obtaining its convex hull
            Polyline cHull = ConvexHull.Lee83(poly);

            //2 - We select an edge at random and search for the vertex of the polygon that is farthest from it
            int firstCorner = 0;
            int firstAntipodal = Antipodal(poly, firstCorner);

            //One would be temped now to assume that we could change the Antipodal algorithm to just go
            //always left (or right) and thus reduce its complexity, but if you run the algorithm you will notice
            //that hopping from one starting corner to the next you need to do a complete 360 turn to find the
            // minimal or maximal distances for each corner of the polygon. But with the rotating callipers you only need
            //to rotate half, since at every step you are doing a minimal rotation.
            //3 - Now we run the rotating calipers algorithm for 180 degrees

            diameter = poly.BoundingBox.Diagonal.Length;//of course this is not true

            return diameter;
        }

        public static double Width2D(Polyline poly)
        {
            double width;

            ///TODO

            width = poly.BoundingBox.Diagonal.Length;//of course this is not true

            return width;
        }

        public static Rectangle3d MBR2D(Polyline poly)
        {
            Rectangle3d mbr = new Rectangle3d();

            return mbr;
        }

        /// <summary>
        /// MinMax Distance: a O(log n) algorithm described by Kurozumi and Davis 1982.
        /// The algorithm described in the paper is more complex in this case we are onlyinterested in the Max distance in a convex polygon.
        /// Houle and Toussaint refer to this algorithm as a binary search in their 1988 paper on rotating caliper
        /// </summary>
        /// <param name="poly">The polygon</param>
        /// <param name="startVertex">the index of the start point of the edge</param>
        /// <returns>The index of the corner of the polygon that is farthest from the edge defined by the start point and the next corner on the polygon</returns>
        public static int Antipodal(Polyline poly, int startVertex)
        {
            int n = poly.Count - 1;// a polygon is a closed polyline. start and end points are repeated
            startVertex %= n; //just making sure the index stays within the bounds of the polygon
            int i = (startVertex + n / 2) % n + 1;
            int prevI = (i == 0) ? n : i - 1;
            int nextI = (i + 1) % n;
            int endVertex = (startVertex + 1) % n;
            int antipodal = startVertex;
            bool done = false;

            double thisDistance = Distance(poly[startVertex], poly[endVertex], poly[i]);
            double prevDistance = Distance(poly[startVertex], poly[endVertex], poly[prevI]);
            double nextDistance = Distance(poly[startVertex], poly[endVertex], poly[nextI]);

            while (!done)
            {
                if (thisDistance >= prevDistance && thisDistance >= nextDistance)
                {
                    antipodal = i;
                    done = true;
                }
                else if (thisDistance < prevDistance)
                {
                    i--;
                    nextDistance = thisDistance;
                    thisDistance = prevDistance;
                    prevI = (i == 0) ? n : i - 1;
                    prevDistance = Distance(poly[startVertex], poly[endVertex], poly[prevI]);
                }
                else if (thisDistance < nextDistance)
                {
                    i++;
                    prevDistance = thisDistance;
                    thisDistance = nextDistance;
                    nextI = (i + 1) % n;
                    nextDistance = Distance(poly[startVertex], poly[endVertex], poly[nextI]);
                }
            }
            return antipodal;
        }

        /// <summary>
        /// given two points on a line an a third point compute the shortest distance from
        /// the point pt to the infinite line defined by the first two points
        /// explanation in http://geomalgorithms.com/a02-_lines.html
        /// </summary>
        /// <param name="from">The start point of the infinite line</param>
        /// <param name="to">A point that the line goes through, ie, that defines the direction of the line</param>
        /// <param name="pt">The point we want to measure the perpendicular distance from</param>
        /// <returns></returns>
        public static double Distance(Point3d from, Point3d to, Point3d pt)
        {
            //distance
            return Math.Abs((from.Y - to.Y) * pt.X + (to.X - from.X) * pt.Y + (from.X * to.Y - to.X * from.Y)) / Math.Sqrt(Math.Pow((to.X - from.X), 2) + Math.Pow((to.Y - from.Y), 2));
        }
    }
}
