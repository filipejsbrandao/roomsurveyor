//
// RandomOgon.cs
//
// Author:
//       Filipe Jorge da Silva Brandao
//
// Copyright (c) 2020 ©2020 Filipe Brandao
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace RoomSurveyor
{
    public class RandomOgon : GH_Component
    {
        public RandomOgon()
          : base("Random Ogon Generator", "Ogon1",
            "Random Ogon Generator - CutPaste. This algorithm generates a random orthogonal polygon by cutting or pasting rectangles to a seed rectangle. Is follows a similiar approach to one described by Tomas and Bajuelos 2004 InflatePaste and InflateCut algorithms, with some improvements to adapt it to floating point and reduce the likelyhood of failure in cuts. It implements one rule to prevent generating features that are smaller than a given dimension. This algorithm may fail if a point is outside the bounding box of the ogon or in certain interior very particular points.",
            "RoomSurveyor", "Polygon Generators")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddRectangleParameter("Region", "R", "A rectangular region within which the ogon will be generated", GH_ParamAccess.item);
            pManager.AddNumberParameter("Minimum Feature", "F", "The length of the smallest feature", GH_ParamAccess.item, 0.10);
            pManager.AddBooleanParameter("Area Bias", "B", "If true the algorithm will favor Cuts or Pastes that minimize the area increase or decrease of the ogon", GH_ParamAccess.item, true);
            pManager.AddIntegerParameter("Number of sides", "N", "The number of sides of the generated Ogon. Only pair numbers are allowed", GH_ParamAccess.item, 6);
            pManager.AddIntegerParameter("Seed", "S", "An integer to be used as a seed for the random generator. If nothing is provided the system clock is used. Provide a seed if you wish to debug your algorithm", GH_ParamAccess.item);
            pManager[4].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddCurveParameter("Ogon", "O", "An orthogonal polygon with n sides", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            int n = 4;
            Rectangle3d r = new Rectangle3d();
            double min_feature = 0.10;
            bool areaBias = true;
            int seed = 0;

            if (!DA.GetData(0, ref r)) return;
            if (!DA.GetData(1, ref min_feature)) return;
            if (!DA.GetData(2, ref areaBias)) return;
            if (!DA.GetData(3, ref n)) return;
            bool user = DA.GetData(4, ref seed);

            //Component.Message = "LATEST VERSION";
            if (!r.IsValid)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The rectangle is invalid");
                return;
            }

            if (double.IsNaN(min_feature))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "min_feature must be a valid number");
                return;
            }

            if (n < 4 || n % 2 == 1)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "n must be pair and larger than 4");
                return;
            }
            else if (n == 4)
            {
                return;
            }

            int rCorners = n / 2 - 2;
            double tol = min_feature; //overall tolerance (we don't want new corners closer than tol to the existing corners

            //first ensure that the provided rectangle is in the xy plane with origin at 0,0
            Plane xy = Plane.WorldXY;
            double maxX = Math.Abs(r.Width);
            double maxY = Math.Abs(r.Height);

            Rectangle3d newR = new Rectangle3d(Plane.WorldXY, maxX, maxY);

            Polyline ogon = newR.ToPolyline();//now make a polyline out of a rectangle (ideally this should be a circular linked list)
            Random rand = new Random();
            Random userRand = new Random(seed);
            for (int i = 0; i < rCorners; i++)
            {

                bool sucess;
                do
                {
                    bool validPt = false;
                    Point3d randomPt = Point3d.Unset;
                    do
                    {
                        double x = (!user) ? maxX * rand.NextDouble() : maxX * userRand.NextDouble();
                        double y = (!user) ? maxY * rand.NextDouble() : maxY * userRand.NextDouble();
                        randomPt = new Point3d(x, y, 0);
                        if (randomPt.DistanceToSquared(ogon.ClosestPoint(randomPt)) < Math.Pow(tol, 2))
                        {
                            //Print("The point is too close to the ogon");
                            //return;
                        }
                        else
                        {
                            validPt = true;
                        }
                    } while (!validPt);

                    sucess = (!user) ? CutPaste(randomPt, ref ogon, maxX, maxY, areaBias, rand.Next(0, 7), tol) : CutPaste(randomPt, ref ogon, maxX, maxY, areaBias, userRand.Next(0, 7), tol);

                } while (!sucess);
            }

            Plane rOrigin = new Plane(r.Corner(0), r.Corner(1), r.Corner(3));
            Transform m = Transform.PlaneToPlane(xy, rOrigin);
            ogon.Transform(m);


            DA.SetData(0, ogon);

        }
        /// <summary>
        /// CutPaste Algorithm for Orthogonal Polygon (ogon) generation. Provided a point and an ogon, the algorithm adds two corners to the ogon.
        /// </summary>
        /// <param name="pt">A point in or out of the ogon. This point will become one of the corners of the ogon</param>
        /// <param name="ogon">A closed planar polyline on the XY plane representing an orthogonal polygon</param>
        /// <param name="maxX">The dimension of the bounding box in the X-axis</param>
        /// <param name="maxY">The dimension of the bounding box in the Y-axis</param>
        /// <param name="maximizeArea">If true the algorithm will select the cut/add candidate that minimizes the area change</param>
        /// <param name="s"></param>
        /// <param name="tol"></param>
        /// <returns></returns>
        private bool CutPaste(Point3d pt, ref Polyline ogon, double maxX, double maxY, bool maximizeArea, int s, double tol)
        {
            bool sucess = false;
            List<Line> sweepLines = new List<Line>();
            tol = (tol < Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance) ? Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance : tol;

            sweepLines.Add(new Line(pt, new Point3d(pt.X, pt.Y + maxY * 2, 0)));//North
            sweepLines.Add(new Line(pt, new Point3d(pt.X - maxX * 2, pt.Y, 0)));//East
            sweepLines.Add(new Line(pt, new Point3d(pt.X, pt.Y - maxY * 2, 0)));//South
            sweepLines.Add(new Line(pt, new Point3d(pt.X + maxX * 2, pt.Y, 0)));//West
            List<Point3d> ccx = new List<Point3d>();

            ccx.AddRange(GetIntersections(sweepLines, ogon, pt));

            List<Rectangle3d> candidates = new List<Rectangle3d>();
            List<double> paramList = new List<double>();
            List<List<Point3d>> cornersToAdd = new List<List<Point3d>>();
            List<List<int>> cornersToRemove = new List<List<int>>();

            foreach (Point3d p in ccx)
                paramList.Add(ogon.ClosestParameter(p));

            if (Math.Abs(PnPoly(pt, ogon.ToArray(), ogon.Count - 1)) > 0)
            {//the point is inside
                //Print("the point is inside");
                //************************************************************************CCW INSIDE***********************************
                for (int i = 0; i < ccx.Count; i++)
                {
                    int prev = (int)paramList[i];
                    int next = (prev + 1) % ogon.SegmentCount;
                    if (ccx[i].DistanceToSquared(ogon[next]) > tol * tol && ccx[i].DistanceToSquared(ogon[prev]) > tol * tol)
                    {
                        Rectangle3d ccw = new Rectangle3d(Plane.WorldXY, pt, ogon[next]);
                        int r_next = NextCorner(pt, ccw, true);
                        int r_prev = NextCorner(pt, ccw, false);
                        //A rule to reject coincident ccx[i] with ogon[next] for ccw or [prev] for cw
                        bool rule1 = (Math.Abs(ccw.Height) > tol && Math.Abs(ccw.Width) > tol) ? true : false;
                        if (!IsCornerInside(ogon, ccw) && CoincidentCorner(ogon, ccw).Count < 3 && rule1)
                        {
                            List<int> corners = CoincidentCorner(ogon, ccw);
                            List<Point3d> newCorners = new List<Point3d>() { ccw.Corner(r_next), pt, ccw.Corner(r_prev) };
                            if (corners.Count == 2)
                            {
                                newCorners.Add(ogon[corners[1]]);
                            }
                            candidates.Add(ccw);
                            cornersToAdd.Add(newCorners);
                            cornersToRemove.Add(corners);
                        }
                    }
                }
                //************************************************************************CW INSIDE************************************
                for (int i = 0; i < ccx.Count; i++)
                {
                    int prev = (int)paramList[i];
                    int next = (prev + 1) % ogon.SegmentCount;
                    int i_prev = (i == 0) ? ccx.Count - 1 : (i - 1) % ccx.Count;
                    if (ccx[i].DistanceToSquared(ogon[next]) > tol * tol && ccx[i].DistanceToSquared(ogon[prev]) > tol * tol)
                    {
                        if (prev != ((int)paramList[i_prev] + 1) % ogon.SegmentCount)
                        {//If the current previous is not equal to the last next (It is necessary to remove the duplicates so that all options have a fair chance of being choosen
                            Rectangle3d cw = new Rectangle3d(Plane.WorldXY, pt, ogon[prev]);
                            int r_next = NextCorner(pt, cw, false);
                            int r_prev = NextCorner(pt, cw, true);
                            //A rule to reject coincident ccx[i] with ogon[next] for ccw or [prev] for cw
                            bool rule1 = (Math.Abs(cw.Height) > tol && Math.Abs(cw.Width) > tol) ? true : false;
                            if (!IsCornerInside(ogon, cw) && CoincidentCorner(ogon, cw).Count < 3 && rule1)
                            {
                                List<int> corners = CoincidentCorner(ogon, cw);
                                List<Point3d> newCorners = new List<Point3d>() { cw.Corner(r_next), pt, cw.Corner(r_prev) };//the current ccx_point (C) and P
                                if (corners.Count == 2)
                                {
                                    newCorners.Add(ogon[corners[0]]);
                                }
                                candidates.Add(cw);
                                newCorners.Reverse();
                                cornersToAdd.Add(newCorners);
                                cornersToRemove.Add(corners);
                            }
                        }
                    }
                }
            }
            else
            {// the point is outside
                //************************************************************************CCW OUTSIDE**********************************
                for (int i = 0; i < ccx.Count; i++)
                {
                    int prev = (int)paramList[i];
                    int next = (prev + 1) % ogon.SegmentCount;
                    if (ccx[i].DistanceToSquared(ogon[next]) > tol * tol && ccx[i].DistanceToSquared(ogon[prev]) > tol * tol)
                    {
                        Rectangle3d ccw = new Rectangle3d(Plane.WorldXY, pt, ogon[prev]);
                        int r_next = NextCorner(pt, ccw, true);
                        int r_prev = NextCorner(pt, ccw, false);
                        //A rule to reject coincident ccx[i] with ogon[next] for ccw or [prev] for cw
                        bool rule1 = (Math.Abs(ccw.Height) > tol && Math.Abs(ccw.Width) > tol) ? true : false;
                        if (!IsCornerInside(ogon, ccw) && CoincidentCorner(ogon, ccw).Count < 3 && rule1)
                        {
                            List<int> corners = CoincidentCorner(ogon, ccw);
                            List<Point3d> newCorners = new List<Point3d>() { ccw.Corner(r_next), pt, ccw.Corner(r_prev) };//the current ccx_point (C) and P
                            if (corners.Count == 2)
                            {
                                newCorners.Add(ogon[corners[0]]);
                            }
                            candidates.Add(ccw);
                            newCorners.Reverse();
                            cornersToAdd.Add(newCorners);
                            cornersToRemove.Add(corners);
                        }
                    }
                }
                //************************************************************************CW OUTSIDE***********************************
                for (int i = 0; i < ccx.Count; i++)
                {
                    int prev = (int)paramList[i];
                    int next = (prev + 1) % ogon.SegmentCount;
                    int i_prev = (i == 0) ? ccx.Count - 1 : (i - 1) % ccx.Count;
                    if (ccx[i].DistanceToSquared(ogon[next]) > tol * tol && ccx[i].DistanceToSquared(ogon[prev]) > tol * tol)
                    {
                        if (prev != ((int)paramList[i_prev] - 1) % ogon.SegmentCount)
                        {// THIS CONDITION MIGHT CAUSE INDEX OUT OF RANGE ERRORS!!
                            Rectangle3d cw = new Rectangle3d(Plane.WorldXY, pt, ogon[next]);
                            int r_next = NextCorner(pt, cw, false);
                            int r_prev = NextCorner(pt, cw, true);
                            //A rule to reject coincident ccx[i] with ogon[next] for ccw or [prev] for cw
                            bool rule1 = (Math.Abs(cw.Height) > tol && Math.Abs(cw.Width) > tol) ? true : false;
                            if (!IsCornerInside(ogon, cw) && CoincidentCorner(ogon, cw).Count < 3 && rule1)
                            {
                                List<int> corners = CoincidentCorner(ogon, cw);
                                List<Point3d> newCorners = new List<Point3d>() { cw.Corner(r_next), pt, cw.Corner(r_prev) };//the current ccx_point (C) and P
                                if (corners.Count == 2)
                                {
                                    newCorners.Add(ogon[corners[1]]);
                                }
                                candidates.Add(cw);
                                cornersToAdd.Add(newCorners);
                                cornersToRemove.Add(corners);
                            }
                        }
                    }
                }
            }

            if (candidates.Count != 0 && !maximizeArea)
            {
                int sel = s % candidates.Count;
                ogon.RemoveAt(ogon.Count - 1);
                ogon.RemoveRange(cornersToRemove[sel][0], cornersToRemove[sel].Count);
                ogon.InsertRange(cornersToRemove[sel][0], cornersToAdd[sel]);
                ogon.Add(ogon[0]);
                sucess = true;
            }
            else if (candidates.Count != 0 && maximizeArea)
            {
                int sel = 0;
                double area = candidates[0].Area;
                for (int i = 1; i < candidates.Count; i++)
                {
                    if (area > candidates[i].Area)
                    {
                        area = candidates[i].Area;
                        sel = i;
                    }
                }
                ogon.RemoveAt(ogon.Count - 1);
                ogon.RemoveRange(cornersToRemove[sel][0], cornersToRemove[sel].Count);
                ogon.InsertRange(cornersToRemove[sel][0], cornersToAdd[sel]);
                ogon.Add(ogon[0]);
                sucess = true;
            }
            return sucess;
        }
        /// <summary>
        /// Rule 2 tests if the proportion of the width (size along X axis) of the rectangle over the distance from P to vPoint is larger than one if the
        /// width is smaller than minConvex and the same for the other axis.
        /// </summary>
        /// <param name="currToP">the squared distance of P to current_ccx point</param>
        /// <param name="perpToP">the squared distance of P to either previous (in CCW) or next (in CW) ccx point</param>
        /// <param name="otherToP">the squared distance of P to the other intersection point, next in CCW, previous in CW</param>
        /// <param name="oppositeToP">the squared distance of P to the opposite to current_ccx intersection point </param>
        /// <param name="minConvex">the minimun allowable distance between opposite internal walls of habitable spaces</param>
        /// <returns>true if the candidate rectangle cut will comply with the rule</returns>
        private bool Rule2(double currToP, double perpToP, double otherToP, double oppositeToP, double minConvex)
        {
            bool rule = true;
            double h_proportion = Math.Round(otherToP / currToP, 3); //We should think a bit more about how to handle tolerance settings
            double v_proportion = Math.Round(oppositeToP / perpToP, 3);
            if (h_proportion < 1.0 && otherToP < minConvex * minConvex || v_proportion < 1.0 && oppositeToP < minConvex * minConvex)
                return false;

            return rule;
        }

        /// <summary>
        /// This method tests if any of the ogon corners lies inside the minimum convex space required by the candidate rectangle cut.
        /// </summary>
        /// <param name="minConvexSpace">A closed polyline around the minimum convex free space.</param>
        /// <param name="ogon">The existing ogon to be cut or pasted.</param>
        /// <param name="n">The number of edges of the minConvexSpace closed polyline.</param>
        /// <returns></returns>
        private bool Rule3(Polyline minConvexSpace, Polyline ogon, int n)
        {
            bool rule = true;
            Point3d[] poly = minConvexSpace.ToArray();

            for (int i = 0; i < ogon.Count - 1; i++)
            {
                rule = (PnPoly(ogon[i], poly, n) == 0) ? true : false;
                if (!rule)
                    break;
            }

            return rule;
        }

        /// <summary>
        /// Generates the minimum convex free space for the candidate rectangle cut as a closed polyline to be tested by Rule2 method.
        /// </summary>
        /// <param name="pts">The list of points to be added to the ogon in CCW.</param>
        /// <param name="minConvex">The smallest dimension of the minimum habitable convex space.</param>
        /// <param name="orientation">true for CCW and false for CW.</param>
        /// <returns>A closed polyline around rhe minimum convex free space.</returns>
        private Polyline MinConvexSpace(List<Point3d> pts, double minConvex, bool orientation)
        {

            Polyline poly = new Polyline(pts);
            Vector3d first_normal = new Vector3d(poly[1] - poly[0]);//1 is P, 0 is C
            Vector3d second_normal = new Vector3d(poly[2] - poly[1]);//1 is P, 2 is C'
            first_normal.Unitize();
            first_normal = first_normal * minConvex;
            second_normal.Unitize();
            second_normal = second_normal * minConvex;
            bool rotate = (orientation) ? first_normal.Rotate(Math.PI / 2, Vector3d.ZAxis) : first_normal.Rotate(-Math.PI / 2, Vector3d.ZAxis);// se a Lista de pontos dos rectangulos CW for fornecida já ordenada CCW podemos rodar sempre na mesma direcção
            rotate = (orientation) ? second_normal.Rotate(Math.PI / 2, Vector3d.ZAxis) : second_normal.Rotate(-Math.PI / 2, Vector3d.ZAxis);
            if (pts.Count == 4)
            {
                Vector3d third_normal = new Vector3d(poly[3] - poly[2]);
                third_normal.Unitize();
                third_normal = third_normal * minConvex;
                rotate = (orientation) ? third_normal.Rotate(Math.PI / 2, Vector3d.ZAxis) : third_normal.Rotate(-Math.PI / 2, Vector3d.ZAxis);
                poly.Add(pts[3] + third_normal + second_normal);
                poly.Add(pts[2] + third_normal + second_normal);
            }
            else
            {
                poly.Add(pts[2] + first_normal + second_normal);
            }
            poly.Add(pts[1] + first_normal + second_normal);
            poly.Add(pts[0] + first_normal + second_normal);
            poly.Add(pts[0]);

            return poly;

        }
        /// <summary>
        /// The objective of this method is to find the closest intersections of the axis oriented semi-segments that start in point pt with the ogon.
        /// THE METHOD NEEDS TO BE EXTENDED FOR IT TO WORK CORRECTLY WITH EXTERNAL POINTS OR THE CASE WHEN NO POINTS ARE FOUND NEEDS TO BE DEALT WITH ELSEWHERE
        /// </summary>
        /// <param name="sweepLines">A list with the four axis-oriented segments that are guaranteed to intersect the ogon</param>
        /// <param name="ogon">a polyline representing the orthogonal polygon (ogon)</param>
        /// <param name="pt">a random point inside or outside the ogon</param>
        /// <returns>A list of points on the ogon which are the projection of pt along the orthogonal axis</returns>
        private List<Point3d> GetIntersections(List<Line> sweepLines, Polyline ogon, Point3d pt)
        {
            Curve c1 = ogon.ToNurbsCurve();
            const double intersection_tolerance = 0.001;
            const double overlap_tolerance = 0.0;
            List<Point3d> ccx = new List<Point3d>();

            foreach (Line ln in sweepLines)
            {
                Curve c2 = ln.ToNurbsCurve();
                var events = Rhino.Geometry.Intersect.Intersection.CurveCurve(c1, c2, intersection_tolerance, overlap_tolerance);

                if (events.Count != 0)
                {
                    Point3d p_ccx = events[0].PointA;
                    for (int i = 0; i < events.Count; i++)
                    {
                        var ccx_event = events[i];
                        if (ccx_event.PointA.DistanceToSquared(pt) < p_ccx.DistanceToSquared(pt))
                            p_ccx = ccx_event.PointA;
                    }
                    ccx.Add(p_ccx);
                }
            }
            return ccx;
        }

        /// <summary>
        /// This methods checks the relation of each corner of the ogon with a given rectangle and returns the indices of the coincident corners
        /// </summary>
        /// <param name="ogon">A polyline representing the ogon</param>
        /// <param name="rectangle">A rectangle</param>
        /// <returns>A list with the indices of each coincident point</returns>
        private List<int> CoincidentCorner(Polyline ogon, Rectangle3d rectangle)
        {
            List<int> relation = new List<int>();

            for (int i = 0; i < ogon.SegmentCount; i++)
            {
                if (rectangle.Contains(ogon[i]) == PointContainment.Coincident)
                    relation.Add(i);
            }
            return relation;
        }

        /// <summary>
        /// This method checks if the rectangle contains any point of the polygon
        /// </summary>
        /// <param name="ogon">A polyline representing the ogon</param>
        /// <param name="rectangle">A rectangle</param>
        /// <returns>true if any corner is Inside the rectangle, false if all the points are either Coincident, Outside or Unset</returns>
        private bool IsCornerInside(Polyline ogon, Rectangle3d rectangle)
        {
            bool inside = false;

            for (int i = 0; i < ogon.SegmentCount; i++)
            {
                if (rectangle.Contains(ogon[i]) == PointContainment.Inside)
                    return true;
            }
            return inside;
        }

        /// <summary>
        /// Winding Number test for a point in a polygon based on http://geomalgorithms.com/a03-_inclusion.html
        /// Copyright 2000 softSurfer, 2012 Dan Sunday
        /// This code may be freely used and modified for any purpose
        /// providing that this copyright notice is included with it.
        /// SoftSurfer makes no warranty for this code, and cannot be held
        /// liable for any real or imagined damage resulting from its use.
        /// Users of this code must verify correctness for their application.
        /// Ported to C# by Filipe Brandão in January 2019
        /// </summary>
        /// <returns> the winding number (=0 only when P is outside)</returns>
        /// <param name="Pt">Point.</param>
        /// <param name="poly">the Poly as an array of points where <paramref name="poly"/>[0] = poly[n]</param>
        /// <param name="n"> the number of edges in the polygon</param>
        private int PnPoly(Point3d Pt, Point3d[] poly, int n)
        {
            int wn = 0; //the winding number counter

            for (int i = 0; i < n; i++)
            { //edge from V[i] to V[i+1]
                if (poly[i].Y <= Pt.Y)
                { // start Y <= Pt.Y
                    if (poly[i + 1].Y > Pt.Y) // an upward crossing
                        if (IsLeft(poly[i], poly[i + 1], Pt) > 0) // Pt is left of edge
                            ++wn; // have a valid up intersect
                }
                else
                { // start Y > Pt.Y (no test needed)
                    if (poly[i + 1].Y <= Pt.Y) // a downward crossing
                        if (IsLeft(poly[i], poly[i + 1], Pt) < 0) // Pt is right of edge
                            --wn; // have a valid down intersect
                }

            }
            return wn;
        }
        /// <summary>
        /// Determines if a third point is Left, On or Right of an infinite 2D line defined by two points
        /// </summary>
        /// <returns>Returns a positive number if it is Left, a negative number if it is Right and 0 if it is on the line</returns>
        /// <param name="P0">P0 - the first point</param>
        /// <param name="P1">P1 - the second point</param>
        /// <param name="P2">P2 - the point to test</param>

        private double IsLeft(Point3d P0, Point3d P1, Point3d P2)
        {
            return ((P1.X - P0.X) * (P2.Y - P0.Y) - (P2.X - P0.X) * (P1.Y - P0.Y));
        }
        /// <summary>
        /// Returns the index of the next corner of a rectangle given a point that corresponds to the current corner
        /// </summary>
        /// <param name="startPt">A point that is close or coincident to one of the corners of the rectangle</param>
        /// <param name="r">The rectangle</param>
        /// <param name="ccw">If true it returns the next corner in the counter-clockwise direction</param>
        /// <returns>the index of the next corner</returns>
        private int NextCorner(Point3d startPt, Rectangle3d r, bool ccw)
        {
            int i = 0;
            double vX = r.Center.X - startPt.X;
            double vY = r.Center.Y - startPt.Y;

            if (vY > 0 && vX < 0)
            {//Caso 4
                i = 1;
            }
            else if (vY < 0 && vX < 0)
            {//Caso 1
                i = 2;
            }
            else if (vY < 0 && vX > 0)
            {//Caso 2
                i = 3;
            }

            i = (ccw) ? (i + 1) % 4 : (i == 0) ? 3 : i - 1;

            return i;
        }

        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return Properties.Resources.RandomOgon_Icon;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("8dd480ad-0f62-4770-9dd1-232f645fbd82"); }
        }
    }
}
