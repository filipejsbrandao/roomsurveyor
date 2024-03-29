﻿//
// IsConvex.cs
//
// Author:
//       Filipe Jorge da Silva Brandao
//
// Copyright (c) 2019 ©2019 Filipe Brandao
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
using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace RoomSurveyor
{
    public class IsConvex : GH_Component
    {

        public IsConvex()
          : base("Is Polygon Convex", "IsCx",
            "Checks if a polygon is Convex. Rory Daulton algorithm ported to Java by Guilherme Campos Hazan and ported to C# by Filipe Jorge da Silva Brandão",
            "RoomSurveyor", "Ises")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Polygon", "P", "A closed polyline or a polygon with 3 or more sides", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddBooleanParameter("Convex", "C", "True if polygon is convex, false if the polygon is either concave or complex, i.e. self-intersecting", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Polyline poly = new Polyline();
            Curve curve = poly.ToNurbsCurve();
            bool convex;

            if (!DA.GetData(0, ref curve)) return;

            if (!curve.IsPlanar())
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The polygon must be a planar polyline");
                return;
            }

            //To ensure this always works we must move the polygon to the plane XY and then return it to it original position
            curve.TryGetPlane(out Plane uPlane);

            if (curve.TryGetPolyline(out poly)) { }
            else
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "P input expects a polygon to be provided");
                return;
            }

            if (poly.IsValid == false)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Closed polylines must have at least 2 segments");
                return;
            }

            if (poly.IsClosed == false)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "A closed polyline must be supplied");
                return;
            }

            Transform transform = Transform.PlaneToPlane(uPlane, Plane.WorldXY);
            poly.Transform(transform);

            convex = IsPolyConvex(poly, false);

            DA.SetData(0, convex);

        }

        /// <summary>
        /// Checks if a polygon represented by a closed polyline is Convex or Not Convex
        /// </summary>
        /// <returns><c>true</c>, if the polygon is convex, <c>false</c> if the polygon is either concave or complex, i.e. self-intersecting.</returns>
        /// <param name="polyline">Polyline.</param>
        public static bool IsPolyConvex(Polyline polyline, bool integer)
        {
            bool convex;

            if (integer) { 
                int[] xcoord = new int[polyline.Count];
                int[] ycoord = new int[polyline.Count];
                int vertex = polyline.Count;

                for (int i = 0; i < xcoord.Length; i++)
                {
                    xcoord[i] = (int)polyline[i].X;
                    ycoord[i] = (int)polyline[i].Y;
                }
                convex = IsPolyConvex(xcoord, ycoord, vertex);
            }
            else
            {
                double[] xcoord = new double[polyline.Count];
                double[] ycoord = new double[polyline.Count];
                int vertex = polyline.Count;

                for(int i= 0; i < xcoord.Length; i++)
                {
                    xcoord[i] = polyline[i].X;
                    ycoord[i] = polyline[i].Y;
                }

                convex = IsPolyConvex(xcoord, ycoord, vertex);
            }


            return convex;
        }

        /// <summary>
        /// Checks if a polygon is Convex.
        /// Rory Daulton algorithm ported to Java by Guilherme Campos Hazan and ported to C# by Filipe Jorge da Silva Brandão
        /// Precision issues with Rhino implementation led to changes from using integer types to double types.
        /// </summary>
        /// <returns><c>true</c>, if convex is convex, <c>false</c> otherwise.</returns>
        /// <param name="x">The x coordinate.</param>
        /// <param name="y">The y coordinate.</param>
        /// <param name="n">N.</param>
        private static bool IsPolyConvex(int[] x, int[] y, int n)
        {
            const double TWO_PI = 2 * Math.PI;

            //in rhino all lists are zero based so the code could be changed to remove the need to explicitly require the user to state that
            int ibase = 0;
            // points is 'strictly convex': points are valid, side lengths non-zero, interior angles are strictly between zero and a straight
            // angle, and the polygon does not intersect itself.
            // NOTES:  1.  Algorithm: the signed changes of the direction angles from one side to the next side must be all positive or
            // all negative, and their sum must equal plus-or-minus one full turn (2 pi radians). Also check for too few,
            // invalid, or repeated points.
            //      2.  No check is explicitly done for zero internal angles(180 degree direction-change angle) as this is covered
            // in other ways, including the `n < 3` check.
            // needed for any bad points or direction changes
            // Check for too few points
            if (n <= 3)
            {
                return true;
            }

            if (x[ibase] == x[n - 1] && y[ibase] == y[n - 1]) // if its a closed polygon, ignore last vertex
                n--;
            // Get starting information
            int old_x = x[n - 2], old_y = y[n - 2];
            int new_x = x[n - 1], new_y = y[n - 1];
            double new_direction = Math.Atan2(new_y - old_y, new_x - old_x), old_direction;
            double angle_sum = 0.0, orientation = 0;
            // Check each point (the side ending there, its angle) and accum. angles for ndx, newpoint in enumerate(polygon):
            for (int i = 0; i < n; i++)
            {
                // Update point coordinates and side directions, check side length
                old_x = new_x; old_y = new_y; old_direction = new_direction;
                int p = ibase++;
                new_x = x[p]; new_y = y[p];
                new_direction = Math.Atan2(new_y - old_y, new_x - old_x);
                if (old_x == new_x && old_y == new_y)
                    return false; // repeated consecutive points
                                  // Calculate & check the normalized direction-change angle
                double angle = new_direction - old_direction;
                if (angle <= -Math.PI)
                    angle += TWO_PI;  // make it in half-open interval (-Pi, Pi]
                else if (angle > Math.PI)
                    angle -= TWO_PI;
                if (i == 0)  // if first time through loop, initialize orientation
                {
                    if (angle == 0.0) return false;
                    orientation = angle > 0 ? 1 : -1;
                }
                else  // if other time through loop, check orientation is stable
                  if (orientation * angle <= 0)  // not both pos. or both neg.
                    return false;
                // Accumulate the direction-change angle
                angle_sum += angle;
                // Check that the total number of full turns is plus-or-minus 1
            }
            return Math.Abs(Math.Round(angle_sum / TWO_PI)) == 1;
        }

        /// <summary>
        /// Checks if a polygon is Convex.
        /// Rory Daulton algorithm ported to Java by Guilherme Campos Hazan and ported to C# by Filipe Jorge da Silva Brandão
        /// Precision issues with Rhino implementation led to changes from using integer types to double types.
        /// </summary>
        /// <returns><c>true</c>, if convex is convex, <c>false</c> otherwise.</returns>
        /// <param name="x">The x coordinate.</param>
        /// <param name="y">The y coordinate.</param>
        /// <param name="n">N.</param>
        private static bool IsPolyConvex(double[] x, double[] y, int n)
        {
            const double TWO_PI = 2 * Math.PI;

            //in rhino all lists are zero based so the code could be changed to remove the need to explicitly require the user to state that
            int ibase = 0;
            // points is 'strictly convex': points are valid, side lengths non-zero, interior angles are strictly between zero and a straight
            // angle, and the polygon does not intersect itself.
            // NOTES:  1.  Algorithm: the signed changes of the direction angles from one side to the next side must be all positive or
            // all negative, and their sum must equal plus-or-minus one full turn (2 pi radians). Also check for too few,
            // invalid, or repeated points.
            //      2.  No check is explicitly done for zero internal angles(180 degree direction-change angle) as this is covered
            // in other ways, including the `n < 3` check.
            // needed for any bad points or direction changes
            // Check for too few points
            if (n <= 3)
            {
                return true;
            }

            if (x[ibase] == x[n - 1] && y[ibase] == y[n - 1]) // if its a closed polygon, ignore last vertex
                n--;
            // Get starting information
            double old_x = x[n - 2], old_y = y[n - 2];
            double new_x = x[n - 1], new_y = y[n - 1];
            double new_direction = Math.Atan2(new_y - old_y, new_x - old_x), old_direction;
            double angle_sum = 0.0, orientation = 0;
            // Check each point (the side ending there, its angle) and accum. angles for ndx, newpoint in enumerate(polygon):
            for (int i = 0; i < n; i++)
            {
                // Update point coordinates and side directions, check side length
                old_x = new_x; old_y = new_y; old_direction = new_direction;
                int p = ibase++;
                new_x = x[p]; new_y = y[p];
                new_direction = Math.Atan2(new_y - old_y, new_x - old_x);
                if (old_x == new_x && old_y == new_y)
                    return false; // repeated consecutive points
                                  // Calculate & check the normalized direction-change angle
                double angle = new_direction - old_direction;
                if (angle <= -Math.PI)
                    angle += TWO_PI;  // make it in half-open interval (-Pi, Pi]
                else if (angle > Math.PI)
                    angle -= TWO_PI;
                if (i == 0)  // if first time through loop, initialize orientation
                {
                    if (angle == 0.0) return false;
                    orientation = angle > 0 ? 1 : -1;
                }
                else  // if other time through loop, check orientation is stable
                  if (orientation * angle <= 0)  // not both pos. or both neg.
                    return false;
                // Accumulate the direction-change angle
                angle_sum += angle;
                // Check that the total number of full turns is plus-or-minus 1
            }
            return Math.Abs(Math.Round(angle_sum / TWO_PI)) == 1;
        }
        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return Properties.Resources.IsConvex_Icon;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("893462c5-f763-4ffa-80b4-dd724756c65b"); }
        }
    }
}
